#include "python_ngcore.hpp"
#include "bitarray.hpp"
#include "taskmanager.hpp"
#include "mpi_wrapper.hpp"

using namespace ngcore;
using namespace std;
using namespace pybind11::literals;

namespace pybind11 { namespace detail {
}} // namespace pybind11::detail




PYBIND11_MODULE(pyngcore, m) // NOLINT
{
  try
    {
      auto numpy = py::module::import("numpy");
      ngcore_have_numpy = !numpy.is_none();
    }
  catch(...) {}
  ExportArray<int>(m);
  ExportArray<unsigned>(m);
  ExportArray<size_t>(m);
  ExportArray<double>(m);
  ExportArray<float>(m);
  ExportArray<signed short>(m);
  ExportArray<signed char>(m);
  ExportArray<unsigned short>(m);
  ExportArray<unsigned char>(m);

  // Compiler dependent implementation of size_t
  if constexpr(!is_same_v<size_t, uint64_t>)
    ExportArray<uint64_t>(m);

  ExportTable<int>(m);

  #ifdef PARALLEL
  py::class_<NG_MPI_Comm> (m, "_NG_MPI_Comm")
          ;
  m.def("InitMPI", &InitMPI, py::arg("mpi_library_path")=nullopt);
  #endif // PARALLEL

  py::class_<BitArray, shared_ptr<BitArray>> (m, "BitArray")
    .def(py::init([] (size_t n) { return make_shared<BitArray>(n); }),py::arg("n"))
    .def(py::init([] (const BitArray& a) { return make_shared<BitArray>(a); } ), py::arg("ba"))
    .def(py::init([] (const vector<bool> & a)
                  {
                    auto ba = make_shared<BitArray>(a.size());
                    ba->Clear();
                    for (size_t i = 0; i < a.size(); i++)
                      if (a[i]) ba->SetBit(i);
                    return ba;
                  } ), py::arg("vec"))
    .def(NGSPickle<BitArray>())
    .def("__str__", &ToString<BitArray>)
    .def("__len__", &BitArray::Size)
    .def("__getitem__", [] (BitArray & self, int i)
                                         {
                                           if (i < 0) i+=self.Size();
                                           if (i < 0 || i >= self.Size())
                                             throw py::index_error();
                                           return self.Test(i);
                                         }, py::arg("pos"), "Returns bit from given position")
    .def("__setitem__", [] (BitArray & self, int i, bool b)
                                         {
                                           if (i < 0) i+=self.Size();
                                           if (i < 0 || i >= self.Size())
                                             throw py::index_error();
                                           if (b) self.SetBit(i); else self.Clear(i);
                                         }, py::arg("pos"), py::arg("value"), "Clear/Set bit at given position")

    .def("__setitem__", [] (BitArray & self, py::slice inds, bool b)
                                         {
                                           size_t start, step, stop, n;
                                           if (!inds.compute(self.Size(), &start, &stop, &step, &n))
                                             throw py::error_already_set();

                                           if (start == 0 && n == self.Size() && step == 1)
                                             { // base branch
                                               if (b)
                                                 self.Set();
                                               else
                                                 self.Clear();
                                             }
                                           else
                                             {
                                               if (b)
                                                 for (size_t i=0; i<n; i++, start+=step)
                                                   self.SetBit(start);
                                               else
                                                 for (size_t i=0; i<n; i++, start+=step)
                                                   self.Clear(start);
                                             }
                                         }, py::arg("inds"), py::arg("value"), "Clear/Set bit at given positions")

    .def("__setitem__", [] (BitArray & self, py::slice inds, BitArray & ba)
                                         {
                                           size_t start, step, stop, n;
                                           if (!inds.compute(self.Size(), &start, &stop, &step, &n))
                                             throw py::error_already_set();

                                           if (start == 0 && n == self.Size() && step == 1)
                                             {
                                               self = ba;
                                             }
                                           else
                                             {
                                               for (size_t i = 0; i < n; i++, start += step)
                                                 {
                                                   bool b = ba.Test(i);
                                                   if (b)
                                                     self.SetBit(start);
                                                   else
                                                     self.Clear(start);
                                                 }
                                             }
                                         }, py::arg("inds"), py::arg("ba"), "copy BitArray")

    .def("__setitem__", [](BitArray & self,  IntRange range, bool b)
      {
        if (b)
          for (size_t i : range)
            self.SetBit(i);
        else
          for (size_t i : range)
            self.Clear(i);
      }, py::arg("range"), py::arg("value"), "Set value for range of indices" )

    .def("NumSet", &BitArray::NumSet)
    .def("Set", [] (BitArray & self) { self.Set(); }, "Set all bits")
    .def("Set", &BitArray::SetBit, py::arg("i"), "Set bit at given position")
    .def("Clear", [] (BitArray & self) { self.Clear(); }, "Clear all bits")
    .def("Clear", [] (BitArray & self, int i)
                                   {
                                       self.Clear(i);
                                   }, py::arg("i"), "Clear bit at given position")


    .def(py::self | py::self)
    .def(py::self & py::self)
    #ifdef __clang__ 
    // see https://github.com/pybind/pybind11/issues/1893
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wself-assign-overloaded"
    #endif
    .def(py::self |= py::self)
    .def(py::self &= py::self)
    #ifdef __clang__
    #pragma GCC diagnostic pop
    #endif
    .def(~py::self)
    ;

  py::class_<Flags>(m, "Flags")
    .def(py::init<>())
    .def("__str__", &ToString<Flags>)
    .def(py::init([](py::dict kwargs) {
          Flags flags;
          for (auto d : kwargs)
            SetFlag(flags, d.first.cast<string>(), d.second.cast<py::object>());
          return flags;
    }), "Create flags from dict")
    .def(py::init([](py::kwargs kwargs) {
          Flags flags;
          for (auto d : kwargs)
            SetFlag(flags, d.first.cast<string>(), d.second.cast<py::object>());
          return flags;
        }), "Create flags from kwargs")
    .def(py::pickle([] (const Flags& self)
        {
          std::stringstream str;
          self.SaveFlags(str);
          return py::make_tuple(py::cast(str.str()));
        },
        [] (py::tuple state)
        {
          string s = state[0].cast<string>();
          std::stringstream str(s);
          Flags flags;
          flags.LoadFlags(str);
          return flags;
        }
    ))
    .def("Set",[](Flags & self,const py::dict & aflags)->Flags&
    {      
      SetFlag(self, "", aflags);
      return self;
    }, py::arg("aflag"), "Set the flags by given dict")

    .def("Set",[](Flags & self, const char * akey, const py::object & value)->Flags&
    {             
        SetFlag(self, akey, value);
        return self;
    }, py::arg("akey"), py::arg("value"), "Set flag by given value.")

    .def("keys", [](Flags & self) -> py::list {
        return CreateDictFromFlags(self).attr("keys")();
    })
    .def("__getitem__", [](Flags & self, const string& name) -> py::object {

	  if(self.NumListFlagDefined(name))
	    return py::cast(self.GetNumListFlag(name));

	  if(self.StringListFlagDefined(name))
	    return py::cast(self.GetStringListFlag(name));
	 
	  if(self.NumFlagDefined(name))
	    return py::cast(*self.GetNumFlagPtr(name));
	  
	  if(self.StringFlagDefined(name))
	    return py::cast(self.GetStringFlag(name));

	  if(self.FlagsFlagDefined(name))
	    return py::cast(self.GetFlagsFlag(name));

	  return py::cast(self.GetDefineFlag(name));
      }, py::arg("name"), "Return flag by given name")
    .def("ToDict", [](const Flags& flags)
    {
      return CreateDictFromFlags(flags);
    })
    .def("items", [](const Flags& flags)
    {
      return CreateDictFromFlags(flags).attr("items")();
    })
  ;
  py::implicitly_convertible<py::dict, Flags>();

  
  py::enum_<level::level_enum>(m, "LOG_LEVEL", "Logging level")
    .value("Trace", level::trace)
    .value("Debug", level::debug)
    .value("Info", level::info)
    .value("Warn", level::warn)
    .value("Error", level::err)
    .value("Critical", level::critical)
    .value("Off", level::off);

  m.def("SetLoggingLevel", &SetLoggingLevel, py::arg("level"), py::arg("logger")="",
        "Set logging level, if name is given only to the specific logger, else set the global logging level");
  m.def("AddFileSink", &AddFileSink, py::arg("filename"), py::arg("level"), py::arg("logger")="",
        "Add File sink, either only to logger specified or globally to all loggers");
  m.def("AddConsoleSink", &AddConsoleSink, py::arg("level"), py::arg("logger")="",
        "Add console output for specific logger or all if none given");
  m.def("ClearLoggingSinks", &ClearLoggingSinks, py::arg("logger")="",
        "Clear sinks of specific logger, or all if none given");
  m.def("FlushOnLoggingLevel", &FlushOnLoggingLevel, py::arg("level"), py::arg("logger")="",
        "Flush every message with level at least `level` for specific logger or all loggers if none given.");

  m.def("RunWithTaskManager",
          [](py::object lam)
                           {
                             GetLogger("TaskManager")->info("running Python function with task-manager");
                             RunWithTaskManager ([&] () { lam(); });
                           }, py::arg("lam"), R"raw_string(
Parameters:

lam : object
  input function

)raw_string")
          ;

  m.def("SetNumThreads", &TaskManager::SetNumThreads, py::arg("threads"), R"raw_string(
Set number of threads

Parameters:

threads : int
  input number of threads

)raw_string");

  // local TaskManager class to be used as context manager in Python
  class ParallelContextManager {
      int num_threads;
    public:
      ParallelContextManager() : num_threads(0) {
        TaskManager::SetPajeTrace(0);
        PajeTrace::SetMaxTracefileSize(0);
      };
      ParallelContextManager(size_t pajesize) : num_threads(0) {
        TaskManager::SetPajeTrace(pajesize > 0);
        PajeTrace::SetMaxTracefileSize(pajesize);
      }
      void Enter() {num_threads = EnterTaskManager(); }
      void Exit(py::object exc_type, py::object exc_value, py::object traceback) {
          ExitTaskManager(num_threads);
      }
    };

  py::class_<ParallelContextManager>(m, "TaskManager")
    .def(py::init<>())
    .def(py::init<size_t>(), "pajetrace"_a, "Run paje-tracer, specify buffersize in bytes")
    .def("__enter__", &ParallelContextManager::Enter)
    .def("__exit__", &ParallelContextManager::Exit)
    .def("__timing__", &TaskManager::Timing)
    ;

  py::class_<PajeTrace>(m, "PajeTrace")
    .def(py::init( [] (string filename, size_t size_mb, bool threads, bool thread_counter, bool memory)
          {
              PajeTrace::SetMaxTracefileSize(size_mb*1014*1024);
              PajeTrace::SetTraceThreads(threads);
              PajeTrace::SetTraceMemory(memory);
              PajeTrace::SetTraceThreadCounter(thread_counter);
              trace = new PajeTrace(TaskManager::GetMaxThreads(), filename);
              return trace;
          }), py::arg("filename")="ng.trace", py::arg("size")=1000,
              py::arg("threads")=true, py::arg("thread_counter")=false,
              py::arg("memory")=true,
              "size in Megabytes"
        )
    .def("__enter__", [](PajeTrace & self) { })
    .def("__exit__", [](PajeTrace & self, py::args) { trace = nullptr; })
    .def_static("SetTraceThreads", &PajeTrace::SetTraceThreads)
    .def_static("SetTraceThreadCounter", &PajeTrace::SetTraceThreadCounter)
    .def_static("SetMaxTracefileSize", &PajeTrace::SetMaxTracefileSize)
#ifdef NETGEN_TRACE_MEMORY
    .def_static("WriteMemoryChart", [](string filename){ if(trace) trace->WriteMemoryChart(filename); }, py::arg("filename")="memory" )
#endif // NETGEN_TRACE_MEMORY
    ;

    m.def("GetTotalMemory", MemoryTracer::GetTotalMemory);

    py::class_<Timer<>> (m, "Timer")
    .def(py::init<const string&>())
    .def("Start", static_cast<void (Timer<>::*)()const>(&Timer<>::Start), "start timer")
    .def("Stop", static_cast<void (Timer<>::*)()const>(&Timer<>::Stop), "stop timer")
    .def_property_readonly("time", &Timer<>::GetTime, "returns time")
    .def("__enter__", static_cast<void (Timer<>::*)()const>(&Timer<>::Start))
    .def("__exit__", [](Timer<>& t, py::object, py::object, py::object)
                     {
                       t.Stop();
                     })
    ;
  
  m.def("Timers",
	  []() 
	   {
	     py::list timers;
	     for (int i = 0; i < NgProfiler::SIZE; i++)
	       if (!NgProfiler::timers[i].name.empty())
               {
                 py::dict timer;
                 timer["name"] = py::str(NgProfiler::timers[i].name);
                 timer["time"] = py::float_(NgProfiler::GetTime(i));
                 timer["counts"] = py::int_(NgProfiler::GetCounts(i));
                 timer["flops"] = py::float_(NgProfiler::GetFlops(i));
                 timer["Gflop/s"] = py::float_(NgProfiler::GetFlops(i)/NgProfiler::GetTime(i)*1e-9);
                 timers.append(timer);
               }
	     return timers;
	   }, "Returns list of timers"
	   );
  m.def("ResetTimers", &NgProfiler::Reset);

  py::class_<NgMPI_Comm> (m, "MPI_Comm")
#ifdef PARALLEL
    .def(py::init([] (mpi4py_comm comm) { return NgMPI_Comm(comm); }))
    .def("WTime", [](NgMPI_Comm  & c) { return NG_MPI_Wtime(); })
    .def_property_readonly ("mpi4py", [](NgMPI_Comm & self) { return NG_MPI_CommToMPI4Py(self); })
#endif  // PARALLEL
    .def_property_readonly ("rank", &NgMPI_Comm::Rank)
    .def_property_readonly ("size", &NgMPI_Comm::Size)
    .def("Barrier", &NgMPI_Comm::Barrier)
    .def("Sum", [](NgMPI_Comm  & c, double x) { return c.AllReduce(x, NG_MPI_SUM); })
    .def("Min", [](NgMPI_Comm  & c, double x) { return c.AllReduce(x, NG_MPI_MIN); })
    .def("Max", [](NgMPI_Comm  & c, double x) { return c.AllReduce(x, NG_MPI_MAX); })
    .def("Sum", [](NgMPI_Comm  & c, int x) { return c.AllReduce(x, NG_MPI_SUM); })
    .def("Min", [](NgMPI_Comm  & c, int x) { return c.AllReduce(x, NG_MPI_MIN); })
    .def("Max", [](NgMPI_Comm  & c, int x) { return c.AllReduce(x, NG_MPI_MAX); })
    .def("Sum", [](NgMPI_Comm  & c, size_t x) { return c.AllReduce(x, NG_MPI_SUM); })
    .def("Min", [](NgMPI_Comm  & c, size_t x) { return c.AllReduce(x, NG_MPI_MIN); })
    .def("Max", [](NgMPI_Comm  & c, size_t x) { return c.AllReduce(x, NG_MPI_MAX); })
    .def("SubComm", [](NgMPI_Comm & c, std::vector<int> proc_list) {
        Array<int> procs(proc_list.size());
        for (int i = 0; i < procs.Size(); i++)
          { procs[i] = proc_list[i]; }
        if (!procs.Contains(c.Rank()))
          { throw Exception("rank "+ToString(c.Rank())+" not in subcomm"); }
	return c.SubCommunicator(procs);
      }, py::arg("procs"));
  ;

    
#ifdef PARALLEL
  py::implicitly_convertible<mpi4py_comm, NgMPI_Comm>();
#endif // PARALLEL
}

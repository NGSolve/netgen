
#include "python_ngcore.hpp"

using namespace ngcore;
using namespace std;

PYBIND11_MODULE(pyngcore, m) // NOLINT
{
  ExportArray<int>(m);
  ExportArray<unsigned>(m);
  ExportArray<size_t>(m);
  ExportArray<double>(m);

  py::class_<Flags>(m, "Flags")
    .def(py::init<>())
    .def("__str__", &ToString<Flags>)
    .def(py::init([](py::object & obj) {
          Flags flags;
          py::dict d(obj);          
          SetFlag (flags, "", d);
          return flags;
        }), py::arg("obj"), "Create Flags by given object")
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
}

#ifndef NETGEN_CORE_PYTHON_NGCORE_HPP
#define NETGEN_CORE_PYTHON_NGCORE_HPP

#include "ngcore_api.hpp" // for operator new
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "array.hpp"
#include "table.hpp"
#include "archive.hpp"
#include "flags.hpp"
#include "ngcore_api.hpp"
#include "ng_mpi.hpp"

namespace py = pybind11;

namespace ngcore
{
#ifdef PARALLEL
  NGCORE_API extern bool (*NG_MPI_CommFromMPI4Py)(py::handle, NG_MPI_Comm &);
  NGCORE_API extern py::handle (*NG_MPI_CommToMPI4Py)(NG_MPI_Comm);
#endif // PARALLEL

  namespace detail
  {
    template<typename T>
    struct HasPyFormat
    {
    private:
      template<typename T2>
      static auto check(T2*) -> std::enable_if_t<std::is_same_v<decltype(std::declval<py::format_descriptor<T2>>().format()), std::string>, std::true_type>;
      static auto check(...) -> std::false_type;
    public:
      static constexpr bool value = decltype(check((T*) nullptr))::value;
    };
  } // namespace detail

#ifdef PARALLEL
  struct mpi4py_comm {
    mpi4py_comm() = default;
    mpi4py_comm(NG_MPI_Comm value) : value(value) {}
    operator NG_MPI_Comm () { return value; }

    NG_MPI_Comm value;
  };
#endif  // PARALLEL
} // namespace ngcore


////////////////////////////////////////////////////////////////////////////////
// automatic conversion of python list to Array<>
namespace pybind11 {
namespace detail {

#ifdef PARALLEL
template <> struct type_caster<ngcore::mpi4py_comm> {
  public:
  PYBIND11_TYPE_CASTER(ngcore::mpi4py_comm, _("mpi4py_comm"));

    // Python -> C++
    bool load(handle src, bool) {
      return ngcore::NG_MPI_CommFromMPI4Py(src, value.value);
    }

    // C++ -> Python
    static handle cast(ngcore::mpi4py_comm src,
                       return_value_policy /* policy */,
                       handle /* parent */)
    {
      // Create an mpi4py handle
      return ngcore::NG_MPI_CommToMPI4Py(src.value);
    }
};
#endif //  PARALLEL

template <typename Type, typename Value> struct ngcore_list_caster {
    using value_conv = make_caster<Value>;

    bool load(handle src, bool convert) {
        if (!isinstance<sequence>(src) || isinstance<str>(src))
            return false;
        auto s = reinterpret_borrow<sequence>(src);
        value.SetSize(s.size());
        value.SetSize0();
        for (auto it : s) {
            value_conv conv;
            if (!conv.load(it, convert))
                return false;
            value.Append(cast_op<Value &&>(std::move(conv)));
        }
        return true;
    }

public:
    template <typename T>
    static handle cast(T &&src, return_value_policy policy, handle parent) {
        if (!std::is_lvalue_reference<T>::value)
            policy = return_value_policy_override<Value>::policy(policy);
        list l(src.Size());
        size_t index = 0;
        for (auto &&value : src) {
            auto value_ = reinterpret_steal<object>(value_conv::cast(forward_like<T>(value), policy, parent));
            if (!value_)
                return handle();
            PyList_SET_ITEM(l.ptr(), (ssize_t) index++, value_.release().ptr()); // steals a reference
        }
        return l.release();
    }

    PYBIND11_TYPE_CASTER(Type, _("Array[") + value_conv::name + _("]"));
};


template <typename Type> struct type_caster<ngcore::Array<Type>, enable_if_t<!ngcore::detail::HasPyFormat<Type>::value>>
 : ngcore_list_caster<ngcore::Array<Type>, Type> { };


  /*
  template <typename Type> struct type_caster<std::shared_ptr<ngcore::Table<Type>>>
  {
    template <typename T>
    static handle cast(T &&src, return_value_policy policy, handle parent)
    {
      std::cout << "handle called with type src = " << typeid(src).name() << std::endl;

      return handle(); // what so ever
    }
    
    PYBIND11_TYPE_CASTER(Type, _("Table[") + make_caster<Type>::name + _("]"));
  };
  */
  
  

} // namespace detail
} // namespace pybind11
////////////////////////////////////////////////////////////////////////////////

namespace ngcore
{
  NGCORE_API extern bool ngcore_have_numpy;
  NGCORE_API extern bool parallel_pickling;
  
  // Python class name type traits
  template <typename T>
  struct PyNameTraits {
    static const std::string & GetName()
    {
      static const std::string name = typeid(T).name();
      return name;
    }
  };

  template <typename T>
  std::string GetPyName(const char *prefix = 0) {
    std::string s;
    if(prefix) s = std::string(prefix);
    s+= PyNameTraits<T>::GetName();
    return s;
  }

  template<>
  struct PyNameTraits<int> {
    static std::string GetName() { return "I"; }
  };

  template<>
  struct PyNameTraits<unsigned> {
    static std::string GetName() { return "U"; }
  };

  template<>
  struct PyNameTraits<float> {
    static std::string GetName() { return "F"; }
  };

  template<>
  struct PyNameTraits<double> {
    static std::string GetName() { return "D"; }
  };

  template<>
  struct PyNameTraits<size_t> {
    static std::string GetName() { return "S"; }
  };

  template<typename T>
  struct PyNameTraits<std::shared_ptr<T>> {
    static std::string GetName()
    { return std::string("sp_")+GetPyName<T>(); }
  };

  template<typename ARCHIVE>
  class NGCORE_API_EXPORT PyArchive : public ARCHIVE
  {
  private:
    pybind11::list lst;
    size_t index = 0;
    std::map<std::string, VersionInfo> version_needed;
  protected:
    using ARCHIVE::stream;
    using ARCHIVE::version_map;
  public:
    PyArchive(const pybind11::object& alst = pybind11::none()) :
      ARCHIVE(std::make_shared<std::stringstream>()),
      lst(alst.is_none() ? pybind11::list() : pybind11::cast<pybind11::list>(alst))
    {
      ARCHIVE::shallow_to_python = true;
      if(Input())
        {
          stream = std::make_shared<std::stringstream>
            (pybind11::cast<pybind11::bytes>(lst[pybind11::len(lst)-1]));
          *this & version_needed;
          for(auto& libversion : version_needed)
            if(libversion.second > GetLibraryVersion(libversion.first))
              throw Exception("Error in unpickling data:\nLibrary " + libversion.first +
                              " must be at least " + libversion.second.to_string());
          stream = std::make_shared<std::stringstream>
            (pybind11::cast<pybind11::bytes>(lst[pybind11::len(lst)-2]));
          *this & version_map;
          stream = std::make_shared<std::stringstream>
            (pybind11::cast<pybind11::bytes>(lst[pybind11::len(lst)-3]));
        }
    }

    void NeedsVersion(const std::string& library, const std::string& version) override
    {
      if(Output())
        {
          version_needed[library] = version_needed[library] > version ? version_needed[library] : version;
        }
    }

    using ARCHIVE::Output;
    using ARCHIVE::Input;
    using ARCHIVE::FlushBuffer;
    using ARCHIVE::operator&;
    using ARCHIVE::operator<<;
    using ARCHIVE::GetVersion;
    void ShallowOutPython(const pybind11::object& val) override { lst.append(val); }
    void ShallowInPython(pybind11::object& val) override { val = lst[index++]; }

    pybind11::list WriteOut()
    {
      auto version_runtime = GetLibraryVersions();
      FlushBuffer();
      lst.append(pybind11::bytes(std::static_pointer_cast<std::stringstream>(stream)->str()));
      stream = std::make_shared<std::stringstream>();
      *this & version_runtime;
      FlushBuffer();
      lst.append(pybind11::bytes(std::static_pointer_cast<std::stringstream>(stream)->str()));
      stream = std::make_shared<std::stringstream>();
      *this & version_needed;
      FlushBuffer();
      lst.append(pybind11::bytes(std::static_pointer_cast<std::stringstream>(stream)->str()));
      return lst;
    }
  };

  template<typename T, typename T_ARCHIVE_OUT=BinaryOutArchive, typename T_ARCHIVE_IN=BinaryInArchive>
  auto NGSPickle()
  {
    return pybind11::pickle([](T* self)
                      {
                        PyArchive<T_ARCHIVE_OUT> ar;
                        ar.SetParallel(parallel_pickling);
                        ar & self;
                        auto output = pybind11::make_tuple(ar.WriteOut());
                        return output;
                      },
                      [](const pybind11::tuple & state)
                      {
                        T* val = nullptr;
                        PyArchive<T_ARCHIVE_IN> ar(state[0]);
                        ar & val;
                        return val;
                      });
  }

  template<typename T>
  Array<T> makeCArray(const py::object& obj)
  {
    Array<T> arr;
    if(py::isinstance<py::list>(obj))
        for(auto& val : py::cast<py::list>(obj))
          arr.Append(py::cast<T>(val));
    else if(py::isinstance<py::tuple>(obj))
      for(auto& val : py::cast<py::tuple>(obj))
        arr.Append(py::cast<T>(val));
    else
      throw py::type_error("Cannot convert Python object to C Array");
    return arr;
  }

  template <typename T>
  // py::object makePyTuple (FlatArray<T> ar)
  py::object makePyTuple (const BaseArrayObject<T> & ar)
  {
    py::tuple res(ar.Size());
    for (auto i : Range(ar))
      res[i] = py::cast(ar[i]);
    return res;
  }

  template <typename T, typename TIND=typename FlatArray<T>::index_type>
  void ExportArray (py::module &m)
  {
      using TFlat = FlatArray<T, TIND>;
      using TArray = Array<T, TIND>;
      std::string suffix = GetPyName<T>() + "_" +
        GetPyName<TIND>();
      std::string fname = std::string("FlatArray_") + suffix;
      auto flatarray_class = py::class_<TFlat>(m, fname.c_str(),
                                               py::buffer_protocol())
        .def ("__len__", [] ( TFlat &self ) { return self.Size(); } )
        .def ("__getitem__",
              [](TFlat & self, TIND i) -> T&
                             {
                               static constexpr int base = IndexBASE<TIND>();
                               if (i < base || i >= self.Size()+base)
                                 throw py::index_error();
                               return self[i]; 
                             },
              py::return_value_policy::reference)
        .def ("__setitem__",
              [](TFlat & self, TIND i, T val) -> T&
                             {
                               static constexpr int base = IndexBASE<TIND>();
                               if (i < base || i >= self.Size()+base)
                                 throw py::index_error();
                               self[i] = val;
                               return self[i];
                             },
              py::return_value_policy::reference)

        .def ("__setitem__",
              [](TFlat & self, py::slice slice, T val)
              {
                size_t start, stop, step, slicelength;
                if (!slice.compute(self.Size(), &start, &stop, &step, &slicelength))
                  throw py::error_already_set();
                static constexpr int base = IndexBASE<TIND>();
                if (start < base || start+(slicelength-1)*step >= self.Size()+base)
                  throw py::index_error();
                for (size_t i = 0; i < slicelength; i++, start+=step)
                  self[start] = val;
              })

        .def("__iter__", [] ( TFlat & self) {
             return py::make_iterator (self.begin(),self.end());
             }, py::keep_alive<0,1>()) // keep array alive while iterator is used

        .def("__str__", [](TFlat& self)
                        {
                          return ToString(self);
                        })
      ;

      if constexpr (detail::HasPyFormat<T>::value)
        {
          if(ngcore_have_numpy && !py::detail::npy_format_descriptor<T>::dtype().is_none())
            {
              flatarray_class
                .def_buffer([](TFlat& self)
                            {
                              return py::buffer_info(
                                self.Addr(0),
                                sizeof(T),
                                py::format_descriptor<T>::format(),
                                1,
                                { self.Size() },
                                { sizeof(T) * (self.Addr(1) - self.Addr(0)) });
                            })
                .def("NumPy", [](py::object self)
                              {
                                return py::module::import("numpy")
                                  .attr("frombuffer")(self, py::detail::npy_format_descriptor<T>::dtype());
                              })
                ;
              }
          }

      std::string aname = std::string("Array_") + suffix;
      auto arr = py::class_<TArray, TFlat> (m, aname.c_str())
        .def(py::init([] (size_t n) { return new TArray(n); }),py::arg("n"), "Makes array of given length")
        .def(py::init([] (std::vector<T> const & x)
                  {
                    size_t s = x.size();
                    TArray tmp(s);
                    for (size_t i : Range(tmp))
                      tmp[TIND(i)] = x[i];
                    return tmp;
                  }), py::arg("vec"), "Makes array with given list of elements")
        ;
      if constexpr(is_archivable<TArray>)
        arr.def(NGSPickle<TArray>());
      py::implicitly_convertible<std::vector<T>, TArray>();
    }

  template <typename T>
  void ExportTable (py::module &m)
  {
    py::class_<ngcore::Table<T>, std::shared_ptr<ngcore::Table<T>>> (m, ("Table_"+GetPyName<T>()).c_str())
      .def(py::init([] (py::list blocks)
                    {
                       size_t size = py::len(blocks);
                       Array<int> cnt(size);
                       size_t i = 0;
                       for (auto block : blocks)
                         cnt[i++] = py::len(block);
                       
                       i = 0;
                       Table<T> blocktable(cnt);
                       for (auto block : blocks)
                         {
                           auto row = blocktable[i++];
                           size_t j = 0;
                           for (auto val : block)
                             row[j++] = val.cast<T>();
                         }
                       // cout << "blocktable = " << *blocktable << endl;
                       return blocktable;
                      
                    }), py::arg("blocks"), "a list of lists")

      .def ("__len__", [] (Table<T> &self ) { return self.Size(); } )
      .def ("__getitem__",
            [](Table<T> & self, size_t i) -> FlatArray<T>
            {
              if (i >= self.Size())
                throw py::index_error();
              return self[i]; 
            })
      .def("__str__", [](Table<T> & self)
           {
             return ToString(self);
           })
      ;
  }

  
  void NGCORE_API SetFlag(Flags &flags, std::string s, py::object value);
  // Parse python kwargs to flags
  Flags NGCORE_API CreateFlagsFromKwArgs(const py::kwargs& kwargs, py::object pyclass = py::none(),
                                         py::list info = py::list());
  // Create python dict from kwargs
  py::dict NGCORE_API CreateDictFromFlags(const Flags& flags);


} // namespace ngcore

#endif // NETGEN_CORE_PYTHON_NGCORE_HPP

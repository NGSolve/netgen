#ifndef NETGEN_CORE_PYTHON_NGCORE_HPP
#define NETGEN_CORE_PYTHON_NGCORE_HPP

#include <pybind11/pybind11.h>

#include "array.hpp"
#include "archive.hpp"
#include "flags.hpp"
#include "ngcore_api.hpp"

namespace ngcore
{
  namespace py = pybind11;

  template<typename T>
  Array<T> makeCArray(const py::object& obj)
  {
    Array<T> arr;
    arr.SetAllocSize(py::len(obj));
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
  
  void NGCORE_API SetFlag(Flags &flags, std::string s, py::object value);
  // Parse python kwargs to flags
  Flags NGCORE_API CreateFlagsFromKwArgs(py::object pyclass, const py::kwargs& kwargs, py::list info = py::list());

  // ***************  Archiving functionality  **************
  
    template<typename T>
    Archive& Archive :: Shallow(T& val)
    {
      static_assert(detail::is_any_pointer<T>, "ShallowArchive must be given pointer type!");
#ifdef NETGEN_PYTHON
      if(shallow_to_python)
        {
          if(is_output)
            ShallowOutPython(pybind11::cast(val));
          else
          {
            pybind11::object obj;
            ShallowInPython(obj);
            val = pybind11::cast<T>(obj);
          }
        }
      else
#endif // NETGEN_PYTHON
        *this & val;
      return *this;
    }

  template<typename ARCHIVE>
  class PyArchive : public ARCHIVE
  {
  private:
    pybind11::list lst;
    size_t index = 0;
    std::map<std::string, VersionInfo> version_needed;
  protected:
    using ARCHIVE::stream;
    using ARCHIVE::version_map;
    using ARCHIVE::logger;
    using ARCHIVE::GetLibraryVersions;
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
          logger->debug("versions needed for unpickling = {}", version_needed);
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
          logger->debug("Need version {} of library {}.", version, library);
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
      FlushBuffer();
      lst.append(pybind11::bytes(std::static_pointer_cast<std::stringstream>(stream)->str()));
      stream = std::make_shared<std::stringstream>();
      *this & GetLibraryVersions();
      FlushBuffer();
      lst.append(pybind11::bytes(std::static_pointer_cast<std::stringstream>(stream)->str()));
      stream = std::make_shared<std::stringstream>();
      logger->debug("Writeout version needed = {}", version_needed);
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
                        ar & self;
                        auto output = pybind11::make_tuple(ar.WriteOut());
                        GetLogger("Archive")->trace("Pickling output for object of type {} = {}",
                                                    Demangle(typeid(T).name()),
                                                    std::string(pybind11::str(output)));
                        return output;
                      },
                      [](const pybind11::tuple & state)
                      {
                        T* val = nullptr;
                        GetLogger("Archive")->trace("State for unpickling of object of type {} = {}",
                                                    Demangle(typeid(T).name()),
                                                    std::string(pybind11::str(state[0])));
                        PyArchive<T_ARCHIVE_IN> ar(state[0]);
                        ar & val;
                        return val;
                      });
  }


} // namespace ngcore

#endif // NETGEN_CORE_PYTHON_NGCORE_HPP

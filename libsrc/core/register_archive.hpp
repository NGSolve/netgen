#ifndef NETGEN_REGISTER_ARCHIVE_HPP
#define NETGEN_REGISTER_ARCHIVE_HPP

#ifdef NETGEN_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/cast.h>
#endif // NETGEN_PYTHON
#include <tuple>

#include "archive.hpp"

namespace ngcore {
  // ***************  Archiving functionality  **************

#ifdef NETGEN_PYTHON
  template<typename T>
  Archive& Archive :: Shallow(T& val)
  {
    static_assert(detail::is_any_pointer<T>, "ShallowArchive must be given pointer type!");
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
      *this & val;
    return *this;
  }
#endif // NETGEN_PYTHON


  template<typename T, typename Bases=std::tuple<>, typename CArgs=std::tuple<>>
  class RegisterClassForArchive
  {
  public:
    std::function<CArgs(T&)> get_cargs;
    RegisterClassForArchive(std::function<CArgs(T&)> _get_cargs =
                            [](T&) -> std::tuple<> { return std::tuple<>{}; }) :
      get_cargs(_get_cargs)
    {
      static_assert(std::is_base_of<Bases, T>::value ||
                    detail::is_base_of_tuple<T, Bases>,
                    "Second argument must be base class or tuple of base classes of T");
      detail::ClassArchiveInfo info {};
      info.creator = [](const std::type_info& ti, Archive& ar) -> void*
      {
        CArgs args;
        ar &args;
        auto nT = detail::constructIfPossible<T>(args);
        return typeid(T) == ti ? nT
          : Archive::Caster<T, Bases>::tryUpcast(ti, nT);
      };
      info.upcaster = [/*this*/](const std::type_info& ti, void* p) -> void*
      { return typeid(T) == ti ? p : Archive::Caster<T, Bases>::tryUpcast(ti, static_cast<T*>(p)); };
      info.downcaster = [/*this*/](const std::type_info& ti, void* p) -> void*
      { return typeid(T) == ti ? p : Archive::Caster<T, Bases>::tryDowncast(ti, p); };
      info.cargs_archiver = [this](Archive &ar, void* p) {
        ar << get_cargs(*static_cast<T*>(p));
      };
#ifdef NETGEN_PYTHON
    info.anyToPyCaster = [](const std::any &a) {
      const T* val = std::any_cast<T>(&a);
      return pybind11::cast(val);
    };
#endif // NETGEN_PYTHON
    Archive::SetArchiveRegister(std::string(Demangle(typeid(T).name())),info);
  }
};
} // namespace ngcore
#endif // NETGEN_REGISTER_ARCHIVE_HPP

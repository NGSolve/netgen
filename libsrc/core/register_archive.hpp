#ifndef NETGEN_REGISTER_ARCHIVE_HPP
#define NETGEN_REGISTER_ARCHIVE_HPP

#ifdef NETGEN_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/cast.h>
#endif // NETGEN_PYTHON

#include "archive.hpp"

namespace ngcore {

  template<typename T, typename ... Bases>
  class RegisterClassForArchive
  {
  public:
    RegisterClassForArchive()
    {
      static_assert(detail::all_of_tmpl<std::is_base_of<Bases,T>::value...>,
                    "Variadic template arguments must be base classes of T");
      detail::ClassArchiveInfo info {};
      info.creator = [](const std::type_info& ti) -> void*
                     { return typeid(T) == ti ? detail::constructIfPossible<T>()
                         : Archive::Caster<T, Bases...>::tryUpcast(ti, detail::constructIfPossible<T>()); };
      info.upcaster = [/*this*/](const std::type_info& ti, void* p) -> void*
                      { return typeid(T) == ti ? p : Archive::Caster<T, Bases...>::tryUpcast(ti, static_cast<T*>(p)); };
      info.downcaster = [/*this*/](const std::type_info& ti, void* p) -> void*
                        { return typeid(T) == ti ? p : Archive::Caster<T, Bases...>::tryDowncast(ti, p); };
#ifdef NETGEN_PYTHON
      info.anyToPyCaster = [](const std::any& a)
      {
        const T* val = std::any_cast<T>(&a);
        return pybind11::cast(val); };
#endif // NETGEN_PYTHON
      Archive::SetArchiveRegister(std::string(Demangle(typeid(T).name())),info);
    }
  };
} // namespace ngcore
#endif // NETGEN_REGISTER_ARCHIVE_HPP

#ifndef NETGEN_REGISTER_ARCHIVE_HPP
#define NETGEN_REGISTER_ARCHIVE_HPP

#ifdef NETGEN_PYTHON
#include <memory>
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

  /*
    // now using has_shared_from_this2 in archive.hpp
  template <typename T>
  struct has_shared_from_this
  {
    template <typename C> static std::true_type check( decltype( sizeof(&C::shared_from_this )) ) { return std::true_type(); }
    template <typename> static std::false_type check(...) { return std::false_type(); }
    typedef decltype( check<T>(sizeof(char)) ) type;
    static constexpr type value = type();
  };
  */
#endif // NETGEN_PYTHON


  template<typename T, typename Bases=std::tuple<>>
  class RegisterClassForArchive
  {
  public:
    RegisterClassForArchive()
    {
      static_assert(std::is_base_of<Bases, T>::value ||
                    detail::is_base_of_tuple<T, Bases>,
                    "Second argument must be base class or tuple of base classes of T");
      detail::ClassArchiveInfo info {};
      info.creator = [](const std::type_info& ti, Archive& ar) -> void*
      {
        detail::TCargs<T> args;
        ar &args;
        auto nT = detail::constructIfPossible<T>(std::move(args));
        return typeid(T) == ti ? nT
          : Archive::Caster<T, Bases>::tryUpcast(ti, nT);
      };
      info.upcaster = [](const std::type_info& ti, void* p) -> void*
      { return typeid(T) == ti ? p : Archive::Caster<T, Bases>::tryUpcast(ti, static_cast<T*>(p)); };
      info.downcaster = [](const std::type_info& ti, void* p) -> void*
      { return typeid(T) == ti ? p : Archive::Caster<T, Bases>::tryDowncast(ti, p); };
      info.cargs_archiver = [](Archive &ar, void* p) {
        if constexpr(detail::has_GetCArgs_v<T>)
          ar << static_cast<T*>(p)->GetCArgs();
      };
#ifdef NETGEN_PYTHON
    info.anyToPyCaster = [](const std::any &a) {
      if constexpr(has_shared_from_this2<T>::value) {
        std::shared_ptr<T> val = std::any_cast<std::shared_ptr<T>>(a);
        return pybind11::cast(val);
      } else {
        const T* val = std::any_cast<T>(&a);
        return pybind11::cast(val);
      }
    };
#endif // NETGEN_PYTHON
    Archive::SetArchiveRegister(std::string(Demangle(typeid(T).name())),info);
  }
};
} // namespace ngcore
#endif // NETGEN_REGISTER_ARCHIVE_HPP

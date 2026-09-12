#ifndef NETGEN_REGISTER_ARCHIVE_HPP
#define NETGEN_REGISTER_ARCHIVE_HPP

#include <tuple>

#include "archive.hpp"

namespace ngcore {
  // ***************  Archiving functionality  **************

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
    Archive::SetArchiveRegister(std::string(Demangle(typeid(T).name())),info);
  }
};
} // namespace ngcore
#endif // NETGEN_REGISTER_ARCHIVE_HPP

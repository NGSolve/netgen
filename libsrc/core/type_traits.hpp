#ifndef NETGEN_CORE_TYPE_TRAITS_HPP
#define NETGEN_CORE_TYPE_TRAITS_HPP

#include <type_traits>

namespace ngcore
{
  namespace detail
  {
    template<bool... b> struct _BoolArray{};

    template<bool ... vals>
    constexpr bool all_of_tmpl = std::is_same<_BoolArray<vals...>, _BoolArray<(vals || true)...>>::value; // NOLINT
  } // namespace detail
} // namespace ngcore

#endif // NETGEN_CORE_TYPE_TRAITS_HPP

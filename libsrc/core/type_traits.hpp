#ifndef NETGEN_CORE_TYPE_TRAITS_HPP
#define NETGEN_CORE_TYPE_TRAITS_HPP

#include <memory>
#include <type_traits>

namespace ngcore
{
  namespace detail
  {
    template<bool... b> struct _BoolArray{};

    template<bool ... vals>
    constexpr bool all_of_tmpl = std::is_same<_BoolArray<vals...>, _BoolArray<(vals || true)...>>::value; // NOLINT

    template<typename T>
    struct is_any_pointer_impl : std::false_type {};

    template<typename T>
    struct is_any_pointer_impl<T*> : std::true_type {};

    template<typename T>
    struct is_any_pointer_impl<std::shared_ptr<T>> : std::true_type {};

    template<typename T>
    struct is_any_pointer_impl<std::unique_ptr<T>> : std::true_type {};

    template<typename T>
    constexpr bool is_any_pointer = is_any_pointer_impl<T>::value;
  } // namespace detail
} // namespace ngcore

#endif // NETGEN_CORE_TYPE_TRAITS_HPP

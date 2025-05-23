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

    // check if second template argument is tuple of base classes to first
    // template argument, return constexpr bool
    template<typename T, typename Tuple>
    constexpr bool is_base_of_tuple = false;

    template<typename T, typename... Ts>
    constexpr bool is_base_of_tuple<T, std::tuple<Ts...>> =
      all_of_tmpl<std::is_base_of<Ts, T>::value...>;

    template<typename T>
    struct is_any_pointer_impl<T*> : std::true_type {};

    template<typename T>
    struct is_any_pointer_impl<std::shared_ptr<T>> : std::true_type {};

    template<typename T>
    struct is_any_pointer_impl<std::unique_ptr<T>> : std::true_type {};

    template<typename T>
    constexpr bool is_any_pointer = is_any_pointer_impl<T>::value;
  } // namespace detail

  
  // Type trait to check if a class implements a 'range_type Range()' function
  namespace detail
  {
    template<typename T>
    struct has_Range
    {
    private:
      template<typename T2>
      static constexpr auto check(T2*) ->
        std::enable_if_t<!std::is_same_v<decltype(std::declval<T2>().Range()), void>, std::true_type>
      { std::true_type(); }
      template<typename>
      static constexpr std::false_type check(...);
      using type = decltype(check<T>(nullptr)); // NOLINT
    public:
      NGCORE_API static constexpr bool value = type::value;
    };
  }
  template<typename T>
  constexpr bool has_range = detail::has_Range<T>::value;

} // namespace ngcore

#endif // NETGEN_CORE_TYPE_TRAITS_HPP

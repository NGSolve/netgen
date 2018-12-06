

namespace ngcore
{
  template<bool... b> struct _BoolArray{};
  template<bool ... T>
  constexpr bool all_of_tmpl = std::is_same<_BoolArray<T...>, _BoolArray<(T || true)...>>::value;
}

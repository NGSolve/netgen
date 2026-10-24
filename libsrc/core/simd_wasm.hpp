#ifndef NETGEN_CORE_SIMD_WASM_HPP
#define NETGEN_CORE_SIMD_WASM_HPP

/**************************************************************************/
/* File:   simd_wasm.hpp                                                  */
/* WebAssembly SIMD128 (-msimd128): native SIMD<double,2>                 */
/* no fused multiply-add in wasm SIMD: FMA(a,b,c) is c+a*b, as in the     */
/* scalar build, so results are identical                                */
/**************************************************************************/

#include <wasm_simd128.h>

namespace ngcore
{
  typedef double ng_f64x2 __attribute__((vector_size(16)));
  typedef int64_t ng_i64x2 __attribute__((vector_size(16)));

  template <>
  class SIMD<mask64,2>
  {
    ng_i64x2 mask;
  public:
    SIMD (int i)
    {
      mask[0] = i > 0 ? -1 : 0;
      mask[1] = i > 1 ? -1 : 0;
    }
    SIMD (bool i0, bool i1) { mask[0] = i0 ? -1 : 0; mask[1] = i1 ? -1 : 0; }
    SIMD (SIMD<mask64,1> i0, SIMD<mask64,1> i1) { mask[0] = i0[0]; mask[1] = i1[0]; }
    SIMD (ng_i64x2 _data) : mask{_data} { }
    auto Data() const { return mask; }
    static constexpr int Size() { return 2; }
    int64_t operator[] (int i) const { return mask[i]; }
    template <int I>
    int64_t Get() const { return mask[I]; }
    auto Lo() const { return mask[0]; }
    auto Hi() const { return mask[1]; }
  };


  // *************************** int64 ***************************

  template<>
  class SIMD<int64_t,2>
  {
    ng_i64x2 data;
  public:
    static constexpr int Size() { return 2; }
    SIMD() {}
    SIMD (int64_t val) : data{val,val} {}
    SIMD (int64_t v0, int64_t v1) : data{v0,v1} { }
    SIMD (std::array<int64_t, 2> arr) : data{arr[0], arr[1]} { }
    SIMD (ng_i64x2 _data) { data = _data; }

    NETGEN_INLINE auto Data() const { return data; }
    NETGEN_INLINE auto & Data() { return data; }

    int64_t Lo() const { return data[0]; }
    int64_t Hi() const { return data[1]; }

    int64_t operator[] (int i) const { return data[i]; }
    int64_t & operator[] (int i)  { return ((int64_t*)&data)[i]; }

    template <int I>
    int64_t Get() const { return data[I]; }
    static SIMD FirstInt(int n0=0) { return { n0+0, n0+1 }; }
  };

  NETGEN_INLINE SIMD<int64_t,2> operator& (SIMD<int64_t,2> a, SIMD<int64_t,2> b)
  { return a.Data() & b.Data(); }

  NETGEN_INLINE SIMD<int64_t,2> operator+ (SIMD<int64_t,2> a, SIMD<int64_t,2> b)
  { return a.Data() + b.Data(); }

  NETGEN_INLINE SIMD<mask64,2> operator== (SIMD<int64_t,2> a, SIMD<int64_t,2> b)
  { return ng_i64x2(a.Data() == b.Data()); }

  NETGEN_INLINE SIMD<mask64,2> operator> (SIMD<int64_t,2> a, SIMD<int64_t,2> b)
  { return ng_i64x2(a.Data() > b.Data()); }

  template <int N>
  SIMD<int64_t,2> operator<< (SIMD<int64_t,2> a, IC<N> n)
  { return a.Data() << N; }


  // *************************** double ***************************

  template<>
  class SIMD<double,2>
  {
    ng_f64x2 data;

  public:
    static constexpr int Size() { return 2; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD (double v0, double v1) : data{v0,v1} { }
    SIMD (SIMD<double,1> v0, SIMD<double,1> v1) : data{v0.Data(), v1.Data()} { }
    SIMD (std::array<double, 2> arr) : data{arr[0], arr[1]} { }

    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) : data{val,val} { }
    SIMD (int val)    : data{double(val),double(val)} { }
    SIMD (size_t val) : data{double(val),double(val)} { }

    SIMD (double const * p) { data = ng_f64x2(wasm_v128_load(p)); }

    SIMD (double const * p, SIMD<mask64,2> mask)
    {
      data[0] = mask[0] ? p[0] : 0;
      data[1] = mask[1] ? p[1] : 0;
    }
    SIMD (ng_f64x2 _data) { data = _data; }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data[0] = func(0);
      data[1] = func(1);
    }

    void Store (double * p) { wasm_v128_store(p, v128_t(data)); }

    void Store (double * p, SIMD<mask64,2> mask)
    {
      if (mask[0]) p[0] = data[0];
      if (mask[1]) p[1] = data[1];
    }

    NETGEN_INLINE double operator[] (int i) const { return data[i]; }
    NETGEN_INLINE double & operator[] (int i)  { return ((double*)&data)[i]; }

    template <int I>
    double Get() const { return data[I]; }

    NETGEN_INLINE auto Data() const { return data; }
    NETGEN_INLINE auto & Data() { return data; }

    operator std::tuple<double&,double&> ()
    {
      auto pdata = (double*)&data;
      return std::tuple<double&,double&>(pdata[0], pdata[1]);
    }

    double Lo() const { return data[0]; }
    double Hi() const { return data[1]; }
  };


  NETGEN_INLINE double HSum (SIMD<double,2> sd)
  { return sd.Lo()+sd.Hi(); }

  // (a0+a1, b0+b1)
  NETGEN_INLINE SIMD<double,2> HSum (SIMD<double,2> a, SIMD<double,2> b)
  {
    ng_f64x2 lo = __builtin_shufflevector(a.Data(), b.Data(), 0, 2);
    ng_f64x2 hi = __builtin_shufflevector(a.Data(), b.Data(), 1, 3);
    return lo + hi;
  }

  NETGEN_INLINE SIMD<double,4> HSum(SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c, SIMD<double,2> d)
  { return SIMD<double,4> (HSum(a,b), HSum(c,d)); }

  NETGEN_INLINE SIMD<double,2> SwapPairs (SIMD<double,2> a)
  { return __builtin_shufflevector(a.Data(), a.Data(), 1, 0); }

  NETGEN_INLINE SIMD<double,2> operator+ (SIMD<double,2> a, SIMD<double,2> b)
  { return a.Data()+b.Data(); }
  NETGEN_INLINE SIMD<double,2> operator- (SIMD<double,2> a, SIMD<double,2> b)
  { return a.Data()-b.Data(); }
  NETGEN_INLINE SIMD<double,2> operator- (SIMD<double,2> a)
  { return -a.Data(); }
  NETGEN_INLINE SIMD<double,2> operator* (SIMD<double,2> a, SIMD<double,2> b)
  { return a.Data()*b.Data(); }
  NETGEN_INLINE SIMD<double,2> operator/ (SIMD<double,2> a, SIMD<double,2> b)
  { return a.Data()/b.Data(); }

  // a*b+c, not fused
  NETGEN_INLINE SIMD<double,2> FMA (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return c.Data() + a.Data()*b.Data(); }
  NETGEN_INLINE SIMD<double,2> FMA (double a, SIMD<double,2> b, SIMD<double,2> c)
  { return FMA(SIMD<double,2> (a), b, c); }
  // -a*b+c
  NETGEN_INLINE SIMD<double,2> FNMA (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return c.Data() - a.Data()*b.Data(); }
  NETGEN_INLINE SIMD<double,2> FNMA (double a, SIMD<double,2> b, SIMD<double,2> c)
  { return FNMA(SIMD<double,2> (a), b, c); }

  NETGEN_INLINE SIMD<double,2> sqrt (SIMD<double,2> x)
  { return ng_f64x2(wasm_f64x2_sqrt(v128_t(x.Data()))); }

  template <>
  NETGEN_INLINE SIMD<double,2> Reinterpret (SIMD<int64_t,2> a)
  { return ng_f64x2(a.Data()); }

  NETGEN_INLINE SIMD<double,2> If (SIMD<mask64,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return ng_f64x2(wasm_v128_bitselect(v128_t(b.Data()), v128_t(c.Data()), v128_t(a.Data()))); }
  NETGEN_INLINE SIMD<int64_t,2> If (SIMD<mask64,2> a, SIMD<int64_t,2> b, SIMD<int64_t,2> c)
  { return ng_i64x2(wasm_v128_bitselect(v128_t(b.Data()), v128_t(c.Data()), v128_t(a.Data()))); }

  NETGEN_INLINE SIMD<mask64,2> operator&& (SIMD<mask64,2> a, SIMD<mask64,2> b)
  { return a.Data() & b.Data(); }
}

#endif // NETGEN_CORE_SIMD_WASM_HPP

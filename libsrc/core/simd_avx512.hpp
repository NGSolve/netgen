#ifndef NETGEN_CORE_SIMD_AVX512_HPP
#define NETGEN_CORE_SIMD_AVX512_HPP

/**************************************************************************/
/* File:   simd_avx512.hpp                                                */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>

namespace ngcore
{

  template <>
  class SIMD<mask64,8>
  {
    __mmask8 mask;
  public:
    SIMD (size_t i)
      : mask(_mm512_cmpgt_epi64_mask(_mm512_set1_epi64(i),
                                     _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)))
    { ; }
    SIMD (int i)
      : mask(_mm512_cmpgt_epi64_mask(_mm512_set1_epi64(i),
                                     _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)))
    { ; }
    SIMD (int64_t i)
      : mask(_mm512_cmpgt_epi64_mask(_mm512_set1_epi64(i),
                                     _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)))
    { ; }
    SIMD (__mmask8 _mask) : mask(_mask) { ; }
    __mmask8 Data() const { return mask; }
    static constexpr int Size() { return 8; }
    static NETGEN_INLINE SIMD<mask64, 8> GetMaskFromBits (unsigned int i)
    {
      return SIMD<mask64, 8>(__mmask8(i));
    }
  };

  template<>
  class alignas(64) SIMD<int64_t,8>
  {
    __m512i data;

  public:
    static constexpr int Size() { return 8; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (int64_t val) { data = _mm512_set1_epi64(val); }
    SIMD (int64_t v0, int64_t v1, int64_t v2, int64_t v3, int64_t v4, int64_t v5, int64_t v6, int64_t v7) { data = _mm512_set_epi64(v7,v6,v5,v4,v3,v2,v1,v0); }
    SIMD (__m512i _data) { data = _data; }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<int64_t(int)>>::value, int>::type = 0>
      SIMD (const T & func)
    {
      data = _mm512_set_epi64(func(7), func(6), func(5), func(4), func(3), func(2), func(1), func(0));
    }


    NETGEN_INLINE auto operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    NETGEN_INLINE __m512i Data() const { return data; }
    NETGEN_INLINE __m512i & Data() { return data; }
    static SIMD FirstInt() { return { 0, 1, 2, 3, 4, 5, 6, 7 }; }
  };

  NETGEN_INLINE SIMD<int64_t,8> operator-(SIMD<int64_t,8> a) { return _mm512_sub_epi64(_mm512_setzero_si512(), a.Data()); }

  NETGEN_INLINE SIMD<int64_t,8> operator+ (SIMD<int64_t,8> a, SIMD<int64_t,8> b) { return _mm512_add_epi64(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<int64_t,8> operator- (SIMD<int64_t,8> a, SIMD<int64_t,8> b) { return _mm512_sub_epi64(a.Data(),b.Data()); }

   NETGEN_INLINE SIMD<int64_t,8> If (SIMD<mask64,8> a, SIMD<int64_t,8> b, SIMD<int64_t,8> c)
  { return _mm512_mask_blend_epi64(a.Data(), c.Data(), b.Data()); }


  template<>
  class alignas(64) SIMD<double,8>
  {
    __m512d data;
  public:
    static constexpr int Size() { return 8; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) { data = _mm512_set1_pd(val); }
    SIMD (int val)    { data = _mm512_set1_pd(val); }
    SIMD (size_t val) { data = _mm512_set1_pd(val); }
    SIMD (double const * p) { data = _mm512_loadu_pd(p); }
    SIMD (double const * p, SIMD<mask64,8> mask)
      { data = _mm512_mask_loadu_pd(_mm512_setzero_pd(), mask.Data(), p); }
    SIMD (__m512d _data) { data = _data; }
    SIMD (SIMD<double,4> v0, SIMD<double,4> v1)
        : data(_mm512_set_pd(v1[3], v1[2], v1[1], v1[0], v0[3], v0[2], v0[1], v0[0]))
    {}
    SIMD (SIMD<double,6> v0, SIMD<double,2> v1)
        : data(_mm512_set_pd(v1[1], v1[0], v0[5], v0[4], v0[3], v0[2], v0[1], v0[0]))
    {}

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = _mm512_set_pd(func(7), func(6), func(5), func(4), func(3), func(2), func(1), func(0));
    }

    void Store (double * p) { _mm512_storeu_pd(p, data); }
    void Store (double * p, SIMD<mask64,8> mask) { _mm512_mask_storeu_pd(p, mask.Data(), data); }

    template <typename Function>
    void SIMD_function (const Function & func, std::true_type)
    {
      data = (__m512d){ func(7), func(6), func(5), func(4),
                       func(3), func(2), func(1), func(0) };
    }

    // not a function
    void SIMD_function (double const * p, std::false_type)
    {
      data = _mm512_loadu_pd(p);
    }

    void SIMD_function (double val, std::false_type)
    {
      data = _mm512_set1_pd(val);
    }

    void SIMD_function (__m512d _data, std::false_type)
    {
      data = _data;
    }

    NETGEN_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    NETGEN_INLINE __m512d Data() const { return data; }
    NETGEN_INLINE __m512d & Data() { return data; }

    SIMD<double,4> Lo() const { return _mm512_extractf64x4_pd(data, 0); }
    SIMD<double,4> Hi() const { return _mm512_extractf64x4_pd(data, 1); }

    template <int I>
    double Get() const
    {
      static_assert(I>=0 && I<8, "Index out of range");
      return (*this)[I];
    }
  };

  NETGEN_INLINE SIMD<double,8> operator- (SIMD<double,8> a) { return _mm512_xor_pd(a.Data(), _mm512_set1_pd(-0.0)); } //{ return -a.Data(); }
  NETGEN_INLINE SIMD<double,8> operator+ (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_add_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,8> operator- (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_sub_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,8> operator* (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_mul_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,8> operator/ (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_div_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,8> operator* (double a, SIMD<double,8> b) { return _mm512_set1_pd(a)*b.Data(); }
  NETGEN_INLINE SIMD<double,8> operator* (SIMD<double,8> b, double a) { return _mm512_set1_pd(a)*b.Data(); }

  NETGEN_INLINE SIMD<double,8> sqrt (SIMD<double,8> a) { return _mm512_sqrt_pd(a.Data()); }
  NETGEN_INLINE SIMD<double,8> floor (SIMD<double,8> a) { return _mm512_floor_pd(a.Data()); }
  NETGEN_INLINE SIMD<double,8> ceil (SIMD<double,8> a) { return _mm512_ceil_pd(a.Data()); }
  NETGEN_INLINE SIMD<double,8> fabs (SIMD<double,8> a) { return _mm512_max_pd(a.Data(), ( - a).Data()); }

  NETGEN_INLINE SIMD<mask64,8> operator<= (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_LE_OQ); }
  NETGEN_INLINE SIMD<mask64,8> operator< (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_LT_OQ); }
  NETGEN_INLINE SIMD<mask64,8> operator>= (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_GE_OQ); }
  NETGEN_INLINE SIMD<mask64,8> operator> (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_GT_OQ); }
  NETGEN_INLINE SIMD<mask64,8> operator== (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_EQ_OQ); }
  NETGEN_INLINE SIMD<mask64,8> operator!= (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_NEQ_OQ); }

  NETGEN_INLINE SIMD<mask64,8> operator<= (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_LE); }
  NETGEN_INLINE SIMD<mask64,8> operator< (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_LT); }
  NETGEN_INLINE SIMD<mask64,8> operator>= (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(),  _MM_CMPINT_NLT); }
  NETGEN_INLINE SIMD<mask64,8> operator> (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_NLE); }
  NETGEN_INLINE SIMD<mask64,8> operator== (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_EQ); }
  NETGEN_INLINE SIMD<mask64,8> operator!= (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_NE); }

  NETGEN_INLINE SIMD<mask64,8> operator&& (SIMD<mask64,8> a, SIMD<mask64,8> b)
  { return (__mmask8)(a.Data() & b.Data()); }
  NETGEN_INLINE SIMD<mask64,8> operator|| (SIMD<mask64,8> a, SIMD<mask64,8> b)
  { return (__mmask8)(a.Data() | b.Data()); }
  NETGEN_INLINE SIMD<mask64,8> operator! (SIMD<mask64,8> a)
  { return (__mmask8)(~a.Data()); }

  NETGEN_INLINE SIMD<double,8> If (SIMD<mask64,8> a, SIMD<double,8> b, SIMD<double,8> c)
  { return _mm512_mask_blend_pd(a.Data(), c.Data(), b.Data()); }

  NETGEN_INLINE SIMD<double,8> IfPos (SIMD<double,8> a, SIMD<double> b, SIMD<double> c)
  {
    auto k = _mm512_cmp_pd_mask(a.Data(),_mm512_setzero_pd(), _CMP_GT_OS);
    return _mm512_mask_blend_pd(k,c.Data(),b.Data());
  }
  NETGEN_INLINE SIMD<double,8> IfZero (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    auto k = _mm512_cmp_pd_mask(a.Data(),_mm512_setzero_pd(), _CMP_EQ_OS);
    return _mm512_mask_blend_pd(k,c.Data(),b.Data());
  }


  NETGEN_INLINE auto Unpack (SIMD<double,8> a, SIMD<double,8> b)
  {
    return std::make_tuple(SIMD<double,8>(_mm512_unpacklo_pd(a.Data(),b.Data())),
                      SIMD<double,8>(_mm512_unpackhi_pd(a.Data(),b.Data())));
  }


  NETGEN_INLINE double HSum (SIMD<double,8> sd)
  {
    SIMD<double,4> low = _mm512_extractf64x4_pd(sd.Data(),0);
    SIMD<double,4> high = _mm512_extractf64x4_pd(sd.Data(),1);
    return HSum(low)+HSum(high);
  }

  NETGEN_INLINE auto HSum (SIMD<double,8> sd1, SIMD<double,8> sd2)
  {
    return SIMD<double,2>(HSum(sd1), HSum(sd2));
  }

  NETGEN_INLINE SIMD<double,4> HSum (SIMD<double,8> v1, SIMD<double,8> v2, SIMD<double,8> v3, SIMD<double,8> v4)
  {
    SIMD<double> lo,hi;
    std::tie(lo,hi) = Unpack(v1, v2);
    SIMD<double> sum01 = lo+hi;
    std::tie(lo,hi) = Unpack(v3, v4);
    SIMD<double> sum23 = lo+hi;
    // sum01  b a b a b a b a
    // sum23  d c d c d c d c
    // __m512 perm = _mm512_permutex2var_pd (sum01.Data(), _mm512_set_epi64(1,2,3,4,5,6,7,8), sum23.Data());
    SIMD<double,4> ab =  _mm512_extractf64x4_pd(sum01.Data(),0) + _mm512_extractf64x4_pd(sum01.Data(),1);
    SIMD<double,4> cd =  _mm512_extractf64x4_pd(sum23.Data(),0) + _mm512_extractf64x4_pd(sum23.Data(),1);
    return _mm256_add_pd (_mm256_permute2f128_pd (ab.Data(), cd.Data(), 1 + 2 * 16), _mm256_blend_pd(ab.Data(), cd.Data(), 12));
  }

  NETGEN_INLINE SIMD<double,8> FMA (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fmadd_pd (a.Data(), b.Data(), c.Data());
  }
  NETGEN_INLINE SIMD<double,8> FMA (const double & a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fmadd_pd (_mm512_set1_pd(a), b.Data(), c.Data());
  }

  NETGEN_INLINE SIMD<double,8> FNMA (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fnmadd_pd (a.Data(), b.Data(), c.Data());
  }
  NETGEN_INLINE SIMD<double,8> FNMA (const double & a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fnmadd_pd (_mm512_set1_pd(a), b.Data(), c.Data());
  }

  NETGEN_INLINE SIMD<double,8> FMAddSub (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fmaddsub_pd(a.Data(), b.Data(), c.Data());
  }

  NETGEN_INLINE SIMD<double,8> SwapPairs (SIMD<double,8> a)
  {
    return _mm512_shuffle_pd (a.Data(), a.Data(), 0b01010101);
  }
  
}

#endif // NETGEN_CORE_SIMD_AVX512_HPP

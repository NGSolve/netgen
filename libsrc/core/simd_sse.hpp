#ifndef NETGEN_CORE_SIMD_SSE_HPP
#define NETGEN_CORE_SIMD_SSE_HPP

/**************************************************************************/
/* File:   simd_sse.hpp                                                   */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>

namespace ngcore
{

  template <>
  class SIMD<mask64,2>
  {
    __m128i mask;
  public:
    SIMD (int i)
      : mask(_mm_cmpgt_epi32(_mm_set1_epi32(i),
                             _mm_set_epi32(1, 1, 0, 0)))
    { ; }
    SIMD (__m128i _mask) : mask(_mask) { ; }
    __m128i Data() const { return mask; }
    static constexpr int Size() { return 2; }
    static NETGEN_INLINE SIMD<mask64, 2> GetMaskFromBits (unsigned int i);
  };

  static SIMD<mask64, 2> masks_from_2bits[4] = {
    _mm_set_epi32 (0,0,0,0), _mm_set_epi32 (0,0,-1,0),
    _mm_set_epi32 (-1,0,0,0), _mm_set_epi32 (-1,0,-1,0),
  };

  NETGEN_INLINE SIMD<mask64, 2> SIMD<mask64, 2> :: GetMaskFromBits (unsigned int i)
  {
    return masks_from_2bits[i & 3];
  }


  template<>
  class alignas(16) SIMD<int64_t,2>
  {
    __m128i data;

  public:
    static constexpr int Size() { return 2; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD (int64_t v0, int64_t v1) { data = _mm_set_epi64x(v1,v0); }
    SIMD (std::array<int64_t, 2> arr)
        : data{_mm_set_epi64x(arr[1],arr[0])}
    {}

    SIMD & operator= (const SIMD &) = default;

    SIMD (int64_t val) { data = _mm_set1_epi64x(val); }
    SIMD (__m128i _data) { data = _data; }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<int64_t(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = _mm_set_epi64(func(1), func(0));
    }

    NETGEN_INLINE auto operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    NETGEN_INLINE __m128i Data() const { return data; }
    NETGEN_INLINE __m128i & Data() { return data; }
    static SIMD FirstInt(int n0=0) { return { n0, n0+1 }; }
  };



NETGEN_INLINE SIMD<int64_t,2> operator-(SIMD<int64_t,2> a) { return _mm_sub_epi64(_mm_setzero_si128(), a.Data()); }
NETGEN_INLINE SIMD<int64_t,2> operator+ (SIMD<int64_t,2> a, SIMD<int64_t,2> b) { return _mm_add_epi64(a.Data(),b.Data()); }
NETGEN_INLINE SIMD<int64_t,2> operator- (SIMD<int64_t,2> a, SIMD<int64_t,2> b) { return _mm_sub_epi64(a.Data(),b.Data()); }


  template<>
  class alignas(16) SIMD<double,2>
  {
    __m128d data;

  public:
    static constexpr int Size() { return 2; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD (double v0, double v1) { data = _mm_set_pd(v1,v0); }
    SIMD (SIMD<double,1> v0, SIMD<double,1> v1)
        : data{_mm_set_pd(v0.Data(), v1.Data())}
    { }
    SIMD (std::array<double, 2> arr)
        : data{_mm_set_pd(arr[1], arr[0])}
    {}

    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) { data = _mm_set1_pd(val); }
    SIMD (int val)    { data = _mm_set1_pd(val); }
    SIMD (size_t val) { data = _mm_set1_pd(val); }

    SIMD (double const * p) { data = _mm_loadu_pd(p); }
    SIMD (double const * p, SIMD<mask64,2> mask)
      {
#ifdef __AVX__
        data = _mm_maskload_pd(p, mask.Data());
#else
        // this versions segfaults if p points to the last allowed element
        // happened on Mac with the new SparseCholesky-factorization
        // data = _mm_and_pd(_mm_castsi128_pd(mask.Data()), _mm_loadu_pd(p));
        auto pmask = (int64_t*)&mask;
        data = _mm_set_pd (pmask[1] ? p[1] : 0.0, pmask[0] ? p[0] : 0.0);
#endif
      }
    SIMD (__m128d _data) { data = _data; }

    void Store (double * p) { _mm_storeu_pd(p, data); }
    void Store (double * p, SIMD<mask64,2> mask)
    {
#ifdef __AVX__
      _mm_maskstore_pd(p, mask.Data(), data);
#else
      /*
      _mm_storeu_pd (p, _mm_or_pd (_mm_and_pd(_mm_castsi128_pd(mask.Data()), data),
                                   _mm_andnot_pd(_mm_castsi128_pd(mask.Data()), _mm_loadu_pd(p))));
      */
      auto pmask = (int64_t*)&mask;
      if (pmask[0]) p[0] = (*this)[0];
      if (pmask[1]) p[1] = (*this)[1];
#endif
    }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = _mm_set_pd(func(1), func(0));
    }

    NETGEN_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    NETGEN_INLINE __m128d Data() const { return data; }
    NETGEN_INLINE __m128d & Data() { return data; }

    template <int I>
    double Get() const
    {
      static_assert(I>=0 && I<2, "Index out of range");
      return (*this)[I];
    }

    double Lo() const { return Get<0>(); } 
    double Hi() const { return Get<1>(); } 
    
    operator std::tuple<double&,double&> ()
    {
      auto pdata = (double*)&data;
      return std::tuple<double&,double&>(pdata[0], pdata[1]);
    }
  };

  NETGEN_INLINE SIMD<double,2> operator- (SIMD<double,2> a) { return _mm_xor_pd(a.Data(), _mm_set1_pd(-0.0)); }
  NETGEN_INLINE SIMD<double,2> operator+ (SIMD<double,2> a, SIMD<double,2> b) { return _mm_add_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,2> operator- (SIMD<double,2> a, SIMD<double,2> b) { return _mm_sub_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,2> operator* (SIMD<double,2> a, SIMD<double,2> b) { return _mm_mul_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,2> operator/ (SIMD<double,2> a, SIMD<double,2> b) { return _mm_div_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,2> operator* (double a, SIMD<double,2> b) { return _mm_set1_pd(a)*b; }
  NETGEN_INLINE SIMD<double,2> operator* (SIMD<double,2> b, double a) { return _mm_set1_pd(a)*b; }

  template<>
  NETGEN_INLINE auto Unpack (SIMD<double,2> a, SIMD<double,2> b)
  {
    return std::make_tuple(SIMD<double,2>(_mm_unpacklo_pd(a.Data(),b.Data())),
                      SIMD<double,2>(_mm_unpackhi_pd(a.Data(),b.Data())));
  }

  NETGEN_INLINE __m128d my_mm_hadd_pd(__m128d a, __m128d b) {
#if defined(__SSE3__) || defined(__AVX__)
    return _mm_hadd_pd(a,b);
#else
    return _mm_add_pd( _mm_unpacklo_pd(a,b), _mm_unpackhi_pd(a,b) );
#endif
  }

#ifndef __AVX__
  NETGEN_INLINE __m128i my_mm_cmpgt_epi64(__m128i a, __m128i b) {
    auto  res_lo = _mm_cvtsi128_si64(a)  > _mm_cvtsi128_si64(b) ? -1:0;
    auto  res_hi = _mm_cvtsi128_si64(_mm_srli_si128(a,8)) > _mm_cvtsi128_si64(_mm_srli_si128(b,8)) ? -1 : 0;
    return _mm_set_epi64x(res_hi,res_lo);
  }
#else
  NETGEN_INLINE __m128i my_mm_cmpgt_epi64(__m128i a, __m128i b) {
    return _mm_cmpgt_epi64(a,b);
  }
#endif


  NETGEN_INLINE SIMD<double,2> sqrt (SIMD<double,2> a) { return _mm_sqrt_pd(a.Data()); }
  NETGEN_INLINE SIMD<double,2> fabs (SIMD<double,2> a) { return _mm_max_pd(a.Data(), (-a).Data()); }
  using std::floor;
  NETGEN_INLINE SIMD<double,2> floor (SIMD<double,2> a)
  { return ngcore::SIMD<double,2>([&](int i)->double { return floor(a[i]); } ); }
  using std::ceil;
  NETGEN_INLINE SIMD<double,2> ceil (SIMD<double,2> a)
  { return ngcore::SIMD<double,2>([&](int i)->double { return ceil(a[i]); } ); }

  NETGEN_INLINE SIMD<mask64,2> operator<= (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmple_pd(a.Data(),b.Data())); }
  NETGEN_INLINE SIMD<mask64,2> operator< (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmplt_pd(a.Data(),b.Data())); }
  NETGEN_INLINE SIMD<mask64,2> operator>= (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpge_pd(a.Data(),b.Data())); }
  NETGEN_INLINE SIMD<mask64,2> operator> (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpgt_pd(a.Data(),b.Data())); }
  NETGEN_INLINE SIMD<mask64,2> operator== (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpeq_pd(a.Data(),b.Data())); }
  NETGEN_INLINE SIMD<mask64,2> operator!= (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpneq_pd(a.Data(),b.Data())); }

  NETGEN_INLINE SIMD<mask64,2> operator<= (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_xor_si128(_mm_cmpgt_epi64(a.Data(),b.Data()),_mm_set1_epi32(-1)); }
  NETGEN_INLINE SIMD<mask64,2> operator< (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  my_mm_cmpgt_epi64(b.Data(),a.Data()); }
  NETGEN_INLINE SIMD<mask64,2> operator>= (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_xor_si128(_mm_cmpgt_epi64(b.Data(),a.Data()),_mm_set1_epi32(-1)); }
  NETGEN_INLINE SIMD<mask64,2> operator> (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  my_mm_cmpgt_epi64(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<mask64,2> operator== (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_cmpeq_epi64(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<mask64,2> operator!= (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_xor_si128(_mm_cmpeq_epi64(a.Data(),b.Data()),_mm_set1_epi32(-1)); }



 NETGEN_INLINE SIMD<mask64,2> operator&& (SIMD<mask64,2> a, SIMD<mask64,2> b)
  { return _mm_castpd_si128(_mm_and_pd (_mm_castsi128_pd(a.Data()),_mm_castsi128_pd( b.Data()))); }
  NETGEN_INLINE SIMD<mask64,2> operator|| (SIMD<mask64,2> a, SIMD<mask64,2> b)
  { return _mm_castpd_si128(_mm_or_pd (_mm_castsi128_pd(a.Data()), _mm_castsi128_pd(b.Data()))); }
  NETGEN_INLINE SIMD<mask64,2> operator! (SIMD<mask64,2> a)
  { return _mm_castpd_si128(_mm_xor_pd (_mm_castsi128_pd(a.Data()),_mm_castsi128_pd( _mm_cmpeq_epi64(a.Data(),a.Data())))); }
#ifdef __SSE4_1__
  NETGEN_INLINE SIMD<double,2> If (SIMD<mask64,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return _mm_blendv_pd(c.Data(), b.Data(), _mm_castsi128_pd(a.Data())); }
#else
  NETGEN_INLINE SIMD<double,2> If (SIMD<mask64,2> a, SIMD<double,2> b, SIMD<double,2> c)
  {
    return _mm_or_pd(
                      _mm_andnot_pd(_mm_castsi128_pd(a.Data()),c.Data()),
                      _mm_and_pd(b.Data(),_mm_castsi128_pd(a.Data()))
                      );}
#endif // __SSE4_1__

  NETGEN_INLINE SIMD<double,2> IfPos (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return ngcore::SIMD<double,2>([&](int i)->double { return a[i]>0 ? b[i] : c[i]; }); }
  NETGEN_INLINE SIMD<double,2> IfZero (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return ngcore::SIMD<double,2>([&](int i)->double { return a[i]==0. ? b[i] : c[i]; }); }


  NETGEN_INLINE double HSum (SIMD<double,2> sd)
  {
    return _mm_cvtsd_f64 (my_mm_hadd_pd (sd.Data(), sd.Data()));
  }

  NETGEN_INLINE auto HSum (SIMD<double,2> sd1, SIMD<double,2> sd2)
  {
    __m128d hv2 = my_mm_hadd_pd(sd1.Data(), sd2.Data());
    return SIMD<double,2> (hv2);
    // return SIMD<double,2>(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  NETGEN_INLINE SIMD<int64_t, 2> If(SIMD<mask64, 2> a, SIMD<int64_t, 2> b,
                             SIMD<int64_t, 2> c) {
    return _mm_or_si128(
                        _mm_andnot_si128(a.Data(),c.Data()),
                        _mm_and_si128(b.Data(),a.Data())
                        );
  }

}

#endif // NETGEN_CORE_SIMD_SSE_HPP

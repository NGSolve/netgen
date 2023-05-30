#ifndef NETGEN_CORE_SIMD_AVX_HPP
#define NETGEN_CORE_SIMD_AVX_HPP

/**************************************************************************/
/* File:   simd_avx.hpp                                                   */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>

namespace ngcore
{

#if defined(__GNUC__) && (__GNUC__ == 7)
  // GCC7 does not have intrinsic _mm256_set_m128i, see
  // https://stackoverflow.com/questions/32630458/setting-m256i-to-the-value-of-two-m128i-values
  NETGEN_INLINE auto _mm256_set_m128i(__m128i v0, __m128i v1) {
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1);
  }
#endif // defined(__GNUC__) && (__GNUC__ == 7)

#if defined(__AVX2__)
  NETGEN_INLINE __m256i my_mm256_cmpgt_epi64 (__m256i a, __m256i b)
  {
    return _mm256_cmpgt_epi64 (a,b);
  }
#else
  NETGEN_INLINE __m256i my_mm256_cmpgt_epi64 (__m256i a, __m256i b)
  {
    __m128i rlo = _mm_cmpgt_epi64(_mm256_extractf128_si256(a, 0),
                                  _mm256_extractf128_si256(b, 0));
    __m128i rhi = _mm_cmpgt_epi64(_mm256_extractf128_si256(a, 1),
                                  _mm256_extractf128_si256(b, 1));
    return _mm256_insertf128_si256 (_mm256_castsi128_si256(rlo), rhi, 1);
  }
#endif


  template <>
  class SIMD<mask64,4>
  {
    __m256i mask;
  public:
    SIMD (int64_t i)
      : mask(my_mm256_cmpgt_epi64(_mm256_set1_epi64x(i),
                                  _mm256_set_epi64x(3, 2, 1, 0)))
    { ; }
    SIMD (__m256i _mask) : mask(_mask) { ; }
    SIMD (__m256d _mask) : mask(_mm256_castpd_si256(_mask)) { ; }
    __m256i Data() const { return mask; }
    static constexpr int Size() { return 4; }
    static SIMD<mask64, 4> GetMaskFromBits (unsigned int i);
  };

  static SIMD<mask64, 4> masks_from_4bits[16] = {
    _mm256_set_epi64x (0,0,0,0), _mm256_set_epi64x (0,0,0,-1),
    _mm256_set_epi64x (0,0,-1,0), _mm256_set_epi64x (0,0,-1,-1),
    _mm256_set_epi64x (0,-1,0,0), _mm256_set_epi64x (0,-1,0,-1),
    _mm256_set_epi64x (0,-1,-1,0), _mm256_set_epi64x (0,-1,-1,-1),
    _mm256_set_epi64x (-1,0,0,0), _mm256_set_epi64x (-1,0,0,-1),
    _mm256_set_epi64x (-1,0,-1,0), _mm256_set_epi64x (-1,0,-1,-1),
    _mm256_set_epi64x (-1,-1,0,0), _mm256_set_epi64x (-1,-1,0,-1),
    _mm256_set_epi64x (-1,-1,-1,0), _mm256_set_epi64x (-1,-1,-1,-1)
  };

  NETGEN_INLINE SIMD<mask64, 4> SIMD<mask64, 4> :: GetMaskFromBits (unsigned int i)
  {
    return masks_from_4bits[i & 15];
  }

  template<>
  class alignas(32) SIMD<int64_t,4>
  {
    __m256i data;

  public:
    static constexpr int Size() { return 4; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (int64_t val) { data = _mm256_set1_epi64x(val); }
    SIMD (int64_t v0, int64_t v1, int64_t v2, int64_t v3) { data = _mm256_set_epi64x(v3,v2,v1,v0); }
    SIMD (std::array<int64_t,4> a)
      : data{_mm256_set_epi64x(a[3],a[2],a[1],a[0])}
    {}
    SIMD (SIMD<int64_t,2> v0, SIMD<int64_t,2> v1)
        : data(_mm256_set_m128i(v0.Data(),v1.Data()))
      {}
    SIMD (__m256i _data) { data = _data; }

    NETGEN_INLINE auto operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    NETGEN_INLINE __m256i Data() const { return data; }
    NETGEN_INLINE __m256i & Data() { return data; }

    SIMD<int64_t,2> Lo() const { return _mm256_extractf128_si256(data, 0); }
    SIMD<int64_t,2> Hi() const { return _mm256_extractf128_si256(data, 1); }
    static SIMD FirstInt(int n0=0) { return { n0+0, n0+1, n0+2, n0+3 }; }
  };


  NETGEN_INLINE SIMD<int64_t,4> operator-(SIMD<int64_t,4> a) { return _mm256_sub_epi64(_mm256_setzero_si256(), a.Data()); }

#ifdef __AVX2__
  NETGEN_INLINE SIMD<int64_t,4> operator+ (SIMD<int64_t,4> a, SIMD<int64_t,4> b) { return _mm256_add_epi64(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<int64_t,4> operator- (SIMD<int64_t,4> a, SIMD<int64_t,4> b) { return _mm256_sub_epi64(a.Data(),b.Data()); }
#endif // __AVX2__

  template<>
  class alignas(32) SIMD<double,4>
  {
    __m256d data;

  public:
    static constexpr int Size() { return 4; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) { data = _mm256_set1_pd(val); }
    SIMD (int val)    { data = _mm256_set1_pd(val); }
    SIMD (size_t val) { data = _mm256_set1_pd(val); }
    SIMD (double v0, double v1, double v2, double v3) { data = _mm256_set_pd(v3,v2,v1,v0); }
    SIMD (SIMD<double,2> v0, SIMD<double,2> v1) : SIMD(v0[0], v0[1], v1[0], v1[1]) { ; }
    SIMD (double const * p) { data = _mm256_loadu_pd(p); }
    SIMD (double const * p, SIMD<mask64,4> mask) { data = _mm256_maskload_pd(p, mask.Data()); }
    SIMD (__m256d _data) { data = _data; }
    SIMD (std::array<double,4> a)
      : data{_mm256_set_pd(a[3],a[2],a[1],a[0])}
    {}

    void Store (double * p) { _mm256_storeu_pd(p, data); }
    void Store (double * p, SIMD<mask64,4> mask) { _mm256_maskstore_pd(p, mask.Data(), data); }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = _mm256_set_pd(func(3), func(2), func(1), func(0));
    }

    NETGEN_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    NETGEN_INLINE double & operator[] (int i) { return ((double*)(&data))[i]; }
    // [[deprecated("don't write to individual elements of SIMD")]]
    // NETGEN_INLINE double & operator[] (int i) { return ((double*)(&data))[i]; }
    NETGEN_INLINE __m256d Data() const { return data; }
    NETGEN_INLINE __m256d & Data() { return data; }

    SIMD<double,2> Lo() const { return _mm256_extractf128_pd(data, 0); }
    SIMD<double,2> Hi() const { return _mm256_extractf128_pd(data, 1); }

    operator std::tuple<double&,double&,double&,double&> ()
    { return std::tuple<double&,double&,double&,double&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }

    template <int I>
    double Get() const
    {
      static_assert(I>=0 && I<4, "Index out of range");
      return (*this)[I];
    }
  };

  NETGEN_INLINE auto Unpack (SIMD<double,4> a, SIMD<double,4> b)
  {
    return std::make_tuple(SIMD<double,4>(_mm256_unpacklo_pd(a.Data(),b.Data())),
                      SIMD<double,4>(_mm256_unpackhi_pd(a.Data(),b.Data())));
  }

  NETGEN_INLINE SIMD<double,4> operator- (SIMD<double,4> a) { return _mm256_xor_pd(a.Data(), _mm256_set1_pd(-0.0)); }
  NETGEN_INLINE SIMD<double,4> operator+ (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_add_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,4> operator- (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_sub_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,4> operator* (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_mul_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,4> operator/ (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_div_pd(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<double,4> operator* (double a, SIMD<double,4> b) { return _mm256_set1_pd(a)*b.Data(); }
  NETGEN_INLINE SIMD<double,4> operator* (SIMD<double,4> b, double a) { return _mm256_set1_pd(a)*b.Data(); }

  NETGEN_INLINE SIMD<double,4> sqrt (SIMD<double,4> a) { return _mm256_sqrt_pd(a.Data()); }
  NETGEN_INLINE SIMD<double,4> floor (SIMD<double,4> a) { return _mm256_floor_pd(a.Data()); }
  NETGEN_INLINE SIMD<double,4> ceil (SIMD<double,4> a) { return _mm256_ceil_pd(a.Data()); }
  NETGEN_INLINE SIMD<double,4> fabs (SIMD<double,4> a) { return _mm256_max_pd(a.Data(), (-a).Data()); }


#ifdef __FMA__
  NETGEN_INLINE SIMD<double,4> FMA (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return _mm256_fmadd_pd (a.Data(), b.Data(), c.Data());
  }
  NETGEN_INLINE SIMD<double,4> FMA (const double & a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return _mm256_fmadd_pd (_mm256_set1_pd(a), b.Data(), c.Data());
  }
  NETGEN_INLINE SIMD<double,4> FNMA (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return _mm256_fnmadd_pd (a.Data(), b.Data(), c.Data());
  }
  NETGEN_INLINE SIMD<double,4> FNMA (const double & a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return _mm256_fnmadd_pd (_mm256_set1_pd(a), b.Data(), c.Data());
  }
#endif

#if defined(__FMA__) && !defined(__AVX512F__)
  // make sure to use the update-version of fma
  // important in matrix kernels using 12 sum-registers, 3 a-values and updated b-value
  // avx512 has enough registers, and gcc seems to use only the first 16 z-regs
  NETGEN_INLINE void FMAasm (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> & sum)
  {
    asm ("vfmadd231pd %[a], %[b], %[sum]"
         : [sum] "+x" (sum.Data())
         : [a] "x" (a.Data()), [b] "x" (b.Data())
         );
  }

  NETGEN_INLINE void FNMAasm (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> & sum)
  {
    asm ("vfnmadd231pd %[a], %[b], %[sum]"
         : [sum] "+x" (sum.Data())
         : [a] "x" (a.Data()), [b] "x" (b.Data())
         );
  }
#endif

#if defined(__FMA__) 
  NETGEN_INLINE SIMD<double,4> FMAddSub (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return _mm256_fmaddsub_pd(a.Data(), b.Data(), c.Data());
  }
#endif  
  
  NETGEN_INLINE SIMD<double,4> SwapPairs (SIMD<double,4> a)
  {
    return _mm256_shuffle_pd (a.Data(), a.Data(), 0b0101);
  }

  
  NETGEN_INLINE SIMD<mask64,4> operator<= (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_LE_OQ); }
  NETGEN_INLINE SIMD<mask64,4> operator< (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_LT_OQ); }
  NETGEN_INLINE SIMD<mask64,4> operator>= (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_GE_OQ); }
  NETGEN_INLINE SIMD<mask64,4> operator> (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_GT_OQ); }
  NETGEN_INLINE SIMD<mask64,4> operator== (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_EQ_OQ); }
  NETGEN_INLINE SIMD<mask64,4> operator!= (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_NEQ_OQ); }

  NETGEN_INLINE SIMD<mask64,4> operator<= (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_xor_si256(_mm256_cmpgt_epi64(a.Data(),b.Data()),_mm256_set1_epi32(-1)); }
  NETGEN_INLINE SIMD<mask64,4> operator< (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  my_mm256_cmpgt_epi64(b.Data(),a.Data()); }
  NETGEN_INLINE SIMD<mask64,4> operator>= (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_xor_si256(_mm256_cmpgt_epi64(b.Data(),a.Data()),_mm256_set1_epi32(-1)); }
  NETGEN_INLINE SIMD<mask64,4> operator> (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  my_mm256_cmpgt_epi64(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<mask64,4> operator== (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_cmpeq_epi64(a.Data(),b.Data()); }
  NETGEN_INLINE SIMD<mask64,4> operator!= (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_xor_si256(_mm256_cmpeq_epi64(a.Data(),b.Data()),_mm256_set1_epi32(-1)); }

#ifdef __AVX2__
  NETGEN_INLINE SIMD<mask64,4> operator&& (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_and_si256 (a.Data(), b.Data()); }
  NETGEN_INLINE SIMD<mask64,4> operator|| (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_or_si256 (a.Data(), b.Data()); }
  NETGEN_INLINE SIMD<mask64,4> operator! (SIMD<mask64,4> a)
  { return _mm256_xor_si256 (a.Data(), _mm256_cmpeq_epi64(a.Data(),a.Data())); }
#else //AVX2 is a superset of AVX. Without it, it is necessary to reinterpret the types
  NETGEN_INLINE SIMD<mask64,4> operator&& (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_castpd_si256(_mm256_and_pd (_mm256_castsi256_pd(a.Data()),_mm256_castsi256_pd( b.Data()))); }
  NETGEN_INLINE SIMD<mask64,4> operator|| (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_castpd_si256(_mm256_or_pd (_mm256_castsi256_pd(a.Data()), _mm256_castsi256_pd(b.Data()))); }
  NETGEN_INLINE SIMD<mask64,4> operator! (SIMD<mask64,4> a)
  { return _mm256_castpd_si256(_mm256_xor_pd (_mm256_castsi256_pd(a.Data()),_mm256_castsi256_pd( _mm256_cmpeq_epi64(a.Data(),a.Data())))); }
#endif
  NETGEN_INLINE SIMD<double,4> If (SIMD<mask64,4> a, SIMD<double,4> b, SIMD<double,4> c)
  { return _mm256_blendv_pd(c.Data(), b.Data(), _mm256_castsi256_pd(a.Data())); }

  NETGEN_INLINE SIMD<double,4> IfPos (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_GT_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  NETGEN_INLINE SIMD<double,4> IfZero (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_EQ_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  NETGEN_INLINE double HSum (SIMD<double,4> sd)
  {
    // __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    __m128d hv = (sd.Lo()+sd.Hi()).Data();
    return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }

  NETGEN_INLINE auto HSum (SIMD<double,4> sd1, SIMD<double,4> sd2)
  {
    __m256d hv = _mm256_hadd_pd(sd1.Data(), sd2.Data());
    __m128d hv2 = _mm_add_pd (_mm256_extractf128_pd(hv,0), _mm256_extractf128_pd(hv,1));
    return SIMD<double,2>(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  NETGEN_INLINE auto HSum (SIMD<double,4> v1, SIMD<double,4> v2, SIMD<double,4> v3, SIMD<double,4> v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1.Data(), v2.Data());
    __m256d hsum2 = _mm256_hadd_pd (v3.Data(), v4.Data());
    SIMD<double,4> hsum = _mm256_add_pd (_mm256_permute2f128_pd (hsum1, hsum2, 1+2*16),
                                         _mm256_blend_pd (hsum1, hsum2, 12));
    return hsum;
    // return make_tuple(hsum[0], hsum[1], hsum[2], hsum[3]);
  }


  NETGEN_INLINE SIMD<int64_t,4> If (SIMD<mask64,4> a, SIMD<int64_t,4> b, SIMD<int64_t,4> c)
  { return _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(c.Data()), _mm256_castsi256_pd(b.Data()),
                                                _mm256_castsi256_pd(a.Data()))); }



}

#endif // NETGEN_CORE_SIMD_AVX_HPP


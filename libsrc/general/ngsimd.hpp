#ifndef FILE_NGSIMD
#define FILE_NGSIMD
/**************************************************************************/
/* File:   ngsimd.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>
#include <tuple>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>

#ifdef WIN32
#ifndef AVX_OPERATORS_DEFINED
#define AVX_OPERATORS_DEFINED
NG_INLINE __m128d operator- (__m128d a) { return _mm_xor_pd(a, _mm_set1_pd(-0.0)); }
NG_INLINE __m128d operator+ (__m128d a, __m128d b) { return _mm_add_pd(a,b); }
NG_INLINE __m128d operator- (__m128d a, __m128d b) { return _mm_sub_pd(a,b); }
NG_INLINE __m128d operator* (__m128d a, __m128d b) { return _mm_mul_pd(a,b); }
NG_INLINE __m128d operator/ (__m128d a, __m128d b) { return _mm_div_pd(a,b); }
NG_INLINE __m128d operator* (double a, __m128d b) { return _mm_set1_pd(a)*b; }
NG_INLINE __m128d operator* (__m128d b, double a) { return _mm_set1_pd(a)*b; }

NG_INLINE __m128d operator+= (__m128d &a, __m128d b) { return a = a+b; }
NG_INLINE __m128d operator-= (__m128d &a, __m128d b) { return a = a-b; }
NG_INLINE __m128d operator*= (__m128d &a, __m128d b) { return a = a*b; }
NG_INLINE __m128d operator/= (__m128d &a, __m128d b) { return a = a/b; }

NG_INLINE __m256d operator- (__m256d a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }
NG_INLINE __m256d operator+ (__m256d a, __m256d b) { return _mm256_add_pd(a,b); }
NG_INLINE __m256d operator- (__m256d a, __m256d b) { return _mm256_sub_pd(a,b); }
NG_INLINE __m256d operator* (__m256d a, __m256d b) { return _mm256_mul_pd(a,b); }
NG_INLINE __m256d operator/ (__m256d a, __m256d b) { return _mm256_div_pd(a,b); }
NG_INLINE __m256d operator* (double a, __m256d b) { return _mm256_set1_pd(a)*b; }
NG_INLINE __m256d operator* (__m256d b, double a) { return _mm256_set1_pd(a)*b; }

NG_INLINE __m256d operator+= (__m256d &a, __m256d b) { return a = a+b; }
NG_INLINE __m256d operator-= (__m256d &a, __m256d b) { return a = a-b; }
NG_INLINE __m256d operator*= (__m256d &a, __m256d b) { return a = a*b; }
NG_INLINE __m256d operator/= (__m256d &a, __m256d b) { return a = a/b; }
#endif
#endif



namespace ngsimd
{

  // MSVC does not define SSE. It's always present on 64bit cpus
#if (defined(_M_AMD64) || defined(_M_X64) || defined(__AVX__))
#ifndef __SSE__
#define __SSE__
#endif
#ifndef __SSE2__
#define __SSE2__
#endif
#endif

  
  constexpr int GetDefaultSIMDSize() {
#if defined __AVX512F__
    return 8;
#elif defined __AVX__
    return 4;
#elif defined __SSE__
    return 2;
#else
    return 1;
#endif
  }

#if defined __AVX512F__
    typedef __m512 tAVX;
    typedef __m512d tAVXd;
#elif defined __AVX__
    typedef __m256 tAVX;
    typedef __m256d tAVXd;
#elif defined __SSE__
    typedef __m128 tAVX;
    typedef __m128d tAVXd;
#endif


  
  template <typename T, int N=GetDefaultSIMDSize()> class SIMD;

  template <typename T>
  struct has_call_operator
  {
      template <typename C> static std::true_type check( decltype( sizeof(&C::operator() )) ) { return std::true_type(); }
      template <typename> static std::false_type check(...) { return std::false_type(); }
      typedef decltype( check<T>(sizeof(char)) ) type;
      static constexpr type value = type();
  };

  template <typename T1, typename T2, typename T3>
  // a*b+c
  NG_INLINE auto FMA(T1 a, T2 b, T3 c)
  {
    return a*b+c;
  }

  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator+ (T a, SIMD<double,N> b) { return SIMD<double,N>(a) + b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator- (T a, SIMD<double,N> b) { return SIMD<double,N>(a) - b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator* (T a, SIMD<double,N> b) { return SIMD<double,N>(a) * b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator/ (T a, SIMD<double,N> b) { return SIMD<double,N>(a) / b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator+ (SIMD<double,N> a, T b) { return a + SIMD<double,N>(b); }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator- (SIMD<double,N> a, T b) { return a - SIMD<double,N>(b); }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator* (SIMD<double,N> a, T b) { return a * SIMD<double,N>(b); }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    NG_INLINE SIMD<double,N> operator/ (SIMD<double,N> a, T b) { return a / SIMD<double,N>(b); }


#ifdef __AVX__
  template <typename T>
  class AlignedAlloc
  {
    protected:
      static void * aligned_malloc(size_t s)
      {
        // Assume 16 byte alignment of standard library
        if(alignof(T)<=16)
            return malloc(s);
        else
            return  _mm_malloc(s, alignof(T));
      }

      static void aligned_free(void *p)
      {
        if(alignof(T)<=16)
            free(p);
        else
            _mm_free(p);
      }

  public:
    void * operator new (size_t s, void *p) { return p; }
    void * operator new (size_t s) { return aligned_malloc(s); }
    void * operator new[] (size_t s) { return aligned_malloc(s); }
    void operator delete (void * p) { aligned_free(p); }
    void operator delete[] (void * p) { aligned_free(p); }
  };
#else
  // it's only a dummy without AVX
  template <typename T>
  class AlignedAlloc { ; };

#endif

using std::sqrt;
using std::fabs;

  class ExceptionNOSIMD : public std::runtime_error
  {
  public:
    using std::runtime_error::runtime_error;
    std::string What() { return what(); }
  };

  using std::exp;
  template<int N> NG_INLINE SIMD<double,N> exp (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return exp(a[i]); } );
  }

  using std::log;
  template<int N> NG_INLINE SIMD<double,N> log (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return log(a[i]); } );
  }

  using std::pow;
  template<int N> NG_INLINE SIMD<double,N> pow (SIMD<double,N> a, double x)
  {
    return SIMD<double,N>([&](int i)->double { return pow(a[i],x); } );
  }

  using std::sin;
  template<int N> NG_INLINE SIMD<double,N> sin (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return sin(a[i]); } );
  }

  using std::cos;
  template<int N> NG_INLINE SIMD<double,N> cos (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return cos(a[i]); } );
  }

  using std::tan;
  template<int N> NG_INLINE SIMD<double,N> tan (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return tan(a[i]); } );
  }

  using std::atan;
  template<int N> NG_INLINE SIMD<double,N> atan (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return atan(a[i]); } );
  }


/////////////////////////////////////////////////////////////////////////////
// SIMD width 1 (in case no AVX support is available)
/////////////////////////////////////////////////////////////////////////////
  template<>
  class SIMD<double,1>
  {
    double data;

  public:
    static constexpr int Size() { return 1; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    // only called if T has a call operator of appropriate type
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = func(0);
    }

    // only called if T is arithmetic (integral or floating point types)
    template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    SIMD (const T & val)
    {
      data = val;
    }

    SIMD (double const * p)
    {
      data = *p;
    }

    NG_INLINE operator double() const { return data; }
    NG_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    NG_INLINE double Data() const { return data; }
    NG_INLINE double & Data() { return data; }

    NG_INLINE SIMD<double,1> &operator+= (SIMD<double,1> b) { data+=b.Data(); return *this; }
    NG_INLINE SIMD<double,1> &operator-= (SIMD<double,1> b) { data-=b.Data(); return *this; }
    NG_INLINE SIMD<double,1> &operator*= (SIMD<double,1> b) { data*=b.Data(); return *this; }
    NG_INLINE SIMD<double,1> &operator/= (SIMD<double,1> b) { data/=b.Data(); return *this; }

  };

  NG_INLINE SIMD<double,1> operator+ (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()+b.Data(); }
  NG_INLINE SIMD<double,1> operator- (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()-b.Data(); }
  NG_INLINE SIMD<double,1> operator- (SIMD<double,1> a) { return -a.Data(); }
  NG_INLINE SIMD<double,1> operator* (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()*b.Data(); }
  NG_INLINE SIMD<double,1> operator/ (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()/b.Data(); }

  NG_INLINE SIMD<double,1> sqrt (SIMD<double,1> a) { return std::sqrt(a.Data()); }
  NG_INLINE SIMD<double,1> floor (SIMD<double,1> a) { return std::floor(a.Data()); }
  NG_INLINE SIMD<double,1> ceil (SIMD<double,1> a) { return std::ceil(a.Data()); }
  NG_INLINE SIMD<double,1> fabs (SIMD<double,1> a) { return std::fabs(a.Data()); }
  NG_INLINE SIMD<double,1> L2Norm2 (SIMD<double,1> a) { return a.Data()*a.Data(); }
  NG_INLINE SIMD<double,1> Trans (SIMD<double,1> a) { return a; }
  NG_INLINE SIMD<double,1> IfPos (SIMD<double,1> a, SIMD<double,1> b, SIMD<double,1> c)
  {
    return (a.Data() > 0) ? b : c;
  }

  NG_INLINE double HSum (SIMD<double,1> sd)
  {
    return sd.Data();
  }

  NG_INLINE auto HSum (SIMD<double,1> sd1, SIMD<double,1> sd2)
  {
    return std::make_tuple(sd1.Data(), sd2.Data());
  }

  NG_INLINE auto HSum (SIMD<double,1> sd1, SIMD<double,1> sd2, SIMD<double,1> sd3, SIMD<double,1> sd4)
  {
    return std::make_tuple(sd1.Data(), sd2.Data(), sd3.Data(), sd4.Data());
  }


/////////////////////////////////////////////////////////////////////////////
// SSE - Simd width 2
/////////////////////////////////////////////////////////////////////////////
#ifdef __SSE__
  template<>
  class alignas(16) SIMD<double,2> //  : public AlignedAlloc<SIMD<double,4>>
  {
    __m128d data;

  public:
    static constexpr int Size() { return 2; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (__m128d adata)
      : data(adata)
      { ; }

    // only called if T has a call operator of appropriate type
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = _mm_set_pd(func(1), func(0));
    }

    // only called if T is arithmetic (integral or floating point types)
    template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    SIMD (const T & val)
    {
      data = _mm_set1_pd(val);
    }

    SIMD (double const * p)
    {
      data = _mm_loadu_pd(p);
    }

    NG_INLINE operator __m128d() const { return data; }
    NG_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    NG_INLINE double& operator[] (int i) { return ((double*)(&data))[i]; }
    NG_INLINE __m128d Data() const { return data; }
    NG_INLINE __m128d & Data() { return data; }

    // NG_INLINE operator std::tuple<double&,double&,double&,double&> ()
    // { return std::tuple<double&,double&,double&,double&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }


    NG_INLINE SIMD<double,2> &operator+= (SIMD<double,2> b) { data+=b.Data(); return *this; }
    NG_INLINE SIMD<double,2> &operator-= (SIMD<double,2> b) { data-=b.Data(); return *this; }
    NG_INLINE SIMD<double,2> &operator*= (SIMD<double,2> b) { data*=b.Data(); return *this; }
    NG_INLINE SIMD<double,2> &operator/= (SIMD<double,2> b) { data/=b.Data(); return *this; }

  };

  NG_INLINE SIMD<double,2> operator+ (SIMD<double,2> a, SIMD<double,2> b) { return a.Data()+b.Data(); }
  NG_INLINE SIMD<double,2> operator- (SIMD<double,2> a, SIMD<double,2> b) { return a.Data()-b.Data(); }
  NG_INLINE SIMD<double,2> operator- (SIMD<double,2> a) { return -a.Data(); }
  NG_INLINE SIMD<double,2> operator* (SIMD<double,2> a, SIMD<double,2> b) { return a.Data()*b.Data(); }
  NG_INLINE SIMD<double,2> operator/ (SIMD<double,2> a, SIMD<double,2> b) { return a.Data()/b.Data(); }

  /*
  NG_INLINE SIMD<double,4> sqrt (SIMD<double,4> a) { return _mm256_sqrt_pd(a.Data()); }
  NG_INLINE SIMD<double,4> floor (SIMD<double,4> a) { return _mm256_floor_pd(a.Data()); }
  NG_INLINE SIMD<double,4> ceil (SIMD<double,4> a) { return _mm256_ceil_pd(a.Data()); }
  NG_INLINE SIMD<double,4> fabs (SIMD<double,4> a) { return _mm256_max_pd(a.Data(), -a.Data()); }
  NG_INLINE SIMD<double,4> L2Norm2 (SIMD<double,4> a) { return a.Data()*a.Data(); }
  NG_INLINE SIMD<double,4> Trans (SIMD<double,4> a) { return a; }
  NG_INLINE SIMD<double,4> IfPos (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_GT_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  NG_INLINE double HSum (SIMD<double,4> sd)
  {
    __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }

  NG_INLINE auto HSum (SIMD<double,4> sd1, SIMD<double,4> sd2)
  {
    __m256d hv = _mm256_hadd_pd(sd1.Data(), sd2.Data());
    __m128d hv2 = _mm_add_pd (_mm256_extractf128_pd(hv,0), _mm256_extractf128_pd(hv,1));
    return std::make_tuple(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  NG_INLINE auto HSum (SIMD<double,4> v1, SIMD<double,4> v2, SIMD<double,4> v3, SIMD<double,4> v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1.Data(), v2.Data());
    __m256d hsum2 = _mm256_hadd_pd (v3.Data(), v4.Data());
    __m256d hsum = _mm256_add_pd (_mm256_permute2f128_pd (hsum1, hsum2, 1+2*16),
                                  _mm256_blend_pd (hsum1, hsum2, 12));
    return SIMD<double,4>(hsum);
  }
  */
#endif // __SSE__

  

/////////////////////////////////////////////////////////////////////////////
// AVX - Simd width 4
/////////////////////////////////////////////////////////////////////////////
#ifdef __AVX__
  template<>
  class alignas(32) SIMD<double,4> : public AlignedAlloc<SIMD<double,4>>
  {
    __m256d data;

  public:
    static constexpr int Size() { return 4; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (__m256d adata)
      : data(adata)
      { ; }

    // only called if T has a call operator of appropriate type
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = _mm256_set_pd(func(3), func(2), func(1), func(0));
    }

    // only called if T is arithmetic (integral or floating point types)
    template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    SIMD (const T & val)
    {
      data = _mm256_set1_pd(val);
    }

    SIMD (double const * p)
    {
      data = _mm256_loadu_pd(p);
    }

    NG_INLINE operator __m256d() const { return data; }
    NG_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    NG_INLINE double& operator[] (int i) { return ((double*)(&data))[i]; }
    NG_INLINE __m256d Data() const { return data; }
    NG_INLINE __m256d & Data() { return data; }

    NG_INLINE operator std::tuple<double&,double&,double&,double&> ()
    { return std::tuple<double&,double&,double&,double&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }


    NG_INLINE SIMD<double,4> &operator+= (SIMD<double,4> b) { data+=b.Data(); return *this; }
    NG_INLINE SIMD<double,4> &operator-= (SIMD<double,4> b) { data-=b.Data(); return *this; }
    NG_INLINE SIMD<double,4> &operator*= (SIMD<double,4> b) { data*=b.Data(); return *this; }
    NG_INLINE SIMD<double,4> &operator/= (SIMD<double,4> b) { data/=b.Data(); return *this; }

  };

  NG_INLINE SIMD<double,4> operator+ (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()+b.Data(); }
  NG_INLINE SIMD<double,4> operator- (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()-b.Data(); }
  NG_INLINE SIMD<double,4> operator- (SIMD<double,4> a) { return -a.Data(); }
  NG_INLINE SIMD<double,4> operator* (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()*b.Data(); }
  NG_INLINE SIMD<double,4> operator/ (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()/b.Data(); }

  NG_INLINE SIMD<double,4> sqrt (SIMD<double,4> a) { return _mm256_sqrt_pd(a.Data()); }
  NG_INLINE SIMD<double,4> floor (SIMD<double,4> a) { return _mm256_floor_pd(a.Data()); }
  NG_INLINE SIMD<double,4> ceil (SIMD<double,4> a) { return _mm256_ceil_pd(a.Data()); }
  NG_INLINE SIMD<double,4> fabs (SIMD<double,4> a) { return _mm256_max_pd(a.Data(), -a.Data()); }
  NG_INLINE SIMD<double,4> L2Norm2 (SIMD<double,4> a) { return a.Data()*a.Data(); }
  NG_INLINE SIMD<double,4> Trans (SIMD<double,4> a) { return a; }
  NG_INLINE SIMD<double,4> IfPos (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_GT_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  NG_INLINE double HSum (SIMD<double,4> sd)
  {
    __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }

  NG_INLINE auto HSum (SIMD<double,4> sd1, SIMD<double,4> sd2)
  {
    __m256d hv = _mm256_hadd_pd(sd1.Data(), sd2.Data());
    __m128d hv2 = _mm_add_pd (_mm256_extractf128_pd(hv,0), _mm256_extractf128_pd(hv,1));
    return std::make_tuple(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  NG_INLINE auto HSum (SIMD<double,4> v1, SIMD<double,4> v2, SIMD<double,4> v3, SIMD<double,4> v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1.Data(), v2.Data());
    __m256d hsum2 = _mm256_hadd_pd (v3.Data(), v4.Data());
    __m256d hsum = _mm256_add_pd (_mm256_permute2f128_pd (hsum1, hsum2, 1+2*16),
                                  _mm256_blend_pd (hsum1, hsum2, 12));
    return SIMD<double,4>(hsum);
  }

#endif // __AVX__

/////////////////////////////////////////////////////////////////////////////
// AVX512 - Simd width 8
/////////////////////////////////////////////////////////////////////////////
#ifdef __AVX512F__
  template<>
  class alignas(64) SIMD<double,8> : public AlignedAlloc<SIMD<double,8>>
  {
    __m512d data;

  public:
    static constexpr int Size() { return 8; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (__m512d adata)
      : data(adata)
      { ; }

    // only called if T has a call operator of appropriate type
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data = _mm512_set_pd(func(7), func(6), func(5), func(4),
                           func(3), func(2), func(1), func(0));
    }

    // only called if T is arithmetic (integral or floating point types)
    template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    SIMD (const T & val)
    {
      data = _mm512_set1_pd(val);
    }

    SIMD (double const * p)
    {
      data = _mm512_loadu_pd(p);
    }

    NG_INLINE operator __m512d() const { return data; }
    NG_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    NG_INLINE __m512d Data() const { return data; }
    NG_INLINE __m512d & Data() { return data; }

    NG_INLINE SIMD<double,8> &operator+= (SIMD<double,8> b) { data+=b.Data(); return *this; }
    NG_INLINE SIMD<double,8> &operator-= (SIMD<double,8> b) { data-=b.Data(); return *this; }
    NG_INLINE SIMD<double,8> &operator*= (SIMD<double,8> b) { data*=b.Data(); return *this; }
    NG_INLINE SIMD<double,8> &operator/= (SIMD<double,8> b) { data/=b.Data(); return *this; }

  };

  NG_INLINE SIMD<double,8> operator- (SIMD<double,8> a) { return _mm512_sub_pd(_mm512_setzero_pd(), a.Data()); }

  NG_INLINE SIMD<double,8> operator+ (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_add_pd(a.Data(),b.Data()); }
  NG_INLINE SIMD<double,8> operator- (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_sub_pd(a.Data(),b.Data()); }
  NG_INLINE SIMD<double,8> operator* (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_mul_pd(a.Data(),b.Data()); }
  NG_INLINE SIMD<double,8> operator/ (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_div_pd(a.Data(),b.Data()); }

  NG_INLINE SIMD<double,8> sqrt (SIMD<double,8> a) { return _mm512_sqrt_pd(a.Data()); }
  NG_INLINE SIMD<double,8> floor (SIMD<double,8> a) { return _mm512_floor_pd(a.Data()); }
  NG_INLINE SIMD<double,8> ceil (SIMD<double,8> a) { return _mm512_ceil_pd(a.Data()); }
  NG_INLINE SIMD<double,8> fabs (SIMD<double,8> a) { return _mm512_max_pd(a.Data(), -a.Data()); }
  NG_INLINE SIMD<double,8> L2Norm2 (SIMD<double,8> a) { return a.Data()*a.Data(); }
  NG_INLINE SIMD<double,8> Trans (SIMD<double,8> a) { return a; }
  NG_INLINE SIMD<double,8> IfPos (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    auto cp = _mm512_cmp_pd_mask (a.Data(), _mm512_setzero_pd(), _MM_CMPINT_GT);
    return _mm512_mask_blend_pd(cp, c.Data(), b.Data());
  }


  template<> NG_INLINE auto FMA (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fmadd_pd (a.Data(), b.Data(), c.Data());
  }

  NG_INLINE double HSum (SIMD<double,8> sd)
  {
    SIMD<double,4> low = _mm512_extractf64x4_pd(sd.Data(),0);
    SIMD<double,4> high = _mm512_extractf64x4_pd(sd.Data(),1);
    return HSum(low)+HSum(high);
  }

  NG_INLINE auto HSum (SIMD<double,8> sd1, SIMD<double,8> sd2)
  {
    return std::make_tuple(HSum(sd1), HSum(sd2));
  }

  NG_INLINE SIMD<double,4> HSum (SIMD<double,8> v1, SIMD<double,8> v2, SIMD<double,8> v3, SIMD<double,8> v4)
  {
    SIMD<double,4> high1 = _mm512_extractf64x4_pd(v1.Data(),1);
    SIMD<double,4> high2 = _mm512_extractf64x4_pd(v2.Data(),1);
    SIMD<double,4> high3 = _mm512_extractf64x4_pd(v3.Data(),1);
    SIMD<double,4> high4 = _mm512_extractf64x4_pd(v4.Data(),1);
    SIMD<double,4> low1 = _mm512_extractf64x4_pd(v1.Data(),0);
    SIMD<double,4> low2 = _mm512_extractf64x4_pd(v2.Data(),0);
    SIMD<double,4> low3 = _mm512_extractf64x4_pd(v3.Data(),0);
    SIMD<double,4> low4 = _mm512_extractf64x4_pd(v4.Data(),0);
    return HSum(low1,low2,low3,low4) + HSum(high1,high2,high3,high4);
  }
#endif // __AVX512F__


////////////////////////////////////////////////////////////////////////////////
// MultiSIMD - Multiple SIMD values in one struct using head-tail implementation
////////////////////////////////////////////////////////////////////////////////
  template <int D, typename T>
  class MultiSIMD : public AlignedAlloc<MultiSIMD<D,T>>
  {
    SIMD<T> head;
    MultiSIMD<D-1,T> tail;
  public:
    MultiSIMD () = default;
    MultiSIMD (const MultiSIMD & ) = default;
    MultiSIMD (T v) : head(v), tail(v) { ; }
    MultiSIMD (SIMD<T> _head, MultiSIMD<D-1,T> _tail)
      : head(_head), tail(_tail) { ; }
    SIMD<T> Head() const { return head; }
    MultiSIMD<D-1,T> Tail() const { return tail; }
    SIMD<T> & Head() { return head; }
    MultiSIMD<D-1,T> & Tail() { return tail; }

    template <int NR>
    SIMD<T> Get() const { return NR==0 ? head : tail.template Get<NR-1>(); }
    template <int NR>
    SIMD<T> & Get() { return NR==0 ? head : tail.template Get<NR-1>(); }
    auto MakeTuple() { return std::tuple_cat(std::tuple<SIMD<T>&> (head), tail.MakeTuple()); }
    // not yet possible for MSVC
    // operator auto () { return MakeTuple(); }
  };

  template <typename T>
  class MultiSIMD<2,T> : public AlignedAlloc<MultiSIMD<2,T>>
  {
    SIMD<T> v0, v1;
  public:
    MultiSIMD () = default;
    MultiSIMD (const MultiSIMD & ) = default;
    MultiSIMD (T v) : v0(v), v1(v) { ; }
    MultiSIMD (SIMD<T> _v0, SIMD<T> _v1) : v0(_v0), v1(_v1) { ; }

    SIMD<T> Head() const { return v0; }
    SIMD<T> Tail() const { return v1; }
    SIMD<T> & Head() { return v0; }
    SIMD<T> & Tail() { return v1; }

    template <int NR>
    SIMD<T> Get() const { return NR==0 ? v0 : v1; }
    template <int NR>
    SIMD<T> & Get() { return NR==0 ? v0 : v1; }
    auto MakeTuple() { return std::tuple<SIMD<T>&, SIMD<T>&> (v0, v1); }
    operator std::tuple<SIMD<T>&, SIMD<T>&>() { return MakeTuple(); }
  };

  template <int D> NG_INLINE MultiSIMD<D,double> operator+ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()+b.Head(), a.Tail()+b.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator+ (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator+ (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }

  template <int D> NG_INLINE MultiSIMD<D,double> operator- (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()-b.Head(), a.Tail()-b.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator- (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a-b.Head(), a-b.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator- (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (b.Head()-a, b.Tail()-a); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator- (MultiSIMD<D,double> a)
  { return MultiSIMD<D,double> (-a.Head(), -a.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator* (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()*b.Head(), a.Tail()*b.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator/ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()/b.Head(), a.Tail()/b.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator* (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }
  template <int D> NG_INLINE MultiSIMD<D,double> operator* (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }

  template <int D> NG_INLINE MultiSIMD<D,double> & operator+= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()+=b.Head(); a.Tail()+=b.Tail(); return a; }
  template <int D> NG_INLINE MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, double b)
  { a.Head()-=b; a.Tail()-=b; return a; }
  template <int D> NG_INLINE MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()-=b.Head(); a.Tail()-=b.Tail(); return a; }
  template <int D> NG_INLINE MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()*=b.Head(); a.Tail()*=b.Tail(); return a; }
  template <int D> NG_INLINE MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, double b)
  { a.Head()*=b; a.Tail()*=b; return a; }
  // NG_INLINE MultiSIMD<double> operator/= (MultiSIMD<double> & a, MultiSIMD<double> b) { return a.Data()/=b.Data(); }

  NG_INLINE SIMD<double> HVSum (SIMD<double> a) { return a; }
  template <int D>
  NG_INLINE SIMD<double> HVSum (MultiSIMD<D,double> a) { return a.Head() + HVSum(a.Tail()); }

  template <int D> NG_INLINE double HSum (MultiSIMD<D,double> a) { return HSum(HVSum(a)); }
  template <int D> NG_INLINE auto HSum (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return HSum(HVSum(a), HVSum(b)); }

  template <int D, typename T>
  std::ostream & operator<< (std::ostream & ost, MultiSIMD<D,T> multi)
  {
    ost << multi.Head() << " " << multi.Tail();
    return ost;
  }

  template <typename T>
  std::ostream & operator<< (std::ostream & ost, SIMD<T> simd)
  {
    ost << simd[0];
    for (int i = 1; i < simd.Size(); i++)
      ost << " " << simd[i];
    return ost;
  }
}

namespace netgen
{
  using namespace ngsimd;
}
#endif

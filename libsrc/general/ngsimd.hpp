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
inline __m128d operator- (__m128d a) { return _mm_xor_pd(a, _mm_set1_pd(-0.0)); }
inline __m128d operator+ (__m128d a, __m128d b) { return _mm_add_pd(a,b); }
inline __m128d operator- (__m128d a, __m128d b) { return _mm_sub_pd(a,b); }
inline __m128d operator* (__m128d a, __m128d b) { return _mm_mul_pd(a,b); }
inline __m128d operator/ (__m128d a, __m128d b) { return _mm_div_pd(a,b); }
inline __m128d operator* (double a, __m128d b) { return _mm_set1_pd(a)*b; }
inline __m128d operator* (__m128d b, double a) { return _mm_set1_pd(a)*b; }

inline __m128d operator+= (__m128d &a, __m128d b) { return a = a+b; }
inline __m128d operator-= (__m128d &a, __m128d b) { return a = a-b; }
inline __m128d operator*= (__m128d &a, __m128d b) { return a = a*b; }
inline __m128d operator/= (__m128d &a, __m128d b) { return a = a/b; }

inline __m256d operator- (__m256d a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }
inline __m256d operator+ (__m256d a, __m256d b) { return _mm256_add_pd(a,b); }
inline __m256d operator- (__m256d a, __m256d b) { return _mm256_sub_pd(a,b); }
inline __m256d operator* (__m256d a, __m256d b) { return _mm256_mul_pd(a,b); }
inline __m256d operator/ (__m256d a, __m256d b) { return _mm256_div_pd(a,b); }
inline __m256d operator* (double a, __m256d b) { return _mm256_set1_pd(a)*b; }
inline __m256d operator* (__m256d b, double a) { return _mm256_set1_pd(a)*b; }

inline __m256d operator+= (__m256d &a, __m256d b) { return a = a+b; }
inline __m256d operator-= (__m256d &a, __m256d b) { return a = a-b; }
inline __m256d operator*= (__m256d &a, __m256d b) { return a = a*b; }
inline __m256d operator/= (__m256d &a, __m256d b) { return a = a/b; }
#endif
#endif



namespace ngsimd
{
  constexpr int GetDefaultSIMDSize() {
#if defined __AVX512F__
    return 8;
#elif defined __AVX__
    return 4;
#else
    return 1;
#endif
  }

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
  inline auto FMA(T1 a, T2 b, T3 c)
  {
    return a*b+c;
  }

  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator+ (T a, SIMD<double,N> b) { return SIMD<double,N>(a) + b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator- (T a, SIMD<double,N> b) { return SIMD<double,N>(a) - b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator* (T a, SIMD<double,N> b) { return SIMD<double,N>(a) * b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator/ (T a, SIMD<double,N> b) { return SIMD<double,N>(a) / b; }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator+ (SIMD<double,N> a, T b) { return a + SIMD<double,N>(b); }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator- (SIMD<double,N> a, T b) { return a - SIMD<double,N>(b); }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator* (SIMD<double,N> a, T b) { return a * SIMD<double,N>(b); }
  template<int N, typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    inline SIMD<double,N> operator/ (SIMD<double,N> a, T b) { return a / SIMD<double,N>(b); }


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
  template<int N> inline SIMD<double,N> exp (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return exp(a[i]); } );
  }

  using std::log;
  template<int N> inline SIMD<double,N> log (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return log(a[i]); } );
  }

  using std::pow;
  template<int N> inline SIMD<double,N> pow (SIMD<double,N> a, double x)
  {
    return SIMD<double,N>([&](int i)->double { return pow(a[i],x); } );
  }

  using std::sin;
  template<int N> inline SIMD<double,N> sin (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return sin(a[i]); } );
  }

  using std::cos;
  template<int N> inline SIMD<double,N> cos (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return cos(a[i]); } );
  }

  using std::tan;
  template<int N> inline SIMD<double,N> tan (SIMD<double,N> a)
  {
    return SIMD<double,N>([&](int i)->double { return tan(a[i]); } );
  }

  using std::atan;
  template<int N> inline SIMD<double,N> atan (SIMD<double,N> a)
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

    inline operator double() const { return data; }
    inline double operator[] (int i) const { return ((double*)(&data))[i]; }
    inline double Data() const { return data; }
    inline double & Data() { return data; }

    inline SIMD<double,1> &operator+= (SIMD<double,1> b) { data+=b.Data(); return *this; }
    inline SIMD<double,1> &operator-= (SIMD<double,1> b) { data-=b.Data(); return *this; }
    inline SIMD<double,1> &operator*= (SIMD<double,1> b) { data*=b.Data(); return *this; }
    inline SIMD<double,1> &operator/= (SIMD<double,1> b) { data/=b.Data(); return *this; }

  };

  inline SIMD<double,1> operator+ (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()+b.Data(); }
  inline SIMD<double,1> operator- (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()-b.Data(); }
  inline SIMD<double,1> operator- (SIMD<double,1> a) { return -a.Data(); }
  inline SIMD<double,1> operator* (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()*b.Data(); }
  inline SIMD<double,1> operator/ (SIMD<double,1> a, SIMD<double,1> b) { return a.Data()/b.Data(); }

  inline SIMD<double,1> sqrt (SIMD<double,1> a) { return std::sqrt(a.Data()); }
  inline SIMD<double,1> fabs (SIMD<double,1> a) { return std::fabs(a.Data()); }
  inline SIMD<double,1> L2Norm2 (SIMD<double,1> a) { return a.Data()*a.Data(); }
  inline SIMD<double,1> Trans (SIMD<double,1> a) { return a; }
  inline SIMD<double,1> IfPos (SIMD<double,1> a, SIMD<double,1> b, SIMD<double,1> c)
  {
    return (a.Data() > 0) ? b : c;
  }

  inline double HSum (SIMD<double,1> sd)
  {
    return sd.Data();
  }

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

    inline operator __m256d() const { return data; }
    inline double operator[] (int i) const { return ((double*)(&data))[i]; }
    inline __m256d Data() const { return data; }
    inline __m256d & Data() { return data; }

    inline SIMD<double,4> &operator+= (SIMD<double,4> b) { data+=b.Data(); return *this; }
    inline SIMD<double,4> &operator-= (SIMD<double,4> b) { data-=b.Data(); return *this; }
    inline SIMD<double,4> &operator*= (SIMD<double,4> b) { data*=b.Data(); return *this; }
    inline SIMD<double,4> &operator/= (SIMD<double,4> b) { data/=b.Data(); return *this; }

  };

  inline SIMD<double,4> operator+ (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()+b.Data(); }
  inline SIMD<double,4> operator- (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()-b.Data(); }
  inline SIMD<double,4> operator- (SIMD<double,4> a) { return -a.Data(); }
  inline SIMD<double,4> operator* (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()*b.Data(); }
  inline SIMD<double,4> operator/ (SIMD<double,4> a, SIMD<double,4> b) { return a.Data()/b.Data(); }

  inline SIMD<double,4> sqrt (SIMD<double,4> a) { return _mm256_sqrt_pd(a.Data()); }
  inline SIMD<double,4> fabs (SIMD<double,4> a) { return _mm256_max_pd(a.Data(), -a.Data()); }
  inline SIMD<double,4> L2Norm2 (SIMD<double,4> a) { return a.Data()*a.Data(); }
  inline SIMD<double,4> Trans (SIMD<double,4> a) { return a; }
  inline SIMD<double,4> IfPos (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_GT_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  inline double HSum (SIMD<double,4> sd)
  {
    __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }

  inline auto HSum (SIMD<double,4> sd1, SIMD<double,4> sd2)
  {
    __m256d hv = _mm256_hadd_pd(sd1.Data(), sd2.Data());
    __m128d hv2 = _mm_add_pd (_mm256_extractf128_pd(hv,0), _mm256_extractf128_pd(hv,1));
    return std::make_tuple(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  inline SIMD<double,4> HSum (SIMD<double,4> v1, SIMD<double,4> v2, SIMD<double,4> v3, SIMD<double,4> v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1.Data(), v2.Data());
    __m256d hsum2 = _mm256_hadd_pd (v3.Data(), v4.Data());
    __m256d hsum = _mm256_add_pd (_mm256_permute2f128_pd (hsum1, hsum2, 1+2*16),
                                  _mm256_blend_pd (hsum1, hsum2, 12));
    return hsum;
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

    inline operator __m512d() const { return data; }
    inline double operator[] (int i) const { return ((double*)(&data))[i]; }
    inline __m512d Data() const { return data; }
    inline __m512d & Data() { return data; }

    inline SIMD<double,8> &operator+= (SIMD<double,8> b) { data+=b.Data(); return *this; }
    inline SIMD<double,8> &operator-= (SIMD<double,8> b) { data-=b.Data(); return *this; }
    inline SIMD<double,8> &operator*= (SIMD<double,8> b) { data*=b.Data(); return *this; }
    inline SIMD<double,8> &operator/= (SIMD<double,8> b) { data/=b.Data(); return *this; }

  };

  inline SIMD<double,8> operator- (SIMD<double,8> a) { return _mm512_sub_pd(_mm512_setzero_pd(), a.Data()); }

  inline SIMD<double,8> operator+ (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_add_pd(a.Data(),b.Data()); }
  inline SIMD<double,8> operator- (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_sub_pd(a.Data(),b.Data()); }
  inline SIMD<double,8> operator* (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_mul_pd(a.Data(),b.Data()); }
  inline SIMD<double,8> operator/ (SIMD<double,8> a, SIMD<double,8> b) { return _mm512_div_pd(a.Data(),b.Data()); }

  inline SIMD<double,8> sqrt (SIMD<double,8> a) { return _mm512_sqrt_pd(a.Data()); }
  inline SIMD<double,8> fabs (SIMD<double,8> a) { return _mm512_max_pd(a.Data(), -a.Data()); }
  inline SIMD<double,8> L2Norm2 (SIMD<double,8> a) { return a.Data()*a.Data(); }
  inline SIMD<double,8> Trans (SIMD<double,8> a) { return a; }
  inline SIMD<double,8> IfPos (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    auto cp = _mm512_cmp_pd_mask (a.Data(), _mm512_setzero_pd(), _MM_CMPINT_GT);
    return _mm512_mask_blend_pd(cp, c.Data(), b.Data());
  }


  template<> inline auto FMA (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fmadd_pd (a.Data(), b.Data(), c.Data());
  }

  inline double HSum (SIMD<double,8> sd)
  {
    SIMD<double,4> low = _mm512_extractf64x4_pd(sd.Data(),0);
    SIMD<double,4> high = _mm512_extractf64x4_pd(sd.Data(),1);
    return HSum(low)+HSum(high);
  }

  inline auto HSum (SIMD<double,8> sd1, SIMD<double,8> sd2)
  {
    return std::make_tuple(HSum(sd1), HSum(sd2));
  }

  inline SIMD<double,4> HSum (SIMD<double,8> v1, SIMD<double,8> v2, SIMD<double,8> v3, SIMD<double,8> v4)
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
  };

  template <int D> inline MultiSIMD<D,double> operator+ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()+b.Head(), a.Tail()+b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator+ (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator+ (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }

  template <int D> inline MultiSIMD<D,double> operator- (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()-b.Head(), a.Tail()-b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator- (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a-b.Head(), a-b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator- (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (b.Head()-a, b.Tail()-a); }
  template <int D> inline MultiSIMD<D,double> operator- (MultiSIMD<D,double> a)
  { return MultiSIMD<D,double> (-a.Head(), -a.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator* (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()*b.Head(), a.Tail()*b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator/ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()/b.Head(), a.Tail()/b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator* (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator* (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }

  template <int D> inline MultiSIMD<D,double> & operator+= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()+=b.Head(); a.Tail()+=b.Tail(); return a; }
  template <int D> inline MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, double b)
  { a.Head()-=b; a.Tail()-=b; return a; }
  template <int D> inline MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()-=b.Head(); a.Tail()-=b.Tail(); return a; }
  template <int D> inline MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()*=b.Head(); a.Tail()*=b.Tail(); return a; }
  template <int D> inline MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, double b)
  { a.Head()*=b; a.Tail()*=b; return a; }
  // inline MultiSIMD<double> operator/= (MultiSIMD<double> & a, MultiSIMD<double> b) { return a.Data()/=b.Data(); }

  inline SIMD<double> HVSum (SIMD<double> a) { return a; }
  template <int D>
  inline SIMD<double> HVSum (MultiSIMD<D,double> a) { return a.Head() + HVSum(a.Tail()); }

  template <int D> inline double HSum (MultiSIMD<D,double> a) { return HSum(HVSum(a)); }
  template <int D> inline auto HSum (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
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

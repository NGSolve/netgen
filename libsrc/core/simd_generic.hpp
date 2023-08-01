#ifndef NETGEN_CORE_SIMD_GENERIC_HPP
#define NETGEN_CORE_SIMD_GENERIC_HPP

/**************************************************************************/
/* File:   simd_base.hpp                                                  */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <type_traits>
#include <functional>
#include <tuple>
#include <cmath>

#include "array.hpp"

namespace ngcore
{
#if defined __AVX512F__
    #define NETGEN_DEFAULT_SIMD_SIZE 8
#elif defined __AVX__
    #define NETGEN_DEFAULT_SIMD_SIZE 4
#else
    #define NETGEN_DEFAULT_SIMD_SIZE 2
#endif

  constexpr int GetDefaultSIMDSize() {
      return NETGEN_DEFAULT_SIMD_SIZE;
  }

  constexpr bool IsNativeSIMDSize(int n) {
    if(n==1) return true;
    if(n==2) return true;
#if defined __AVX__
    if(n==4) return true;
#endif
#if defined __AVX512F__
    if(n==8) return true;
#endif
    return false;
  }

  // split n = k+l such that k is the largest natively supported simd size < n
  constexpr int GetLargestNativeSIMDPart(int n) {
      int k = n-1;
      while(!IsNativeSIMDSize(k))
          k--;
      return k;
  }


  template <typename T, int N=GetDefaultSIMDSize()> class SIMD;

  class mask64;

  ////////////////////////////////////////////////////////////////////////////
  namespace detail {
    template <typename T, size_t N, size_t... I>
    auto array_range_impl(std::array<T, N> const& arr,
                   size_t first,
                   std::index_sequence<I...>)
    -> std::array<T, sizeof...(I)> {
        return {arr[first + I]...};
    }

    template <size_t S, typename T, size_t N>
    auto array_range(std::array<T, N> const& arr, size_t first) {
      return array_range_impl(arr, first, std::make_index_sequence<S>{});
    }
  
  } // namespace detail

  ////////////////////////////////////////////////////////////////////////////
  // mask

  template <>
  class SIMD<mask64,1>
  {
    int64_t mask;
  public:
    SIMD (int64_t i)
      : mask(i > 0 ? -1 : 0) { ; }
    bool Data() const { return mask; }
    static constexpr int Size() { return 1; }
    auto operator[] (int /* i */) const { return mask; }
  };


  template <int N>
  class alignas(GetLargestNativeSIMDPart(N)*sizeof(int64_t)) SIMD<mask64,N>
  {
    static constexpr int N1 = GetLargestNativeSIMDPart(N);
    static constexpr int N2 = N-N1;

    SIMD<mask64,N1> lo;
    SIMD<mask64,N2> hi;
  public:

    SIMD (int64_t i) : lo(i), hi(i-N1 ) { ; }
    SIMD (SIMD<mask64,N1> lo_, SIMD<mask64,N2> hi_) : lo(lo_), hi(hi_) { ; }
    SIMD<mask64,N1> Lo() const { return lo; }
    SIMD<mask64,N2> Hi() const { return hi; }
    static constexpr int Size() { return N; }
  };

  template<int N>
  NETGEN_INLINE SIMD<mask64,N> operator&& (SIMD<mask64,N> a, SIMD<mask64,N> b)
    {
      if constexpr(N==1) return a.Data() && b.Data();
      else               return { a.Lo() && b.Lo(), a.Hi() && b.Hi() };
    }


  ////////////////////////////////////////////////////////////////////////////
  // int64

  template<>
  class SIMD<int64_t,1>
  {
    int64_t data;

  public:
    static constexpr int Size() { return 1; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;
    SIMD (int val) : data{val} {}
    SIMD (int64_t val) : data{val} {}
    SIMD (size_t val) : data(val) {}
    explicit SIMD (std::array<int64_t, 1> arr)
        : data{arr[0]}
    {}

    int64_t operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    auto Data() const { return data; }
    static SIMD FirstInt(int64_t n0=0) { return {n0}; }
    template <int I>
    int64_t Get()
    {
      static_assert(I==0);
      return data;
    }
  };

  template<int N>
  class alignas(GetLargestNativeSIMDPart(N)*sizeof(int64_t)) SIMD<int64_t,N>
  {
    static constexpr int N1 = GetLargestNativeSIMDPart(N);
    static constexpr int N2 = N-N1;

    SIMD<int64_t,N1> lo;
    SIMD<int64_t,N2> high;

  public:
    static constexpr int Size() { return N; }

    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (int val) : lo{val}, high{val} { ; }
    SIMD (int64_t val) : lo{val}, high{val} { ; }
    SIMD (size_t val) : lo{val}, high{val} { ; }
    SIMD (SIMD<int64_t,N1> lo_, SIMD<int64_t,N2> high_) : lo(lo_), high(high_) { ; }

    explicit SIMD( std::array<int64_t, N> arr )
        : lo(detail::array_range<N1>(arr, 0)),
          high(detail::array_range<N2>(arr, N1))
      {}

    template<typename ...T>
    explicit SIMD(const T... vals)
    : lo(detail::array_range<N1>(std::array<int64_t, N>{vals...}, 0)),
      high(detail::array_range<N2>(std::array<int64_t, N>{vals...}, N1))
      {
        static_assert(sizeof...(vals)==N, "wrong number of arguments");
      }


    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<int64_t(int)>>::value, int>::type = 0>
      SIMD (const T & func)
    {
      for(auto i : IntRange(N1))
          lo[i] = func(i);
      for(auto i : IntRange(N2))
          high[i] = func(N1+i);
    }

    auto Lo() const { return lo; }
    auto Hi() const { return high; }

    int64_t operator[] (int i) const { return ((int64_t*)(&lo))[i]; }

    /*
    operator tuple<int64_t&,int64_t&,int64_t&,int64_t&> ()
    { return tuple<int64_t&,int64_t&,int64_t&,int64_t&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }
    */

    /*
    static SIMD FirstInt() { return { 0, 1, 2, 3 }; }
    */
    static SIMD FirstInt(int64_t n0=0) { return {SIMD<int64_t,N1>::FirstInt(n0), SIMD<int64_t,N2>::FirstInt(n0+N1)}; }
    template <int I>
    int64_t Get()
    {
      static_assert(I>=0 && I<N, "Index out of range");
      if constexpr(I<N1) return lo.template Get<I>();
      else               return high.template Get<I-N1>();
    }
  };


  ////////////////////////////////////////////////////////////////////////////
  // double

  template<>
  class SIMD<double,1>
  {
    double data;

  public:
    static constexpr int Size() { return 1; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;
    SIMD (double val) { data = val; }
    SIMD (int val)    { data = val; }
    SIMD (size_t val) { data = val; }
    SIMD (double const * p) { data = *p; }
    SIMD (double const * p, SIMD<mask64,1> mask) { data = mask.Data() ? *p : 0.0; }
    explicit SIMD (std::array<double, 1> arr)
        : data{arr[0]}
    {}

    template <typename T, typename std::enable_if<std::is_convertible<T,std::function<double(int)>>::value,int>::type = 0>
    SIMD (const T & func)
    {
      data = func(0);
    }

    template <typename T, typename std::enable_if<std::is_convertible<T,std::function<double(int)>>::value,int>::type = 0>
    SIMD & operator= (const T & func)
    {
      data = func(0);
      return *this;
    }

    void Store (double * p) { *p = data; }
    void Store (double * p, SIMD<mask64,1> mask) { if (mask.Data()) *p = data; }

    double operator[] (int i) const { return ((double*)(&data))[i]; }
    double Data() const { return data; }
    template <int I>
    double Get()
    {
      static_assert(I==0);
      return data;
    }
  };


  template<int N>
  class alignas(GetLargestNativeSIMDPart(N)*sizeof(double)) SIMD<double, N>
  {
    static constexpr int N1 = GetLargestNativeSIMDPart(N);
    static constexpr int N2 = N-N1;

    SIMD<double, N1> lo;
    SIMD<double, N2> high;

  public:
    static constexpr int Size() { return N; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD (SIMD<double,N1> lo_, SIMD<double,N2> hi_) : lo(lo_), high(hi_) { ; }

    template <typename T, typename std::enable_if<std::is_convertible<T,std::function<double(int)>>::value,int>::type = 0>
    SIMD (const T & func)
    {
      double  *p = (double*)this;
      for(auto i : IntRange(N))
          p[i] = func(i);
    }

    template <typename T, typename std::enable_if<std::is_convertible<T,std::function<double(int)>>::value,int>::type = 0>
    SIMD & operator= (const T & func)
    {
      double  *p = (double*)this;
      for(auto i : IntRange(N))
          p[i] = func(i);
      return *this;
    }


    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) : lo{val}, high{val} { ; }
    SIMD (int val)    : lo{val}, high{val} { ; }
    SIMD (size_t val) : lo{val}, high{val} { ; }

    SIMD (double const * p) : lo{p}, high{p+N1} { ; }
    SIMD (double const * p, SIMD<mask64,N> mask)
        : lo{p, mask.Lo()}, high{p+N1, mask.Hi()}
      { }
    SIMD (double * p) : lo{p}, high{p+N1} { ; }
    SIMD (double * p, SIMD<mask64,N> mask)
        : lo{p, mask.Lo()}, high{p+N1, mask.Hi()}
      { }

    explicit SIMD( std::array<double, N> arr )
        : lo(detail::array_range<N1>(arr, 0)),
          high(detail::array_range<N2>(arr, N1))
      {}

    template<typename ...T>
    explicit SIMD(const T... vals)
    : lo(detail::array_range<N1>(std::array<double, N>{vals...}, 0)),
      high(detail::array_range<N2>(std::array<double, N>{vals...}, N1))
      {
        static_assert(sizeof...(vals)==N, "wrong number of arguments");
      }

    void Store (double * p) { lo.Store(p); high.Store(p+N1); }
    void Store (double * p, SIMD<mask64,N> mask)
    {
      lo.Store(p, mask.Lo());
      high.Store(p+N1, mask.Hi());
    }

    auto Lo() const { return lo; }
    auto Hi() const { return high; }

    double operator[] (int i) const { return ((double*)(&lo))[i]; }

    template<typename=std::enable_if<N==2>>
    operator std::tuple<double&,double&> ()
    { 
	double *p = (double*)this;
	return std::tuple<double&,double&>(p[0], p[1]); 
    }

    template<typename=std::enable_if<N==4>>
    operator std::tuple<double&,double&,double&,double&> ()
    { return std::tuple<double&,double&,double&,double&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }

    template <int I>
    double Get()
    {
      static_assert(I>=0 && I<N, "Index out of range");
      if constexpr(I<N1) return lo.template Get<I>();
      else               return high.template Get<I-N1>();
    }
    auto Data() const { return *this; }
  };


  // Generic operators for any arithmetic type/simd width
  template <typename T, int N>
  NETGEN_INLINE SIMD<T,N> operator+ (SIMD<T,N> a, SIMD<T,N> b) {
      if constexpr(N==1) return a.Data()+b.Data();
      else               return { a.Lo()+b.Lo(), a.Hi()+b.Hi() };
  }

  template <typename T, int N>
  NETGEN_INLINE SIMD<T,N> operator- (SIMD<T,N> a, SIMD<T,N> b) {
      if constexpr(N==1) return a.Data()-b.Data();
      else               return { a.Lo()-b.Lo(), a.Hi()-b.Hi() };
  }
  template <typename T, int N>
  NETGEN_INLINE SIMD<T,N> operator- (SIMD<T,N> a) {
      if constexpr(N==1) return -a.Data();
      else               return { -a.Lo(), -a.Hi() };
  }

  template <typename T, int N>
  NETGEN_INLINE SIMD<T,N> operator* (SIMD<T,N> a, SIMD<T,N> b) {
      if constexpr(N==1) return a.Data()*b.Data();
      else               return { a.Lo()*b.Lo(), a.Hi()*b.Hi() };
  }

  template <typename T, int N>
  NETGEN_INLINE SIMD<T,N> operator/ (SIMD<T,N> a, SIMD<T,N> b) {
      if constexpr(N==1) return a.Data()/b.Data();
      else               return { a.Lo()/b.Lo(), a.Hi()/b.Hi() };
  }

  template <typename T, int N>
  NETGEN_INLINE SIMD<mask64,N> operator< (SIMD<T,N> a, SIMD<T,N> b)
    {
      if constexpr(N==1) return a.Data() < b.Data();
      else               return { a.Lo()<b.Lo(), a.Hi()<b.Hi() };
    }

  template <typename T, int N>
  NETGEN_INLINE SIMD<mask64,N> operator<= (SIMD<T,N> a, SIMD<T,N> b)
    {
      if constexpr(N==1) return a.Data() <= b.Data();
      else               return { a.Lo()<=b.Lo(), a.Hi()<=b.Hi() };
    }

  template <typename T, int N>
  NETGEN_INLINE SIMD<mask64,N> operator> (SIMD<T,N> a, SIMD<T,N> b)
    {
      if constexpr(N==1) return a.Data() > b.Data();
      else               return { a.Lo()>b.Lo(), a.Hi()>b.Hi() };
    }

  template <typename T, int N>
  NETGEN_INLINE SIMD<mask64,N> operator>= (SIMD<T,N> a, SIMD<T,N> b)
    {
      if constexpr(N==1) return a.Data() >= b.Data();
      else               return { a.Lo()>=b.Lo(), a.Hi()>=b.Hi() };
    }

  template <typename T, int N>
  NETGEN_INLINE SIMD<mask64,N> operator== (SIMD<T,N> a, SIMD<T,N> b)
    {
      if constexpr(N==1) return a.Data() == b.Data();
      else               return { a.Lo()==b.Lo(), a.Hi()==b.Hi() };
    }

  template <typename T, int N>
  NETGEN_INLINE SIMD<mask64,N> operator!= (SIMD<T,N> a, SIMD<T,N> b)
    {
      if constexpr(N==1) return a.Data() != b.Data();
      else               return { a.Lo()!=b.Lo(), a.Hi()!=b.Hi() };
    }

  // int64_t operators with scalar operand (implement overloads to allow implicit casts for second operand)
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator+ (SIMD<int64_t,N> a, int64_t b) { return a+SIMD<int64_t,N>(b); }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator+ (int64_t a, SIMD<int64_t,N> b) { return SIMD<int64_t,N>(a)+b; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator- (int64_t a, SIMD<int64_t,N> b) { return SIMD<int64_t,N>(a)-b; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator- (SIMD<int64_t,N> a, int64_t b) { return a-SIMD<int64_t,N>(b); }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator* (int64_t a, SIMD<int64_t,N> b) { return SIMD<int64_t,N>(a)*b; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator* (SIMD<int64_t,N> b, int64_t a) { return SIMD<int64_t,N>(a)*b; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator/ (SIMD<int64_t,N> a, int64_t b) { return a/SIMD<int64_t,N>(b); }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> operator/ (int64_t a, SIMD<int64_t,N> b) { return SIMD<int64_t,N>(a)/b; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> & operator+= (SIMD<int64_t,N> & a, SIMD<int64_t,N> b) { a=a+b; return a; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> & operator+= (SIMD<int64_t,N> & a, int64_t b) { a+=SIMD<int64_t,N>(b); return a; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> & operator-= (SIMD<int64_t,N> & a, SIMD<int64_t,N> b) { a = a-b; return a; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> & operator-= (SIMD<int64_t,N> & a, int64_t b) { a-=SIMD<int64_t,N>(b); return a; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> & operator*= (SIMD<int64_t,N> & a, SIMD<int64_t,N> b) { a=a*b; return a; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> & operator*= (SIMD<int64_t,N> & a, int64_t b) { a*=SIMD<int64_t,N>(b); return a; }
  template <int N>
  NETGEN_INLINE SIMD<int64_t,N> & operator/= (SIMD<int64_t,N> & a, SIMD<int64_t,N> b) { a = a/b; return a; }

  // double operators with scalar operand (implement overloads to allow implicit casts for second operand)
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator+ (SIMD<double,N> a, double b) { return a+SIMD<double,N>(b); }
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator+ (double a, SIMD<double,N> b) { return SIMD<double,N>(a)+b; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator- (double a, SIMD<double,N> b) { return SIMD<double,N>(a)-b; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator- (SIMD<double,N> a, double b) { return a-SIMD<double,N>(b); }
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator* (double a, SIMD<double,N> b) { return SIMD<double,N>(a)*b; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator* (SIMD<double,N> b, double a) { return SIMD<double,N>(a)*b; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator/ (SIMD<double,N> a, double b) { return a/SIMD<double,N>(b); }
  template <int N>
  NETGEN_INLINE SIMD<double,N> operator/ (double a, SIMD<double,N> b) { return SIMD<double,N>(a)/b; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> & operator+= (SIMD<double,N> & a, SIMD<double,N> b) { a=a+b; return a; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> & operator+= (SIMD<double,N> & a, double b) { a+=SIMD<double,N>(b); return a; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> & operator-= (SIMD<double,N> & a, SIMD<double,N> b) { a = a-b; return a; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> & operator-= (SIMD<double,N> & a, double b) { a-=SIMD<double,N>(b); return a; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> & operator*= (SIMD<double,N> & a, SIMD<double,N> b) { a=a*b; return a; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> & operator*= (SIMD<double,N> & a, double b) { a*=SIMD<double,N>(b); return a; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> & operator/= (SIMD<double,N> & a, SIMD<double,N> b) { a = a/b; return a; }

  // double functions

  template <int N>
  NETGEN_INLINE SIMD<double,N> L2Norm2 (SIMD<double,N> a) { return a*a; }
  template <int N>
  NETGEN_INLINE SIMD<double,N> Trans (SIMD<double,N> a) { return a; }

  template <int N>
  NETGEN_INLINE double HSum (SIMD<double,N> a)
    {
      if constexpr(N==1)
          return a.Data();
      else
          return HSum(a.Lo()) + HSum(a.Hi());
    }


  template<typename T, int N>
  NETGEN_INLINE SIMD<T,N> IfPos (SIMD<T,N> a, SIMD<T,N> b, SIMD<T,N> c)
    {
      if constexpr(N==1) return a.Data()>0.0 ? b : c;
      else               return { IfPos(a.Lo(), b.Lo(), c.Lo()), IfPos(a.Hi(), b.Hi(), c.Hi())};

    }

  template<typename T, int N>
  NETGEN_INLINE SIMD<T,N> IfZero (SIMD<T,N> a, SIMD<T,N> b, SIMD<T,N> c)
    {
      if constexpr(N==1) return a.Data()==0.0 ? b : c;
      else               return { IfZero(a.Lo(), b.Lo(), c.Lo()), IfZero(a.Hi(), b.Hi(), c.Hi())};

    }

  template<typename T, int N>
  NETGEN_INLINE SIMD<T,N> If (SIMD<mask64,N> a, SIMD<T,N> b, SIMD<T,N> c)
    {
      if constexpr(N==1) return a.Data() ? b : c;
      else               return { If(a.Lo(), b.Lo(), c.Lo()), If(a.Hi(), b.Hi(), c.Hi())};

    }

  // a*b+c
  template <typename T1, typename T2, typename T3>
  NETGEN_INLINE auto FMA(T1 a, T2 b, T3 c)
  {
    return c+a*b;
  }

  template <typename T1, typename T2, typename T3>
  NETGEN_INLINE auto FNMA(T1 a, T2 b, T3 c)
  {
    return c-a*b;
  }

  // update form of fma
  template <int N>
  void FMAasm (SIMD<double,N> a, SIMD<double,N> b, SIMD<double,N> & sum)
  {
    sum = FMA(a,b,sum);
  }

  // update form of fms
  template <int N>
  void FNMAasm (SIMD<double,N> a, SIMD<double,N> b, SIMD<double,N> & sum)
  {
    // sum -= a*b;
    sum = FNMA(a,b,sum);
  }

  // c += a*b    (a0re, a0im, a1re, a1im, ...), 
  template <int N>
  void FMAComplex (SIMD<double,N> a, SIMD<double,N> b, SIMD<double,N> & c)
  {
    auto [are, aim] = Unpack(a, a);
    SIMD<double,N> bswap = SwapPairs(b);
    SIMD<double,N> aim_bswap = aim*bswap;
    c += FMAddSub (are, b, aim_bswap);
  }
  
  template <int i, typename T, int N>
  T get(SIMD<T,N> a) { return a.template Get<i>(); }

  template <int NUM, typename FUNC>
  NETGEN_INLINE void Iterate2 (FUNC f)
  {
    if constexpr (NUM > 1) Iterate2<NUM-1> (f);
    if constexpr (NUM >= 1) f(std::integral_constant<int,NUM-1>());
  }


  template <typename T, int N>
  ostream & operator<< (ostream & ost, SIMD<T,N> simd)
  {
    /*
    ost << simd[0];
    for (int i = 1; i < simd.Size(); i++)
      ost << " " << simd[i];
    */
    Iterate2<simd.Size()> ([&] (auto I) {
        if (I.value != 0) ost << " ";
        ost << get<I.value>(simd);
      });
    return ost;
  }

  using std::sqrt;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> sqrt (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return sqrt(a[i]); } );
  }

  using std::fabs;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> fabs (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return fabs(a[i]); } );
  }

  using std::floor;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> floor (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return floor(a[i]); } );
  }

  using std::ceil;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> ceil (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return ceil(a[i]); } );
  }

  using std::exp;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> exp (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return exp(a[i]); } );
  }

  using std::log;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> log (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return log(a[i]); } );
  }

  using std::erf;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> erf (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return erf(a[i]); } );
  }

  using std::pow;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> pow (ngcore::SIMD<double,N> a, double x) {
    return ngcore::SIMD<double,N>([a,x](int i)->double { return pow(a[i],x); } );
  }

  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> pow (ngcore::SIMD<double,N> a, ngcore::SIMD<double,N> b) {
    return ngcore::SIMD<double,N>([a,b](int i)->double { return pow(a[i],b[i]); } );
  }

  using std::sin;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> sin (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return sin(a[i]); } );
  }

  using std::cos;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> cos (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return cos(a[i]); } );
  }

  using std::tan;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> tan (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return tan(a[i]); } );
  }

  using std::atan;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> atan (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return atan(a[i]); } );
  }

  using std::atan2;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> atan2 (ngcore::SIMD<double,N> y, ngcore::SIMD<double,N> x) {
    return ngcore::SIMD<double,N>([y,x](int i)->double { return atan2(y[i], x[i]); } );
  }

  using std::acos;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> acos (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return acos(a[i]); } );
  }

  using std::asin;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> asin (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return asin(a[i]); } );
  }

  using std::sinh;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> sinh (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return sinh(a[i]); } );
  }

  using std::cosh;
  template <int N>
  NETGEN_INLINE ngcore::SIMD<double,N> cosh (ngcore::SIMD<double,N> a) {
    return ngcore::SIMD<double,N>([a](int i)->double { return cosh(a[i]); } );
  }

  template<int N, typename T>
  using MultiSIMD = SIMD<T, N*GetDefaultSIMDSize()>;

  template<int N>
  NETGEN_INLINE auto Unpack (SIMD<double,N> a, SIMD<double,N> b)
  {
    if constexpr(N==1)
      {
        return std::make_tuple(SIMD<double,N>{a.Data()}, SIMD<double,N>{b.Data()} );
      }
    else if constexpr(N==2)
      {
        return std::make_tuple(SIMD<double,N>{ a.Lo(), b.Lo() },
            SIMD<double,N>{ a.Hi(), b.Hi() });
      }
    else
      {
        auto [a1,b1] = Unpack(a.Lo(), b.Lo());
        auto [a2,b2] = Unpack(a.Hi(), b.Hi());
        return std::make_tuple(SIMD<double,N>{ a1, a2 },
            SIMD<double,N>{ b1, b2 });
      }
  }

  // TODO: specialize for AVX, ... 
  template<int N>
  NETGEN_INLINE auto SwapPairs (SIMD<double,N> a)
  {
    if constexpr(N==1) {
        // static_assert(false);
        return a;
      }
    else if constexpr(N==2) {
        return SIMD<double,N> (a.Hi(), a.Lo());
      }
    else {
      return SIMD<double,N> (SwapPairs(a.Lo()), SwapPairs(a.Hi()));
    }
  }


  template<int N>
  NETGEN_INLINE auto HSum128 (SIMD<double,N> a)
  {
    if constexpr(N==1) {
        // static_assert(false);
        return a;
      }
    else if constexpr(N==2) {
        return a;
      }
    else {
      return HSum128(a.Lo()) + HSum128(a.Hi());
    }
  }

  
  // TODO: specialize for AVX, ... 
  // a*b+-c   (even: -, odd: +)
  template<int N>
  NETGEN_INLINE auto FMAddSub (SIMD<double,N> a, SIMD<double,N> b, SIMD<double,N> c)
  {
    if constexpr(N==1) {
        // static_assert(false);
        return a*b-c;
      }
    else if constexpr(N==2) {
        return SIMD<double,N> (a.Lo()*b.Lo()-c.Lo(),
                               a.Hi()*b.Hi()+c.Hi());
      }
    else {
      return SIMD<double,N> (FMAddSub(a.Lo(), b.Lo(), c.Lo()),
                             FMAddSub(a.Hi(), b.Hi(), c.Hi()));
    }
  }
}


namespace std
{
  // structured binding support
  template <typename T, int N >
  struct tuple_size<ngcore::SIMD<T,N>> : std::integral_constant<std::size_t, N> {};
  template<size_t N, typename T, int M> struct tuple_element<N,ngcore::SIMD<T,M>> { using type = T; };
}

#endif // NETGEN_CORE_SIMD_GENERIC_HPP

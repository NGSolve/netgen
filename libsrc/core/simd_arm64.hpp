#include "arm_neon.h"

namespace ngcore
{

  template <>
  class SIMD<mask64,2>
  {
    int64x2_t mask;
  public:
    SIMD (int i)
    {
      mask[0] = i > 0 ? -1 : 0;
      mask[1] = i > 1 ? -1 : 0;
    }

    SIMD (bool i0, bool i1) { mask[0] = i0 ? -1:0; mask[1] = i1 ? -1 : 0; }
    SIMD (SIMD<mask64,1> i0, SIMD<mask64,1> i1) { mask[0] = i0[0]; mask[1] = i1[0]; }
    SIMD (float64x2_t _data) : mask{_data} { }
    auto Data() const { return mask; }
    static constexpr int Size() { return 2; }
    // static NETGEN_INLINE SIMD<mask64, 2> GetMaskFromBits (unsigned int i);
    int64_t operator[] (int i) const { return mask[i]; }

    template <int I>
    int64_t Get() const { return mask[I]; }
    
    auto Lo() const { return mask[0]; }
    auto Hi() const { return mask[1]; }
  };

  
  template<>
  class SIMD<double,2>
  {
    float64x2_t data;

  public:
    static constexpr int Size() { return 2; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    // SIMD (double v0, double v1) : data{v0,v1} { }
    SIMD (double v0, double v1) : data{vcombine_f64(float64x1_t{v0}, float64x1_t{v1})} { }
    SIMD (std::array<double, 2> arr) : data{arr[0], arr[1]} { } 

    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) : data{val,val} { }
    SIMD (int val)    : data{double(val),double(val)} { }
    SIMD (size_t val) : data{double(val),double(val)} { }

    SIMD (double const * p)
    {
      data = vld1q_f64(p);      
      // data[0] = p[0];
      // data[1] = p[1];
    }
    
    SIMD (double const * p, SIMD<mask64,2> mask)
      {
	data[0] = mask[0] ? p[0] : 0;
	data[1] = mask[1] ? p[1] : 0;
      }
    SIMD (float64x2_t _data) { data = _data; }
    
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {
      data[0] = func(0);
      data[1] = func(1);
    }

    void Store (double * p)
    {
      vst1q_f64(p, data);
      /*
      p[0] = data[0];
      p[1] = data[1];
      */
    }
    
    void Store (double * p, SIMD<mask64,2> mask)
    {
      if (mask[0]) p[0] = data[0];
      if (mask[1]) p[1] = data[1];
    }
    
    // NETGEN_INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
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
    
    double Lo() const { return Get<0>(); } // data[0]; }
    double Hi() const { return Get<1>(); } // data[1]; }
    // double Hi() const { return vget_high_f64(data)[0]; }
  };



  NETGEN_INLINE double HSum (SIMD<double,2> sd)
  {
    return sd.Lo()+sd.Hi();  // sd[0]+sd[1];
  }

  NETGEN_INLINE SIMD<double,2> HSum (SIMD<double,2> a, SIMD<double,2> b)
  {
    // return SIMD<double,2> (a[0]+a[1], b[0]+b[1]);
    return vpaddq_f64(a.Data(), b.Data());
  }

  NETGEN_INLINE SIMD<double,4> HSum(SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c, SIMD<double,2> d)
  {
    return SIMD<double,4> (HSum(a,b), HSum(c,d));
  }

  

  // a*b+c
  NETGEN_INLINE SIMD<double,2> FMA (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  {
    return vmlaq_f64(c.Data(), a.Data(), b.Data());
  }
  NETGEN_INLINE SIMD<double,2> FMA (const double & a, SIMD<double,2> b, SIMD<double,2> c)
  {
    return FMA(SIMD<double,2> (a), b, c);
  }
  // -a*b+c
  NETGEN_INLINE SIMD<double,2> FNMA (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  {
    return vmlsq_f64(c.Data(), a.Data(), b.Data());    
    // return c-a*b;
  }
  NETGEN_INLINE SIMD<double,2> FNMA (const double & a, SIMD<double,2> b, SIMD<double,2> c)
  {
    return FNMA(SIMD<double,2> (a), b, c);
  }

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


  
  NETGEN_INLINE SIMD<double,2> If (SIMD<mask64,2> a, SIMD<double,2> b, SIMD<double,2> c)
  {
    // return { a[0] ? b[0] : c[0], a[1] ? b[1] : c[1] };
    return vbslq_f64(a.Data(), b.Data(), c.Data());
  }
  NETGEN_INLINE SIMD<int64_t,2> If (SIMD<mask64,2> a, SIMD<int64_t,2> b, SIMD<int64_t,2> c)
  {
    return SIMD<int64_t,2> (a[0] ? b[0] : c[0], a[1] ? b[1] : c[1]);
  }

  NETGEN_INLINE SIMD<mask64,2> operator&& (SIMD<mask64,2> a, SIMD<mask64,2> b)
  {
    return vandq_u64 (a.Data(), b.Data());
  }
  
}


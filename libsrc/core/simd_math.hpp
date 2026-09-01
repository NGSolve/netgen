#ifndef NETGEN_CORE_SIMD_MATH_HPP
#define NETGEN_CORE_SIMD_MATH_HPP

#include <tuple>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace ngcore
{

  /*
    based on:
    Stephen L. Moshier: Methods and Programs For Mathematical Functions
    https://www.moshier.net/methprog.pdf
    
    CEPHES MATHEMATICAL FUNCTION LIBRARY
    https://www.netlib.org/cephes/
  */

  static constexpr double sincof[] = {
    1.58962301576546568060E-10,
    -2.50507477628578072866E-8,
    2.75573136213857245213E-6,
    -1.98412698295895385996E-4,
    8.33333333332211858878E-3,
    -1.66666666666666307295E-1,
  };

  static constexpr double coscof[6] = {
    -1.13585365213876817300E-11,
    2.08757008419747316778E-9,
    -2.75573141792967388112E-7,
    2.48015872888517045348E-5,
    -1.38888888888730564116E-3,
    4.16666666666665929218E-2,
  };


  // highly accurate on [-pi/4, pi/4]
  template <int N>
  auto sincos_reduced (SIMD<double,N> x)
  {
    auto x2 = x*x;
  
    auto s = ((((( sincof[0]*x2 + sincof[1]) * x2 + sincof[2]) * x2 + sincof[3]) * x2 + sincof[4]) * x2 + sincof[5]);
    s = x + x*x*x * s;

    auto c = ((((( coscof[0]*x2 + coscof[1]) * x2 + coscof[2]) * x2 + coscof[3]) * x2 + coscof[4]) * x2 + coscof[5]);
    c = 1.0 - 0.5*x2 + x2*x2*c;

    return std::tuple{ s, c };
  }

  template <int N>
  auto sincos (SIMD<double,N> x)
  {
    // Cody-Waite: pi/2 split into a 33-bit part and correction (fdlibm)
    static constexpr double pio2_hi = 1.57079632673412561417E0;
    static constexpr double pio2_lo = 6.07710050650619224932E-11;

    auto y = round((2/M_PI) * x);
    auto q = lround(y);

    auto [s1,c1] = sincos_reduced(x - y * pio2_hi - y * pio2_lo);

    auto s2 = If((q & SIMD<int64_t,N>(1)) == SIMD<int64_t,N>(0), s1,  c1);
    auto s  = If((q & SIMD<int64_t,N>(2)) == SIMD<int64_t,N>(0), s2, -s2);
  
    auto c2 = If((q & SIMD<int64_t,N>(1)) == SIMD<int64_t,N>(0), c1, -s1);
    auto c  = If((q & SIMD<int64_t,N>(2)) == SIMD<int64_t,N>(0), c2, -c2);
  
    return std::tuple{ s, c };
  }

  template <int N>
  auto sincos_reduced (SIMD<float,N> x)
  {
    auto x2 = x*x;

    auto s = ((((( float(sincof[0])*x2 + float(sincof[1])) * x2 + float(sincof[2])) * x2 + float(sincof[3])) * x2 + float(sincof[4])) * x2 + float(sincof[5]));
    s = x + x*x*x * s;

    auto c = ((((( float(coscof[0])*x2 + float(coscof[1])) * x2 + float(coscof[2])) * x2 + float(coscof[3])) * x2 + float(coscof[4])) * x2 + float(coscof[5]));
    c = 1.0f - 0.5f*x2 + x2*x2*c;

    return std::tuple{ s, c };
  }

  template <int N>
  auto sincos (SIMD<float,N> x)
  {
    // Cody-Waite split of pi/2, from cephes sinf
    static constexpr float DP1 = 2*0.78515625f;
    static constexpr float DP2 = 2*2.4187564849853515625e-4f;
    static constexpr float DP3 = 2*3.77489497744594108e-8f;

    auto y = round(float(2/M_PI) * x);
    auto q = lround(y);

    auto [s1,c1] = sincos_reduced(((x - y*DP1) - y*DP2) - y*DP3);

    auto s2 = If((q & SIMD<int32_t,N>(1)) == SIMD<int32_t,N>(0), s1,  c1);
    auto s  = If((q & SIMD<int32_t,N>(2)) == SIMD<int32_t,N>(0), s2, -s2);

    auto c2 = If((q & SIMD<int32_t,N>(1)) == SIMD<int32_t,N>(0), c1, -s1);
    auto c  = If((q & SIMD<int32_t,N>(2)) == SIMD<int32_t,N>(0), c2, -c2);

    return std::tuple{ s, c };
  }






  
  template <int N>
  SIMD<double,N> exp_reduced (SIMD<double,N> x)
  {
    static constexpr double P[] = {
      1.26177193074810590878E-4,
      3.02994407707441961300E-2,
      9.99999999999999999910E-1,
    };
  
    static constexpr double Q[] = {
      3.00198505138664455042E-6,
      2.52448340349684104192E-3,
      2.27265548208155028766E-1,
      2.00000000000000000009E0,
    };
  
    /*
    // from:  https://www.netlib.org/cephes/
    rational approximation for exponential
    * of the fractional part:
    * e**x = 1 + 2x P(x**2)/( Q(x**2) - x P(x**2) )

    xx = x * x;
    px = x * polevl( xx, P, 2 );
    x =  px/( polevl( xx, Q, 3 ) - px );
    x = 1.0 + 2.0 * x;
    */

    auto xx = x*x;
    auto px = (P[0]*xx + P[1]) * xx + P[2];
    auto qx = ((Q[0]*xx+Q[1])*xx+Q[2])*xx+Q[3];
    return 1.0 + 2.0*x * px / (qx- x * px);
  }  


  template <int N>
  SIMD<double,N> pow2_int64_to_float64(SIMD<int64_t,N> n)
  {
    // thx to deepseek
    
    // Step 1: Clamp the input to valid exponent range [-1022, 1023]
    // (We use saturated operations to handle out-of-range values)
    SIMD<int64_t,N> max_exp(1023);
    SIMD<int64_t,N> min_exp(-1022);
    n = If(n > max_exp, max_exp, n);
    n = If(min_exp > n, min_exp, n);

    // Step 2: Add exponent bias (1023)
    n = n + SIMD<int64_t,N>(1023);

    // Step 3: Shift to exponent bit position (bit 52)
    auto shifted_exp = (n << IC<52>());
  
    // Step 4: Reinterpret as double
    return Reinterpret<double> (shifted_exp);
  }


  template <int N>
  SIMD<double,N> myexp (SIMD<double,N> x)
  {
    constexpr double log2 = 0.693147180559945286;  //  log(2.0);
    // Cody-Waite split of log(2), from cephes exp.c
    constexpr double C1 = 6.93145751953125E-1;
    constexpr double C2 = 1.42860682030941723212E-6;

    // keep the reduction arithmetic in range
    x = If(x > SIMD<double,N>(1000.0), SIMD<double,N>(1000.0), x);
    x = If(SIMD<double,N>(-1000.0) > x, SIMD<double,N>(-1000.0), x);

    auto r = round(1/log2 * x);

    // 2^r in two factors, to reach inf and denormals without clamping artifacts
    auto r1 = round(0.5*r);
    SIMD<double,N> pow2_1 = pow2_int64_to_float64 (lround(r1));
    SIMD<double,N> pow2_2 = pow2_int64_to_float64 (lround(r-r1));

    return exp_reduced(x - r*C1 - r*C2) * pow2_1 * pow2_2;

    // maybe better:
    // x = ldexp( x, n );
  }


  // *************************** float versions ***************************

  template <int N>
  SIMD<float,N> exp_reduced (SIMD<float,N> x)
  {
    // from cephes expf
    static constexpr float P[] = {
      1.9875691500E-4f,
      1.3981999507E-3f,
      8.3334519073E-3f,
      4.1665795894E-2f,
      1.6666665459E-1f,
      5.0000001201E-1f,
    };

    auto z = ((((P[0]*x + P[1])*x + P[2])*x + P[3])*x + P[4])*x + P[5];
    return 1.0f + x + x*x*z;
  }

  template <int N>
  SIMD<float,N> pow2_int32_to_float32 (SIMD<int32_t,N> n)
  {
    // biased exponent 0 gives 0.0, 255 gives inf
    SIMD<int32_t,N> max_exp(128);
    SIMD<int32_t,N> min_exp(-127);
    n = If(n > max_exp, max_exp, n);
    n = If(min_exp > n, min_exp, n);

    n = n + SIMD<int32_t,N>(127);
    auto shifted_exp = (n << IC<23>());
    return Reinterpret<float> (shifted_exp);
  }

  template <int N>
  SIMD<float,N> myexp (SIMD<float,N> x)
  {
    constexpr float log2e = 1.44269504088896341f;
    // Cody-Waite split of log(2), from cephes expf
    constexpr float C1 = 0.693359375f;
    constexpr float C2 = -2.12194440e-4f;

    // keep the reduction arithmetic in range (inf above 88.73, 0 below -103.98)
    x = Min(Max(x, SIMD<float,N>(-105.0f)), SIMD<float,N>(90.0f));

    auto r = round(log2e * x);

    // 2^r in two factors, to reach inf and denormals without clamping artifacts
    auto r1 = round(0.5f*r);
    SIMD<float,N> pow2_1 = pow2_int32_to_float32 (lround(r1));
    SIMD<float,N> pow2_2 = pow2_int32_to_float32 (lround(r-r1));

    return exp_reduced((x - r*C1) - r*C2) * pow2_1 * pow2_2;
  }

  /*
  inline auto Test1 (SIMD<double> x)
  {
    return myexp(x);
  }

  inline auto Test2 (SIMD<double> x)
  {
    return sincos(x);
  }

  inline auto Test3 (SIMD<double,4> x)
  {
    return myexp(x);
  }

  inline auto Test4 (SIMD<double,4> x)
  {
    return sincos(x);
  }
  */
  
}

#endif

#ifndef NETGEN_CORE_SIMD_MATH_HPP
#define NETGEN_CORE_SIMD_MATH_HPP

#include <tuple>
#include <limits>

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
  NETGEN_INLINE auto sincos_reduced (SIMD<double,N> x)
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
  NETGEN_INLINE auto sincos_reduced (SIMD<float,N> x)
  {
    auto x2 = x*x;

    auto s = ((((( float(sincof[0])*x2 + float(sincof[1])) * x2 + float(sincof[2])) * x2 + float(sincof[3])) * x2 + float(sincof[4])) * x2 + float(sincof[5]));
    s = x + x*x*x * s;

    auto c = ((((( float(coscof[0])*x2 + float(coscof[1])) * x2 + float(coscof[2])) * x2 + float(coscof[3])) * x2 + float(coscof[4])) * x2 + float(coscof[5]));
    c = 1.0f - 0.5f*x2 + x2*x2*c;

    return std::tuple{ s, c };
  }

  template <int N>
  NETGEN_INLINE auto sincos (SIMD<float,N> x)
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
  NETGEN_INLINE SIMD<double,N> exp_reduced (SIMD<double,N> x)
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
  NETGEN_INLINE SIMD<double,N> pow2_int64_to_float64(SIMD<int64_t,N> n)
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

  // atan: reduction to [0, tan(pi/8)] and rational approximation, Moshier Sec. 4.9
  template <int N>
  NETGEN_INLINE SIMD<double,N> atan2 (SIMD<double,N> y, SIMD<double,N> x)
  {
    constexpr double T3P8 = 2.41421356237309504880;   // tan(3pi/8)
    constexpr double TP8  = 0.41421356237309504880;   // tan(pi/8)
    constexpr double pio2_hi = 1.57079632679489655800e+00, pio2_lo = 6.12323399573676603587e-17;
    constexpr double pi_hi = 3.14159265358979311600e+00, pi_lo = 1.22464679914735317720e-16;
    SIMD<double,N> zero(0.0), one(1.0), inf(std::numeric_limits<double>::infinity());

    auto ysign = SIMD<int64_t,N>(0) > Reinterpret<int64_t>(y);
    auto xsign = SIMD<int64_t,N>(0) > Reinterpret<int64_t>(x);
    // inf/inf: use +-1/+-1, 0/0: use y/1 to keep the sign of zero
    auto bothinf = (fabs(x) == inf) && (fabs(y) == inf);
    auto bothzero = (x == zero) && (y == zero);
    auto num = If(bothinf, If(ysign, -one, one), y);
    auto den = If(bothinf, If(xsign, -one, one), If(bothzero, one, x));
    auto t = num/den;

    auto at = fabs(t);
    auto big = at > SIMD<double,N>(T3P8);
    auto mid = at > SIMD<double,N>(TP8);
    auto xr = If(big, -one/at, If(mid, (at-one)/(at+one), at));
    auto z = xr*xr;
    auto p = ((-0.8409808780644997716001*z - 8.83860837023772394279)*z - 21.8476213081316705724)*z - 14.8307050340438946993;
    auto q = (((z + 15.4974124675307267552)*z + 62.7906555762653017263)*z + 92.2381329856214406485)*z + 44.4921151021319438465;
    auto r = xr + xr*z*(p/q);
    r = If(big, SIMD<double,N>(pio2_hi), If(mid, SIMD<double,N>(0.5*pio2_hi), zero))
      + (r + If(big, SIMD<double,N>(pio2_lo), If(mid, SIMD<double,N>(0.5*pio2_lo), zero)));
    r = If(SIMD<int64_t,N>(0) > Reinterpret<int64_t>(t), -r, r);

    // x < 0 (including -0): shift by +-pi
    auto rx = If(ysign, (r - pi_hi) - pi_lo, (r + pi_hi) + pi_lo);
    return If(xsign, rx, r);
  }

  template <int N>
  NETGEN_INLINE SIMD<double,N> acos (SIMD<double,N> x)
  {
    x = If(x > SIMD<double,N>(1.0), SIMD<double,N>(1.0), x);
    x = If(SIMD<double,N>(-1.0) > x, SIMD<double,N>(-1.0), x);
    return atan2(sqrt((1.0-x)*(1.0+x)), x);
  }

  template <int N>
  NETGEN_INLINE SIMD<double,N> cbrt (SIMD<double,N> x)
  {
    double mantissa[N], scale[N];
    for (int i = 0; i < N; i++)
      {
        double a = x[i] < 0.0 ? -x[i] : x[i];
        int correction = 0;
        if (a < 0x1p-1022) { a *= 0x1p54; correction = -54; }
        uint64_t bits = BitCast<uint64_t>(a);
        int e = int((bits >> 52) & 2047) - 1023 + correction;
        int q = (e + 1074) / 3 - 358;
        mantissa[i] = BitCast<double>((bits & 0xfffffffffffffull)
                                     | (uint64_t(1023 + e - 3*q) << 52));
        scale[i] = BitCast<double>(uint64_t(1023 + q) << 52);
      }
    SIMD<double,N> a(mantissa), factor(scale);
    auto r = (((((SIMD<double,N>(6.067507661319302e-05)*a + -0.0016718798884351356)*a + 0.018694798810849338)*a + -0.11290585995314265)*a + 0.4880841295041025)*a + 0.6105363149353815);
    for (int i = 0; i < 3; i++)
      r = r + (a/(r*r)-r)/3.0;
    r = r*factor;
    r = If(SIMD<double,N>(0.0) > x, -r, r);
    r = If(fabs(x) == SIMD<double,N>(std::numeric_limits<double>::infinity()), x, r);
    r = If(x == x, r, x);
    return If(x == SIMD<double,N>(0.0), x, r);
  }

  namespace math
  {
    inline double sin (double x) { return std::get<0>(sincos(SIMD<double,1>(x)))[0]; }
    inline double cos (double x) { return std::get<1>(sincos(SIMD<double,1>(x)))[0]; }
    inline double atan2 (double y, double x) { return ngcore::atan2(SIMD<double,1>(y), SIMD<double,1>(x))[0]; }
    inline double acos (double x) { return ngcore::acos(SIMD<double,1>(x))[0]; }
    inline double cbrt (double x) { return ngcore::cbrt(SIMD<double,1>(x))[0]; }
  }
}

#endif

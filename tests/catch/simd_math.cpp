#include <catch2/catch.hpp>

#include <cmath>
#include <random>

#include <core/simd.hpp>

using namespace ngcore;

namespace
{

  // error of got relative to the spacing of doubles at ref
  double UlpError (double got, double ref)
  {
    if (std::isnan(got) || std::isnan(ref)) return 1e300;
    double u = std::nextafter(std::fabs(ref), 1e300) - std::fabs(ref);
    return std::fabs(got - ref) / u;
  }

  template <int N>
  void TestSinCosExp (double lo, double hi, double tol_ulp_exp,
                      double tol_abs_sincos, int samples = 20000)
  {
    CAPTURE(N, lo, hi);
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(lo, hi);

    double max_exp_ulp = 0, max_sin_abs = 0, max_cos_abs = 0;
    double worst_exp_x = 0, worst_sin_x = 0, worst_cos_x = 0;

    for (int it = 0; it < samples; it += N)
      {
        double xs[N];
        for (int i = 0; i < N; i++) xs[i] = dist(rng);
        SIMD<double,N> x(&xs[0]);

        auto [s,c] = sincos(x);
        auto e = myexp(x);

        for (int i = 0; i < N; i++)
          {
            double es = std::fabs(s[i] - std::sin(xs[i]));
            double ec = std::fabs(c[i] - std::cos(xs[i]));
            double ee = UlpError(e[i], std::exp(xs[i]));
            if (es > max_sin_abs) { max_sin_abs = es; worst_sin_x = xs[i]; }
            if (ec > max_cos_abs) { max_cos_abs = ec; worst_cos_x = xs[i]; }
            if (ee > max_exp_ulp) { max_exp_ulp = ee; worst_exp_x = xs[i]; }
          }
      }

    CAPTURE(worst_sin_x, worst_cos_x, worst_exp_x);
    CHECK(max_sin_abs <= tol_abs_sincos);
    CHECK(max_cos_abs <= tol_abs_sincos);
    CHECK(max_exp_ulp <= tol_ulp_exp);
  }

  template <int N>
  void TestLanesIndependent ()
  {
    // each lane must be computed from its own input
    double xs[N];
    for (int i = 0; i < N; i++) xs[i] = 0.3 + 1.7*i;
    SIMD<double,N> x(&xs[0]);

    auto [s,c] = sincos(x);
    auto e = myexp(x);
    for (int i = 0; i < N; i++)
      {
        CHECK(s[i] == Approx(std::sin(xs[i])).epsilon(1e-14));
        CHECK(c[i] == Approx(std::cos(xs[i])).epsilon(1e-14));
        CHECK(e[i] == Approx(std::exp(xs[i])).epsilon(1e-14));
      }
  }

  template <int N>
  void TestSpecialValues ()
  {
    auto [s0,c0] = sincos(SIMD<double,N>(0.0));
    for (int i = 0; i < N; i++)
      {
        CHECK(s0[i] == 0.0);
        CHECK(c0[i] == 1.0);
      }
    CHECK(myexp(SIMD<double,N>(0.0))[0] == 1.0);
    CHECK(myexp(SIMD<double,N>(1.0))[0] == Approx(std::exp(1.0)).epsilon(1e-15));

    // quadrant boundaries
    for (double x : { M_PI/2, M_PI, 3*M_PI/2, 2*M_PI, -M_PI/2, -M_PI })
      {
        auto [s,c] = sincos(SIMD<double,N>(x));
        CHECK(std::fabs(s[0] - std::sin(x)) < 1e-15);
        CHECK(std::fabs(c[0] - std::cos(x)) < 1e-15);
      }
  }

} // namespace


TEST_CASE("sincos/myexp SIMD<double> accuracy", "[simd_math]")
{
  // near-zero reduced range: everything ~1 ulp
  TestSinCosExp<1>(-0.785, 0.785, 3, 3e-16);
  TestSinCosExp<2>(-0.785, 0.785, 3, 3e-16);
  TestSinCosExp<4>(-0.785, 0.785, 3, 3e-16);
  TestSinCosExp<8>(-0.785, 0.785, 3, 3e-16);

  // moderate range
  TestSinCosExp<1>(-100, 100, 5, 5e-16);
  TestSinCosExp<2>(-100, 100, 5, 5e-16);
  TestSinCosExp<3>(-100, 100, 5, 5e-16);   // composed, non-power-of-two
  TestSinCosExp<4>(-100, 100, 5, 5e-16);
  TestSinCosExp<8>(-100, 100, 5, 5e-16);

  // full double exp range; myexp clamps instead of returning inf/0 outside +-708
  TestSinCosExp<2>(-700, 700, 5, 5e-16);
  TestSinCosExp<4>(-700, 700, 5, 5e-16);
}

TEST_CASE("sincos/myexp SIMD<double> lanes", "[simd_math]")
{
  TestLanesIndependent<1>();
  TestLanesIndependent<2>();
  TestLanesIndependent<4>();
  TestLanesIndependent<8>();
}

TEST_CASE("sincos/myexp SIMD<double> special values", "[simd_math]")
{
  TestSpecialValues<1>();
  TestSpecialValues<2>();
  TestSpecialValues<4>();
  TestSpecialValues<8>();

  // overflow / underflow of exp
  CHECK(std::isinf(myexp(SIMD<double,2>(710.0))[0]));
  CHECK(myexp(SIMD<double,2>(709.7))[0] == Approx(std::exp(709.7)).epsilon(1e-13));
  CHECK(myexp(SIMD<double,2>(-710.0))[0] == Approx(std::exp(-710.0)).epsilon(1e-12));
  CHECK(myexp(SIMD<double,2>(-746.0))[0] == 0.0);
}


////////////////////////////////////////////////////////////////////////////
// float

namespace
{

  double Ulp32Error (float got, double ref)
  {
    if (std::isnan(got) || std::isnan(ref)) return 1e300;
    float rf = float(ref);
    float u = std::nextafterf(std::fabs(rf), 3e38f) - std::fabs(rf);
    return std::fabs(double(got) - ref) / u;
  }

  template <int N>
  void TestSinCosExp32 (double lo, double hi, double tol_ulp, int samples = 20000)
  {
    CAPTURE(N, lo, hi);
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(lo, hi);

    double max_sin = 0, max_cos = 0, max_exp = 0;
    double worst_sin_x = 0, worst_cos_x = 0, worst_exp_x = 0;

    for (int it = 0; it < samples; it += N)
      {
        float xs[N];
        for (int i = 0; i < N; i++) xs[i] = float(dist(rng));
        SIMD<float,N> x(&xs[0]);

        auto [s,c] = sincos(x);
        auto e = myexp(x);

        for (int i = 0; i < N; i++)
          {
            double es = Ulp32Error(s[i], std::sin(double(xs[i])));
            double ec = Ulp32Error(c[i], std::cos(double(xs[i])));
            double ee = Ulp32Error(e[i], std::exp(double(xs[i])));
            if (es > max_sin) { max_sin = es; worst_sin_x = xs[i]; }
            if (ec > max_cos) { max_cos = ec; worst_cos_x = xs[i]; }
            if (ee > max_exp) { max_exp = ee; worst_exp_x = xs[i]; }
          }
      }

    CAPTURE(worst_sin_x, worst_cos_x, worst_exp_x);
    CHECK(max_sin <= tol_ulp);
    CHECK(max_cos <= tol_ulp);
    CHECK(max_exp <= tol_ulp);
  }

} // namespace


TEST_CASE("sincos/myexp SIMD<float> accuracy", "[simd_math]")
{
  TestSinCosExp32<1>(-0.785, 0.785, 2);
  TestSinCosExp32<2>(-0.785, 0.785, 2);
  TestSinCosExp32<4>(-0.785, 0.785, 2);
  TestSinCosExp32<8>(-0.785, 0.785, 2);

  TestSinCosExp32<1>(-87, 87, 3);
  TestSinCosExp32<2>(-87, 87, 3);
  TestSinCosExp32<3>(-87, 87, 3);   // composed, non-power-of-two
  TestSinCosExp32<4>(-87, 87, 3);
  TestSinCosExp32<8>(-87, 87, 3);
  TestSinCosExp32<16>(-87, 87, 3);

  // exp edges: still finite just below overflow, denormal results below -87
  TestSinCosExp32<4>(88.3, 88.72, 3);
  TestSinCosExp32<4>(-103, -87, 3);
}

TEST_CASE("sincos/myexp SIMD<float> lanes and special values", "[simd_math]")
{
  constexpr int N = 4;
  float xs[N];
  for (int i = 0; i < N; i++) xs[i] = 0.3f + 1.7f*i;
  SIMD<float,N> x(&xs[0]);
  auto [s,c] = sincos(x);
  auto e = myexp(x);
  for (int i = 0; i < N; i++)
    {
      CHECK(s[i] == Approx(std::sin(double(xs[i]))).epsilon(1e-6));
      CHECK(c[i] == Approx(std::cos(double(xs[i]))).epsilon(1e-6));
      CHECK(e[i] == Approx(std::exp(double(xs[i]))).epsilon(1e-6));
    }

  auto [s0,c0] = sincos(SIMD<float,N>(0.0f));
  CHECK(s0[0] == 0.0f);
  CHECK(c0[0] == 1.0f);
  CHECK(myexp(SIMD<float,N>(0.0f))[0] == 1.0f);

  for (double xq : { M_PI/2, M_PI, 3*M_PI/2, 2*M_PI, -M_PI/2, -M_PI })
    {
      auto [sq,cq] = sincos(SIMD<float,N>(float(xq)));
      CHECK(std::fabs(double(sq[0]) - std::sin(double(float(xq)))) < 1e-7);
      CHECK(std::fabs(double(cq[0]) - std::cos(double(float(xq)))) < 1e-7);
    }

  // overflow / underflow of exp
  CHECK(std::isinf(myexp(SIMD<float,N>(89.0f))[0]));
  CHECK(myexp(SIMD<float,N>(88.5f))[0] == Approx(std::exp(88.5)).epsilon(1e-6));
  CHECK(myexp(SIMD<float,N>(-100.0f))[0] == Approx(std::exp(-100.0)).epsilon(0.1)); // denormal, ~4% spacing
  CHECK(myexp(SIMD<float,N>(-110.0f))[0] == 0.0f);
}

namespace
{
  uint64_t ResultBits (double x) { return BitCast<uint64_t>(x); }

  struct MathRandom
  {
    uint64_t state = 1;
    uint64_t Bits ()
    {
      uint64_t z = (state += 0x9e3779b97f4a7c15ull);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
      return z ^ (z >> 31);
    }
    double Unit () { return double(Bits() >> 11)*0x1p-53; }
  };

  template <int N>
  void TestDeterministicMath ()
  {
    CAPTURE(N);
    MathRandom rng;
    std::array<uint64_t,6> fingerprints;
    fingerprints.fill(0xcbf29ce484222325ull);
    double errors[6] = {};
    for (int j = 0; j < 99984; j += N)
      {
        double xs[N], ys[N], positive[N];
        for (int i = 0; i < N; i++)
          {
            xs[i] = 2*rng.Unit()-1;
            ys[i] = 2*rng.Unit()-1;
            uint64_t mantissa = rng.Bits() & 0xfffffffffffffull;
            uint64_t exponent = 1+rng.Bits()%2046;
            positive[i] = BitCast<double>(mantissa | (exponent << 52));
          }
        SIMD<double,N> x(xs), y(ys), p(positive);
        auto a = ngcore::atan2(y,x);
        auto c = ngcore::acos(x);
        auto r = ngcore::cbrt(p);
        auto [s,co] = sincos(100*x);
        auto e = myexp(100*x);
        for (int i = 0; i < N; i++)
          {
            double got[] = {a[i], c[i], r[i], s[i], co[i], e[i]};
            double ref[] = {std::atan2(ys[i],xs[i]), std::acos(xs[i]), std::cbrt(positive[i]),
                            std::sin(100*xs[i]), std::cos(100*xs[i]), std::exp(100*xs[i])};
            for (int k = 0; k < 6; k++)
              {
                errors[k] = std::max(errors[k], UlpError(got[k], ref[k]));
                uint64_t bits = ResultBits(got[k]);
                for (int b = 0; b < 8; b++)
                  {
                    fingerprints[k] ^= (bits >> (8*b)) & 255;
                    fingerprints[k] *= 0x100000001b3ull;
                  }
              }
            if (j < 96)
              {
                CHECK(ResultBits(a[i]) == ResultBits(math::atan2(ys[i],xs[i])));
                CHECK(ResultBits(c[i]) == ResultBits(math::acos(xs[i])));
                CHECK(ResultBits(r[i]) == ResultBits(math::cbrt(positive[i])));
                CHECK(ResultBits(s[i]) == ResultBits(math::sin(100*xs[i])));
                CHECK(ResultBits(co[i]) == ResultBits(math::cos(100*xs[i])));
                CHECK(ResultBits(e[i]) == ResultBits(myexp(SIMD<double,1>(100*xs[i]))[0]));
              }
          }
      }
    const double tolerances[] = {3, 3, 4, 4, 4, 5};
    const uint64_t expected[] = {0x2512244d33384504ull, 0xa2e297791dd116d3ull,
                                 0x68510e590e11a5efull, 0x53d3fa72e12dec13ull,
                                 0x1c84ebe136f7855dull, 0xa2cbfe21da809e0dull};
    for (int k = 0; k < 6; k++)
      {
        CAPTURE(k, errors[k], fingerprints[k]);
        CHECK(errors[k] <= tolerances[k]);
        CHECK(fingerprints[k] == expected[k]);
      }
  }

  template <int N>
  void TestMathEdges ()
  {
    CAPTURE(N);
    constexpr double inf = std::numeric_limits<double>::infinity();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    const double values[] = {0.0, -0.0, 0x1p-1074, -0x1p-1074, 0x1p-1022, -0x1p-1022,
                             1.0, -1.0, 0x1.fffffffffffffp1023, -0x1.fffffffffffffp1023,
                             inf, -inf};
    for (double x : values)
      {
        auto r = ngcore::cbrt(SIMD<double,N>(x));
        for (int i = 0; i < N; i++)
          if (x == 0.0 || std::isinf(x)) CHECK(ResultBits(r[i]) == ResultBits(x));
          else CHECK(UlpError(r[i], std::cbrt(x)) <= 4);
        for (double y : values)
          {
            auto a = ngcore::atan2(SIMD<double,N>(y), SIMD<double,N>(x));
            double ref = std::atan2(y,x);
            for (int i = 0; i < N; i++)
              if (ref == 0.0) CHECK(ResultBits(a[i]) == ResultBits(ref));
              else CHECK(UlpError(a[i], ref) <= 3);
          }
      }
    for (double x : {-1.0, std::nextafter(-1.0,0.0), -0.5, 0.0, 0.5, std::nextafter(1.0,0.0), 1.0})
      CHECK(UlpError(ngcore::acos(SIMD<double,N>(x))[0], std::acos(x)) <= 3);
    CHECK(ngcore::acos(SIMD<double,N>(1.01))[0] == 0.0);
    CHECK(ngcore::acos(SIMD<double,N>(-1.01))[0] == M_PI);
    CHECK(std::isnan(ngcore::acos(SIMD<double,N>(nan))[0]));
    CHECK(std::isnan(ngcore::atan2(SIMD<double,N>(nan),SIMD<double,N>(1.0))[0]));
    CHECK(std::isnan(ngcore::atan2(SIMD<double,N>(1.0),SIMD<double,N>(nan))[0]));
    CHECK(std::isnan(ngcore::cbrt(SIMD<double,N>(nan))[0]));
    CHECK(std::isnan(myexp(SIMD<double,N>(nan))[0]));
    CHECK(myexp(SIMD<double,N>(inf))[0] == inf);
    CHECK(myexp(SIMD<double,N>(-inf))[0] == 0.0);
    CHECK(std::isnan(std::get<0>(sincos(SIMD<double,N>(inf)))[0]));
  }

  // round/lround at ties (x.5) must agree between SIMD widths: libm round is
  // half away from zero, the x86/arm vector instructions are half to even.
  // Otherwise the argument reduction of sincos/exp differs between widths.
  template <int N>
  void TestRoundTies ()
  {
    CAPTURE(N);
    constexpr double inf = std::numeric_limits<double>::infinity();
    for (int j = -100; j <= 100; j++)
      {
        double ties[N];
        for (int i = 0; i < N; i++) ties[i] = double(j+i)+0.5;
        SIMD<double,N> t(ties);
        SIMD<float,N> tf = SIMD<float,N>([&] (int i) { return float(ties[i]); });
        for (int i = 0; i < N; i++)
          {
            CHECK(ngcore::round(t)[i] == ngcore::round(SIMD<double,1>(ties[i]))[0]);
            CHECK(ngcore::lround(t)[i] == ngcore::lround(SIMD<double,1>(ties[i]))[0]);
            CHECK(ngcore::round(tf)[i] == ngcore::round(SIMD<float,1>(float(ties[i])))[0]);
            CHECK(ngcore::lround(tf)[i] == ngcore::lround(SIMD<float,1>(float(ties[i])))[0]);
          }
        // not a tie: must round to nearest, not truncate
        CHECK(ngcore::lround(SIMD<double,N>(ties[0]+0.25))[0] == int64_t(j+1));
        CHECK(ngcore::lround(SIMD<double,N>(ties[0]-0.25))[0] == int64_t(j));

        for (double x : {ties[0]*(M_PI/2), ties[0]*0.693147180559945286})
          for (double near : {std::nextafter(x,-inf), x, std::nextafter(x,inf)})
            {
              auto [s,c] = sincos(SIMD<double,N>(near));
              auto e = myexp(SIMD<double,N>(near));
              auto [s1,c1] = sincos(SIMD<double,1>(near));
              auto e1 = myexp(SIMD<double,1>(near));
              for (int i = 0; i < N; i++)
                {
                  CHECK(ResultBits(s[i]) == ResultBits(s1[0]));
                  CHECK(ResultBits(c[i]) == ResultBits(c1[0]));
                  CHECK(ResultBits(e[i]) == ResultBits(e1[0]));
                }
            }
      }
  }
}

TEST_CASE("deterministic math accuracy and fingerprints", "[simd_math]")
{
  TestDeterministicMath<1>();
  TestDeterministicMath<2>();
  TestDeterministicMath<3>();
  TestDeterministicMath<4>();
  TestDeterministicMath<8>();
}

TEST_CASE("deterministic math edges", "[simd_math]")
{
  TestMathEdges<1>();
  TestMathEdges<2>();
  TestMathEdges<3>();
  TestMathEdges<4>();
  TestMathEdges<8>();
}

/*
TEST_CASE("round ties to even, consistent across SIMD widths", "[simd_math]")
{
  TestRoundTies<1>();
  TestRoundTies<2>();
  TestRoundTies<3>();
  TestRoundTies<4>();
  TestRoundTies<8>();
}
*/

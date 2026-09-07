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

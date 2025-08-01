#ifndef NETGEN_CORE_SIMD_MATH_HPP
#define NETGEN_CORE_SIMD_MATH_HPP

#include <tuple>

namespace ngcore
{

  /*
    base on:
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
    auto y = Round(x / (M_PI/2));
    auto q = RoundI(y);
  
    auto [s1,c1] = sincos_reduced(x - y * (M_PI/2));

    auto s2 = If((q & SIMD<int64_t,N>(1)) == SIMD<int64_t,N>(0), s1,  c1);
    auto s  = If((q & SIMD<int64_t,N>(2)) == SIMD<int64_t,N>(0), s2, -s2);
  
    auto c2 = If((q & SIMD<int64_t,N>(1)) == SIMD<int64_t,N>(0), c1, -s1);
    auto c  = If((q & SIMD<int64_t,N>(2)) == SIMD<int64_t,N>(0), c2, -c2);
  
    return std::tuple{ s, c };
  }

}

#endif

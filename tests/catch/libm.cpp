#include <catch2/catch.hpp>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

// Last-bit differences of libm functions between platforms (glibc versions,
// UCRT, Apple libm) change mesh point coordinates, e.g. atan2 in
// Cylinder::ToPlane. The fingerprints below are from the machine generating
// tests/pytest/results.json.

namespace
{
  struct SplitMix64
  {
    uint64_t state;
    uint64_t operator() ()
    {
      uint64_t z = (state += 0x9e3779b97f4a7c15ull);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
      return z ^ (z >> 31);
    }
    // uniform in [lo, lo + 2^e), exact without rounding
    double operator() (double lo, int e)
    {
      return lo + std::ldexp(double((*this)() >> 11), e - 53);
    }
  };

  struct Fingerprint
  {
    uint64_t h = 0xcbf29ce484222325ull;
    void Add (double x)
    {
      uint64_t bits;
      std::memcpy(&bits, &x, sizeof(bits));
      for (int i = 0; i < 8; i++)
        {
          h ^= (bits >> (8*i)) & 0xff;
          h *= 0x100000001b3ull;
        }
    }
    std::string Hex () const
    {
      std::stringstream s;
      s << std::hex << std::setw(16) << std::setfill('0') << h;
      return s.str();
    }
  };

  const int samples = 100000;

  std::string Fingerprint1 (std::function<double(double)> f, double lo, int e)
  {
    SplitMix64 rng{1};
    Fingerprint fp;
    for (int i = 0; i < samples; i++)
      fp.Add(f(rng(lo, e)));
    return fp.Hex();
  }

  std::string Fingerprint2 (std::function<double(double,double)> f,
                            double lo0, int e0, double lo1, int e1)
  {
    SplitMix64 rng{2};
    Fingerprint fp;
    for (int i = 0; i < samples; i++)
      {
        double x = rng(lo0, e0);
        double y = rng(lo1, e1);
        fp.Add(f(x, y));
      }
    return fp.Hex();
  }
} // namespace


TEST_CASE("libm bitwise reproducibility", "[libm][!mayfail]")
{
  volatile double half = 0.5;
  volatile double third = 1./3.;

  // argument ranges as recorded while meshing twocyl.geo and cylsphere.geo
  struct { const char * name; std::string fp; const char * ref; } cases[] =
    {
      { "sqrt",  Fingerprint1([](double x) { return std::sqrt(x); }, 0, 4), "61bde8d69fe962d4" },
      { "sin",   Fingerprint1([](double x) { return std::sin(x); }, 0, 1), "cd307cd50b841459" },
      { "cos",   Fingerprint1([](double x) { return std::cos(x); }, -4, 3), "570ae882b46fdccd" },
      { "acos",  Fingerprint1([](double x) { return std::acos(x); }, -1, 1), "b14e812a4398ac1e" },
      { "exp",   Fingerprint1([](double x) { return std::exp(x); }, -16, 5), "9b32ed99535a9851" },
      { "log",   Fingerprint1([](double x) { return std::log(x); }, 0, 4), "97de66c9efd09fbd" },
      { "pow(x,1/2)", Fingerprint1([&](double x) { return std::pow(x, double(half)); }, 1, 3), "3e1b2c6b38b4685f" },
      { "pow(x,1/3)", Fingerprint1([&](double x) { return std::pow(x, double(third)); }, 0, 13), "c22603d4442fe8ac" },
      { "pow(x,y)",   Fingerprint2([](double x, double y) { return std::pow(x, y); }, 0, 4, -4, 3), "d719ce35d116156e" },
      { "atan2", Fingerprint2([](double y, double x) { return std::atan2(y, x); }, -1, 1, -1, 1), "b8f9d7b8af967373" },
      { "atan2 small ratio", Fingerprint2([](double y, double x) { return std::atan2(y, x); }, -0.0625, -3, 0.25, -2), "37ae2be52ebc55cf" },
    };

  for (auto & c : cases)
    {
      WARN(c.name << ": " << c.fp);
      CAPTURE(c.name);
      CHECK(c.fp == c.ref);
    }
}

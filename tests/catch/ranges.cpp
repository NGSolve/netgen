
#ifdef CATCH2_v3
#include "catch2/catch_all.hpp"
#else
#include <catch2/catch.hpp>
#endif

#include <core/array.hpp>
#include <core/ranges.hpp>

using namespace ngcore;

TEST_CASE("ranges")
{
  Array<int> a { 3, -1, 10, -5 };
  Array<int> positive { 3, 10 };
  int i = 0;
  for(auto pos_val : a | filter([](auto val) { return val >= 0; }))
    {
      CHECK(pos_val == positive[i++]);
    }
}


#include "catch.hpp"
#include <array.hpp>
using namespace ngcore;
using namespace std;


TEST_CASE("Array")
{
  Array<int> array;
#ifdef DEBUG
  CHECK_THROWS_AS(array[0], RangeException);
  CHECK_THROWS_AS(array.DeleteLast(), RangeException);
  CHECK_THROWS_AS(array.Last(), RangeException);
#endif // DEBUG
  Array<double> a_initlst = { 1., 2., 3.};
  CHECK(a_initlst[1] == 2.);
  CHECK(a_initlst.size() == 3);
  FlatArray fa_a = a_initlst;
  CHECK(typeid(fa_a) == typeid(FlatArray<double>));
  CHECK(fa_a.size() == 3);
  CHECK(a.Last() == 3.);
  a.DeleteLast();
  CHECK(a.Last() == 2. && a.Size() == 2);
#ifdef DEBUG
  CHECK_THROWS_AS(fa_a[5], RangeException);
#endif // DEBUG
  Array b = Array<int>(4);
  b = 2;
  int count = 0;
  for(auto val : b)
  {
    count++;
    CHECK(val == 2);
  }
  CHECK(count == 4);
}

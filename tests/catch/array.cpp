
#include "catch.hpp"
#include <core/array.hpp>
using namespace ngcore;
using namespace std;

#include "meshing.hpp"


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
  CHECK(a_initlst.Size() == 3);
  FlatArray fa_a = a_initlst;
  CHECK(typeid(fa_a) == typeid(FlatArray<double>));
  CHECK(fa_a.Size() == 3);
  CHECK(fa_a.Last() == 3.);
  a_initlst.DeleteLast();
  CHECK(a_initlst.Last() == 2.);
  CHECK(a_initlst.Size() == 2);
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

  // range tests
  CHECK(typeid(array.Range()) == typeid(T_Range<size_t>));
  Array<int, int> intarray;
  CHECK(typeid(intarray.Range()) == typeid(T_Range<int>));
  CHECK(typeid(Range(intarray)) == typeid(T_Range<int>));
  int i = 0;
  for(auto j : Range(b))
    CHECK(j == i++);
  i = 0;
  for(auto j : b.Range())
    CHECK(j == i++);

  // pointindex is still 1 based
  Array<double, netgen::PointIndex> piarray(2);
  i = 1;
  for(auto j : Range(piarray))
    CHECK(j == i++);
  i = 1;
  for(auto j : piarray.Range())
    CHECK(j == i++);
}

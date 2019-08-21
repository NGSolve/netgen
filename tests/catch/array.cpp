
#include "catch.hpp"
#include <core/array.hpp>
using namespace ngcore;
using namespace std;

#include "meshing.hpp"

template<typename TIND>
class ClsWithIndexType
{
  size_t size;
public:
  ClsWithIndexType(size_t asize) : size(asize) {}
  using index_type = TIND;
  size_t Size() const { return size; }
};

template<typename TIND>
class ClsWithRange : public ClsWithIndexType<TIND>
{
public:
  ClsWithRange(size_t size) : ClsWithIndexType<TIND>(size) {}
  T_Range<size_t> Range() const { return {1, 1+this->Size()}; }
};


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
  // a class can implement index_type and Size as well.
  ClsWithIndexType<int> clsi(3);
  CHECK(typeid(Range(clsi)) == typeid(T_Range<int>));
  i = 0;
  for(auto j : Range(clsi))
    CHECK(j == i++);
  // if the class has a Range function prefer that one
  ClsWithRange<int> clsr(3);
  CHECK(typeid(Range(clsr)) == typeid(T_Range<size_t>));
  i=1;
  for(auto j : Range(clsr))
    CHECK(j == i++);
  CHECK(typeid(Range(size_t(4))) == typeid(T_Range<size_t>));
  CHECK(typeid(Range(4)) == typeid(T_Range<int>));
}

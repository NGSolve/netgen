
#include <catch2/catch.hpp>
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
  netgen::PointIndex pi = IndexBASE<netgen::PointIndex>();
  for(auto j : Range(piarray))
    CHECK(j == pi++);
  pi = IndexBASE<netgen::PointIndex>();
  for(auto j : piarray.Range())
    CHECK(j == pi++);
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

TEST_CASE("Array constructors with index type")
{
  // these used to initialize FlatArray<T> instead of FlatArray<T,IndexType>,
  // so they did not compile for an array with a non-default index type
  Array<double, netgen::PointIndex> a { 1.0, 2.0, 3.0 };
  CHECK(a.Size() == 3);
  auto pi = IndexBASE<netgen::PointIndex>();
  CHECK(a[pi] == 1.0);
  CHECK(a[pi+2] == 3.0);
  for (auto i : a.Range())
    CHECK(a[i] == double(i-IndexBASE<netgen::PointIndex>()+1));

  Array<double> b { 1.0, 2.0 }, c { 3.0 };
  Array<double, netgen::PointIndex> m (b, c);   // merge-copy
  CHECK(m.Size() == 3);
  CHECK(m[pi] == 1.0);
  CHECK(m[pi+2] == 3.0);

  LocalHeap lh(10000, "test");
  Array<double, netgen::PointIndex> h (5, lh);
  CHECK(h.Size() == 5);
}

TEST_CASE("Array brace-init from another array")
{
  // IVec has a converting ctor from any array-like. Unconstrained, it made
  // Array<IVec<..>,..> x { other_array } deduce a one-element initializer_list
  // instead of copying. Constrained on value_type, the copy ctor wins again.
  using PI = netgen::PointIndex;
  Array<IVec<2,PI>, PI> src(3);
  auto b = IndexBASE<PI>();
  for (int k = 0; k < 3; k++)
    src[b+k] = IVec<2,PI>(b+k, b+k+1);

  Array<IVec<2,PI>, PI> dst { src };
  CHECK(dst.Size() == 3);          // was 1 with the unconstrained ctor
  CHECK(dst[b][0] == b);
  CHECK(dst[b+2][1] == b+3);
}

TEST_CASE("IVec from array-like")
{
  // ngsolve builds IVec<3> from an AOWrapper, which has no value_type -
  // the ctor constraint must key on indexing, not on that typedef
  Array<int> verts { 10, 20, 30 };
  IVec<3> fromwrapper (ArrayObject(verts));
  CHECK(fromwrapper[0] == 10);
  CHECK(fromwrapper[2] == 30);

  IVec<3> fromarray (verts);
  CHECK(fromarray[1] == 20);

  IVec<2> shorter (verts);         // reading fewer than available is fine
  CHECK(shorter[1] == 20);
}

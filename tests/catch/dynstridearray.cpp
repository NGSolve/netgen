#include <catch2/catch.hpp>
#include <core/ngcore.hpp>
using namespace ngcore;

struct Hdr { int typ = 7; int8_t np = 0; float badness = 1.5f; int index = -1; };
struct GI { int trignum = -2; double u = 0.5, v = 0.25; };
struct P { int i = 42; };   // like PointIndex

// base-1 index like netgen Index
struct Idx {
  int i = 1;
  Idx () = default;
  constexpr Idx (int ai) : i(ai) { }
  static constexpr Idx Base () { return Idx(1); }
  Idx operator+ (size_t n) const { return Idx(i+int(n)); }
  Idx operator- (Idx b) const { return Idx(i-b.i); }
  Idx & operator++ () { i++; return *this; }
  operator size_t () const { return size_t(i); }   // for stride arithmetic
  bool operator!= (Idx b) const { return i != b.i; }
  bool operator< (Idx b) const { return i < b.i; }
};

TEST_CASE("DynStrideArray")
{
  // layout: header 16 B, GI 8-aligned, P 4-aligned
  using L = DynStrideLayout<Hdr, GI, P>;
  L l3 (3);
  CHECK (l3.offset[0] == 16);
  CHECK (l3.offset[1] == 16 + 3*24);
  CHECK (l3.stride == 104);      // 16 + 3*24 + 3*4 = 100, rounded to 8
  CHECK (l3.stride % 8 == 0);
  DynStrideLayout<Hdr, GI, P> l0;
  CHECK (l0.width == 0);
  CHECK (l0.stride == 16);
  DynStrideLayout<void, int> lv (6);
  CHECK (lv.offset[0] == 0);
  CHECK (lv.stride == 24);

  using A2 = DynStrideArray<Hdr, TailList<GI, P>>;
  A2 a (3);
  CHECK (a.Size() == 0);
  CHECK (a.Width() == 3);
  auto i0 = a.Append();
  CHECK (i0 == 0);
  CHECK (a.Size() == 1);
  auto v = a[0];
  CHECK (v.Head().typ == 7);   // default-initialized header
  CHECK (v.Head().badness == 1.5f);
  CHECK (v.TailArray<0>()[2].trignum == -2);
  CHECK (v.TailArray<1>()[1].i == 42);
  v.Head().np = 3;
  for (int k = 0; k < 3; k++) { v.TailArray<1>()[k].i = 10+k; v.TailArray<0>()[k].u = k; }

  // growth keeps contents
  for (int e = 1; e < 100; e++) { auto i = a.Append(); a[i].Head().index = e; a[i].TailArray<1>(3)[0].i = 1000+e; }
  CHECK (a.Size() == 100);
  CHECK (a[0].TailArray<1>()[2].i == 12);
  CHECK (a[0].TailArray<0>()[2].u == 2);
  CHECK (a[57].Head().index == 57);
  CHECK (a[57].TailArray<1>()[0].i == 1057);
  CHECK (a[57].TailArray<1>()[1].i == 42);   // untouched entry default

  // widen: old entries kept, new ones default
  a.SetWidth (8);
  CHECK (a.Width() == 8);
  CHECK (a.Stride() == ((16 + 8*24 + 8*4 + 7)/8)*8);
  CHECK (a[0].Head().np == 3);
  CHECK (a[0].TailArray<1>()[2].i == 12);
  CHECK (a[0].TailArray<0>()[2].u == 2);
  CHECK (a[0].TailArray<1>()[7].i == 42);
  CHECK (a[0].TailArray<0>()[5].trignum == -2);
  CHECK (a[99].Head().index == 99);
  CHECK (a[99].TailArray<1>()[0].i == 1099);
  a[99].TailArray<1>()[7].i = 777;

  // narrow: truncates
  a.SetWidth (4);
  CHECK (a[0].TailArray<1>()[2].i == 12);
  CHECK (a[99].Head().index == 99);
  CHECK (a[99].TailArray<1>().Size() == 4);

  // append a view from another array with smaller width
  A2 b (2);
  auto ib = b.Append();
  b[ib].Head().index = 5555; b[ib].TailArray<1>()[1].i = -9;
  auto ia = a.Append (b[ib]);
  CHECK (a.Width() == 4);
  CHECK (a[ia].Head().index == 5555);
  CHECK (a[ia].TailArray<1>()[1].i == -9);
  CHECK (a[ia].TailArray<1>()[3].i == 42);
  // append a wider view grows the array
  auto ic = b.Append (a[0]);
  CHECK (b.Width() == 4);
  CHECK (b[ic].TailArray<1>()[2].i == 12);
  CHECK (b[ib].Head().index == 5555);
  CHECK (b[ib].TailArray<1>()[3].i == 42);

  // remove
  size_t n = a.Size();
  a.RemoveElement (0);
  CHECK (a.Size() == n-1);
  CHECK (a[0].Head().index == 1);
  a.RemoveElementIf ([] (auto el) { return el.Head().index % 2 == 0; });
  for (auto el : a) CHECK (el.Head().index % 2 == 1);
  CHECK (a.Last().Head().index == 5555);

  // const access & iteration
  const A2 & ca = a;
  int cnt = 0;
  for (auto el : ca) { cnt++; (void) el.TailArray<0>(); }
  CHECK (cnt == int(a.Size()));
  A2::ConstView cv = a[0];
  CHECK (cv.Head().index == 1);

  // copy / move
  A2 c = a;
  CHECK (c.Size() == a.Size());
  CHECK (c[0].TailArray<1>()[0].i == a[0].TailArray<1>()[0].i);
  c[0].TailArray<1>()[0].i = -1;
  CHECK (a[0].TailArray<1>()[0].i != -1);
  A2 d = std::move (c);
  CHECK (c.Size() == 0);
  CHECK (d.Size() == a.Size());

  // header-less table with base-1 index
  DynStrideArray<void, TailList<int>, Idx> t (6);
  t.SetSize (10);
  for (auto i : t.Range()) for (int k = 0; k < 6; k++) t[i].TailArray()[k] = i.i*10+k;
  CHECK (t[Idx(1)].TailArray()[0] == 10);
  CHECK (t[Idx(10)].TailArray()[5] == 105);
  CHECK (t.Range().First() == 1);
  CHECK (t.Range().Next() == 11);
  auto ti = t.Append();
  CHECK (ti.i == 11);
  CHECK (t[ti].TailArray()[3] == 0);
  t.RemoveElement (Idx(1));
  CHECK (t[Idx(1)].TailArray()[0] == 20);

  // default-constructed array is usable
  A2 z;
  z.SetSize (3);
  CHECK (z.Stride() == 16);
  CHECK (z[2].Head().typ == 7);
  z.Append (a[0]);
  CHECK (z.Width() == a.Width());

  // sized ctor
  A2 e (5, 2);
  CHECK (e.Size() == 5);
  CHECK (e[4].Head().typ == 7);
  CHECK (e[4].TailArray<1>()[1].i == 42);

}

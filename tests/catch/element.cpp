#include <catch2/catch.hpp>
#include "meshing.hpp"
using namespace ngcore;
using namespace netgen;

TEST_CASE("ElementRef")
{
  T_VOLELEMENTS els(4);
  CHECK (els.Width() == 4);
  CHECK (els.Stride() == sizeof(ElementHeader) + 4*sizeof(PointIndex));

  Element tet(TET);
  for (int i = 0; i < 4; i++) tet[i] = PointIndex::FromNr0(10+i);
  tet.SetIndex(3);
  CHECK (tet.MaxNP() == ELEMENT_MAXPOINTS);

  // a value copy is independent, a handle copy aliases
  Element tet2 = tet;
  tet2[0] = PointIndex::FromNr0(77);
  CHECK (tet[0] == PointIndex::FromNr0(10));
  ElementRef href = tet;
  href[0] = PointIndex::FromNr0(78);
  CHECK (tet[0] == PointIndex::FromNr0(78));
  tet[0] = PointIndex::FromNr0(10);

  // store in the strided array
  auto ei = els.Append (tet);
  ElementRef v = els[ei];
  CHECK (v.GetNP() == 4);
  CHECK (v.GetType() == TET);
  CHECK (v.GetIndex().Nr1() == 3);
  CHECK (v[2] == PointIndex::FromNr0(12));
  CHECK (v.PNum(1) == PointIndex::FromNr0(10));
  CHECK (v.PNums().Size() == 4);
  CHECK (v.Vertices().Size() == 4);
  CHECK (v.MaxNP() == 4);

  // mutate through the handle, read back as value
  v[0] = PointIndex::FromNr0(99);
  v.SetIndex(7);
  v.SetRefinementFlag(false);
  v.SetOrder(2,3,4);
  Element back (els[ei]);
  CHECK (back[0] == PointIndex::FromNr0(99));
  CHECK (back.GetIndex().Nr1() == 7);
  CHECK (!back.TestRefinementFlag());
  int ox, oy, oz;
  back.GetOrder(ox, oy, oz);
  CHECK (ox == 2); CHECK (oy == 3); CHECK (oz == 4);
  CHECK (Copy(els[ei]).GetIndex().Nr1() == 7);

  const T_VOLELEMENTS & cels = els;
  const ElementRef cv = cels[ei];
  CHECK (cv.GetNP() == 4);
  CHECK (cv[0] == PointIndex::FromNr0(99));

  // appending a wider element re-strides, keeps the old one, invalidates spare slots
  Element prism(PRISM);
  for (int i = 0; i < 6; i++) prism[i] = PointIndex::FromNr0(20+i);
  auto ei2 = els.Append (prism);
  CHECK (els.Width() == 6);
  CHECK (els.Size() == 2);
  ElementRef v2 = els[ei2];
  CHECK (v2.GetType() == PRISM);
  CHECK (v2[5] == PointIndex::FromNr0(25));
  CHECK (!els[ei][4].IsValid());
  CHECK (v2.GetNV() == 6);
  CHECK (v2.GetNFaces() == 5);
  ElementRef v1 = els[ei];
  CHECK (v1[0] == PointIndex::FromNr0(99));
  CHECK (v1.GetIndex().Nr1() == 7);
  CHECK (v1.GetNP() == 4);

  // handle assignment copies contents
  els[ei] = els[ei2];
  CHECK (els[ei].GetType() == PRISM);
  CHECK (els[ei][4] == PointIndex::FromNr0(24));
  els[ei] = tet;
  CHECK (els[ei].GetNP() == 4);

  // too wide for the slot
  T_VOLELEMENTS narrow(4);
  auto ni = narrow.Append();
  CHECK_THROWS (narrow[ni] = prism);

  v1.Invert();
  CHECK (v1.PNum(3) == tet.PNum(4));
  CHECK (v1.PNum(4) == tet.PNum(3));

  int cnt = 0;
  for (auto el : els) cnt += el.GetNP();
  for (auto el : cels) cnt += el.GetNP();
  for (auto el : els.Range(els.Range())) cnt += el.GetNP();
  CHECK (cnt == 30);
}

TEST_CASE("ElementRef geometry")
{
  T_POINTS pts;
  Element tet(TET);
  tet[0] = pts.Append (MeshPoint(Point<3>(0,0,0)));
  tet[1] = pts.Append (MeshPoint(Point<3>(1,0,0)));
  tet[2] = pts.Append (MeshPoint(Point<3>(0,1,0)));
  tet[3] = pts.Append (MeshPoint(Point<3>(0,0,1)));

  T_VOLELEMENTS els(4);
  auto ei = els.Append (tet);
  const ElementRef v = els[ei];
  CHECK (v.Volume(pts) == Approx(tet.Volume(pts)));
  CHECK (fabs(v.Volume(pts)) == Approx(1.0/6));
  CHECK (v.GetNIP() == tet.GetNIP());
  CHECK (v.CalcJacobianBadness(pts) == Approx(tet.CalcJacobianBadness(pts)));
  Element2d face, face2;
  v.GetFace (1, face);
  tet.GetFace (1, face2);
  CHECK (face[0] == face2[0]);
  CHECK (face[2] == face2[2]);
  Box3d box;
  v.GetBox (pts, box);
  CHECK (box.PMax()(0) == Approx(1.0));
}

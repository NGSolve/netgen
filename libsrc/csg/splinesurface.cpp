
#include <csg.hpp>

namespace netgen
{
void SplineSurface :: AppendPoint(const Point<3> & p, const double reffac, const bool hpref)
{
  auto pp = p;
  Project(pp);
  geompoints.Append(GeomPoint<3>(pp,reffac));
  geompoints.Last().hpref = hpref;
}
  
  void SplineSurface :: AppendSegment(shared_ptr<SplineSeg<3>> sp, string & bcname, double amaxh)
  {
    splines.Append(sp);
    bcnames.Append(bcname);
    maxh.Append(amaxh);
  }

  string SplineSurface :: GetBCNameOf (Point<3> p1, Point<3> p2) const
  {
    
    for(int i=0; i<splines.Size(); i++)
      {
	auto pp1 = Point<3>(splines[i]->GetPoint(0));
	Project(pp1);
	auto pp2 = Point<3>(splines[i]->GetPoint(1));
	Project(pp2);
        double eps = (p1-p2).Length() * 1e-4;
	if (((pp1-p1).Length()<eps && (pp2-p2).Length() < eps) || ((pp1-p2).Length() < eps && (pp2-p1).Length() < eps))
	  {
	    return bcnames[i];
	  }
      }
    return "default";
  }

  const shared_ptr<NgArray<shared_ptr<OneSurfacePrimitive>>> SplineSurface :: CreateCuttingSurfaces()
  {
    if(all_cuts)
      return all_cuts;
    auto cuttings = make_shared<NgArray<shared_ptr<OneSurfacePrimitive>>>();
    for (auto cut : *cuts)
      cuttings->Append(cut);
    for(int i = 0; i<splines.Size(); i++)
      {
	auto spline = splines[i];
	auto lineseg = dynamic_cast<LineSeg<3>*>(spline.get());
        if(lineseg)
          {
            auto p1 = Point<3>(spline->GetPoint(0));
            Project(p1);
            auto p2 = Point<3>(spline->GetPoint(1));
            Project(p2);
            auto vec = Vec<3>(p2)-Vec<3>(p1);
            auto plane = make_shared<Plane>(p1,-Cross(vec,baseprimitive->GetNormalVector(p1)));
            if(maxh[i]>0)
              {
                plane->SetMaxH(maxh[i]);
              }
            cuttings->Append(plane);
          }
        else
          {
            auto spline3 = dynamic_cast<SplineSeg3<3>*>(spline.get());
            if(spline3)
              {
                auto p1 = Point<3>(spline3->StartPI());
                Project(p1);
                auto p2 = Point<3>(spline3->TangentPoint());
                Project(p2);
                auto p3 = Point<3>(spline3->EndPI());
                Project(p3);
                Vec<3> v1 = p2-p1;
                Vec<3> v2 = p2-p3;
                Point<3> mid = p1 - v2;
                cout << "mid point = " << mid << endl;
                cout << "v1 = " << v1 << endl;
                cout << "v2 = " << v2 << endl;
                auto cyl = make_shared<EllipticCylinder>(mid, v1, v2);
                if(maxh[i] > 0)
                  cyl->SetMaxH(maxh[i]);
                cuttings->Append(cyl);
              }
            else
              throw NgException("Spline type not implemented for SplineSurface!");
          }
      }
    all_cuts = cuttings;
    return cuttings;
  }
  
  void SplineSurface :: Print(ostream & str) const
{
  str << "SplineSurface with base " << *baseprimitive << endl;
}

  static RegisterClassForArchive<SplineSurface, OneSurfacePrimitive> regss;
}

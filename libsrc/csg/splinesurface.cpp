
#include <csg.hpp>

namespace netgen
{
void SplineSurface :: AppendPoint(const Point<3> & p, const double reffac, const bool hpref)
{
  auto pp = Point<3>(p);
  geompoints.push_back(GeomPoint<3>(pp,reffac));
  geompoints.back().hpref = hpref;
}
  
  void SplineSurface :: AppendSegment(shared_ptr<SplineSeg<3>> sp, string & bcname, double amaxh)
  {
    splines.push_back(sp);
    bcnames.push_back(bcname);
    maxh.Append(amaxh);
  }

  string SplineSurface :: GetBCNameOf (Point<3> p1, Point<3> p2) const
  {
    
    double eps = 1e-5;
    for(int i=0; i<splines.size(); i++)
      {
	auto pp1 = Point<3>(splines[i]->GetPoint(0));
	Project(pp1);
	auto pp2 = Point<3>(splines[i]->GetPoint(1));
	Project(pp2);
	if (((pp1-p1).Length()<eps && (pp2-p2).Length() < eps) || ((pp1-p2).Length() < eps && (pp2-p1).Length() < eps))
	  {
	    return bcnames[i];
	  }
      }
    return "default";
  }

  const shared_ptr<std::vector<shared_ptr<OneSurfacePrimitive>>> SplineSurface :: CreateCuttingSurfaces()
  {
    if(all_cuts)
      return all_cuts;
    auto cuttings = make_shared<std::vector<shared_ptr<OneSurfacePrimitive>>>();
    for (auto cut : *cuts)
      cuttings->push_back(cut);
    for(int i = 0; i<splines.size(); i++)
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
            cuttings->push_back(plane);
          }
        else
          throw NgException("Spline type not implemented for SplineSurface!");
      }
    all_cuts = cuttings;
    return cuttings;
  }
  
  void SplineSurface :: Print(ostream & str) const
{
  str << "SplineSurface with base " << *baseprimitive << endl;
}

}

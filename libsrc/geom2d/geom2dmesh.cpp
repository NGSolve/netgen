#include <meshing.hpp>
#include <geometry2d.hpp>

namespace netgen
{

  Refinement2d :: Refinement2d (const SplineGeometry2d & ageometry)
    : Refinement(), geometry(ageometry)
  {
    ;
  }

  Refinement2d :: ~Refinement2d ()
  {
    ;
  }
  

  void Refinement2d :: 
  PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
		int surfi, 
		const PointGeomInfo & gi1, 
		const PointGeomInfo & gi2,
		Point<3> & newp, PointGeomInfo & newgi) const
  {
    newp = p1+secpoint*(p2-p1);
    newgi.trignum = 1;
  }



  void Refinement2d :: 
  PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint, 
		int surfi1, int surfi2, 
		const EdgePointGeomInfo & ap1, 
		const EdgePointGeomInfo & ap2,
		Point<3> & newp, EdgePointGeomInfo & newgi) const
  {
    Point<2> p2d;
    double newdist;
    auto spline = geometry.GetSplines().Get(ap1.edgenr);
    if( (ap1.dist == 0.0) && (ap2.dist == 0.0) )
      {
        // used for manually generated meshes
        const SplineSeg3<2> * ss3;
        const LineSeg<2> * ls;
        auto ext = dynamic_cast<const SplineSegExt *>(spline);
        if(ext)
          {
            ss3 = dynamic_cast<const SplineSeg3<2> *>(ext->seg);
            ls = dynamic_cast<const LineSeg<2> *>(ext->seg);
          }
        else
          {
            ss3 = dynamic_cast<const SplineSeg3<2> *>(spline);
            ls = dynamic_cast<const LineSeg<2> *>(spline);
          }
        Point<2> p12d(p1(0),p1(1)), p22d(p2(0),p2(1));
        Point<2> p1_proj(0.0,0.0), p2_proj(0.0,0.0);
        double t1_proj = 0.0;
        double t2_proj = 0.0;
        if(ss3)
          {
            ss3->Project(p12d,p1_proj,t1_proj);
            ss3->Project(p22d,p2_proj,t2_proj);
          }
        else if(ls)
          {
            ls->Project(p12d,p1_proj,t1_proj);
            ls->Project(p22d,p2_proj,t2_proj);
          }
        p2d = spline->GetPoint (((1-secpoint)*t1_proj+secpoint*t2_proj));
        newdist = (1-secpoint)*t1_proj+secpoint*t2_proj;
      }
    else
      {
        p2d = spline->GetPoint (((1-secpoint)*ap1.dist+secpoint*ap2.dist));
        newdist = (1-secpoint)*ap1.dist+secpoint*ap2.dist;
      }
  
    //  (*testout) << "refine 2d line, ap1.dist, ap2.dist = " << ap1.dist << ", " << ap2.dist << endl;
    //  (*testout) << "p1, p2 = " << p1 << p2 << ", newp = " << p2d << endl;

    newp = Point3d (p2d(0), p2d(1), 0);
    newgi.edgenr = ap1.edgenr;
    newgi.dist = newdist;
  };



  Vec<3> Refinement2d :: GetTangent (const Point<3> & p, int surfi1, int surfi2,
                                     const EdgePointGeomInfo & ap1) const
  {
    Vec<2> t2d = geometry.GetSplines().Get(ap1.edgenr) -> GetTangent(ap1.dist);
    return Vec<3> (t2d(0), t2d(1), 0);
  }

  Vec<3> Refinement2d :: GetNormal (const Point<3> & p, int surfi1, 
                                    const PointGeomInfo & gi) const
  {
    return Vec<3> (0,0,1);
  }


  void Refinement2d :: ProjectToSurface (Point<3> & p, int surfi, const PointGeomInfo & /* gi */) const
  {
    p(2) = 0;
  }


  void Refinement2d :: ProjectToEdge (Point<3> & p, int surfi1, int surfi2, 
                                      const EdgePointGeomInfo & egi) const
  {
    Point<2> p2d (p(0), p(1)), pp;
    double t;
    geometry.GetSplines().Get(egi.edgenr) -> Project (p2d, pp, t);
    p = Point<3> (pp(0), pp(1), 0);
  }
}

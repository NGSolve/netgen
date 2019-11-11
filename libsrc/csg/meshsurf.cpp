#include <mystdlib.h>

#include <csg.hpp>
#include <meshing.hpp>



namespace netgen
{
  /*
Meshing2Surfaces :: Meshing2Surfaces (const Surface & asurface)
  : surface(asurface)
{
  ;
}
  */
  Meshing2Surfaces :: Meshing2Surfaces (const CSGeometry& geo,
                                        const Surface & asurf,
                                        const MeshingParameters & mp,
                                        const Box<3> & abb)
    : Meshing2(geo, mp, abb), surface(asurf), mparam (mp)
  {
    ;
  }


void Meshing2Surfaces :: DefineTransformation (const Point<3> & p1, const Point<3> & p2,
					       const PointGeomInfo * geominfo1,
					       const PointGeomInfo * geominfo2)
{
  ((Surface&)surface).DefineTangentialPlane (p1, p2);
}

void Meshing2Surfaces :: TransformToPlain (const Point<3> & locpoint, 
					   const MultiPointGeomInfo & geominfo,
					   Point<2> & planepoint, 
					   double h, int & zone)
{
  surface.ToPlane (locpoint, planepoint, h, zone);
}

int Meshing2Surfaces :: TransformFromPlain (const Point<2> & planepoint,
                                            Point<3> & locpoint, 
                                            PointGeomInfo & gi,
                                            double h)
{
  surface.FromPlane (planepoint, locpoint, h);
  gi.trignum = 1;
  return 0;
}



double Meshing2Surfaces :: CalcLocalH (const Point<3> & p, double gh) const
{
  return surface.LocH (p, 3, 1, mparam, gh);
  /*
    double loch = mesh.lochfunc->GetH(p);
    if (gh < loch) loch = gh;
    return loch;
    */
}
}

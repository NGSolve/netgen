#ifdef OCCGEOMETRY

#ifndef FILE_OCCMESHSURF
#define FILE_OCCMESHSURF

#include "occgeom.hpp"
#include "mydefs.hpp"

#include <TopoDS_Face.hxx>
#include <Geom_Surface.hxx>
#include <ShapeAnalysis.hxx>

#define PARAMETERSPACE -1
#define PLANESPACE     1

namespace netgen
{
class OCCGeometry;

class SingularMatrixException
{};

class UVBoundsException
{};

class OCCSurface
{
public:
  TopoDS_Face topods_face;
  Handle(Geom_Surface) occface;
  TopAbs_Orientation orient;
  int projecttype;

protected:
  Point<3> p1;
  Point<3> p2;

  /// in plane, directed p1->p2
  Vec<3> ex;
  /// in plane
  Vec<3> ey;
  /// outer normal direction
  Vec<3> ez;

  /// normal vector in p2
  Vec<3> n2;

  /// average normal vector
  Vec<3> nmid;

  // for transformation to parameter space
  Point<2> psp1;
  Point<2> psp2;
  Vec<2> psex;
  Vec<2> psey;
  Mat<2,2> Amat, Amatinv;

  // UV Bounds
  double umin, umax, vmin, vmax;

public:
  OCCSurface (const TopoDS_Face & aface, int aprojecttype)
  {
    static Timer t("occurface ctor"); RegionTimer r(t);
    topods_face = aface;
    occface = BRep_Tool::Surface(topods_face);
    orient = topods_face.Orientation();
    projecttype = aprojecttype;
    ShapeAnalysis::GetFaceUVBounds (topods_face, umin, umax, vmin, vmax);
    umin -= fabs(umax-umin)/100.0;
    vmin -= fabs(vmax-vmin)/100.0;
    umax += fabs(umax-umin)/100.0;
    vmax += fabs(vmax-vmin)/100.0;
    // projecttype = PLANESPACE;
    /*
    TopExp_Explorer exp1;
    exp1.Init (topods_face, TopAbs_WIRE);
    orient = TopAbs::Compose (orient, exp1.Current().Orientation());
    */
  };
  
  ~OCCSurface()
  {};

  void Project (Point<3> & p, PointGeomInfo & gi);

  void GetNormalVector (const Point<3> & p,
			const PointGeomInfo & geominfo,
			Vec<3> & n) const;

  /**
    Defines tangential plane in ap1.
    The local x-coordinate axis point to the direction of ap2 */
  void DefineTangentialPlane (const Point<3> & ap1, 
			      const PointGeomInfo & geominfo1,
			      const Point<3> & ap2,
			      const PointGeomInfo & geominfo2);


  /// Transforms 3d point p3d to local coordinates pplane
  void ToPlane (const Point<3> & p3d, const PointGeomInfo & geominfo,
		Point<2> & pplane, double h, int & zone) const;
  
  /// Transforms point pplane in local coordinates to 3d point
  void FromPlane (const Point<2> & pplane, 
		  Point<3> & p3d,
		  PointGeomInfo & gi,
		  double h);
};



///
class Meshing2OCCSurfaces : public Meshing2
{
  ///
  OCCSurface surface;
 

public:
  ///
  Meshing2OCCSurfaces (const NetgenGeometry& geo,
                       const TopoDS_Shape & asurf, const Box<3> & aboundingbox,
                       int aprojecttype, const MeshingParameters & mparam);

  ///
  int GetProjectionType ()
  { return surface.projecttype; }

protected:
  ///
  void DefineTransformation (const Point<3> & p1, const Point<3> & p2,
                             const PointGeomInfo * geominfo1,
                             const PointGeomInfo * geominfo2) override;
  ///
  void TransformToPlain (const Point<3> & locpoint, 
                         const MultiPointGeomInfo & geominfo,
                         Point<2> & plainpoint, 
                         double h, int & zone) override;
  ///

  int TransformFromPlain (const Point<2> & plainpoint,
                          Point<3> & locpoint,
                          PointGeomInfo & gi,
                          double h) override;
  ///
  double CalcLocalH (const Point<3> & p, double gh) const override;
  
};

class OCCGeometry;

} // namespace netgen

#endif

#endif

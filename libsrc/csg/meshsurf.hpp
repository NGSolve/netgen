#ifndef FILE_MESHSURF
#define FILE_MESHSURF

namespace netgen
{

  ///
  class Meshing2Surfaces : public Meshing2
  {
    ///
    const Surface & surface;
    
    /// should be movec to base ... 
    const MeshingParameters & mparam;
  public:
    ///
    //  Meshing2Surfaces (const Surface & asurf);
    ///
    Meshing2Surfaces (const Surface & asurf, const MeshingParameters & mp, 
		      const Box<3> & aboundingbox);

  protected:
    ///
    void DefineTransformation(const Point<3> & p1,
                              const Point<3> & p2,
                              const PointGeomInfo * geominfo1,
                              const PointGeomInfo * geominfo2) override;
    ///
    void TransformToPlain(const Point<3> & locpoint, 
                          const MultiPointGeomInfo & geominfo,
                          Point<2> & plainpoint, 
                          double h, int & zone) override;
    ///
    int TransformFromPlain(const Point<2>& plainpoint,
                           Point<3>& locpoint, 
                           PointGeomInfo & gi,
                           double h) override;
    ///
    double CalcLocalH(const Point3d & p, double gh) const override;
  };



  ///
  class MeshOptimize2dSurfaces : public MeshOptimize2d
  {
    ///
    const CSGeometry & geometry;

  public:
    ///
    MeshOptimize2dSurfaces (const CSGeometry & ageometry); 
   
    ///
    virtual void ProjectPoint (INDEX surfind, Point<3> & p) const override;
    ///
    virtual void ProjectPoint2 (INDEX surfind, INDEX surfind2, Point<3> & p) const override;
    ///
    virtual void GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const override;
  };



  class RefinementSurfaces : public Refinement
  {
    const CSGeometry & geometry;

  public:
    RefinementSurfaces (const CSGeometry & ageometry);
    virtual ~RefinementSurfaces ();
  
    virtual void PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
			       int surfi, 
			       const PointGeomInfo & gi1, 
			       const PointGeomInfo & gi2,
			       Point<3> & newp, PointGeomInfo & newgi) const override;

    virtual void PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
			       int surfi1, int surfi2, 
			       const EdgePointGeomInfo & ap1, 
			       const EdgePointGeomInfo & ap2,
			       Point<3> & newp, EdgePointGeomInfo & newgi) const override;

    virtual Vec<3> GetTangent (const Point<3> & p, int surfi1, int surfi2,
			       const EdgePointGeomInfo & ap1) const override;

    virtual Vec<3> GetNormal (const Point<3> & p, int surfi1, 
			      const PointGeomInfo & gi) const override;


    virtual void ProjectToSurface (Point<3> & p, int surfi) const override;

    virtual void ProjectToEdge (Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & egi) const override;

  };

}

#endif


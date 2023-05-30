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
    Meshing2Surfaces (const CSGeometry& geo,
                      const Surface & asurf,
                      const MeshingParameters & mp,
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
    double CalcLocalH(const Point<3> & p, double gh) const override;
  };
}

#endif


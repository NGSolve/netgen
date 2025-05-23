#ifndef FILE_SURFACEGEOM
#define FILE_SURFACEGEOM

/* *************************************************************************/
/* File:   surfacegeom.hpp                                                 */
/* Author: Michael Neunteufel                                              */
/* Date:   Jun. 2020                                                       */
/* *************************************************************************/


#include <functional>


namespace netgen
{

  class DLL_HEADER SurfaceGeometry : public NetgenGeometry
  {
    function<Vec<3>(Point<2>)> func;
    double eps=1e-4;

  private:

    void CalcHesse(double u, double v, Vec<3>& f_uu, Vec<3>& f_vv, Vec<3>& f_uv) const;
  public:

    SurfaceGeometry();
    SurfaceGeometry(function<Vec<3>(Point<2>)> func);
    SurfaceGeometry(const SurfaceGeometry& geom);
    SurfaceGeometry& operator =(const SurfaceGeometry& geom)
    {
      func = geom.func;
      eps = geom.eps;
      return *this;
    }

    Array<Vec<3>> GetTangentVectors(double u, double v) const;

    void GetTangentVectors(double u, double v, Array<Vec<3>>& tang) const;


    virtual Vec<3> GetNormal(int surfind, const Point<3> & p, const PointGeomInfo* gi) const override;

    virtual PointGeomInfo ProjectPoint(int surfind, Point<3> & p) const override;
    
    virtual void ProjectPointEdge (int surfind, int surfind2, Point<3> & p,
                           EdgePointGeomInfo* gi = nullptr) const override;
    
    virtual bool ProjectPointGI (int surfind, Point<3> & p, PointGeomInfo & gi) const override;

    virtual bool CalcPointGeomInfo(int surfind, PointGeomInfo& gi, const Point<3> & p3) const override;

    virtual void PointBetweenEdge(const Point<3> & p1, const Point<3> & p2, double secpoint,
                          int surfi1, int surfi2, 
                          const EdgePointGeomInfo & ap1, 
                          const EdgePointGeomInfo & ap2,
                          Point<3> & newp, EdgePointGeomInfo & newgi) const override;
    
    virtual void PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint,
                      int surfi, 
                      const PointGeomInfo & gi1, 
                      const PointGeomInfo & gi2,
                      Point<3> & newp, PointGeomInfo & newgi) const override;

    int GenerateStructuredMesh(shared_ptr<Mesh> & mesh, bool quads, int nx, int ny, bool flip_triangles, const Array<Point<3>>& bbbpts, const Array<string>& bbbnames, const Array<Point<3>>& hppoints, const Array<float>& hpptsfac, const Array<string>& hpbnd, const Array<float>& hpbndfac, Array<double> layer_thickness[4], bool layer_quad);

  };
  
}
   


#endif //SURFACEGEOM

#ifndef FILE_CALLBACKGEOM
#define FILE_CALLBACKGEOM

/* *************************************************************************/
/* File:   callbackgeom.hpp                                                */
/* Author: K. Sugahara                                                     */
/* Date:   Mar. 2026                                                       */
/*                                                                         */
/* NetgenGeometry subclass that delegates surface projection to Python      */
/* callbacks. Enables mesh.Curve(order) with any CAD kernel (e.g. Cubit    */
/* ACIS) without requiring OCC or STEP files.                              */
/* *************************************************************************/

#include <functional>

namespace netgen
{

  class DLL_HEADER CallbackGeometry : public NetgenGeometry
  {
  public:
    // Callback types:
    //   project_point: (surfnr, x, y, z, u_hint, v_hint) -> (x_proj, y_proj, z_proj, u, v)
    //   get_normal:    (surfnr, x, y, z) -> (nx, ny, nz)
    using ProjectFunc = std::function<std::tuple<double,double,double,double,double>
                                      (int surfnr, double x, double y, double z,
                                       double u_hint, double v_hint)>;
    using NormalFunc = std::function<std::tuple<double,double,double>
                                    (int surfnr, double x, double y, double z)>;

  private:
    ProjectFunc project_func;
    NormalFunc normal_func;
    int num_surfaces;

  public:
    CallbackGeometry() : num_surfaces(0) {}

    CallbackGeometry(ProjectFunc _project, NormalFunc _normal, int _num_surfaces)
      : project_func(_project), normal_func(_normal), num_surfaces(_num_surfaces) {}

    virtual Vec<3> GetNormal(int surfind, const Point<3> & p,
                             const PointGeomInfo* gi) const override;

    virtual PointGeomInfo ProjectPoint(int surfind, Point<3> & p) const override;

    virtual void ProjectPointEdge(int surfind, int surfind2, Point<3> & p,
                                  EdgePointGeomInfo* gi = nullptr) const override;

    virtual bool ProjectPointGI(int surfind, Point<3> & p,
                                PointGeomInfo & gi) const override;

    virtual bool CalcPointGeomInfo(int surfind, PointGeomInfo& gi,
                                   const Point<3> & p3) const override;

    virtual void PointBetweenEdge(const Point<3> & p1, const Point<3> & p2,
                                  double secpoint,
                                  int surfi1, int surfi2,
                                  const EdgePointGeomInfo & ap1,
                                  const EdgePointGeomInfo & ap2,
                                  Point<3> & newp,
                                  EdgePointGeomInfo & newgi) const override;

    virtual void PointBetween(const Point<3> & p1, const Point<3> & p2,
                              double secpoint,
                              int surfi,
                              const PointGeomInfo & gi1,
                              const PointGeomInfo & gi2,
                              Point<3> & newp,
                              PointGeomInfo & newgi) const override;
  };

}

#endif // FILE_CALLBACKGEOM

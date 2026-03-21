/* *************************************************************************/
/* File:   callbackgeom.cpp                                                */
/* Author: K. Sugahara                                                     */
/* Date:   Mar. 2026                                                       */
/*                                                                         */
/* NetgenGeometry with Python callback-based surface projection.           */
/* Enables mesh.Curve(order) with Cubit ACIS or any external CAD kernel.   */
/* *************************************************************************/

#include <meshing.hpp>
#include "callbackgeom.hpp"

namespace netgen
{

  Vec<3> CallbackGeometry::GetNormal(int surfind, const Point<3> & p,
                                     const PointGeomInfo* gi) const
  {
    if (!normal_func)
      return Vec<3>(0, 0, 1);  // fallback

    auto [nx, ny, nz] = normal_func(surfind, p[0], p[1], p[2]);
    Vec<3> n(nx, ny, nz);
    double len = n.Length();
    if (len > 1e-15)
      n /= len;
    return n;
  }

  PointGeomInfo CallbackGeometry::ProjectPoint(int surfind, Point<3> & p) const
  {
    PointGeomInfo gi;
    gi.trignum = surfind;
    gi.u = 0.5;
    gi.v = 0.5;
    ProjectPointGI(surfind, p, gi);
    return gi;
  }

  void CallbackGeometry::ProjectPointEdge(int surfind, int surfind2, Point<3> & p,
                                           EdgePointGeomInfo* gi) const
  {
    // For edges shared between two surfaces, project onto first surface
    PointGeomInfo pgi;
    pgi.trignum = surfind;
    pgi.u = (gi != nullptr) ? gi->u : 0.5;
    pgi.v = (gi != nullptr) ? gi->v : 0.5;
    ProjectPointGI(surfind, p, pgi);
    if (gi)
    {
      gi->u = pgi.u;
      gi->v = pgi.v;
    }
  }

  bool CallbackGeometry::ProjectPointGI(int surfind, Point<3> & p,
                                         PointGeomInfo & gi) const
  {
    if (!project_func)
      return false;

    auto [xp, yp, zp, u, v] = project_func(surfind, p[0], p[1], p[2],
                                             gi.u, gi.v);
    p = Point<3>(xp, yp, zp);
    gi.u = u;
    gi.v = v;
    gi.trignum = surfind;
    return true;
  }

  bool CallbackGeometry::CalcPointGeomInfo(int surfind, PointGeomInfo& gi,
                                            const Point<3> & p3) const
  {
    if (!project_func)
      return false;

    Point<3> pp = p3;
    gi.trignum = surfind;
    gi.u = 0.5;
    gi.v = 0.5;
    return ProjectPointGI(surfind, pp, gi);
  }

  void CallbackGeometry::PointBetweenEdge(const Point<3> & p1, const Point<3> & p2,
                                           double secpoint,
                                           int surfi1, int surfi2,
                                           const EdgePointGeomInfo & ap1,
                                           const EdgePointGeomInfo & ap2,
                                           Point<3> & newp,
                                           EdgePointGeomInfo & newgi) const
  {
    // Interpolate UV, then project
    newgi.u = ap1.u + secpoint * (ap2.u - ap1.u);
    newgi.v = ap1.v + secpoint * (ap2.v - ap1.v);
    newgi.edgenr = ap1.edgenr;
    newgi.body = -1;
    newgi.dist = -1.0;

    // Linear interpolation as initial guess
    newp = p1 + secpoint * (p2 - p1);

    // Project onto surface using callback
    if (project_func)
    {
      auto [xp, yp, zp, u, v] = project_func(surfi1, newp[0], newp[1], newp[2],
                                               newgi.u, newgi.v);
      newp = Point<3>(xp, yp, zp);
      newgi.u = u;
      newgi.v = v;
    }
  }

  void CallbackGeometry::PointBetween(const Point<3> & p1, const Point<3> & p2,
                                       double secpoint,
                                       int surfi,
                                       const PointGeomInfo & gi1,
                                       const PointGeomInfo & gi2,
                                       Point<3> & newp,
                                       PointGeomInfo & newgi) const
  {
    // Interpolate UV hints
    newgi.u = gi1.u + secpoint * (gi2.u - gi1.u);
    newgi.v = gi1.v + secpoint * (gi2.v - gi1.v);
    newgi.trignum = surfi;

    // Linear interpolation as initial guess
    newp = p1 + secpoint * (p2 - p1);

    // Project onto surface using callback
    if (project_func)
    {
      auto [xp, yp, zp, u, v] = project_func(surfi, newp[0], newp[1], newp[2],
                                               newgi.u, newgi.v);
      newp = Point<3>(xp, yp, zp);
      newgi.u = u;
      newgi.v = v;
    }
  }

} // namespace netgen

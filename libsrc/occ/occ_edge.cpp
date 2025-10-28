#include <StdFail_NotDone.hxx>
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GCPnts_AbscissaPoint.hxx>

#pragma clang diagnostic pop

#include "occ_edge.hpp"
#include "occgeom.hpp"

namespace netgen
{
    OCCEdge::OCCEdge(TopoDS_Shape edge_, GeometryVertex & start_, GeometryVertex & end_)
        : GeometryEdge(start_, end_),
          edge(TopoDS::Edge(edge_))
    {
        curve = BRep_Tool::Curve(edge, s0, s1);
        BRepGProp::LinearProperties(edge, props);

        auto verts = GetVertices(edge);
        if(verts.size() != 2)
            throw Exception("OCC edge does not have 2 vertices");

        if(start != end)
        {
            // swap start/end if necessary
            double d00 = Dist(GetPoint(0), start->GetPoint());
            double d01 = Dist(GetPoint(0), end->GetPoint());
            if(d01 < d00)
                swap(start, end);
          static_cast<OCCVertex*>(start)->SetPoint(GetPoint(0));
          static_cast<OCCVertex*>(end)->SetPoint(GetPoint(1));
        }
    }

    double OCCEdge::GetLength() const
    {
        return props.Mass();
    }

    Point<3> OCCEdge::GetCenter() const
    {
        return occ2ng( props.CentreOfMass() );
    }

    Point<3> OCCEdge::GetPoint(double t) const
    {
        return occ2ng( curve->Value(s0+t*(s1-s0)) );
    }

    double OCCEdge::CalcStep(double t, double sag) const
    {
        throw Exception(ToString("not implemented") + __FILE__ + ":" + ToString(__LINE__));
    }

    void OCCEdge::ProjectPoint(Point<3>& p, EdgePointGeomInfo* gi) const
    {
        auto pnt = ng2occ(p);
        // extend the projection parameter range, else projection might fail
        // for an endpoint
        // see discussion here: https://forum.ngsolve.org/t/how-to-apply-occidentification-correctly/2555
        // I do not see a better way using occ tolerances?
        double eps = 1e-7 * (s1-s0);
        GeomAPI_ProjectPointOnCurve proj(pnt, curve, s0-eps, s1+eps);
        pnt = proj.NearestPoint();
        if(gi)
            gi->dist = (proj.LowerDistanceParameter() - s0)/(s1-s0);
        p = occ2ng(pnt);
    }

    void OCCEdge::PointBetween(const Point<3>& p1,
                      const Point<3>& p2,
                      double secpoint,
                      const EdgePointGeomInfo& gi1,
                      const EdgePointGeomInfo& gi2,
                      Point<3>& newp,
                      EdgePointGeomInfo& newgi) const
    {
      static Timer tim("OCCEdge::PointBetween");
      RegionTimer rtim(tim);
      // try {
      //   GeometryEdge::PointBetween(p1, p2, secpoint, gi1, gi2, newp, newgi);
      // }
      // catch (const StdFail_NotDone &)
      {
        static Timer tim("OCCEdge::PointBetween-Fallback");
        RegionTimer rtim(tim);

        double t0 = gi1.dist;
        double t1 = gi2.dist;
        if (t0 > t1)
        {
          swap(t0, t1);
          secpoint = 1.0 - secpoint;
        }

        GeomAdaptor_Curve adaptor(curve, t0, t1);
        Standard_Real length = GCPnts_AbscissaPoint::Length(adaptor, t0, t1);
        newgi.dist = GCPnts_AbscissaPoint(adaptor, secpoint * length, t0).Parameter();
        newp = GetPoint(newgi.dist);
      }
    }

    Vec<3> OCCEdge::GetTangent(double t) const
    {
        t = s0 + t*(s1-s0);
        gp_Pnt p;
        gp_Vec v;
        curve->D1(t, p, v);
        return occ2ng(v) * (s1-s0);
    }

}

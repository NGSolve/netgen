#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

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
        if(proj.NbPoints() > 0)
        {
            if(gi)
                gi->dist = (proj.LowerDistanceParameter() - s0)/(s1-s0);
            p = occ2ng(proj.NearestPoint());
        }
        else
        {
            double bests = s0, bestd = 1e99;
            const int N = 100;
            for(int i = 0; i <= N; i++)
            {
                double s = s0 + (s1-s0)*double(i)/N;
                double d = curve->Value(s).SquareDistance(pnt);
                if(d < bestd) { bestd = d; bests = s; }
            }
            for(int it = 0; it < 20; it++)
            {
                gp_Pnt c; gp_Vec d1, d2;
                curve->D2(bests, c, d1, d2);
                gp_Vec r(pnt, c);
                double f = r * d1;
                double fp = d1 * d1 + r * d2;
                if(fabs(fp) < 1e-30) break;
                double snew = min(s1, max(s0, bests - f/fp));
                if(fabs(snew-bests) < 1e-12*(s1-s0)) { bests = snew; break; }
                bests = snew;
            }
            if(gi)
                gi->dist = (bests - s0)/(s1-s0);
            p = occ2ng(curve->Value(bests));
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

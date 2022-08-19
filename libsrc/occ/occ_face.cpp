#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

#include "occ_edge.hpp"
#include "occ_face.hpp"
#include "occgeom.hpp"

namespace netgen
{
    OCCFace::OCCFace(TopoDS_Shape dshape)
        : face(TopoDS::Face(dshape))
    {
        BRepGProp::SurfaceProperties (dshape, props);
        bbox = ::netgen::GetBoundingBox(face);

        surface = BRep_Tool::Surface(face);
        shape_analysis = new ShapeAnalysis_Surface( surface );
        tolerance = BRep_Tool::Tolerance( face );
    }

    size_t OCCFace::GetNBoundaries() const
    {
        return 0;
    }

    size_t OCCFace::GetHash() const
    {
      return face.HashCode(std::numeric_limits<Standard_Integer>::max());
    }

    Point<3> OCCFace::GetCenter() const
    {
        return occ2ng( props.CentreOfMass() );
    }

    Array<Segment> OCCFace::GetBoundary(const Mesh& mesh) const
    {
        auto & geom = dynamic_cast<OCCGeometry&>(*mesh.GetGeometry());

        auto n_edges = geom.edge_map.size();
        constexpr int UNUSED = 0;
        constexpr int FORWARD = 1;
        constexpr int REVERSED = 2;
        constexpr int BOTH = 3;

        Array<int> edge_orientation(n_edges);
        edge_orientation = UNUSED;

        Array<Handle(Geom2d_Curve)> curve_on_face[BOTH];
        curve_on_face[FORWARD].SetSize(n_edges);
        curve_on_face[REVERSED].SetSize(n_edges);

        Array<TopoDS_Edge> edge_on_face[BOTH];
        edge_on_face[FORWARD].SetSize(n_edges);
        edge_on_face[REVERSED].SetSize(n_edges);

        for(auto edge_ : GetEdges(face))
        {
            auto edge = TopoDS::Edge(edge_);
            if(geom.edge_map.count(edge)==0)
                continue;
            auto edgenr = geom.edge_map[edge];
            auto & orientation = edge_orientation[edgenr];
            double s0, s1;
            auto cof = BRep_Tool::CurveOnSurface (edge, face, s0, s1);
            if(edge.Orientation() == TopAbs_FORWARD)
            {
                curve_on_face[FORWARD][edgenr] = cof;
                orientation += FORWARD;
                edge_on_face[FORWARD][edgenr] = edge;
            }
            if(edge.Orientation() == TopAbs_REVERSED)
            {
                curve_on_face[REVERSED][edgenr] = cof;
                orientation += REVERSED;
                edge_on_face[REVERSED][edgenr] = edge;
            }

            if(orientation > BOTH)
                throw Exception("have edge more than twice in face " + ToString(nr) + " " + properties.GetName() + ", orientation: " + ToString(orientation));
        }

        Array<Segment> boundary;
        for (auto seg : mesh.LineSegments())
        {
            auto edgenr = seg.epgeominfo[0].edgenr;
            auto orientation = edge_orientation[edgenr];

            if(orientation == UNUSED)
                continue;

            for(const auto ORIENTATION : {FORWARD, REVERSED})
            {
                if((orientation & ORIENTATION) == 0)
                    continue;

                // auto cof = curve_on_face[ORIENTATION][edgenr];
                auto edge = edge_on_face[ORIENTATION][edgenr];
                double s0, s1;
                auto cof = BRep_Tool::CurveOnSurface (edge, face, s0, s1);

                double s[2] = { seg.epgeominfo[0].dist, seg.epgeominfo[1].dist };

                // dist is in [0,1], map parametrization to [s0, s1]
                s[0] = s0 + s[0]*(s1-s0);
                s[1] = s0 + s[1]*(s1-s0);

                // fixes normal-vector roundoff problem when endpoint is cone-tip
                double delta = s[1]-s[0];
                s[0] += 1e-10*delta;
                s[1] -= 1e-10*delta;

                for(auto i : Range(2))
                {
                    auto uv = cof->Value(s[i]);
                    seg.epgeominfo[i].u = uv.X();
                    seg.epgeominfo[i].v = uv.Y();
                }

                if(ORIENTATION == REVERSED)
                {
                    swap(seg[0], seg[1]);
                    swap(seg.epgeominfo[0].dist, seg.epgeominfo[1].dist);
                    swap(seg.epgeominfo[0].u, seg.epgeominfo[1].u);
                    swap(seg.epgeominfo[0].v, seg.epgeominfo[1].v);
                }

                boundary.Append(seg);
            }
        }
        return boundary;
    }

    PointGeomInfo OCCFace::Project(Point<3>& p) const
    {
        auto suval = shape_analysis->ValueOfUV(ng2occ(p), tolerance);
        double u,v;
        suval.Coord(u, v);
        p = occ2ng(surface->Value( u, v ));

        PointGeomInfo gi;
        gi.trignum = nr+1;
        gi.u = u;
        gi.v = v;
        return gi;
    }

    bool OCCFace::ProjectPointGI(Point<3>& p_, PointGeomInfo& gi) const
    {
        double u = gi.u;
        double v = gi.v;
        auto p = ng2occ(p_);
        auto x = surface->Value (u,v);
      
        if (p.SquareDistance(x) <= sqr(PROJECTION_TOLERANCE)) return true;
      
        gp_Vec du, dv;
        surface->D1(u,v,x,du,dv);
      
        int count = 0;
        gp_Pnt xold;
        gp_Vec n;
        double det, lambda, mu;
      
        do {
           count++;

           n = du^dv;

           det = Det3 (n.X(), du.X(), dv.X(),
              n.Y(), du.Y(), dv.Y(),
              n.Z(), du.Z(), dv.Z());

           if (det < 1e-15) return false;

           lambda = Det3 (n.X(), p.X()-x.X(), dv.X(),
              n.Y(), p.Y()-x.Y(), dv.Y(),
              n.Z(), p.Z()-x.Z(), dv.Z())/det;

           mu     = Det3 (n.X(), du.X(), p.X()-x.X(),
              n.Y(), du.Y(), p.Y()-x.Y(),
              n.Z(), du.Z(), p.Z()-x.Z())/det;

           u += lambda;
           v += mu;

           xold = x;
           surface->D1(u,v,x,du,dv);

        } while (xold.SquareDistance(x) > sqr(PROJECTION_TOLERANCE) && count < 50);

        //    (*testout) << "FastProject count: " << count << endl;

        if (count == 50) return false;

        p_ = occ2ng(x);

        return true;
    }

    Point<3> OCCFace::GetPoint(const PointGeomInfo& gi) const
    {
        return occ2ng(surface->Value( gi.u, gi.v ));
    }

    void OCCFace::CalcEdgePointGI(const GeometryEdge& edge,
            double t,
            EdgePointGeomInfo& egi) const
    {
        throw Exception(ToString("not implemented") + __FILE__ + ":" + ToString(__LINE__));
    }

    Box<3> OCCFace::GetBoundingBox() const
    {
        return bbox;
    }


    double OCCFace::GetCurvature(const PointGeomInfo& gi) const
    {
        throw Exception(ToString("not implemented") + __FILE__ + ":" + ToString(__LINE__));
    }

    void OCCFace::RestrictH(Mesh& mesh, const MeshingParameters& mparam) const
    {
        throw Exception(ToString("not implemented") + __FILE__ + ":" + ToString(__LINE__));
    }

    Vec<3> OCCFace::GetNormal(const Point<3>& p, const PointGeomInfo* gi) const
    {
        PointGeomInfo gi_;
        if(gi==nullptr)
        {
            auto p_ = p;
            gi_ = Project(p_);
            gi = &gi_;
        }

        gp_Pnt pnt;
        gp_Vec du, dv;
        surface->D1(gi->u,gi->v,pnt,du,dv);
        auto n = Cross (occ2ng(du), occ2ng(dv));
        n.Normalize();
        if (face.Orientation() == TopAbs_REVERSED)
            n *= -1;
        return n;
    }


}

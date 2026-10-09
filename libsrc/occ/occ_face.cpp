#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <BRepLProp_SLProps.hxx>
#include <ShapeAnalysis.hxx>

#pragma clang diagnostic pop

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
        if (surface->IsUPeriodic()) uperiod = surface->UPeriod();
        if (surface->IsVPeriodic()) vperiod = surface->VPeriod();
        double umin, vmin;
        BRepTools::UVBounds (face, umin, umax, vmin, vmax);
        surface->Bounds (sumin, sumax, svmin, svmax);
    }

    void OCCFace::AlignGeomInfo(PointGeomInfo& gi1, PointGeomInfo& gi2) const
    {
        // at a singular point (pole, apex) the parameter is arbitrary
        auto singular = [&] (const PointGeomInfo & gi, bool & su, bool & sv)
          {
            gp_Pnt p;
            gp_Vec du, dv;
            surface->D1 (gi.u, gi.v, p, du, dv);
            double dum = du.Magnitude(), dvm = dv.Magnitude();
            su = dum < 1e-8 * dvm;
            sv = dvm < 1e-8 * dum;
          };
        bool su1, sv1, su2, sv2;
        singular (gi1, su1, sv1);
        singular (gi2, su2, sv2);
        if (su1 && !su2) gi1.u = gi2.u;
        if (su2 && !su1) gi2.u = gi1.u;
        if (sv1 && !sv2) gi1.v = gi2.v;
        if (sv2 && !sv1) gi2.v = gi1.v;
        auto periodic = [] (double & a, double & b, double period, double pmax)
          {
            if (period == 0 || fabs (a-b) <= 0.5*period) return;
            double & lo = (a < b) ? a : b;
            double & hi = (a < b) ? b : a;
            if (lo + period <= pmax + 1e-8*period) lo += period;
            else hi -= period;
          };
        periodic (gi1.u, gi2.u, uperiod, umax);
        periodic (gi1.v, gi2.v, vperiod, vmax);
    }

    size_t OCCFace::GetNBoundaries() const
    {
        return 0;
    }

    Point<3> OCCFace::GetCenter() const
    {
        return occ2ng( props.CentreOfMass() );
    }

    Array<Segment> OCCFace::GetBoundary(const Mesh& mesh) const
    {
        auto & geom = dynamic_cast<OCCGeometry&>(*mesh.GetGeometry());

        auto n_edges = geom.GetNEdges();
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

        // In case the face is INTERNAL, we need to orient it to FORWARD to get proper orientation for the edges
        // (relative to the face) otherwise, all edges are also INTERNAL
        auto oriented_face = TopoDS_Face(face);
        if(oriented_face.Orientation() == TopAbs_INTERNAL)
          oriented_face.Orientation(TopAbs_FORWARD);

        for(auto edge_ : GetEdges(oriented_face))
        {
            auto edge = TopoDS::Edge(edge_);
            auto edgenr = geom.GetEdge(edge).nr;
            auto & orientation = edge_orientation[edgenr];
            double s0, s1;
            auto cof = BRep_Tool::CurveOnSurface (edge, oriented_face, s0, s1);
            if(edge.Orientation() == TopAbs_FORWARD || edge.Orientation() == TopAbs_INTERNAL)
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
            if(edge.Orientation() == TopAbs_INTERNAL)
            {
              // add reversed edge
              auto r_edge = TopoDS::Edge(edge.Reversed());
              auto cof = BRep_Tool::CurveOnSurface (r_edge, oriented_face, s0, s1);
              curve_on_face[REVERSED][edgenr] = cof;
              orientation += REVERSED;
              edge_on_face[REVERSED][edgenr] = r_edge;
            }

            if(orientation > BOTH)
                throw Exception("have edge more than twice in face " + ToString(nr) + " " + properties.GetName() + ", orientation: " + ToString(orientation));
        }

        auto n_faces = static_cast<int>(geom.GetNFaces());

        double umin, umax, vmin, vmax;
        ShapeAnalysis::GetFaceUVBounds (face, umin, umax, vmin, vmax);
        double du = 0.01*(umax-umin), dv = 0.01*(vmax-vmin);

        Array<Segment> boundary;
        for (auto seg : mesh.LineSegments())
        {
            const auto & ed = mesh.GetEdgeDescriptor(seg.GetIndex());
            auto edgenr = ed.EdgeNr() - 1;
            if(edgenr < 0 || edgenr >= n_edges)
                continue;  // not on an edge of the geometry, handled below

            if((ed.SurfNr(0) > n_faces || ed.SurfNr(1) > n_faces) &&
               ed.SurfNr(0) != nr+1 && ed.SurfNr(1) != nr+1)
                continue;

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

                double s[2] = { seg.EPGeomInfo(0).dist, seg.EPGeomInfo(1).dist };

                // dist is in [0,1], map parametrization to [s0, s1]
                s[0] = s0 + s[0]*(s1-s0);
                s[1] = s0 + s[1]*(s1-s0);

                // fixes normal-vector roundoff problem when endpoint is cone-tip
                double delta = s[1]-s[0];
                s[0] += 1e-10*delta;
                s[1] -= 1e-10*delta;

                for(auto i : Range(2))
                {
                    // take uv from CurveOnSurface as start value but project again for better accuracy
                    // (cof->Value yields wrong values (outside of surface) for complicated faces
                    auto uv = cof->Value(s[i]);
                    PointGeomInfo gi;
                    gi.u = uv.X();
                    gi.v = uv.Y();
                    Point<3> pproject = mesh[seg[i]];
                    ProjectPointGI(pproject, gi);
                    // points off the surface (large vertex tolerance) may project outside the face
                    if(gi.u < umin-du || gi.u > umax+du || gi.v < vmin-dv || gi.v > vmax+dv)
                      {
                        gi.u = uv.X();
                        gi.v = uv.Y();
                      }
                    seg.GeomInfo(i).u = gi.u;
                    seg.GeomInfo(i).v = gi.v;
                }

                bool do_swap = ORIENTATION == REVERSED;
                if(seg.EPGeomInfo(1).dist < seg.EPGeomInfo(0).dist)
                  do_swap = !do_swap;

                if(do_swap)
                {
                    swap(seg[0], seg[1]);
                    swap(seg.EPGeomInfo(0).dist, seg.EPGeomInfo(1).dist);
                    swap(seg.GeomInfo(0).u, seg.GeomInfo(1).u);
                    swap(seg.GeomInfo(0).v, seg.GeomInfo(1).v);
                }

                boundary.Append(seg);
            }
        }

        for (auto seg : mesh.LineSegments())
        {
            const auto & ed = mesh.GetEdgeDescriptor(seg.GetIndex());
            auto edgenr = ed.EdgeNr() - 1;
            if(edgenr >= 0 && edgenr < n_edges)
                continue;

            bool forward = ed.SurfNr(0) == nr+1;
            bool reversed = ed.SurfNr(1) == nr+1;
            if(forward == reversed)
                continue;  // not adjacent to this face, or interior to it

            if(reversed)
            {
                swap(seg[0], seg[1]);
                swap(seg.EPGeomInfo(0), seg.EPGeomInfo(1));
            }
            for(auto i : Range(2))
            {
                Point<3> p = mesh[seg[i]];
                seg.GeomInfo(i) = Project(p);
            }
            boundary.Append(seg);
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
      /*
        static Timer t("OCCFace::ProjectPointGI");
        RegionTimer rt(t);
        // *testout << "input, uv = " << gi.u << ", " << gi.v << endl;
        auto suval = shape_analysis->NextValueOfUV({gi.u, gi.v}, ng2occ(p_), tolerance);
        gi.trignum = nr+1;
        suval.Coord(gi.u, gi.v);
        // *testout << "result, uv = " << gi.u << ", " << gi.v << endl;
        p_ = occ2ng(surface->Value( gi.u, gi.v ));        
        return true;
      */
        // Old code: do newton iterations manually
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
        bool out = false;
      
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

           // a non-uniform parametrisation can make a step leave the surface
           // domain, where OCC extrapolates: clamp once, give up if it wants
           // out again
           double uc = uperiod ? u : std::clamp (u, sumin, sumax);
           double vc = vperiod ? v : std::clamp (v, svmin, svmax);
           if (!std::isfinite (u) || !std::isfinite (v)) return false;
           if (uc != u || vc != v)
             {
               if (out) return false;
               out = true;
               u = uc; v = vc;
             }
           else
             out = false;

           xold = x;
           surface->D1(u,v,x,du,dv);

        } while (xold.SquareDistance(x) > sqr(PROJECTION_TOLERANCE) && count < 50);

        //    (*testout) << "FastProject count: " << count << endl;

        if (count == 50) return false;

        p_ = occ2ng(x);
        gi.u = u; gi.v = v;

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
        BRepAdaptor_Surface sf(face, Standard_True);
        BRepLProp_SLProps prop2(sf, 2, 1e-5);
        prop2.SetParameters (gi.u, gi.v);
        return max(fabs(prop2.MinCurvature()),
                   fabs(prop2.MaxCurvature()));
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

#ifndef FILE_OCC_EDGE_INCLUDED
#define FILE_OCC_EDGE_INCLUDED


// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <Geom_Curve.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>

// #pragma clang diagnostic pop

#include "occ_vertex.hpp"
#include "meshing.hpp"

namespace netgen
{
    class OCCEdge : public GeometryEdge 
    {
        public:
        TopoDS_Edge edge;
        Handle(Geom_Curve) curve;
        double s0, s1;
        GProp_GProps props;

        public:
        OCCEdge(TopoDS_Shape edge_, GeometryVertex & start_, GeometryVertex & end_);

        auto Shape() const { return edge; }

        double GetLength() const override;
        Point<3> GetCenter() const override;
        Point<3> GetPoint(double t) const override;
        double CalcStep(double t, double sag) const override;
        void ProjectPoint(Point<3>& p, EdgePointGeomInfo* gi) const override;
        Vec<3> GetTangent(double t) const override;
        bool IsDegenerated(double) const override {
          return BRep_Tool::Degenerated(edge);
        }
    };
}

#endif // FILE_OCCEDGE_INCLUDED

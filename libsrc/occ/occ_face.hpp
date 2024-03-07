#ifndef FILE_OCC_FACE_INCLUDED
#define FILE_OCC_FACE_INCLUDED

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <ShapeAnalysis_Surface.hxx>

#pragma clang diagnostic pop

#include "occ_vertex.hpp"
#include "meshing.hpp"

namespace netgen
{
    class OCCFace : public GeometryFace 
    {
        TopoDS_Face face;
        GProp_GProps props;
        Box<3> bbox;

        Handle( Geom_Surface ) surface;
        Handle( ShapeAnalysis_Surface ) shape_analysis;
        double tolerance;

        public:
        OCCFace(TopoDS_Shape dshape);

        const TopoDS_Face Shape() const { return face; }

        Point<3> GetCenter() const override;
        virtual size_t GetNBoundaries() const override;
        virtual Array<Segment> GetBoundary(const Mesh& mesh) const override;
        virtual PointGeomInfo Project(Point<3>& p) const override;
        virtual bool ProjectPointGI(Point<3>& p, PointGeomInfo& gi) const override;
        virtual Point<3> GetPoint(const PointGeomInfo& gi) const override;
        virtual void CalcEdgePointGI(const GeometryEdge& edge,
                double t,
                EdgePointGeomInfo& egi) const override;
        virtual Box<3> GetBoundingBox() const override;

        virtual double GetCurvature(const PointGeomInfo& gi) const override;

        virtual void RestrictH(Mesh& mesh, const MeshingParameters& mparam) const override;
        virtual Vec<3> GetNormal(const Point<3>& p, const PointGeomInfo* gi = nullptr) const override;
    };
}

#endif // FILE_OCC_FACE_INCLUDED

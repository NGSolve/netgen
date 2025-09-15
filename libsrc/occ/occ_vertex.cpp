#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>

#include "occ_vertex.hpp"

namespace netgen
{

    OCCVertex::OCCVertex( TopoDS_Shape s )
        : vertex(TopoDS::Vertex(s))
    {
        p = occ2ng(vertex);
    }

    Point<3> OCCVertex::GetPoint() const
    {
        return p;
    }

    void OCCVertex::SetPoint(Point<3> newp) {
      p = newp;
    }
}

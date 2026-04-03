#ifndef FILE_OCC_VERTEX_INCLUDED
#define FILE_OCC_VERTEX_INCLUDED

// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <TopoDS.hxx>
#include <BRep_TVertex.hxx>

// #pragma clang diagnostic pop

#include "meshing.hpp"
#include "occ_utils.hpp"

namespace netgen
{
    class OCCVertex : public GeometryVertex
    {
        TopoDS_Vertex vertex;
        Point<3> p;

        public:
        OCCVertex( ) = default;
        OCCVertex( TopoDS_Shape s );
        ~OCCVertex() {}
        Point<3> GetPoint() const override;
    };
}

#endif // FILE_OCC_VERTEX_INCLUDED

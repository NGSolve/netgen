#ifndef FILE_OCC_VERTEX_INCLUDED
#define FILE_OCC_VERTEX_INCLUDED

#include <TopoDS.hxx>
#include <BRep_TVertex.hxx>

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
        size_t GetHash() const override;
    };
}

#endif // FILE_OCC_VERTEX_INCLUDED

#ifndef FILE_OCC_SOLID_INCLUDED
#define FILE_OCC_SOLID_INCLUDED

#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>

#include "meshing.hpp"

namespace netgen
{
    class OCCSolid : public GeometrySolid 
    {
        TopoDS_Solid solid;

        public:
        OCCSolid(TopoDS_Shape dshape)
            : solid(TopoDS::Solid(dshape))
        { }

        size_t GetHash() const override { return solid.HashCode(std::numeric_limits<Standard_Integer>::max()); }
    };
}

#endif // FILE_OCC_SOLID_INCLUDED

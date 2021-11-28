#ifndef FILE_OCC_SOLID_INCLUDED
#define FILE_OCC_SOLID_INCLUDED

#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>

#include "meshing.hpp"

namespace netgen
{
    class OCCSolid : public GeometrySolid 
    {
        T_Shape tsolid;
        TopoDS_Solid solid;

        public:
        OCCSolid(TopoDS_Shape dshape)
            : tsolid(dshape.TShape()),
              solid(TopoDS::Solid(dshape))
        { }

        T_Shape TShape() { return tsolid; }
        size_t GetHash() const override { return reinterpret_cast<size_t>(tsolid.get()); }
    };
}

#endif // FILE_OCC_SOLID_INCLUDED

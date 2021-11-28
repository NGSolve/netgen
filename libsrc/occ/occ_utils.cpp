#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRep_TVertex.hxx>

#include "occ_utils.hpp"

namespace netgen
{
    Point<3> occ2ng (Handle(TopoDS_TShape) shape)
    {
        return occ2ng( Handle(BRep_TVertex)::DownCast(shape)->Pnt() );
    }

    Box<3> GetBoundingBox( const TopoDS_Shape & shape )
    {
        Bnd_Box bb;
#if OCC_VERSION_HEX < 0x070000
        BRepBndLib::Add (shape, bb);
#else
        BRepBndLib::Add (shape, bb, true);
#endif
        return {occ2ng(bb.CornerMin()), occ2ng(bb.CornerMax())};
    }
}

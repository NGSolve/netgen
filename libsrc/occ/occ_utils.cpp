#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRep_TVertex.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>

#include "occ_utils.hpp"
#include "occgeom.hpp"

namespace netgen
{
    Point<3> occ2ng (const TopoDS_Shape& shape)
    {
      if(shape.ShapeType() != TopAbs_VERTEX)
        throw Exception("Try to convert non vertex to point!");
      return occ2ng( BRep_Tool::Pnt(TopoDS::Vertex(shape)) );
    }

    Transformation<3> occ2ng (const gp_Trsf & occ_trafo)
    {
        Transformation<3> trafo;
        auto v = occ_trafo.TranslationPart();
        auto m = occ_trafo.VectorialPart();
        auto & tv = trafo.GetVector();
        auto & tm = trafo.GetMatrix();
        for(auto i : Range(3))
        {
            tv[i] = v.Coord(i+1);
            for(auto k : Range(3))
                tm(i,k) = m(i+1,k+1);
        }
        return trafo;
    }

  Transformation<3> occ2ng (const gp_GTrsf & occ_trafo)
    {
        Transformation<3> trafo;
        auto v = occ_trafo.TranslationPart();
        auto m = occ_trafo.VectorialPart();
        auto & tv = trafo.GetVector();
        auto & tm = trafo.GetMatrix();
        for(auto i : Range(3))
        {
            tv[i] = v.Coord(i+1);
            for(auto k : Range(3))
                tm(i,k) = m(i+1,k+1);
        }
        return trafo;
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

    Standard_Integer BuildTriangulation( const TopoDS_Shape & shape )
    {
       BRepTools::Clean (shape);
       // double deflection = 0.01;

       // https://dev.opencascade.org/doc/overview/html/occt_user_guides__mesh.html
       // from Standard_Boolean meshing_imeshtools_parameters()
       IMeshTools_Parameters aMeshParams;
       aMeshParams.Deflection               = 0.01;
       aMeshParams.Angle                    = 0.5;
       aMeshParams.Relative                 = Standard_True;
       aMeshParams.InParallel               = Standard_True;
       aMeshParams.MinSize                  = Precision::Confusion();
       aMeshParams.InternalVerticesMode     = Standard_True;
       aMeshParams.ControlSurfaceDeflection = Standard_True;

       BRepMesh_IncrementalMesh aMesher (shape, aMeshParams);
       return  aMesher.GetStatusFlags();
    }
}

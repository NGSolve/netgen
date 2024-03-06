
#ifdef OCCGEOMETRY

#include <cstdio>
#include <set>

#include <mystdlib.h>
#include <core/register_archive.hpp>

#include "occ_vertex.hpp"
#include "occ_edge.hpp"
#include "occ_face.hpp"
#include "occ_solid.hpp"
#include "occgeom.hpp"
#include "Partition_Spliter.hxx"

#include <BOPAlgo_Builder.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepOffsetAPI_Sewing.hxx>
#include <BRepTools.hxx>
#include <IGESCAFControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <Interface_InterfaceModel.hxx>
#include <Interface_Static.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPConstruct.hxx>
#include <STEPControl_Writer.hxx>
#include <ShapeAnalysis_CheckSmallFace.hxx>
#include <ShapeAnalysis_DataMapOfShapeListOfReal.hxx>
#include <ShapeAnalysis_ShapeContents.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_FixSmallFace.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <StepBasic_ProductDefinitionRelationship.hxx>
#include <StepRepr_RepresentationItem.hxx>
#include <StlAPI_Writer.hxx>
#include <TopoDS_Shape.hxx>
#include <TransferBRep.hxx>
#include <Transfer_FinderProcess.hxx>
#include <Transfer_TransientProcess.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFPrs.hxx>
#include <XCAFPrs_Style.hxx>
#include <XSControl_TransferReader.hxx>
#include <XSControl_WorkSession.hxx>

#if OCC_VERSION_HEX < 0x070000
// pass
#elif OCC_VERSION_HEX < 0x070200
   #include <StlTransfer.hxx>
   #include <TopoDS_Iterator.hxx>
#else
   #include <TopoDS_Iterator.hxx>
#endif

namespace netgen
{
  void LoadOCCInto(OCCGeometry* occgeo, const filesystem::path & filename);
  void PrintContents (OCCGeometry * geom);

  // Utility function to apply builder and propagate properties
  template <typename T>
  static TopoDS_Shape Apply(T & builder, TopoDS_Shape & shape) {
    auto newshape = builder->Apply(shape);
    PropagateProperties(*builder, newshape);
    return newshape;
  };


  TopTools_IndexedMapOfShape OCCGeometry::global_shape_property_indices;
  std::vector<ShapeProperties> OCCGeometry::global_shape_properties;
  TopTools_IndexedMapOfShape OCCGeometry::global_identification_indices;
  std::vector<std::vector<OCCIdentification>> OCCGeometry::global_identifications;

  TopoDS_Shape ListOfShapes::Max(gp_Vec dir)
  {
    double maxval = -1e99;
    TopoDS_Shape maxshape;
    for (auto shape : *this)
      {
        GProp_GProps props;
        gp_Pnt center;

        switch (shape.ShapeType())
          {
          case TopAbs_VERTEX:
            center = BRep_Tool::Pnt (TopoDS::Vertex(shape)); break;
          case TopAbs_FACE:
            BRepGProp::SurfaceProperties (shape, props);
            center = props.CentreOfMass();
            break;
          default:
            BRepGProp::LinearProperties(shape, props);
            center = props.CentreOfMass();
          }

        double val = center.X()*dir.X() + center.Y()*dir.Y() + center.Z() * dir.Z();
        if (val > maxval)
          {
            maxval = val;
            maxshape = shape;
          }
      }
    return maxshape;
  }
  TopoDS_Shape ListOfShapes::Nearest(gp_Pnt pnt)
  {
    double mindist = 1e99;
    TopoDS_Shape nearestshape;
    auto vertex = BRepBuilderAPI_MakeVertex (pnt).Vertex();
    
    for (auto shape : *this)
      {
        double dist = BRepExtrema_DistShapeShape(shape, vertex).Value();
        if (dist < mindist)
          {
            nearestshape = shape;
            mindist = dist;
          }
      }
    return nearestshape;
  }

  ListOfShapes ListOfShapes::SubShapes(TopAbs_ShapeEnum type) const
  {
    TopTools_MapOfShape check_unique;
    ListOfShapes sub;
    for(const auto& shape : *this)
      for(TopExp_Explorer e(shape, type); e.More(); e.Next())
        if(const auto& s = e.Current(); !check_unique.Contains(s))
          {
            check_unique.Add(s);
            sub.push_back(s);
          }
    return sub;
  }

  OCCGeometry::OCCGeometry(const TopoDS_Shape& _shape, int aoccdim, bool copy)
  {
    if(copy)
      {
        auto filename = GetTempFilename();
        step_utils::WriteSTEP(_shape, filename);
        LoadOCCInto(this, filename);
        dimension = aoccdim;
        filesystem::remove(filename);
      }
    else
      {
        shape = _shape;
        changed = 1;
        dimension = aoccdim;
        BuildFMap();
        CalcBoundingBox();
        PrintContents (this);
      }
  }

  const GeometryShape & OCCGeometry :: GetShape(const TopoDS_Shape & shape) const
  {
    switch (shape.ShapeType())
    {
      case TopAbs_VERTEX:
        return GetVertex(shape);
      case TopAbs_EDGE:
        return GetEdge(shape);
      case TopAbs_FACE:
        return GetFace(shape);
      default:
        throw Exception("unknown shape type");
    }
  }

  const GeometryVertex & OCCGeometry :: GetVertex(const TopoDS_Shape & shape) const
  {
      return *vertices[vmap.FindIndex(shape)-1];
  }

  const GeometryEdge & OCCGeometry :: GetEdge(const TopoDS_Shape & shape) const
  {
      return *edges[emap.FindIndex(shape)-1];
  }

  const GeometryFace & OCCGeometry :: GetFace(const TopoDS_Shape & shape) const
  {
      return *faces[fmap.FindIndex(shape)-1];
  }


  string STEP_GetEntityName(const TopoDS_Shape & theShape, STEPCAFControl_Reader * aReader)
  {
    const Handle(XSControl_WorkSession)& theSession = aReader->Reader().WS();
    const Handle(XSControl_TransferReader)& aTransferReader =
        theSession->TransferReader();

    Handle(Standard_Transient) anEntity =
        aTransferReader->EntityFromShapeResult(theShape, 1);

    if (anEntity.IsNull()) // as just mapped
        anEntity = aTransferReader->EntityFromShapeResult (theShape,-1);

    if (anEntity.IsNull()) // as anything
        anEntity = aTransferReader->EntityFromShapeResult (theShape,4);

    if (anEntity.IsNull())
      {
        cout<<"Warning: cannot get entity from shape" <<endl;
        return "none";
      }

    auto aReprItem = Handle(StepRepr_RepresentationItem)::DownCast(anEntity);
    if(!aReprItem.IsNull())
        return aReprItem->Name()->ToCString();;

    auto bReprItem = Handle(StepBasic_ProductDefinitionRelationship)::DownCast(anEntity);
    if (!bReprItem.IsNull())
        return bReprItem->Description()->ToCString();

    cout<<"Warning: unknown entity type " << anEntity->DynamicType() << endl;
    return "none";
  }

  void OCCGeometry :: Analyse(Mesh& mesh,
                              const MeshingParameters& mparam) const
  {
    OCCSetLocalMeshSize(*this, mesh, mparam, occparam);
  }

  bool OCCGeometry :: MeshFace(Mesh& mesh,
                                  const MeshingParameters& mparam, int nr, FlatArray<int, PointIndex> glob2loc) const
  {
    MeshingParameters local_mp = mparam;
    auto face = TopoDS::Face(fmap(nr+1));
    if(auto quad_dominated = OCCGeometry::GetProperties(face).quad_dominated; quad_dominated.has_value())
      local_mp.quad = *quad_dominated;

    bool failed = OCCMeshFace(*this, mesh, glob2loc, local_mp, nr, PARAMETERSPACE, true);
    if(failed)
        failed = OCCMeshFace(*this, mesh, glob2loc, local_mp, nr, PLANESPACE, false);

    if(failed)
    {
        facemeshstatus[nr] = -1;
        PrintError ("Problem in Surface mesh generation");
    }
    else
    {
        facemeshstatus[nr] = 1;
    }
    return failed;
  }

   void OCCGeometry :: PrintNrShapes ()
   {
      TopExp_Explorer e;
      int count = 0;
      for (e.Init(shape, TopAbs_COMPSOLID); e.More(); e.Next()) count++;
      cout << "CompSolids: " << count << endl;

      cout << "Solids    : " << somap.Extent() << endl;
      cout << "Shells    : " << shmap.Extent() << endl;
      cout << "Faces     : " << fmap.Extent() << endl;
      cout << "Edges     : " << emap.Extent() << endl;
      cout << "Vertices  : " << vmap.Extent() << endl;
   }




   void PrintContents (OCCGeometry * geom)
   {
      ShapeAnalysis_ShapeContents cont;
      cont.Clear();
      cont.Perform(geom->shape);

      (*testout) << "OCC CONTENTS" << endl;
      (*testout) << "============" << endl;
      (*testout) << "SOLIDS   : " << cont.NbSolids() << endl;
      (*testout) << "SHELLS   : " << cont.NbShells() << endl;
      (*testout) << "FACES    : " << cont.NbFaces() << endl;
      (*testout) << "WIRES    : " << cont.NbWires() << endl;
      (*testout) << "EDGES    : " << cont.NbEdges() << endl;
      (*testout) << "VERTICES : " << cont.NbVertices() << endl;

      TopExp_Explorer e;
      int count = 0;
      for (e.Init(geom->shape, TopAbs_COMPOUND); e.More(); e.Next())
         count++;
      (*testout) << "Compounds: " << count << endl;

      count = 0;
      for (e.Init(geom->shape, TopAbs_COMPSOLID); e.More(); e.Next())
         count++;
      (*testout) << "CompSolids: " << count << endl;

      (*testout) << endl;

      cout << IM(3) << "Highest entry in topology hierarchy: " << endl;
      if (count)
         cout << IM(3) << count << " composite solid(s)" << endl;
      else
         if (geom->somap.Extent())
            cout << IM(3) << geom->somap.Extent() << " solid(s)" << endl;
         else
            if (geom->shmap.Extent())
               cout << IM(3) << geom->shmap.Extent() << " shells(s)" << endl;
            else
               if (geom->fmap.Extent())
                  cout << IM(3) << geom->fmap.Extent() << " face(s)" << endl;
               else
                  if (geom->wmap.Extent())
                     cout << IM(3) << geom->wmap.Extent() << " wire(s)" << endl;
                  else
                     if (geom->emap.Extent())
                        cout << IM(3) << geom->emap.Extent() << " edge(s)" << endl;
                     else
                        if (geom->vmap.Extent())
                           cout << IM(3) << geom->vmap.Extent() << " vertices(s)" << endl;
                        else
                           cout << IM(3) << "no entities" << endl;

   }

  void OCCGeometry :: GlueGeometry()
  {
    PrintMessage(1, "OCC Glue Geometry");
    /*
      // 
    BRep_Builder builder;
    TopoDS_Shape my_fuse;
    int cnt = 0;
    for (TopExp_Explorer exp_solid(shape, TopAbs_SOLID); exp_solid.More(); exp_solid.Next())
      {
        cout << "cnt = " << cnt << endl;
	if (cnt == 0)
	  my_fuse = exp_solid.Current();
	else
          // my_fuse = BRepAlgoAPI_Fuse (my_fuse, exp_solid.Current());
          my_fuse = QANewModTopOpe_Glue::QANewModTopOpe_Glue(my_fuse, exp_solid.Current());
	cnt++;
      }
    cout << "remove" << endl;
    // for (int i = 1; i <= somap.Size(); i++)
    // builder.Remove (shape, somap(i));
    cout << "now add" << endl;
    // builder.Add (shape, my_fuse);
    shape = my_fuse;
    cout << "build fmap" << endl;
    BuildFMap();
    */


    // from 
    // https://www.opencascade.com/doc/occt-7.4.0/overview/html/occt_user_guides__boolean_operations.html
    BOPAlgo_Builder aBuilder;

    // Setting arguments
    TopTools_ListOfShape aLSObjects; 
    for (TopExp_Explorer exp_solid(shape, TopAbs_SOLID); exp_solid.More(); exp_solid.Next())
      aLSObjects.Append (exp_solid.Current());
    aBuilder.SetArguments(aLSObjects);

    // Setting options for GF
    // Set parallel processing mode (default is false)
    // Standard_Boolean bRunParallel = Standard_True;
    // aBuilder.SetRunParallel(bRunParallel);
    
    // Set Fuzzy value (default is Precision::Confusion())
    // Standard_Real aFuzzyValue = 1.e-5;
    // aBuilder.SetFuzzyValue(aFuzzyValue);
    
    // Set safe processing mode (default is false)
    // Standard_Boolean bSafeMode = Standard_True;
    // aBuilder.SetNonDestructive(bSafeMode);
    
    // Set Gluing mode for coinciding arguments (default is off)
    // BOPAlgo_GlueEnum aGlue = BOPAlgo_GlueShift;
    // aBuilder.SetGlue(aGlue);
    
    // Disabling/Enabling the check for inverted solids (default is true)
    // Standard Boolean bCheckInverted = Standard_False;
    // aBuilder.SetCheckInverted(bCheckInverted);
    
    // Set OBB usage (default is false)
    // Standard_Boolean bUseOBB = Standard_True;
    // aBuilder.SetUseOBB(buseobb);
    
    // Perform the operation
    aBuilder.Perform();
    // Check for the errors
#if OCC_VERSION_HEX >= 0x070200
    if (aBuilder.HasErrors())
      {
        cout << "builder has errors" << endl;
        return;
      }
    // Check for the warnings
    if (aBuilder.HasWarnings())
      {
        // treatment of the warnings
        ;
      }
#endif

#ifdef OCC_HAVE_HISTORY    
    Handle(BRepTools_History) history = aBuilder.History ();
    
    for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next())
      {
        if (auto name = OCCGeometry::GetProperties(e.Current()).name)
          for (auto mods : history->Modified(e.Current()))
            OCCGeometry::GetProperties(mods).name = *name;
      }
#endif // OCC_HAVE_HISTORY    
    
    // result of the operation
    shape = aBuilder.Shape();
    BuildFMap();
  }

   void OCCGeometry :: HealGeometry ()
   {
      int nrc = 0, nrcs = 0,
         nrso = somap.Extent(),
         nrsh = shmap.Extent(),
         nrf = fmap.Extent(),
         nrw = wmap.Extent(),
         nre = emap.Extent(),
         nrv = vmap.Extent();

      TopExp_Explorer exp0;
      TopExp_Explorer exp1;


      for (exp0.Init(shape, TopAbs_COMPOUND); exp0.More(); exp0.Next()) nrc++;
      for (exp0.Init(shape, TopAbs_COMPSOLID); exp0.More(); exp0.Next()) nrcs++;

      double surfacecont = 0;
      
      {
         Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
         Apply(rebuild, shape);
         for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
         {
            TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
            if ( BRep_Tool::Degenerated(edge) )
               rebuild->Remove(edge);
         }
         shape = Apply(rebuild, shape);
      }

      BuildFMap();


      for (exp0.Init (shape, TopAbs_FACE); exp0.More(); exp0.Next())
      {
         TopoDS_Face face = TopoDS::Face(exp0.Current());

         GProp_GProps system;
         BRepGProp::SurfaceProperties(face, system);
         surfacecont += system.Mass();
      }


      cout << "Starting geometry healing procedure (tolerance: " << tolerance << ")" << endl
         << "-----------------------------------" << endl;

      {
         cout << endl << "- repairing faces" << endl;

         Handle(ShapeFix_Face) sff;
         Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
         Apply(rebuild, shape);

         for (exp0.Init (shape, TopAbs_FACE); exp0.More(); exp0.Next())
         {
            TopoDS_Face face = TopoDS::Face (exp0.Current());
            auto props = GetProperties(face);

            sff = new ShapeFix_Face (face);
            sff->FixAddNaturalBoundMode() = Standard_True;
            sff->FixSmallAreaWireMode() = Standard_True;
            sff->Perform();

            if(sff->Status(ShapeExtend_DONE1) ||
               sff->Status(ShapeExtend_DONE2) ||
               sff->Status(ShapeExtend_DONE3) ||
               sff->Status(ShapeExtend_DONE4) ||
               sff->Status(ShapeExtend_DONE5))
            {
               cout << "repaired face " << fmap.FindIndex(face) << " ";
               if(sff->Status(ShapeExtend_DONE1))
                  cout << "(some wires are fixed)" <<endl;
               else if(sff->Status(ShapeExtend_DONE2))
                  cout << "(orientation of wires fixed)" <<endl;
               else if(sff->Status(ShapeExtend_DONE3))
                  cout << "(missing seam added)" <<endl;
               else if(sff->Status(ShapeExtend_DONE4))
                  cout << "(small area wire removed)" <<endl;
               else if(sff->Status(ShapeExtend_DONE5))
                  cout << "(natural bounds added)" <<endl;
               TopoDS_Face newface = sff->Face();

               rebuild->Replace(face, newface);
            }

            // Set the original properties of the face to the newly created 
            // face (after the healing process)
            // GetProperties(face);
         }
         shape = Apply(rebuild, shape);
      }


      {
         Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
         Apply(rebuild, shape);
         for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
         {
            TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
            if ( BRep_Tool::Degenerated(edge) )
               rebuild->Remove(edge);
         }
         shape = Apply(rebuild, shape);
      }


      if (fixsmalledges)
      {
         cout << endl << "- fixing small edges" << endl;

         Handle(ShapeFix_Wire) sfw;
         Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
         Apply(rebuild, shape);


         for (exp0.Init (shape, TopAbs_FACE); exp0.More(); exp0.Next())
         {
            TopoDS_Face face = TopoDS::Face(exp0.Current());

            for (exp1.Init (face, TopAbs_WIRE); exp1.More(); exp1.Next())
            {
               TopoDS_Wire oldwire = TopoDS::Wire(exp1.Current());
               sfw = new ShapeFix_Wire (oldwire, face ,tolerance);
               sfw->ModifyTopologyMode() = Standard_True;

               sfw->ClosedWireMode() = Standard_True;

               bool replace = false;

               replace = sfw->FixReorder() || replace;

               replace = sfw->FixConnected() || replace;



               if (sfw->FixSmall (Standard_False, tolerance) && ! (sfw->StatusSmall(ShapeExtend_FAIL1) ||
                  sfw->StatusSmall(ShapeExtend_FAIL2) ||
                  sfw->StatusSmall(ShapeExtend_FAIL3)))
               {
                  cout << "Fixed small edge in wire " << wmap.FindIndex (oldwire) << endl;
                  replace = true;

               }
               else if (sfw->StatusSmall(ShapeExtend_FAIL1))
                  cerr << "Failed to fix small edge in wire " << wmap.FindIndex (oldwire)
                  << ", edge cannot be checked (no 3d curve and no pcurve)" << endl;
               else if (sfw->StatusSmall(ShapeExtend_FAIL2))
                  cerr << "Failed to fix small edge in wire " << wmap.FindIndex (oldwire)
                  << ", edge is null-length and has different vertives at begin and end, and lockvtx is True or ModifiyTopologyMode is False" << endl;
               else if (sfw->StatusSmall(ShapeExtend_FAIL3))
                  cerr << "Failed to fix small edge in wire " << wmap.FindIndex (oldwire)
                  << ", CheckConnected has failed" << endl;

               replace = sfw->FixEdgeCurves() || replace;

               replace = sfw->FixDegenerated() || replace;

               replace = sfw->FixSelfIntersection() || replace;

               replace = sfw->FixLacking(Standard_True) || replace;

               if(replace)
               {
                  TopoDS_Wire newwire = sfw->Wire();
                  rebuild->Replace(oldwire, newwire);
               }

               //delete sfw; sfw = NULL;

            }
         }

         shape = Apply(rebuild, shape);



         {
            BuildFMap();
            Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
            Apply(rebuild, shape);

            for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
            {
               TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
               if (vmap.FindIndex(TopExp::FirstVertex (edge)) ==
                  vmap.FindIndex(TopExp::LastVertex (edge)))
               {
                  GProp_GProps system;
                  BRepGProp::LinearProperties(edge, system);
                  if (system.Mass() < tolerance)
                  {
                     cout << "removing degenerated edge " << emap.FindIndex(edge)
                        << " from vertex " << vmap.FindIndex(TopExp::FirstVertex (edge))
                        << " to vertex " << vmap.FindIndex(TopExp::LastVertex (edge)) << endl;
                     rebuild->Remove(edge);
                  }
               }
            }
            shape = Apply(rebuild, shape);

            //delete rebuild; rebuild = NULL;
         }



         {
            Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
            Apply(rebuild, shape);
            for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
            {
               TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
               if ( BRep_Tool::Degenerated(edge) )
                  rebuild->Remove(edge);
            }
            shape = Apply(rebuild, shape);
         }




         Handle(ShapeFix_Wireframe) sfwf = new ShapeFix_Wireframe;
         sfwf->SetPrecision(tolerance);
         sfwf->Load (shape);
         sfwf->ModeDropSmallEdges() = Standard_True;

         sfwf->SetPrecision(boundingbox.Diam());

         if (sfwf->FixWireGaps())
         {
            cout << endl << "- fixing wire gaps" << endl;
            if (sfwf->StatusWireGaps(ShapeExtend_OK)) cout << "no gaps found" << endl;
            if (sfwf->StatusWireGaps(ShapeExtend_DONE1)) cout << "some 2D gaps fixed" << endl;
            if (sfwf->StatusWireGaps(ShapeExtend_DONE2)) cout << "some 3D gaps fixed" << endl;
            if (sfwf->StatusWireGaps(ShapeExtend_FAIL1)) cout << "failed to fix some 2D gaps" << endl;
            if (sfwf->StatusWireGaps(ShapeExtend_FAIL2)) cout << "failed to fix some 3D gaps" << endl;
         }

         sfwf->SetPrecision(tolerance);


         {
            for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
            {
               TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
               if ( BRep_Tool::Degenerated(edge) )
                  cout << "degenerated edge at position 4" << endl;
            }
         }



         if (sfwf->FixSmallEdges())
         {
            cout << endl << "- fixing wire frames" << endl;
            if (sfwf->StatusSmallEdges(ShapeExtend_OK)) cout << "no small edges found" << endl;
            if (sfwf->StatusSmallEdges(ShapeExtend_DONE1)) cout << "some small edges fixed" << endl;
            if (sfwf->StatusSmallEdges(ShapeExtend_FAIL1)) cout << "failed to fix some small edges" << endl;
         }



         auto newshape = sfwf->Shape();
         PropagateProperties(*sfwf->Context(), newshape);
         shape = newshape;

         //delete sfwf; sfwf = NULL;
         //delete rebuild; rebuild = NULL;

      }





      {
         for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
         {
            TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
            if ( BRep_Tool::Degenerated(edge) )
               cout << "degenerated edge at position 5" << endl;
         }
      }




      if (fixspotstripfaces)
      {

         cout << endl << "- fixing spot and strip faces" << endl;
         Handle(ShapeFix_FixSmallFace) sffsm = new ShapeFix_FixSmallFace();
         sffsm -> Init (shape);
         sffsm -> SetPrecision (tolerance);
         sffsm -> Perform();

         auto newshape = sffsm -> FixShape();
         PropagateProperties(*sffsm->Context(), newshape);
         shape = newshape;
         //delete sffsm; sffsm = NULL;
      }


      {
         for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
         {
            TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
            if ( BRep_Tool::Degenerated(edge) )
               cout << "degenerated edge at position 6" << endl;
         }
      }



      if (sewfaces)
      {
         cout << endl << "- sewing faces" << endl;

         BRepOffsetAPI_Sewing sewedObj(tolerance);

         for (exp0.Init (shape, TopAbs_FACE); exp0.More(); exp0.Next())
         {
            TopoDS_Face face = TopoDS::Face (exp0.Current());
            sewedObj.Add (face);
         }

         sewedObj.Perform();
         PropagateProperties(sewedObj, shape);

         if (!sewedObj.SewedShape().IsNull())
            shape = sewedObj.SewedShape();
         else
            cout << " not possible";
      }



      {
         Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
         rebuild->Apply(shape);
         for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
         {
            TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
            if ( BRep_Tool::Degenerated(edge) )
               rebuild->Remove(edge);
         }
         shape = Apply(rebuild, shape);
      }


      if (makesolids)
      {
         cout << endl << "- making solids" << endl;

         BRepBuilderAPI_MakeSolid ms;
         int count = 0;
         for (exp0.Init(shape, TopAbs_SHELL); exp0.More(); exp0.Next())
         {
            count++;
            ms.Add (TopoDS::Shell(exp0.Current()));
         }

         if (!count)
         {
            cout << " not possible (no shells)" << endl;
         }
         else
         {
            BRepCheck_Analyzer ba(ms);
            if (ba.IsValid ())
            {
               Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
               sfs->Init (ms);
               sfs->SetPrecision(tolerance);
               sfs->SetMaxTolerance(tolerance);
               sfs->Perform();
               shape = sfs->Shape();

               for (exp0.Init(shape, TopAbs_SOLID); exp0.More(); exp0.Next())
               {
                  TopoDS_Solid solid = TopoDS::Solid(exp0.Current());
                  TopoDS_Solid newsolid = solid;
                  BRepLib::OrientClosedSolid (newsolid);
                  Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
                  //		  rebuild->Apply(shape);
                  rebuild->Replace(solid, newsolid);
                  TopoDS_Shape newshape = rebuild->Apply(shape, TopAbs_COMPSOLID);//, 1);
                  //		  TopoDS_Shape newshape = rebuild->Apply(shape);
                  shape = newshape;
               }

               //delete sfs; sfs = NULL;
            }
            else
               cout << " not possible" << endl;
         }
      }



      if (splitpartitions)
      {
         cout << "- running SALOME partition splitter" << endl;

         TopExp_Explorer e2;
         Partition_Spliter ps;
         int count = 0;

         for (e2.Init (shape, TopAbs_SOLID);
            e2.More(); e2.Next())
         {
            count++;
            ps.AddShape (e2.Current());
         }

         ps.Compute();
         shape = ps.Shape();

         cout << " before: " << count << " solids" << endl;

         count = 0;
         for (e2.Init (shape, TopAbs_SOLID);
            e2.More(); e2.Next()) count++;

            cout << " after : " << count << " solids" << endl;
      }

      BuildFMap();



      {
         for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
         {
            TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
            if ( BRep_Tool::Degenerated(edge) )
               cout << "degenerated edge at position 8" << endl;
         }
      }


      double newsurfacecont = 0;


      for (exp0.Init (shape, TopAbs_FACE); exp0.More(); exp0.Next())
      {
         TopoDS_Face face = TopoDS::Face(exp0.Current());
         GProp_GProps system;
         BRepGProp::SurfaceProperties(face, system);
         newsurfacecont += system.Mass();
      }


      int nnrc = 0, nnrcs = 0,
         nnrso = somap.Extent(),
         nnrsh = shmap.Extent(),
         nnrf = fmap.Extent(),
         nnrw = wmap.Extent(),
         nnre = emap.Extent(),
         nnrv = vmap.Extent();

      for (exp0.Init(shape, TopAbs_COMPOUND); exp0.More(); exp0.Next()) nnrc++;
      for (exp0.Init(shape, TopAbs_COMPSOLID); exp0.More(); exp0.Next()) nnrcs++;

      cout << "-----------------------------------" << endl;
      cout << "Compounds       : " << nnrc << " (" << nrc << ")" << endl;
      cout << "Composite solids: " << nnrcs << " (" << nrcs << ")" << endl;
      cout << "Solids          : " << nnrso << " (" << nrso << ")" << endl;
      cout << "Shells          : " << nnrsh << " (" << nrsh << ")" << endl;
      cout << "Wires           : " << nnrw << " (" << nrw << ")" << endl;
      cout << "Faces           : " << nnrf << " (" << nrf << ")" << endl;
      cout << "Edges           : " << nnre << " (" << nre << ")" << endl;
      cout << "Vertices        : " << nnrv << " (" << nrv << ")" << endl;
      cout << endl;
      cout << "Total surface area : " << newsurfacecont << " (" << surfacecont << ")" << endl;
      cout << endl;
   }


   // For 2d geometries, make sure all faces have a normal vector with positive z-component
   void OCCGeometry :: FixFaceOrientation()
   {
     if(dimension!=2) return;

     bool needs_fix = false;
     Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
     for (auto face : GetFaces(shape))
     {
       auto occface = OCCFace(face);
       auto normal = occface.GetNormal(occ2ng(GetVertices(face)[0]));
       if(normal[2] < 0) {
         needs_fix = true;
         // Need do copy the face, otherwise replace is ignored
         BRepBuilderAPI_Copy copy(face);
         auto newface = copy.Shape().Reversed();
         GetProperties(newface).Merge(GetProperties(face));
         rebuild->Replace(face, newface);
       }
     }

     if(needs_fix )
       shape = Apply(rebuild, shape);
   }

   void OCCGeometry :: BuildFMap()
   {
      somap.Clear();
      shmap.Clear();
      fmap.Clear();
      wmap.Clear();
      emap.Clear();
      vmap.Clear();

      TopExp_Explorer exp0, exp1, exp2, exp3, exp4, exp5;

      // Check face orientation in 2d geometries
      FixFaceOrientation();

      for (exp0.Init(shape, TopAbs_COMPOUND);
         exp0.More(); exp0.Next())
      {
         TopoDS_Compound compound = TopoDS::Compound (exp0.Current());
         (*testout) << "compound" << endl;
         int i = 0;
         for (exp1.Init(compound, TopAbs_SHELL);
            exp1.More(); exp1.Next())
         {
            (*testout) << "shell " << ++i << endl;
         }
      }

      for (exp0.Init(shape, TopAbs_SOLID);
         exp0.More(); exp0.Next())
      {
         TopoDS_Solid solid = TopoDS::Solid (exp0.Current());

         if (somap.FindIndex(solid) < 1)
         {
            somap.Add (solid);

            for (exp1.Init(solid, TopAbs_SHELL);
               exp1.More(); exp1.Next())
            {
               TopoDS_Shell shell = TopoDS::Shell (exp1.Current());
               if (shmap.FindIndex(shell) < 1)
               {
                  shmap.Add (shell);

                  for (exp2.Init(shell, TopAbs_FACE);
                     exp2.More(); exp2.Next())
                  {
                     TopoDS_Face face = TopoDS::Face(exp2.Current());
                     if (fmap.FindIndex(face) < 1)
                     {
                        fmap.Add (face);
                        (*testout) << "face " << fmap.FindIndex(face) << " ";
                        (*testout) << ((face.Orientation() == TopAbs_REVERSED) ? "-" : "+") << ", ";
                        (*testout) << ((exp2.Current().Orientation() == TopAbs_REVERSED) ? "-" : "+") << endl;
                        for (exp3.Init(exp2.Current(), TopAbs_WIRE);
                           exp3.More(); exp3.Next())
                        {
                           TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
                           if (wmap.FindIndex(wire) < 1)
                           {
                              wmap.Add (wire);

                              for (exp4.Init(exp3.Current(), TopAbs_EDGE);
                                 exp4.More(); exp4.Next())
                              {
                                 TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
                                 if (emap.FindIndex(edge) < 1)
                                 {
                                    emap.Add (edge);
                                    for (exp5.Init(exp4.Current(), TopAbs_VERTEX);
                                       exp5.More(); exp5.Next())
                                    {
                                       TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                                       if (vmap.FindIndex(vertex) < 1)
                                          vmap.Add (vertex);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }

      // Free Shells
      for (exp1.Init(shape, TopAbs_SHELL, TopAbs_SOLID); exp1.More(); exp1.Next())
      {
         TopoDS_Shell shell = TopoDS::Shell(exp1.Current());
         if (shmap.FindIndex(shell) < 1)
         {
            shmap.Add (shell);

            (*testout) << "shell " << shmap.FindIndex(shell) << " ";
            (*testout) << ((shell.Orientation() == TopAbs_REVERSED) ? "-" : "+") << ", ";
            (*testout) << ((exp1.Current().Orientation() == TopAbs_REVERSED) ? "-" : "+") << endl;

            for (exp2.Init(shell, TopAbs_FACE); exp2.More(); exp2.Next())
            {
               TopoDS_Face face = TopoDS::Face(exp2.Current());
               if (fmap.FindIndex(face) < 1)
               {
                 fmap.Add (face);

                  for (exp3.Init(face, TopAbs_WIRE); exp3.More(); exp3.Next())
                  {
                     TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
                     if (wmap.FindIndex(wire) < 1)
                     {
                        wmap.Add (wire);

                        for (exp4.Init(wire, TopAbs_EDGE); exp4.More(); exp4.Next())
                        {
                           TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
                           if (emap.FindIndex(edge) < 1)
                           {
                              emap.Add (edge);
                              for (exp5.Init(edge, TopAbs_VERTEX); exp5.More(); exp5.Next())
                              {
                                 TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                                 if (vmap.FindIndex(vertex) < 1)
                                    vmap.Add (vertex);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }


      // Free Faces
      for (auto face : MyExplorer(shape, TopAbs_FACE, TopAbs_SHELL))
        if (!fmap.Contains(face))
          {
            fmap.Add (face);
            for (auto wire : MyExplorer(face, TopAbs_WIRE))
              if (!wmap.Contains(wire))
                {
                  wmap.Add (wire);
                  for (auto edge : MyExplorer(wire, TopAbs_EDGE))
                    if (!emap.Contains(edge))
                      {
                        emap.Add (edge);
                        for (auto vertex : MyExplorer(edge, TopAbs_VERTEX))
                          if (!vmap.Contains(vertex))
                            vmap.Add (vertex);
                      }
                }
          }


      // Free Wires

      for (exp3.Init(shape, TopAbs_WIRE, TopAbs_FACE); exp3.More(); exp3.Next())
      {
         TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
         if (wmap.FindIndex(wire) < 1)
         {
            wmap.Add (wire);

            for (exp4.Init(exp3.Current(), TopAbs_EDGE); exp4.More(); exp4.Next())
            {
               TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
               if (emap.FindIndex(edge) < 1)
               {
                  emap.Add (edge);
                  for (exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next())
                  {
                     TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                     if (vmap.FindIndex(vertex) < 1)
                        vmap.Add (vertex);
                  }
               }
            }
         }
      }


      // Free Edges
      /*
      for (exp4.Init(shape, TopAbs_EDGE, TopAbs_WIRE); exp4.More(); exp4.Next())
      {
         TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
         if (emap.FindIndex(edge) < 1)
         {
            emap.Add (edge);
            for (exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next())
            {
               TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
               if (vmap.FindIndex(vertex) < 1)
                  vmap.Add (vertex);
            }
         }
      }
      */
      for (auto edge : MyExplorer(shape, TopAbs_EDGE, TopAbs_WIRE))
        if (!emap.Contains(edge))
          {
            emap.Add (edge);
            for (auto vertex : MyExplorer(edge, TopAbs_VERTEX))
              if (!vmap.Contains(vertex))
                vmap.Add (vertex);
          }

      
      // Free Vertices
      /*
      for (exp5.Init(shape, TopAbs_VERTEX, TopAbs_EDGE); exp5.More(); exp5.Next())
      {
         TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
         if (vmap.FindIndex(vertex) < 1)
            vmap.Add (vertex);
      }
      */
      for (auto vertex : MyExplorer(shape, TopAbs_VERTEX, TopAbs_EDGE))
        if (!vmap.Contains(vertex))
          vmap.Add (vertex);

      facemeshstatus.DeleteAll();
      facemeshstatus.SetSize (fmap.Extent());
      facemeshstatus = 0;

      // Philippose - 15/01/2009
      face_maxh.DeleteAll();
      face_maxh.SetSize (fmap.Extent());
      face_maxh = 1e99; // mparam.maxh;

      // Philippose - 15/01/2010      
      face_maxh_modified.DeleteAll();      
      face_maxh_modified.SetSize(fmap.Extent());      
      face_maxh_modified = 0;
      

      // Philippose - 17/01/2009
      face_sel_status.DeleteAll();
      face_sel_status.SetSize (fmap.Extent());
      face_sel_status = 0;

      fvispar.SetSize (fmap.Extent());
      evispar.SetSize (emap.Extent());
      vvispar.SetSize (vmap.Extent());

      fsingular.SetSize (fmap.Extent());
      esingular.SetSize (emap.Extent());
      vsingular.SetSize (vmap.Extent());

      fsingular = esingular = vsingular = false;

      NetgenGeometry::Clear();

      // Add shapes
      for(auto i1 : Range(1, vmap.Extent()+1))
      {
          auto v = vmap(i1);
          auto occ_vertex = make_unique<OCCVertex>(TopoDS::Vertex(v));
          occ_vertex->nr = vertices.Size();

          if(HaveProperties(v))
            occ_vertex->properties = GetProperties(v);
          vertices.Append(std::move(occ_vertex));
      }

      for(auto i1 : Range(1, emap.Extent()+1))
      {
          auto e = emap(i1);
          auto edge = TopoDS::Edge(e);
          auto verts = GetVertices(e);
          if(verts.size() == 0)
            continue;
          auto occ_edge = make_unique<OCCEdge>(edge, GetVertex(verts[0]), GetVertex(verts[1]) );
          occ_edge->properties = GetProperties(e);
          edges.Append(std::move(occ_edge));
      }

      for(auto i1 : Range(1, fmap.Extent()+1))
      {
          auto f = fmap(i1);

          auto k = faces.Size();
          auto occ_face = make_unique<OCCFace>(f);

          for(auto e : GetEdges(f))
              occ_face->edges.Append( &GetEdge(e) );

          if(HaveProperties(f))
            occ_face->properties = GetProperties(f);
          faces.Append(std::move(occ_face));

          if(dimension==2)
              for(auto e : GetEdges(f))
              {
                  auto & edge = GetEdge(e);
                  if(e.Orientation() == TopAbs_REVERSED)
                      edge.domout = k;
                  else
                      edge.domin = k;
              }
      }


      for(auto i1 : Range(1, somap.Extent()+1))
      {
          auto s = somap(i1);
          int k = solids.Size();
          auto occ_solid = make_unique<OCCSolid>(s);
          if(HaveProperties(s))
            occ_solid->properties = GetProperties(s);
          solids.Append(std::move(occ_solid));

          for(auto f : GetFaces(s))
          {
              auto & face = static_cast<OCCFace&>(GetFace(f));
              if(face.domin==-1)
                  face.domin = k;
              else
                  face.domout = k;
              if(face.Shape().Orientation() == TopAbs_INTERNAL)
                  face.domout = k;
          }
      }

      // Add identifications
      auto add_identifications = [&](auto & shapes, auto & shape_map)
      {
          for(auto i1 : Range(1, shape_map.Extent()+1))
          {
            auto shape = shape_map(i1);
            if(HaveIdentifications(shape))
              for(auto & ident : GetIdentifications(shape))
                {
                    if(!shape_map.Contains(ident.from) || !shape_map.Contains(ident.to))
                        continue;
                    ShapeIdentification si{
                        &GetShape(ident.from),
                        &GetShape(ident.to),
                        ident.trafo,
                        ident.type,
                        ident.name
                    };
                    shapes[i1-1]->identifications.Append(si);
                }
          }
      };
      add_identifications( vertices, vmap );
      add_identifications( edges, emap );
      add_identifications( faces, fmap );

      bounding_box = ::netgen::GetBoundingBox( shape );
      ProcessIdentifications();
   }



   void OCCGeometry :: SewFaces ()
   {
      (*testout) << "Trying to sew faces ..." << endl;
      cout << "Trying to sew faces ..." << flush;

      BRepOffsetAPI_Sewing sewedObj(1);
 
      for (int i = 1; i <= fmap.Extent(); i++)
      {
         TopoDS_Face face = TopoDS::Face (fmap(i));
         sewedObj.Add (face);
      }

      sewedObj.Perform();

      if (!sewedObj.SewedShape().IsNull())
      {
         shape = sewedObj.SewedShape();
         cout << " done" << endl;
      }
      else
         cout << " not possible";
   }





   void OCCGeometry :: MakeSolid ()
   {
      TopExp_Explorer exp0;

      (*testout) << "Trying to build solids ..." << endl;
      cout << "Trying to build solids ..." << flush;

      BRepBuilderAPI_MakeSolid ms;
      int count = 0;
      for (exp0.Init(shape, TopAbs_SHELL); exp0.More(); exp0.Next())
      {
         count++;
         ms.Add (TopoDS::Shell(exp0.Current()));
      }

      if (!count)
      {
         cout << " not possible (no shells)" << endl;
         return;
      }

      BRepCheck_Analyzer ba(ms);
      if (ba.IsValid ())
      {
         Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
         sfs->Init (ms);

         sfs->SetPrecision(1e-5);
         sfs->SetMaxTolerance(1e-5);

         sfs->Perform();

         shape = sfs->Shape();

         for (exp0.Init(shape, TopAbs_SOLID); exp0.More(); exp0.Next())
         {
            TopoDS_Solid solid = TopoDS::Solid(exp0.Current());
            TopoDS_Solid newsolid = solid;
            BRepLib::OrientClosedSolid (newsolid);
            Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
            rebuild->Replace(solid, newsolid);

            TopoDS_Shape newshape = rebuild->Apply(shape, TopAbs_SHAPE, 1);
            shape = newshape;
         }

         cout << " done" << endl;
      }
      else
         cout << " not possible" << endl;
   }


  Array<const GeometryVertex*> OCCGeometry :: GetFaceVertices(const GeometryFace& face) const
  {
    Array<const GeometryVertex*> verts;
    const auto& occface = dynamic_cast<const OCCFace&>(face);
    for(auto& vert : GetVertices(occface.Shape()))
      verts.Append(&GetVertex(vert));
    return verts;
  }


   void OCCGeometry :: BuildVisualizationMesh (double deflection)
   {
      // cout << IM(5) << "Preparing visualization (deflection = " << deflection << ") ... " << flush;
      BuildTriangulation(shape);
      // cout << IM(5) << "done" << endl;
   }




   void OCCGeometry :: CalcBoundingBox ()
   {
      boundingbox = ::netgen::GetBoundingBox(shape);
      (*testout) << "Bounding Box = [" << boundingbox.PMin() << " - " << boundingbox.PMax() << "]" << endl;
      SetCenter();
   }


//    void OCCGeometry :: WriteOCC_STL(char * filename)
//    {
//       cout << "writing stl..."; cout.flush();
//       StlAPI_Writer writer;
//       writer.RelativeMode() = Standard_False;
// 
//       writer.SetDeflection(0.02);
//       writer.Write(shape,filename);
// 
//       cout << "done" << endl;
//    }


  void LoadOCCInto(OCCGeometry* occgeo, const filesystem::path & filename)
  {
      static Timer timer_all("LoadOCC"); RegionTimer rtall(timer_all);
      static Timer timer_readfile("LoadOCC-ReadFile");
      static Timer timer_transfer("LoadOCC-Transfer");
      static Timer timer_getnames("LoadOCC-get names");

      // Initiate a dummy XCAF Application to handle the STEP XCAF Document
      static Handle(XCAFApp_Application) dummy_app = XCAFApp_Application::GetApplication();

      // Create an XCAF Document to contain the STEP file itself
      Handle(TDocStd_Document) step_doc;

      // Check if a STEP File is already open under this handle, if so, close it to prevent
      // Segmentation Faults when trying to create a new document
      if(dummy_app->NbDocuments() > 0)
      {
         dummy_app->GetDocument(1,step_doc);
         dummy_app->Close(step_doc);
      }
      dummy_app->NewDocument ("STEP-XCAF",step_doc);

      timer_readfile.Start();
      STEPCAFControl_Reader reader;

      // Enable transfer of colours
      reader.SetColorMode(Standard_True);
      reader.SetNameMode(Standard_True);
      Standard_Integer stat = reader.ReadFile(filename.string().c_str());
      timer_readfile.Stop();

      timer_transfer.Start();
      if(stat != IFSelect_RetDone)
      {
        throw NgException("Couldn't load OCC geometry");
      }

      reader.Transfer(step_doc);
      timer_transfer.Stop();

      // Read in the shape(s) and the colours present in the STEP File
      auto step_shape_contents = XCAFDoc_DocumentTool::ShapeTool(step_doc->Main());

      TDF_LabelSequence step_shapes;
      step_shape_contents->GetShapes(step_shapes);

      // For the STEP File Reader in OCC, the 1st Shape contains the entire 
      // compound geometry as one shape
      auto main_shape = step_shape_contents->GetShape(step_shapes.Value(1)); 

      step_utils::LoadProperties(main_shape, reader, step_doc);

      occgeo->shape = main_shape;
      occgeo->changed = 1;
      occgeo->BuildFMap();
      occgeo->CalcBoundingBox();
      PrintContents (occgeo);
  }

   // Philippose - 23/02/2009
   /* Special IGES File load function including the ability
   to extract individual surface colours via the extended
   OpenCascade XDE and XCAF Feature set.
   */
   OCCGeometry *LoadOCC_IGES(const filesystem::path & filename)
   {
      OCCGeometry *occgeo;
      occgeo = new OCCGeometry;
      // Initiate a dummy XCAF Application to handle the IGES XCAF Document
      static Handle(XCAFApp_Application) dummy_app = XCAFApp_Application::GetApplication();

      // Create an XCAF Document to contain the IGES file itself
      Handle(TDocStd_Document) iges_doc;

      // Check if a IGES File is already open under this handle, if so, close it to prevent
      // Segmentation Faults when trying to create a new document
      if(dummy_app->NbDocuments() > 0)
      {
         dummy_app->GetDocument(1,iges_doc);
         dummy_app->Close(iges_doc);
      }
      dummy_app->NewDocument ("IGES-XCAF",iges_doc);

      IGESCAFControl_Reader reader;

      Standard_Integer stat = reader.ReadFile(filename.string().c_str());

      if(stat != IFSelect_RetDone)
      {
        throw NgException("Couldn't load occ");
      }

      // Enable transfer of colours
      reader.SetColorMode(Standard_True);
      reader.SetNameMode(Standard_True);

      reader.Transfer(iges_doc);

      // Read in the shape(s) and the colours present in the IGES File
      Handle(XCAFDoc_ShapeTool) iges_shape_contents = XCAFDoc_DocumentTool::ShapeTool(iges_doc->Main());
      Handle(XCAFDoc_ColorTool) iges_colour_contents = XCAFDoc_DocumentTool::ColorTool(iges_doc->Main());

      TDF_LabelSequence iges_shapes;
      iges_shape_contents->GetShapes(iges_shapes);

      // List out the available colours in the IGES File as Colour Names
      // TDF_LabelSequence all_colours;
      // iges_colour_contents->GetColors(all_colours);
      // PrintMessage(1,"Number of colours in IGES File: ",all_colours.Length());
      // for(int i = 1; i <= all_colours.Length(); i++)
      // {
      //    Quantity_Color col;
      //    stringstream col_rgb;
      //    iges_colour_contents->GetColor(all_colours.Value(i),col);
      //    col_rgb << " : (" << col.Red() << "," << col.Green() << "," << col.Blue() << ")";
      //    PrintMessage(1, "Colour [", i, "] = ",col.StringName(col.Name()),col_rgb.str());
      // }


      // For the IGES Reader, all the shapes can be exported as one compound shape
      // using the "OneShape" member
      auto shape = reader.OneShape();
      auto shapeTool = XCAFDoc_DocumentTool::ShapeTool(iges_doc->Main());
      // load colors
      for (auto typ : {TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
        for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
          {
            TDF_Label label;
            shapeTool->Search(e.Current(), label);

            if(label.IsNull())
              continue;

            XCAFPrs_IndexedDataMapOfShapeStyle set;
            TopLoc_Location loc;
            XCAFPrs::CollectStyleSettings(label, loc, set);
            XCAFPrs_Style aStyle;
            set.FindFromKey(e.Current(), aStyle);

            auto & prop = OCCGeometry::GetProperties(e.Current());
            if(aStyle.IsSetColorSurf())
              prop.col = step_utils::ReadColor(aStyle.GetColorSurfRGBA());
          }

      // load names
      auto workSession = reader.WS();
      auto model = workSession->Model();
      auto transProc = workSession->TransferReader()->TransientProcess();
      Standard_Integer nb = model->NbEntities();
      for (Standard_Integer i = 1; i <= nb; i++)
        {
          Handle(Standard_Transient) entity = model->Value(i);
          auto item = Handle(StepRepr_RepresentationItem)::DownCast(entity);

          if(item.IsNull())
            continue;

          TopoDS_Shape shape = TransferBRep::ShapeResult(transProc->Find(item));
          string name = item->Name()->ToCString();
          if (!transProc->IsBound(item))
            continue;

          OCCGeometry::GetProperties(shape).name = name;
        }

      occgeo->shape = shape;
      occgeo->changed = 1;
      occgeo->BuildFMap();

      occgeo->CalcBoundingBox();
      PrintContents (occgeo);
      return occgeo;
   }





   // Philippose - 29/01/2009
   /* Special STEP File load function including the ability
   to extract individual surface colours via the extended
   OpenCascade XDE and XCAF Feature set.
   */
   OCCGeometry * LoadOCC_STEP (const filesystem::path & filename)
   {
      OCCGeometry * occgeo;
      occgeo = new OCCGeometry;

      LoadOCCInto(occgeo, filename);
      return occgeo;
   }




   OCCGeometry *LoadOCC_BREP (const filesystem::path & filename)
   {
      OCCGeometry * occgeo;
      occgeo = new OCCGeometry;

      BRep_Builder aBuilder;
      Standard_Boolean result = BRepTools::Read(occgeo->shape, filename.string().c_str(), aBuilder);

      if(!result)
      {
         delete occgeo;
         return NULL;
      }

      occgeo->changed = 1;
      occgeo->BuildFMap();

      occgeo->CalcBoundingBox();
      PrintContents (occgeo);

      return occgeo;
   }


  void OCCGeometry :: Save (const filesystem::path & filename) const
  {
    string ext = ToLower(filename.extension());
    auto s_filename = filename.string();
    auto c_filename = s_filename.c_str();

    if (ext == ".igs")
      {
	IGESControl_Writer writer("millimeters", 1);
	writer.AddShape (shape);
	writer.Write (c_filename);
      }
    else if (ext == ".stp")
      {
          step_utils::WriteSTEP(*this, filename);
      }
    else if (ext == ".stl")
      {
	StlAPI_Writer writer;
	writer.ASCIIMode() = Standard_True;
	writer.Write (shape, c_filename);
      }
    else if (ext == ".stlb")
      {
	StlAPI_Writer writer;
	writer.ASCIIMode() = Standard_False;
	writer.Write (shape, c_filename);
      }

    throw NgException ("Unknown target format: " + filename);
  }

  void OCCGeometry :: SaveToMeshFile (ostream & ost) const
  {
    auto ss = make_shared<stringstream>();
    TextOutArchive out(ss);
    NetgenGeometry *geo = const_cast<OCCGeometry*>(this);
    out & geo;

    ost << "TextOutArchive" << endl;
    ost << ss->str().size() << endl;
    ost << ss->str();
  }

  void OCCGeometry :: DoArchive(Archive& ar)
  {
    constexpr int current_format_version = 0;

    int format_version = current_format_version;
    auto netgen_version = GetLibraryVersion("netgen");
    ar & netgen_version & format_version;

    if(ar.Output())
      {
        std::stringstream ss;
#if OCC_VERSION_HEX < 0x070600
        BRepTools::Write(shape, ss);
#else
        BRepTools::Write(shape, ss, false, false, TopTools_FormatVersion_VERSION_1);
#endif
        ar << ss.str();
      }
    else
      {
        if(format_version>current_format_version)
            throw Exception("Loading OCCGeometry from archive: unknown format version "
                    + ToString(format_version)
                    + ", written by netgen version "
                    + ToString(netgen_version));
        std::string str;
        ar & str;
        stringstream ss(str);
        BRep_Builder builder;
        BRepTools::Read(shape, ss, builder);
      }

    TopTools_IndexedMapOfShape shape_map;
    Array<TopoDS_Shape> shape_list;

    ar & dimension;
    for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
      for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
        {
          auto ds = e.Current();
          if(shape_map.FindIndex(ds) == 0)
            {
              shape_map.Add(ds);
              shape_list.Append(ds);
            }
        }

    for (auto s : shape_list)
      {
        bool has_properties = HaveProperties(s);
        ar & has_properties;
        if(has_properties)
          ar & GetProperties(s);

        bool has_identifications = HaveIdentifications(s);
        ar & has_identifications;
        if(has_identifications)
          {
            auto & idents = GetIdentifications(s);
            auto n_idents = idents.size();
            ar & n_idents;
            idents.resize(n_idents);
            for(auto i : Range(n_idents))
              {
                auto & id = idents[i];
                int id_from, id_to;
                if(ar.Output())
                {
                  id_from = shape_map.FindIndex(id.from)-1;
                  id_to = shape_map.FindIndex(id.to)-1;
                }
                ar & id_from & id_to & id.trafo & id.name;
                if(ar.Input())
                {
                    id.from = shape_list[id_from];
                    id.to = shape_list[id_to];
                }
              }
          }
      }

    if(ar.Input())
      {
        changed = 1;
        BuildFMap();
        CalcBoundingBox();
      }
  }
  
  const char * shapesname[] =
   {" ", "CompSolids", "Solids", "Shells",

   "Faces", "Wires", "Edges", "Vertices"};

  const char * shapename[] =
   {" ", "CompSolid", "Solid", "Shell",
   "Face", "Wire", "Edge", "Vertex"};

  const char * orientationstring[] =
    {"+", "-", "i", "e"};




   void OCCGeometry :: RecursiveTopologyTree (const TopoDS_Shape & sh,
      stringstream & str,
      TopAbs_ShapeEnum l,
      bool isfree,
      const char * lname)
   {
      if (l > TopAbs_VERTEX) return;

      TopExp_Explorer e;
      int count = 0;
      int count2 = 0;

      if (isfree)
         e.Init(sh, l, TopAbs_ShapeEnum(l-1));
      else
         e.Init(sh, l);

      for (; e.More(); e.Next())
      {
         count++;

         stringstream lname2;
         lname2 << lname << "/" << shapename[l] << count;
         str << lname2.str() << " ";

         switch (e.Current().ShapeType())
	   {
	   case TopAbs_SOLID:
	     count2 = somap.FindIndex(TopoDS::Solid(e.Current())); break;
	   case TopAbs_SHELL:
	     count2 = shmap.FindIndex(TopoDS::Shell(e.Current())); break;
	   case TopAbs_FACE:
	     count2 = fmap.FindIndex(TopoDS::Face(e.Current())); break;
	   case TopAbs_WIRE:
	     count2 = wmap.FindIndex(TopoDS::Wire(e.Current())); break;
	   case TopAbs_EDGE:
	     count2 = emap.FindIndex(TopoDS::Edge(e.Current())); break;
	   case TopAbs_VERTEX:
	     count2 = vmap.FindIndex(TopoDS::Vertex(e.Current())); break;
	   default:
	     cout << "RecursiveTopologyTree: Case " << e.Current().ShapeType() << " not handled" << endl;
         }

         int nrsubshapes = 0;

         if (l <= TopAbs_WIRE)
         {
            TopExp_Explorer e2;
            for (e2.Init (e.Current(), TopAbs_ShapeEnum (l+1));
               e2.More(); e2.Next())
               nrsubshapes++;
         }

         str << "{" << shapename[l] << " " << count2;
         if(HaveProperties(e.Current()))
           {
             const auto& props = GetProperties(e.Current());
             if(props.name || props.maxh < 1e99)
               str << " - ";
             if(props.name)
               str << props.GetName();
             if(props.maxh < 1e99)
               str << " maxh(" << props.maxh << ")";
           }

         if (l <= TopAbs_EDGE)
         {
            str << " (" << orientationstring[e.Current().Orientation()];
            if (nrsubshapes != 0) str << ", " << nrsubshapes;
            str << ") } ";
         }
         else
            str << " } ";

         RecursiveTopologyTree (e.Current(), str, TopAbs_ShapeEnum (l+1),
            false, (char*)lname2.str().c_str());

      }
   }




   void OCCGeometry :: GetTopologyTree (stringstream & str)
   {
      cout << "Building topology tree ... " << flush;
      RecursiveTopologyTree (shape, str, TopAbs_COMPSOLID, false, "CompSolids");
      RecursiveTopologyTree (shape, str, TopAbs_SOLID, true, "FreeSolids");
      RecursiveTopologyTree (shape, str, TopAbs_SHELL, true, "FreeShells");
      RecursiveTopologyTree (shape, str, TopAbs_FACE, true, "FreeFaces");
      RecursiveTopologyTree (shape, str, TopAbs_WIRE, true, "FreeWires");
      RecursiveTopologyTree (shape, str, TopAbs_EDGE, true, "FreeEdges");
      RecursiveTopologyTree (shape, str, TopAbs_VERTEX, true, "FreeVertices");
      str << flush;
      //  cout << "done" << endl;
   }




   void OCCGeometry :: CheckIrregularEntities(stringstream & str)
   {
      ShapeAnalysis_CheckSmallFace csm;

      csm.SetTolerance (1e-6);

      TopTools_DataMapOfShapeListOfShape mapEdges;
      ShapeAnalysis_DataMapOfShapeListOfReal mapParam;
      TopoDS_Compound theAllVert;

      int spotfaces = 0;
      int stripsupportfaces = 0;
      int singlestripfaces = 0;
      int stripfaces = 0;
      int facessplitbyvertices = 0;
      int stretchedpinfaces = 0;
      int smoothpinfaces = 0;
      int twistedfaces = 0;
      // int edgessamebutnotidentified = 0;

      cout << "checking faces ... " << flush;

      int i;
      for (i = 1; i <= fmap.Extent(); i++)
      {
         TopoDS_Face face = TopoDS::Face (fmap(i));
         TopoDS_Edge e1, e2;

         if (csm.CheckSpotFace (face))
         {
            if (!spotfaces++)
               str << "SpotFace {Spot face} ";

            (*testout) << "Face " << i << " is a spot face" << endl;
            str << "SpotFace/Face" << i << " ";
            str << "{Face " << i << " } ";
         }

         if (csm.IsStripSupport (face))
         {
            if (!stripsupportfaces++)
               str << "StripSupportFace {Strip support face} ";

            (*testout) << "Face " << i << " has strip support" << endl;
            str << "StripSupportFace/Face" << i << " ";
            str << "{Face " << i << " } ";
         }

         if (csm.CheckSingleStrip(face, e1, e2))
         {
            if (!singlestripfaces++)
               str << "SingleStripFace {Single strip face} ";

            (*testout) << "Face " << i << " is a single strip (edge " << emap.FindIndex(e1)
               << " and edge " << emap.FindIndex(e2) << " are identical)" << endl;
            str << "SingleStripFace/Face" << i << " ";
            str << "{Face " << i << " (edge " << emap.FindIndex(e1)
               << " and edge " << emap.FindIndex(e2) << " are identical)} ";
         }

         if (csm.CheckStripFace(face, e1, e2))
         {
            if (!stripfaces++)
               str << "StripFace {Strip face} ";

            (*testout) << "Face " << i << " is a strip (edge " << emap.FindIndex(e1)
               << " and edge " << emap.FindIndex(e2)
               << " are identical)" << endl;
            str << "StripFace/Face" << i << " ";
            str << "{Face " << i << " (edge " << emap.FindIndex(e1)
               << " and edge " << emap.FindIndex(e2) << " are identical)} ";
         }

         if (int count = csm.CheckSplittingVertices(face, mapEdges, mapParam, theAllVert))
         {
            if (!facessplitbyvertices++)
               str << "FaceSplitByVertices {Face split by vertices} ";

            (*testout) << "Face " << i << " is split by " << count
               << " vertex/vertices " << endl;
            str << "FaceSplitByVertices/Face" << i << " ";
            str << "{Face " << i << " (split by " << count << "vertex/vertices)} ";
         }

         int whatrow, sens;
         if (int type = csm.CheckPin (face, whatrow, sens))
         {
            if (type == 1)
            {
               if (!smoothpinfaces++)
                  str << "SmoothPinFace {Smooth pin face} ";

               (*testout) << "Face " << i << " is a smooth pin" << endl;
               str << "SmoothPinFace/Face" << i << " ";
               str << "{Face " << i << " } ";
            }
            else
            {
               if (!stretchedpinfaces++)
                  str << "StretchedPinFace {Stretched pin face} ";

               (*testout) << "Face " << i << " is a stretched pin" << endl;
               str << "StretchedPinFace/Face" << i << " ";
               str << "{Face " << i << " } ";
            }
         }

         double paramu, paramv;
         if (csm.CheckTwisted (face, paramu, paramv))
         {
            if (!twistedfaces++)
               str << "TwistedFace {Twisted face} ";

            (*testout) << "Face " << i << " is twisted" << endl;
            str << "TwistedFace/Face" << i << " ";
            str << "{Face " << i << " } ";
         }
      }

      cout << "done" << endl;
      cout << "checking edges ... " << flush;

      // double dmax;
      // int cnt = 0;
      NgArray <double> edgeLengths;
      NgArray <int> order;
      edgeLengths.SetSize (emap.Extent());
      order.SetSize (emap.Extent());

      for (i = 1; i <= emap.Extent(); i++)
      {
         TopoDS_Edge edge1 = TopoDS::Edge (emap(i));
         GProp_GProps system;
         BRepGProp::LinearProperties(edge1, system);
         edgeLengths[i-1] = system.Mass();
      }

      Sort (edgeLengths, order);

      str << "ShortestEdges {Shortest edges} ";
      for (i = 1; i <= min(20, emap.Extent()); i++)
      {
         str << "ShortestEdges/Edge" << i;
         str << " {Edge " << order[i-1] << " (L=" << edgeLengths[order[i-1]-1] << ")} ";
      }

      str << flush;

      cout << "done" << endl;
   }




   void OCCGeometry :: GetUnmeshedFaceInfo (stringstream & str)
   {
      for (int i = 1; i <= fmap.Extent(); i++)
      {
         if (facemeshstatus[i-1] == -1)
            str << "Face" << i << " {Face " << i << " } ";
      }
      str << flush;
   }




   void OCCGeometry :: GetNotDrawableFaces (stringstream & str)
   {
      for (int i = 1; i <= fmap.Extent(); i++)
      {
         if (!fvispar[i-1].IsDrawable())
            str << "Face" << i << " {Face " << i << " } ";
      }
      str << flush;
   }




   bool OCCGeometry :: ErrorInSurfaceMeshing ()
   {
      for (int i = 1; i <= fmap.Extent(); i++)
         if (facemeshstatus[i-1] == -1)
            return true;

      return false;
   }

  bool IsMappedShape(const Transformation<3> & trafo, const TopoDS_Shape & me, const TopoDS_Shape & you)
  {
      if(me.ShapeType() != you.ShapeType()) return false;

      Bnd_Box bbox;
      BRepBndLib::Add(me, bbox);
      BRepBndLib::Add(you, bbox);
      BoxTree<3> tree( occ2ng(bbox.CornerMin()), occ2ng(bbox.CornerMax()) );

      Point<3> c_me = occ2ng(Center(me));
      Point<3> c_you = occ2ng(Center(you));
      if(tree.GetTolerance() < Dist(trafo(c_me), c_you))
          return false;

      TopTools_IndexedMapOfShape vmap;
      std::vector<optional<TopoDS_Shape>> verts;

      auto verts_me = GetVertices(me);
      auto verts_you = GetVertices(you);

      if(verts_me.size() != verts_you.size())
          return false;

      for (auto i : Range(verts_me.size()))
      {
          auto s = verts_me[i];
          if(vmap.FindIndex(s) > 0)
              continue;
          auto p = trafo(occ2ng(s));
          tree.Insert( p, i );
          vmap.Add(s);
          verts.push_back(nullopt);
      }
          
      for (auto vert : verts_you)
      {
          auto s = vert;
          auto p = occ2ng(s);
          bool vert_mapped = false;
          tree.GetFirstIntersecting( p, p, [&](size_t i ) {
            vmap.Add(verts_me[i]);
            verts[vmap.FindIndex(verts_me[i])-1] = s;
            vert_mapped = true;
            return true;
          });
          if(!vert_mapped)
              return false;
      }
      return true;
  }

  void Identify(const TopoDS_Shape & me, const TopoDS_Shape & you, string name, Identifications::ID_TYPE type, std::optional<std::variant<gp_Trsf, gp_GTrsf>> opt_trafo) 
  {
    Transformation<3> trafo;
    if(opt_trafo)
    {
        trafo = occ2ng(*opt_trafo);
    }
    else
    {
        auto v = occ2ng(Center(you)) - occ2ng(Center(me));
        trafo = Transformation<3>(v);
    }

    ListOfShapes list_me, list_you;
    list_me.push_back(me);
    list_you.push_back(you);
    Identify(list_me, list_you, name, type, trafo);
  }

  void Identify(const ListOfShapes & me, const ListOfShapes & you, string name, Identifications::ID_TYPE type, Transformation<3> trafo) 
  {
    ListOfShapes id_me;
    ListOfShapes id_you;

    if(auto faces_me = me.Faces(); faces_me.size()>0)
    {
        id_me = faces_me;
        id_you = you.Faces();
    }
    else if(auto edges_me = me.Edges(); edges_me.size()>0)
    {
        id_me = edges_me;
        id_you = you.Edges();
    }
    else
    {
        id_me = me.Vertices();
        id_you = you.Vertices();
    }

    for(auto shape_me : id_me)
        for(auto shape_you : id_you)
        {
            if(!IsMappedShape(trafo, shape_me, shape_you))
                continue;

            OCCGeometry::GetIdentifications(shape_me).push_back
                (OCCIdentification { shape_me, shape_you, trafo, name, type });
        }
  }

  void OCCParameters :: Print(ostream & ost) const
   {
      ost << "OCC Parameters:" << endl
		 << "minimum edge length: " << resthminedgelenenable
		 << ", min len = " << resthminedgelen << endl;
   }

  DLL_HEADER extern OCCParameters occparam;
  OCCParameters occparam;













  // int OCCGeometry :: GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam)
  //  {
  //    return OCCGenerateMesh (*this, mesh, mparam, occparam);
  //  }
  static RegisterClassForArchive<OCCGeometry, NetgenGeometry> regnggeo;

  namespace step_utils
  {
      void LoadProperties(const TopoDS_Shape & shape,
                          const STEPCAFControl_Reader & reader,
                          const Handle(TDocStd_Document) step_doc)
      {
        static Timer t("step_utils::LoadProperties"); RegionTimer rt(t);

        auto workSession = reader.Reader().WS();
        auto model = workSession->Model();
        auto transferReader = workSession->TransferReader();
        auto transProc = transferReader->TransientProcess();
        auto shapeTool = XCAFDoc_DocumentTool::ShapeTool(step_doc->Main());

        // load colors
        for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
          for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
          {
            TDF_Label label;
            shapeTool->Search(e.Current(), label);

            if(label.IsNull())
                continue;

            XCAFPrs_IndexedDataMapOfShapeStyle set;
            TopLoc_Location loc;
            XCAFPrs::CollectStyleSettings(label, loc, set);
            XCAFPrs_Style aStyle;
            set.FindFromKey(e.Current(), aStyle);

            auto & prop = OCCGeometry::GetProperties(e.Current());
            if(aStyle.IsSetColorSurf())
                prop.col = step_utils::ReadColor(aStyle.GetColorSurfRGBA());
          }

        // load names
        Standard_Integer nb = model->NbEntities();
        for (Standard_Integer i = 1; i <= nb; i++)
          {
            Handle(Standard_Transient) entity = model->Value(i);
            auto item = Handle(StepRepr_RepresentationItem)::DownCast(entity);

            if(item.IsNull())
                continue;

            TopoDS_Shape shape = TransferBRep::ShapeResult(transProc->Find(item));
            string name = item->Name()->ToCString();
            if (!transProc->IsBound(item))
              continue;

            OCCGeometry::GetProperties(shape).name = name;
          }


        // load custom data (maxh etc.)
        for (Standard_Integer i = 1; i <= nb; i++)
          {
            Handle(Standard_Transient) entity = model->Value(i);

            auto item = Handle(StepRepr_CompoundRepresentationItem)::DownCast(entity);

            if(item.IsNull())
                continue;

            auto shape_item = item->ItemElementValue(1);
            TopoDS_Shape shape = TransferBRep::ShapeResult(transProc->Find(shape_item));
            string name = item->Name()->ToCString();

            if(name == "netgen_geometry_identification")
                ReadIdentifications(item, transProc);

            if(name != "netgen_geometry_properties")
                continue;

            auto & prop = OCCGeometry::GetProperties(shape);

            auto nprops = item->NbItemElement();

            for(auto i : Range(2, nprops+1))
            {
                auto prop_item = item->ItemElementValue(i);
                string prop_name = prop_item->Name()->ToCString();

                if(prop_name=="maxh")
                    prop.maxh = Handle(StepRepr_ValueRepresentationItem)::DownCast(prop_item)
                        ->ValueComponentMember()->Real();

                if(prop_name=="hpref")
                    prop.hpref = Handle(StepRepr_ValueRepresentationItem)::DownCast(prop_item)
                        ->ValueComponentMember()->Real();
            }
          }
      }

      void WriteProperties(const Handle(Interface_InterfaceModel) model, const Handle(Transfer_FinderProcess) finder, const TopoDS_Shape & shape)
      {
          static const ShapeProperties default_props;
          Handle(StepRepr_RepresentationItem) item = STEPConstruct::FindEntity(finder, shape);
          if(!item)
              return;
          auto prop = OCCGeometry::GetProperties(shape);

          if(auto n = prop.name)
              item->SetName(MakeName(*n));

          Array<Handle(StepRepr_RepresentationItem)> props;
          props.Append(item);

          if(auto maxh = prop.maxh; maxh != default_props.maxh)
              props.Append( MakeReal(maxh, "maxh") );

          if(auto hpref = prop.hpref; hpref != default_props.hpref)
              props.Append( MakeReal(hpref, "hpref") );

          if(props.Size()>1)
          {
              for(auto & item : props.Range(1, props.Size()))
                  model->AddEntity(item);

              auto compound = MakeCompound(props, "netgen_geometry_properties");
              model->AddEntity(compound);
          }

          WriteIdentifications(model, shape, finder);
      }

      void WriteIdentifications(const Handle(Interface_InterfaceModel) model, const TopoDS_Shape & shape, const Handle(Transfer_FinderProcess) finder)
      {
          Handle(StepRepr_RepresentationItem) item = STEPConstruct::FindEntity(finder, shape);
          if(!OCCGeometry::HaveIdentifications(shape))
            return;
          auto & identifications = OCCGeometry::GetIdentifications(shape);
          if(identifications.size()==0)
              return;
          // auto n = identifications.size();
          Array<Handle(StepRepr_RepresentationItem)> ident_items;
          ident_items.Append(item);

          for(auto & ident : identifications)
          {
              const auto& to = STEPConstruct::FindEntity(finder, ident.from == shape ? ident.to : ident.from);
              if(to.IsNull())
                  continue;
              Array<Handle(StepRepr_RepresentationItem)> items;
              items.Append(MakeReal(ident.from == shape ? 1 : 0));
              items.Append(to);
              auto & m = ident.trafo.GetMatrix();
              for(auto i : Range(9))
                  items.Append(MakeReal(m(i)));
              auto & v = ident.trafo.GetVector();
              for(auto i : Range(3))
                  items.Append(MakeReal(v(i)));
              items.Append(MakeInt(ident.type));
              for(auto & item : items.Range(0, items.Size()))
                  model->AddEntity(item);
              ident_items.Append(MakeCompound(items, ident.name));
          }
          for(auto & item : ident_items.Range(1, ident_items.Size()))
            model->AddEntity(item);
          auto comp = MakeCompound(ident_items, "netgen_geometry_identification");
          model->AddEntity(comp);
      }

      void ReadIdentifications(Handle(StepRepr_RepresentationItem) item, Handle(Transfer_TransientProcess) transProc)
      {
          auto idents = Handle(StepRepr_CompoundRepresentationItem)::DownCast(item);
          auto n = idents->NbItemElement();
          std::vector<OCCIdentification> result;
          auto shape_origin = TransferBRep::ShapeResult(transProc->Find(idents->ItemElementValue(1)));

          for(auto i : Range(2,n+1))
          {
              auto id_item = Handle(StepRepr_CompoundRepresentationItem)::DownCast(idents->ItemElementValue(i));
              OCCIdentification ident;
              ident.name = id_item->Name()->ToCString();
              auto is_from = ReadReal(id_item->ItemElementValue(1));
              if(is_from)
                {
                  ident.from = shape_origin;
                  ident.to = TransferBRep::ShapeResult(transProc->Find(id_item->ItemElementValue(2)));
                }
              else
                {
                  ident.from = TransferBRep::ShapeResult(
                      transProc->Find(id_item->ItemElementValue(2)));
                  ident.to = shape_origin;
                }

              auto & m = ident.trafo.GetMatrix();
              for(auto i : Range(9))
                  m(i) = ReadReal(id_item->ItemElementValue(3+i));
              auto & v = ident.trafo.GetVector();
              for(auto i : Range(3))
                  v(i) = ReadReal(id_item->ItemElementValue(12+i));
              ident.type = Identifications::ID_TYPE(ReadInt(id_item->ItemElementValue(15)));
              result.push_back(ident);
          }
          OCCGeometry::GetIdentifications(shape_origin) = result;
      }

      void WriteSTEP(const TopoDS_Shape & shape, const filesystem::path & filename)
      {
          Interface_Static::SetCVal("write.step.schema", "AP242IS");
          Interface_Static::SetIVal("write.step.assembly",1);
          Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication();
          Handle(TDocStd_Document) doc;

          app->NewDocument("STEP-XCAF", doc);
          Handle(XCAFDoc_ShapeTool) shapetool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
          Handle(XCAFDoc_ColorTool) colortool = XCAFDoc_DocumentTool::ColorTool(doc->Main());
          TDF_Label label = shapetool->NewShape();
          shapetool->SetShape(label, shape);

          Handle(XSControl_WorkSession) session = new XSControl_WorkSession;
          STEPCAFControl_Writer writer(session);
          const Handle(Interface_InterfaceModel) model = session->Model();

          // Set colors (BEFORE transferring shape into step data structures)
          for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
            for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
              {
                auto prop = OCCGeometry::GetProperties(e.Current());
                if(auto col = prop.col)
                    colortool->SetColor(e.Current(), step_utils::MakeColor(*col), XCAFDoc_ColorGen);
              }

          // Transfer shape into step data structures -> now we can manipulate/add step representation items
          writer.Transfer(doc, STEPControl_AsIs);

          // Write all other properties
          auto finder = session->TransferWriter()->FinderProcess();

          for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
            for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
                WriteProperties(model, finder, e.Current());

          writer.Write(filename.string().c_str());
      }

  } // namespace step_utils
}


#endif

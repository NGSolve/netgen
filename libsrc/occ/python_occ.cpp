#ifdef NG_PYTHON
#ifdef OCCGEOMETRY

#include <memory>

#include <general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include <meshing/python_mesh.hpp>
#include <meshing.hpp>

#include "occgeom.hpp"

#include <BOPAlgo_Builder.hxx>
#include <BRepLProp_SLProps.hxx>
#include <Message.hxx>
#include <Standard_GUID.hxx>
#include <Standard_Version.hxx>
#include <TDF_Attribute.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_MaterialTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>

using namespace netgen;

namespace netgen
{
  extern std::shared_ptr<NetgenGeometry> ng_geometry;
  extern std::shared_ptr<Mesh> mesh;
}

static string occparameter_description = R"delimiter(
OCC Specific Meshing Parameters
-------------------------------

closeedgefac: Optional[float] = 2.
  Factor for meshing close edges, if None it is disabled.

minedgelen: Optional[float] = 0.001
  Minimum edge length to be used for dividing edges to mesh points. If
  None this is disabled.

)delimiter";

void CreateOCCParametersFromKwargs(OCCParameters& occparam, py::dict kwargs)
{
  if(kwargs.contains("minedgelen"))
    {
      auto val = kwargs.attr("pop")("minedgelen");
      if(val.is_none())
        occparam.resthminedgelenenable = false;
      else
        {
          occparam.resthminedgelen = py::cast<double>(val);
          occparam.resthminedgelenenable = true;
        }
    }
}

extern py::object CastShape(const TopoDS_Shape & s);

DLL_HEADER void ExportNgOCCBasic(py::module &m);
DLL_HEADER void ExportNgOCCShapes(py::module &m);



DLL_HEADER void ExportNgOCC(py::module &m) 
{
  m.attr("occ_version") = OCC_VERSION_COMPLETE;

  // suppress info messages from occ (like statistics on Transfer)
  Message_Gravity aGravity = Message_Alarm;
  for (Message_SequenceOfPrinters::Iterator aPrinterIter (Message::DefaultMessenger()->Printers());
       aPrinterIter.More(); aPrinterIter.Next())
  {
    aPrinterIter.Value()->SetTraceLevel (aGravity);
  }

  ExportNgOCCBasic(m);
  ExportNgOCCShapes(m);

  static py::exception<Standard_Failure> exc(m, "OCCException");
  py::register_exception_translator([](std::exception_ptr p)
  {
    try {
      if(p) std::rethrow_exception(p);
    } catch (const Standard_Failure& e) {
      exc((string(e.DynamicType()->Name()) + ": " + e.GetMessageString()).c_str());
    }
  });
  
  py::class_<OCCGeometry, shared_ptr<OCCGeometry>, NetgenGeometry> (m, "OCCGeometry", R"raw_string(Use LoadOCCGeometry to load the geometry from a *.step file.)raw_string")
    /*
    .def(py::init<const TopoDS_Shape&>(), py::arg("shape"),
         "Create Netgen OCCGeometry from existing TopoDS_Shape")
    */
    .def(py::init([] (const TopoDS_Shape& shape, int occdim, bool copy)
                  {
                    auto geo = make_shared<OCCGeometry> (shape, occdim);
                    // ng_geometry = geo;
                    
                    // geo->BuildFMap();
                    // geo->CalcBoundingBox();
                    return geo;
                  }), py::arg("shape"), py::arg("dim")=3, py::arg("copy")=false,
         "Create Netgen OCCGeometry from existing TopoDS_Shape")
    
    .def(py::init([] (const std::vector<TopoDS_Shape> shapes)
                  {
                    BOPAlgo_Builder builder;
                    for (auto & s : shapes)
                      builder.AddArgument(s);                    
                    builder.Perform();
                    cout << "glued together" << endl;
                    
#ifdef OCC_HAVE_HISTORY
                    Handle(BRepTools_History) history = builder.History ();
                    
                    for (auto & s : shapes)
                      for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next())
                        if (auto name = OCCGeometry::global_shape_properties[e.Current()].name)
                          {
                            TopTools_ListOfShape modlist = history->Modified(e.Current());
                            for (auto mods : modlist)
                              OCCGeometry::global_shape_properties[mods].name = *name;
                          }
#endif // OCC_HAVE_HISTORY

                    auto geo = make_shared<OCCGeometry> (builder.Shape());
                    ng_geometry = geo;
                    // geo->BuildFMap();
                    // geo->CalcBoundingBox();
                    return geo;
                  }), py::arg("shape"),
         "Create Netgen OCCGeometry from existing TopoDS_Shape")
    
    .def(py::init([] (const string& filename)
                  {
                    shared_ptr<OCCGeometry> geo;
                    if(EndsWith(filename, ".step") || EndsWith(filename, ".stp"))
                      geo.reset(LoadOCC_STEP(filename));
                    else if(EndsWith(filename, ".brep"))
                      geo.reset(LoadOCC_BREP(filename));
                    else if(EndsWith(filename, ".iges"))
                      geo.reset(LoadOCC_IGES(filename));
                    else
                      throw Exception("Cannot load file " + filename + "\nValid formats are: step, stp, brep, iges");
                    ng_geometry = geo;
                    return geo;
                  }), py::arg("filename"),
        "Load OCC geometry from step, brep or iges file")
    .def(NGSPickle<OCCGeometry>())
    .def("Glue", &OCCGeometry::GlueGeometry)
    .def("Heal",[](OCCGeometry & self, double tolerance, bool fixsmalledges, bool fixspotstripfaces, bool sewfaces, bool makesolids, bool splitpartitions)
         {
           self.tolerance = tolerance;
           self.fixsmalledges = fixsmalledges;
           self.fixspotstripfaces = fixspotstripfaces;
           self.sewfaces = sewfaces;
           self.makesolids = makesolids;
           self.splitpartitions = splitpartitions;

           self.HealGeometry();
           self.BuildFMap();
         },py::arg("tolerance")=1e-3, py::arg("fixsmalledges")=true, py::arg("fixspotstripfaces")=true, py::arg("sewfaces")=true, py::arg("makesolids")=true, py::arg("splitpartitions")=false,R"raw_string(Heal the OCCGeometry.)raw_string",py::call_guard<py::gil_scoped_release>())
    .def("SetFaceMeshsize", [](OCCGeometry& self, size_t fnr, double meshsize)
                            {
                              self.SetFaceMaxH(fnr, meshsize);
                            }, "Set maximum meshsize for face fnr. Face numbers are 0 based.")
    .def("Draw", [](shared_ptr<OCCGeometry> geo)
    {
      ng_geometry = geo;
    })
    .def("_visualizationData", [] (shared_ptr<OCCGeometry> occ_geo)
         {
           std::vector<float> vertices;
           std::vector<int> trigs;
           std::vector<float> normals;
           std::vector<float> min = {std::numeric_limits<float>::max(),
                               std::numeric_limits<float>::max(),
                               std::numeric_limits<float>::max()};
           std::vector<float> max = {std::numeric_limits<float>::lowest(),
                               std::numeric_limits<float>::lowest(),
                               std::numeric_limits<float>::lowest()};
           std::vector<string> surfnames;
           auto box = occ_geo->GetBoundingBox();
           for(int i = 0; i < 3; i++)
             {
               min[i] = box.PMin()[i];
               max[i] = box.PMax()[i];
             }
           occ_geo->BuildVisualizationMesh(0.01);
           gp_Pnt2d uv;
           gp_Pnt pnt;
           gp_Vec n;
           gp_Pnt p[3];
           int count = 0;
           for (int i = 1; i <= occ_geo->fmap.Extent(); i++)
             {
               surfnames.push_back("occ_surface" + to_string(i));
               auto face = TopoDS::Face(occ_geo->fmap(i));
               auto surf = BRep_Tool::Surface(face);
               TopLoc_Location loc;
               BRepAdaptor_Surface sf(face, Standard_False);
               BRepLProp_SLProps prop(sf, 1, 1e-5);
               Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);
               if (triangulation.IsNull())
                 cout << "cannot visualize face " << i << endl;
               trigs.reserve(trigs.size() + triangulation->NbTriangles()*4);
               vertices.reserve(vertices.size() + triangulation->NbTriangles()*3*3);
               normals.reserve(normals.size() + triangulation->NbTriangles()*3*3);
               for (int j = 1; j < triangulation->NbTriangles()+1; j++)
                 {
                   auto triangle = triangulation->Triangle(j);
                   for (int k = 1; k < 4; k++)
                     p[k-1] = triangulation->Node(triangle(k)).Transformed(loc);
                   for (int k = 1; k < 4; k++)
                     {
                       vertices.insert(vertices.end(),{float(p[k-1].X()), float(p[k-1].Y()), float(p[k-1].Z())});
                       trigs.insert(trigs.end(),{count, count+1, count+2,i});
                       count += 3;
                       uv = triangulation->UVNode(triangle(k));
                       prop.SetParameters(uv.X(), uv.Y());
                       if (prop.IsNormalDefined())
                         n = prop.Normal();
                       else
                         {
                           gp_Vec a(p[0], p[1]);
                           gp_Vec b(p[0], p[2]);
                           n = b^a;
                         }
                       if (face.Orientation() == TopAbs_REVERSED) n*= -1;
                       normals.insert(normals.end(),{float(n.X()), float(n.Y()), float(n.Z())});
                     }
                 }
             }
            py::gil_scoped_acquire ac;
            py::dict res;
            py::list snames;
            for(auto name : surfnames)
              snames.append(py::cast(name));
            res["vertices"] = MoveToNumpy(vertices);
            res["triangles"] = MoveToNumpy(trigs);
            res["normals"] = MoveToNumpy(normals);
            res["surfnames"] = snames;
            res["min"] = MoveToNumpy(min);
            res["max"] = MoveToNumpy(max);
            return res;
         }, py::call_guard<py::gil_scoped_release>())
    .def("GenerateMesh", [](shared_ptr<OCCGeometry> geo,
                            MeshingParameters* pars, NgMPI_Comm comm, py::kwargs kwargs)
                         {
                           MeshingParameters mp;
                           OCCParameters occparam;
                           {
                             py::gil_scoped_acquire aq;
                             if(pars)
                               {
                                 auto mp_kwargs = CreateDictFromFlags(pars->geometrySpecificParameters);
                                 CreateOCCParametersFromKwargs(occparam, mp_kwargs);
                                 mp = *pars;
                               }
                             CreateOCCParametersFromKwargs(occparam, kwargs);
                             CreateMPfromKwargs(mp, kwargs);
                           }
                           geo->SetOCCParameters(occparam);
                           auto mesh = make_shared<Mesh>();
                           mesh->SetCommunicator(comm);
                           mesh->SetGeometry(geo);

                           if (comm.Rank()==0)
                             {
                               SetGlobalMesh(mesh);
                               auto result = geo->GenerateMesh(mesh, mp);
                               if(result != 0)
                                 {
                                   netgen::mesh = mesh;   // keep mesh for debugging
                                   throw Exception("Meshing failed!");
                                 }
                               ng_geometry = geo;
                               if (comm.Size() > 1)
                                 mesh->Distribute();
                             }
                           else
                             {
                               mesh->SendRecvMesh();
                             }
                           return mesh;
                         }, py::arg("mp") = nullptr, py::arg("comm")=NgMPI_Comm{},
      py::call_guard<py::gil_scoped_release>(),
         (meshingparameter_description + occparameter_description).c_str())
    .def_property_readonly("shape", [](const OCCGeometry & self) { return self.GetShape(); })
    ;

  
  
  m.def("LoadOCCGeometry",[] (filesystem::path filename)
        {
          cout << "WARNING: LoadOCCGeometry is deprecated! Just use the OCCGeometry(filename) constructor. It is able to read brep and iges files as well!" << endl;
          ifstream ist(filename);
          OCCGeometry * instance = new OCCGeometry();
          instance = LoadOCC_STEP(filename.c_str());
          ng_geometry = shared_ptr<OCCGeometry>(instance, NOOP_Deleter);
          return ng_geometry;
        },py::call_guard<py::gil_scoped_release>());


  m.def("TestXCAF", [] (TopoDS_Shape shape) {

      /*static*/ Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication();
      cout << endl << endl << endl;
      cout << "app = " << *reinterpret_cast<void**>(&app) << endl;
      Handle(TDocStd_Document) doc;
      cout << "nbdocs = " << app->NbDocuments() << endl;
      if(app->NbDocuments() > 0)
      {
         app->GetDocument(1,doc);
         // app->Close(doc);
      }
      else
        app->NewDocument ("STEP-XCAF",doc);
      Handle(XCAFDoc_ShapeTool) shape_tool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
      Handle(XCAFDoc_MaterialTool) material_tool = XCAFDoc_DocumentTool::MaterialTool(doc->Main());
      // Handle(XCAFDoc_VisMaterialTool) vismaterial_tool = XCAFDoc_DocumentTool::VisMaterialTool(doc->Main());

      // TDF_LabelSequence doc_shapes;
      // shape_tool->GetShapes(doc_shapes);
      // cout << "shape tool nbentities: " << doc_shapes.Size() << endl;
      TDF_Label label = shape_tool -> FindShape(shape);
      cout << "shape label = " << endl << label << endl;
      if (label.IsNull()) return;
      cout << "nbattr = " << label.NbAttributes() << endl;
                                                     
                                                     
      if (!label.IsNull())
        {
          Handle(TDF_Attribute) attribute;
          cout << "create guid" << endl;
          // Standard_GUID guid("c4ef4200-568f-11d1-8940-080009dc3333");
          Standard_GUID guid("2a96b608-ec8b-11d0-bee7-080009dc3333");      
          cout << "have guid" << endl;
          cout << "find attrib " << label.FindAttribute(guid, attribute) << endl;
          cout << "attrib = " << attribute << endl;
          cout << "tag = " << label.Tag() << endl;
          cout << "father.tag = " << label.Father().Tag() << endl;
          cout << "Data = " << label.Data() << endl;
          
          cout << "nbchild = " << label.NbChildren() << endl;
          for (auto i : Range(label.NbChildren()))
            {
              TDF_Label child = label.FindChild(i+1);
              cout << "child[" << i << "] = " << child << endl;
              cout << "find attrib " << child.FindAttribute(guid, attribute) << endl;
              cout << "attrib = " << attribute << endl;
            }
          
          // cout << "findshape = " << shape_tool -> FindShape(shape) << endl;
          cout << "IsMaterial = " << material_tool->IsMaterial(label) << endl;
          // cout << "IsVisMaterial = " << vismaterial_tool->IsMaterial(label) << endl;
        }
    }, py::arg("shape")=TopoDS_Shape());
        
}

PYBIND11_MODULE(libNgOCC, m) {
  ExportNgOCC(m);
}

#endif // OCCGEOMETRY
#endif // NG_PYTHON

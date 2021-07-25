#ifdef NG_PYTHON
#ifdef OCCGEOMETRY

#include <../general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include "../meshing/python_mesh.hpp"

#include <meshing.hpp>
#include <occgeom.hpp>
#include <Standard_Version.hxx>
#include <gp_Ax2.hxx>

#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <XCAFDoc_VisMaterialTool.hxx>
#include <TDF_Attribute.hxx>
#include <Standard_GUID.hxx>

using namespace netgen;

namespace netgen
{
  extern std::shared_ptr<NetgenGeometry> ng_geometry;
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


DLL_HEADER void ExportNgOCC(py::module &m) 
{
  m.attr("occ_version") = OCC_VERSION_COMPLETE;
  py::class_<OCCGeometry, shared_ptr<OCCGeometry>, NetgenGeometry> (m, "OCCGeometry", R"raw_string(Use LoadOCCGeometry to load the geometry from a *.step file.)raw_string")
    .def(py::init<>())
    /*
    .def(py::init<const TopoDS_Shape&>(), py::arg("shape"),
         "Create Netgen OCCGeometry from existing TopoDS_Shape")
    */
    .def(py::init([] (const TopoDS_Shape& shape)
                  {
                    auto geo = make_shared<OCCGeometry> (shape);
                    ng_geometry = geo;
                    
                    geo->BuildFMap();
                    geo->CalcBoundingBox();
                    return geo;
                  }), py::arg("shape"),
         "Create Netgen OCCGeometry from existing TopoDS_Shape")
    
    .def(py::init([] (const std::vector<TopoDS_Shape> shapes)
                  {
                    BOPAlgo_Builder aBuilder;
                    
                    // Setting arguments
                    TopTools_ListOfShape aLSObjects;
                    /*
                    for (TopExp_Explorer exp_solid(shape, TopAbs_SOLID); exp_solid.More(); exp_solid.Next())
                      aLSObjects.Append (exp_solid.Current());
                    */
                    for (auto & s : shapes)
                      aLSObjects.Append (s);
                    aBuilder.SetArguments(aLSObjects);
                    aBuilder.Perform();
    
                    auto geo = make_shared<OCCGeometry> (aBuilder.Shape());
                    ng_geometry = geo;
                    geo->BuildFMap();
                    geo->CalcBoundingBox();
                    return geo;
                  }), py::arg("shape"),
         "Create Netgen OCCGeometry from existing TopoDS_Shape")
    
    .def(py::init([] (const string& filename)
                  {
                    shared_ptr<OCCGeometry> geo;
                    if(EndsWith(filename, ".step") || EndsWith(filename, ".stp"))
                      geo.reset(LoadOCC_STEP(filename.c_str()));
                    else if(EndsWith(filename, ".brep"))
                      geo.reset(LoadOCC_BREP(filename.c_str()));
                    else if(EndsWith(filename, ".iges"))
                      geo.reset(LoadOCC_IGES(filename.c_str()));
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
                   auto triangle = (triangulation->Triangles())(j);
                   for (int k = 1; k < 4; k++)
                     p[k-1] = (triangulation->Nodes())(triangle(k)).Transformed(loc);
                   for (int k = 1; k < 4; k++)
                     {
                       vertices.insert(vertices.end(),{float(p[k-1].X()), float(p[k-1].Y()), float(p[k-1].Z())});
                       trigs.insert(trigs.end(),{count, count+1, count+2,i});
                       count += 3;
                       uv = (triangulation->UVNodes())(triangle(k));
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
                            MeshingParameters* pars, py::kwargs kwargs)
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
                           mesh->SetGeometry(geo);
                           auto result = geo->GenerateMesh(mesh, mp);
                           if(result != 0)
                             throw Exception("Meshing failed!");
                           SetGlobalMesh(mesh);
                           ng_geometry = geo;
                           return mesh;
                         }, py::arg("mp") = nullptr,
      py::call_guard<py::gil_scoped_release>(),
         (meshingparameter_description + occparameter_description).c_str())
    .def_property_readonly("shape", [](const OCCGeometry & self) { return self.GetShape(); })
    ;

  py::enum_<TopAbs_ShapeEnum>(m, "TopAbs_ShapeEnum", "Enumeration of all supported TopoDS_Shapes")
    .value("COMPOUND", TopAbs_COMPOUND)   .value("COMPSOLID", TopAbs_COMPSOLID)
    .value("SOLID", TopAbs_SOLID)       .value("SHELL", TopAbs_SHELL)
    .value("FACE", TopAbs_FACE)         .value("WIRE", TopAbs_WIRE)
    .value("EDGE", TopAbs_EDGE) .value("VERTEX", TopAbs_VERTEX)
    .value("SHAPE", TopAbs_SHAPE)
    .export_values()
    ;

  py::class_<gp_Pnt>(m, "gp_Pnt")
    .def(py::init([] (py::tuple pnt)
                  {
                    return gp_Pnt(py::cast<double>(pnt[0]),
                                  py::cast<double>(pnt[1]),
                                  py::cast<double>(pnt[2]));
                  }))
    ;

  py::class_<gp_Dir>(m, "gp_Dir")
    .def(py::init([] (py::tuple dir)
                  {
                    return gp_Dir(py::cast<double>(dir[0]),
                                  py::cast<double>(dir[1]),
                                  py::cast<double>(dir[2]));
                  }))
    ;
  py::implicitly_convertible<py::tuple, gp_Pnt>();
  py::implicitly_convertible<py::tuple, gp_Dir>();
  
  
  py::class_<TopoDS_Shape> (m, "TopoDS_Shape")
    .def("__str__", [] (const TopoDS_Shape & shape)
         {
           stringstream str;
           shape.DumpJson(str);
           return str.str();
         })
    
    .def("ShapeType", [] (const TopoDS_Shape & shape)
         { return shape.ShapeType(); })
    
    .def("SubShapes", [] (const TopoDS_Shape & shape, TopAbs_ShapeEnum & type)
         {
           py::list sub;
           TopExp_Explorer e;
           for (e.Init(shape, type); e.More(); e.Next())
             sub.append(e.Current());
           return sub;
         })

    .def("bc", [](const TopoDS_Shape & shape, const string & name)
         {
           TopExp_Explorer e;
           for (e.Init(shape, TopAbs_FACE); e.More(); e.Next())
             OCCGeometry::global_shape_names[e.Current().TShape()] = name;
           return shape;
         })
    
    .def("__add__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        return BRepAlgoAPI_Fuse(shape1, shape2).Shape();
      })
    
    .def("__mul__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        // return BRepAlgoAPI_Common(shape1, shape2).Shape();
        
        BRepAlgoAPI_Common builder(shape1, shape2);
        Handle(BRepTools_History) history = builder.History ();

        /*
          // work in progress ...
        TopTools_ListOfShape modlist = history->Modified(shape1);
        for (auto s : modlist)
          cout << "modified from list el: " << s.ShapeType() << endl;
        */
        
        for (TopExp_Explorer e(shape1, TopAbs_FACE); e.More(); e.Next())
          {
            const string & name = OCCGeometry::global_shape_names[e.Current().TShape()];
            TopTools_ListOfShape modlist = history->Modified(e.Current());
            for (auto s : modlist)
              OCCGeometry::global_shape_names[s.TShape()] = name;
          }
        
        for (TopExp_Explorer e(shape2, TopAbs_FACE); e.More(); e.Next())
          {
            const string & name = OCCGeometry::global_shape_names[e.Current().TShape()];
            TopTools_ListOfShape modlist = history->Modified(e.Current());
            for (auto s : modlist)
              OCCGeometry::global_shape_names[s.TShape()] = name;
          }        
        
        return builder.Shape();
      })
    
    .def("__sub__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        return BRepAlgoAPI_Cut(shape1, shape2).Shape();
      })
    ;


  m.def("Sphere", [] (gp_Pnt cc, double r) {
      return BRepPrimAPI_MakeSphere (cc, r).Shape();
    });
  
  m.def("Cylinder", [] (gp_Pnt cpnt, gp_Dir cdir, double r, double h) {
      return BRepPrimAPI_MakeCylinder (gp_Ax2(cpnt, cdir), r, h).Shape();
    });
  
  m.def("Box", [] (gp_Pnt cp1, gp_Pnt cp2) {
      return BRepPrimAPI_MakeBox (cp1, cp2).Shape();
    });
  
  m.def("LoadOCCGeometry",[] (const string & filename)
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
      Handle(XCAFDoc_VisMaterialTool) vismaterial_tool = XCAFDoc_DocumentTool::VisMaterialTool(doc->Main());

      cout << "handle(shape) = " << *(void**)(void*)(&(shape.TShape())) << endl;
      
      TDF_LabelSequence doc_shapes;
      shape_tool->GetShapes(doc_shapes);
      cout << "shape tool nbentities: " << doc_shapes.Size() << endl;
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
          cout << "IsVisMaterial = " << vismaterial_tool->IsMaterial(label) << endl;
        }
    }, py::arg("shape")=TopoDS_Shape());
        
}

PYBIND11_MODULE(libNgOCC, m) {
  ExportNgOCC(m);
}

#endif // OCCGEOMETRY
#endif // NG_PYTHON

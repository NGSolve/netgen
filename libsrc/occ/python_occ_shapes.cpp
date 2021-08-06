#ifdef NG_PYTHON
#ifdef OCCGEOMETRY

#include <../general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include "../meshing/python_mesh.hpp"

#include <meshing.hpp>
#include <occgeom.hpp>

#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax2d.hxx>
#include <gp_Trsf.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
// #include <XCAFDoc_VisMaterialTool.hxx>
#include <TDF_Attribute.hxx>
#include <Standard_GUID.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GC_MakeSegment.hxx>
#include <GC_MakeCircle.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>

#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepLib.hxx>

#include <Geom2d_Curve.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <GCE2d_MakeCircle.hxx>

#if OCC_VERSION_MAJOR>=7 && OCC_VERSION_MINOR>=4
#define OCC_HAVE_DUMP_JSON
#endif

using namespace netgen;


py::object CastShape(const TopoDS_Shape & s)
{
  switch (s.ShapeType())
    {
    case TopAbs_VERTEX:
      return py::cast(TopoDS::Vertex(s));
    case TopAbs_FACE:
      return py::cast(TopoDS::Face(s));      
    case TopAbs_EDGE:
      return py::cast(TopoDS::Edge(s));      
    case TopAbs_WIRE:
      return py::cast(TopoDS::Wire(s));      

    case TopAbs_COMPOUND:
    case TopAbs_COMPSOLID:
    case TopAbs_SOLID:
    case TopAbs_SHELL:
    case TopAbs_SHAPE:
      return py::cast(s);
    }
};


DLL_HEADER void ExportNgOCCShapes(py::module &m) 
{
  cout << "export shapes" << endl;
  py::enum_<TopAbs_ShapeEnum>(m, "TopAbs_ShapeEnum", "Enumeration of all supported TopoDS_Shapes")
    .value("COMPOUND", TopAbs_COMPOUND)   .value("COMPSOLID", TopAbs_COMPSOLID)
    .value("SOLID", TopAbs_SOLID)       .value("SHELL", TopAbs_SHELL)
    .value("FACE", TopAbs_FACE)         .value("WIRE", TopAbs_WIRE)
    .value("EDGE", TopAbs_EDGE) .value("VERTEX", TopAbs_VERTEX)
    .value("SHAPE", TopAbs_SHAPE)
    .export_values()
    ;


  class ListOfShapes : public std::vector<TopoDS_Shape> { };
  
  
  py::class_<TopoDS_Shape> (m, "TopoDS_Shape")
    .def("__str__", [] (const TopoDS_Shape & shape)
         {
           stringstream str;
#ifdef OCC_HAVE_DUMP_JSON
           shape.DumpJson(str);
#endif // OCC_HAVE_DUMP_JSON
           return str.str();
         })
    
    .def("ShapeType", [] (const TopoDS_Shape & shape)
         {
           cout << "WARNING: pls use 'shape' instead of 'ShapeType()'" << endl;
           return shape.ShapeType();
         })
    .def_property_readonly("type", [](const TopoDS_Shape & shape)
                           { return shape.ShapeType(); })    
    
    .def("SubShapes", [] (const TopoDS_Shape & shape, TopAbs_ShapeEnum & type)
         {
           /*
           py::list sub;
           TopExp_Explorer e;
           for (e.Init(shape, type); e.More(); e.Next())
             {
               switch (type)
                 {
                 case TopAbs_FACE:
                   sub.append(TopoDS::Face(e.Current())); break;
                 default:
                   sub.append(e.Current());
                 }
             }
           return sub;
           */
           ListOfShapes sub;
           for (TopExp_Explorer e(shape, type); e.More(); e.Next())
             sub.push_back(e.Current());
           return sub;
         })
    
    .def_property_readonly("faces", [] (const TopoDS_Shape & shape)
         {
           ListOfShapes sub;
           for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
             sub.push_back(e.Current());
           return sub;
         })
    .def_property_readonly("edges", [] (const TopoDS_Shape & shape)
         {
           ListOfShapes sub;
           for (TopExp_Explorer e(shape, TopAbs_EDGE); e.More(); e.Next())
             sub.push_back(e.Current());
           return sub;
         })
    .def_property_readonly("vertices", [] (const TopoDS_Shape & shape)
         {
           ListOfShapes sub;
           for (TopExp_Explorer e(shape, TopAbs_VERTEX); e.More(); e.Next())
             sub.push_back(e.Current());
           return sub;
         })

    .def("Properties", [] (const TopoDS_Shape & shape)
         {
           GProp_GProps props;
           switch (shape.ShapeType())
             {
             case TopAbs_FACE:
               BRepGProp::SurfaceProperties (shape, props); break;
             default:
               BRepGProp::LinearProperties(shape, props);
               // throw Exception("Properties implemented only for FACE");
             }
           double mass = props.Mass();
           gp_Pnt center = props.CentreOfMass();
           return tuple( py::cast(mass), py::cast(center) );
         })
    .def_property_readonly("center", [](const TopoDS_Shape & shape) {
           GProp_GProps props;
           switch (shape.ShapeType())
             {
             case TopAbs_FACE:
               BRepGProp::SurfaceProperties (shape, props); break;
             default:
               BRepGProp::LinearProperties(shape, props);
             }
           return props.CentreOfMass();
      })
    
    .def("bc", [](const TopoDS_Shape & shape, const string & name)
         {
           for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
             OCCGeometry::global_shape_properties[e.Current().TShape()].name = name;
           return shape;
         })

    .def("mat", [](const TopoDS_Shape & shape, const string & name)
         {
           for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next())
             OCCGeometry::global_shape_properties[e.Current().TShape()].name = name;
           return shape;
         })
    
    .def_property("name", [](const TopoDS_Shape & self) {
        if (auto name = OCCGeometry::global_shape_properties[self.TShape()].name)
          return *name;
        else
          return string();
      }, [](const TopoDS_Shape & self, string name) {
        OCCGeometry::global_shape_properties[self.TShape()].name = name;            
      })

    .def_property("col", [](const TopoDS_Shape & self) {
        auto it = OCCGeometry::global_shape_properties.find(self.TShape());
        Vec<3> col(0.2, 0.2, 0.2);
        if (it != OCCGeometry::global_shape_properties.end() && it->second.col)
          col = *it->second.col; // .value();
        return std::vector<double> ( { col(0), col(1), col(2) } );
      }, [](const TopoDS_Shape & self, std::vector<double> c) {
        Vec<3> col(c[0], c[1], c[2]);
        OCCGeometry::global_shape_properties[self.TShape()].col = col;    
      })
    
    .def_property("location",
                  [](const TopoDS_Shape & shape) { return shape.Location(); },
                  [](TopoDS_Shape & shape, const TopLoc_Location & loc)
                  { shape.Location(loc); })
    .def("Located", [](const TopoDS_Shape & shape, const TopLoc_Location & loc)
                  { return shape.Located(loc); })

    .def("__add__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        return BRepAlgoAPI_Fuse(shape1, shape2).Shape();
      })
    
    .def("__mul__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        // return BRepAlgoAPI_Common(shape1, shape2).Shape();
        
        BRepAlgoAPI_Common builder(shape1, shape2);
#ifdef OCC_HAVE_HISTORY
        Handle(BRepTools_History) history = builder.History ();

        /*
          // work in progress ...
        TopTools_ListOfShape modlist = history->Modified(shape1);
        for (auto s : modlist)
          cout << "modified from list el: " << s.ShapeType() << endl;
        */

        for (auto & s : { shape1, shape2 })
          for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
            {
              auto & prop = OCCGeometry::global_shape_properties[e.Current().TShape()];
              for (auto smod : history->Modified(e.Current()))            
                OCCGeometry::global_shape_properties[smod.TShape()].Merge(prop);
            }        
#endif // OCC_HAVE_HISTORY
        
        return builder.Shape();
      })
    
    .def("__sub__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        // return BRepAlgoAPI_Cut(shape1, shape2).Shape();
        
        BRepAlgoAPI_Cut builder(shape1, shape2);
#ifdef OCC_HAVE_HISTORY        
        Handle(BRepTools_History) history = builder.History ();

        for (auto s : { shape1, shape2 })
          for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
            {
              /*
              const string & name = OCCGeometry::global_shape_names[e.Current().TShape()];
              for (auto s : history->Modified(e.Current()))            
                OCCGeometry::global_shape_names[s.TShape()] = name;
              */
              /*
              auto it = OCCGeometry::global_shape_cols.find(e.Current().TShape());
              if (it != OCCGeometry::global_shape_cols.end())
                for (auto s : history->Modified(e.Current()))
                  OCCGeometry::global_shape_cols[s.TShape()] = it->second;
              */
              auto propit = OCCGeometry::global_shape_properties.find(e.Current().TShape());
              if (propit != OCCGeometry::global_shape_properties.end())
                for (auto s : history->Modified(e.Current()))
                  OCCGeometry::global_shape_properties[s.TShape()].Merge(propit->second);
            }

        /*
        for (TopExp_Explorer e(shape2, TopAbs_FACE); e.More(); e.Next())
          {
            auto it = OCCGeometry::global_shape_cols[e.Current().TShape()];
            if (it != OCCGeometry::global_shape_cols.end())
              for (auto s : history->Modified(e.Current()))
                OCCGeometry::global_shape_cols[s.TShape()] = it->second;
          }        
        */
#endif // OCC_HAVE_HISTORY

        
        return builder.Shape();        
      })

    .def("Reversed", [](const TopoDS_Shape & shape) {
        return CastShape(shape.Reversed()); })

    .def("Find", [](const TopoDS_Shape & shape, gp_Pnt p)
         {
           // find sub-shape contianing point
           // BRepClass_FaceClassifier::Perform  (p);
         })
    
    .def("MakeFillet", [](const TopoDS_Shape & shape, std::vector<TopoDS_Shape> edges, double r) {
        BRepFilletAPI_MakeFillet mkFillet(shape);
        for (auto e : edges)
          mkFillet.Add (r, TopoDS::Edge(e));
        return mkFillet.Shape();
      })
  
    .def("MakeThickSolid", [](const TopoDS_Shape & body, std::vector<TopoDS_Shape> facestoremove,
                              double offset, double tol) {
           TopTools_ListOfShape faces;
           for (auto f : facestoremove)
             faces.Append(f);
           
           BRepOffsetAPI_MakeThickSolid maker;
           maker.MakeThickSolidByJoin(body, faces, offset, tol);
           return maker.Shape();
         })
    
    .def("MakeTriangulation", [](const TopoDS_Shape & shape)
         {
           BRepTools::Clean (shape);
           double deflection = 0.01;
           BRepMesh_IncrementalMesh (shape, deflection, true);
         })
    
    .def("Triangulation", [](const TopoDS_Shape & shape)
         {
           // extracted from vsocc.cpp
           TopoDS_Face face;
           try
             {
               face = TopoDS::Face(shape);
             }
           catch (Standard_Failure & e)
             {
               e.Print (cout);
               throw NgException ("Triangulation: shape is not a face");
             }

           /*
           BRepTools::Clean (shape);
           double deflection = 0.01;
           BRepMesh_IncrementalMesh (shape, deflection, true);
           */

           Handle(Geom_Surface) surf = BRep_Tool::Surface (face);

           TopLoc_Location loc;
           Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);
           
           if (triangulation.IsNull())
             {
               BRepTools::Clean (shape);
               double deflection = 0.01;
               BRepMesh_IncrementalMesh (shape, deflection, true);
               triangulation = BRep_Tool::Triangulation (face, loc);               
             }
           // throw Exception("Don't have a triangulation, call 'MakeTriangulation' first");

           int ntriangles = triangulation -> NbTriangles();
           Array< std::array<Point<3>,3> > triangles;
           for (int j = 1; j <= ntriangles; j++)
             {
               Poly_Triangle triangle = (triangulation -> Triangles())(j);
               std::array<Point<3>,3> pts;
               for (int k = 0; k < 3; k++)
                 pts[k] = occ2ng( (triangulation -> Nodes())(triangle(k+1)).Transformed(loc) );
               
               triangles.Append ( pts );
             }
           
           // return MoveToNumpyArray(triangles);
           return triangles;
         })
    ;
  
  py::class_<TopoDS_Vertex, TopoDS_Shape> (m, "TopoDS_Vertex")
    .def(py::init([] (const TopoDS_Shape & shape) {
          return TopoDS::Vertex(shape);
        }))
    .def_property_readonly("p", [] (const TopoDS_Vertex & v) -> gp_Pnt {
        return BRep_Tool::Pnt (v); })
    ;
  
  py::class_<TopoDS_Edge, TopoDS_Shape> (m, "TopoDS_Edge")
    .def(py::init([] (const TopoDS_Shape & shape) {
          return TopoDS::Edge(shape);
        }))
    .def_property_readonly("start",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           return curve->Value(s0);
                           })
    .def_property_readonly("end",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           return curve->Value(s1);
                           })
    ;
  py::class_<TopoDS_Wire, TopoDS_Shape> (m, "TopoDS_Wire");
  py::class_<TopoDS_Face, TopoDS_Shape> (m, "TopoDS_Face")
    .def(py::init([] (const TopoDS_Shape & shape) {
          return TopoDS::Face(shape);
        }))
    /*
    .def("surf", [] (TopoDS_Face face) -> Handle(Geom_Surface)
         {
           Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
           return surf;
         })
    */
    ;
  py::class_<TopoDS_Solid, TopoDS_Shape> (m, "TopoDS_Solid");

  
  
  py::implicitly_convertible<TopoDS_Shape, TopoDS_Face>();
  


  class ListOfShapesIterator 
  {
    TopoDS_Shape * ptr;
  public:
    ListOfShapesIterator (TopoDS_Shape * aptr) : ptr(aptr) { }
    ListOfShapesIterator operator++ () { return ListOfShapesIterator(++ptr); }
    auto operator*() const { return CastShape(*ptr); }
    bool operator!=(ListOfShapesIterator it2) const { return ptr != it2.ptr; }
    bool operator==(ListOfShapesIterator it2) const { return ptr == it2.ptr; }
    
  };
  
  py::class_<ListOfShapes> (m, "ListOfShapes")
    .def("__iter__", [](ListOfShapes &s) {
        return py::make_iterator(ListOfShapesIterator(&*s.begin()),
                                 ListOfShapesIterator(&*s.end()));
      },
      py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("__getitem__", [](const ListOfShapes & list, size_t i) {
        return CastShape(list[i]); })
    
    .def("__getitem__", [](const ListOfShapes & self, py::slice inds) {
        size_t start, step, n, stop;
        if (!inds.compute(self.size(), &start, &stop, &step, &n))                                          
          throw py::error_already_set();
        ListOfShapes sub;
        sub.reserve(n);
        for (size_t i = 0; i < n; i++)
          sub.push_back (self[start+i*step]);
        return sub;
      })
    
    .def("__add__", [](const ListOfShapes & l1, const ListOfShapes & l2) {
        ListOfShapes l = l1;
        for (auto s : l2) l.push_back(s);
        return l;
      } )
    .def("__add__", [](const ListOfShapes & l1, py::list l2) {
        ListOfShapes l = l1;
        for (auto s : l2) l.push_back(py::cast<TopoDS_Shape>(s));
        return l;
      } )
    .def("__len__", [](const ListOfShapes & l) { return l.size(); })
    .def("Max", [] (ListOfShapes & shapes, gp_Vec dir)
         {
           double maxval = -1e99;
           TopoDS_Shape maxshape;
           for (auto shape : shapes)
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
           return CastShape(maxshape);
         })
    
    ;
         










  
  py::class_<Handle(Geom2d_Curve)> (m, "Geom2d_Curve")
    .def("Trim", [](Handle(Geom2d_Curve) curve, double u1, double u2) -> Handle(Geom2d_Curve)
         {
           return new Geom2d_TrimmedCurve (curve, u1, u2);
         })
    .def("Value", [](Handle(Geom2d_Curve) curve, double s) {
        return curve->Value(s);
      })
    .def_property_readonly("start", [](Handle(Geom2d_Curve) curve) {
        return curve->Value(curve->FirstParameter());
      })
    .def_property_readonly("end", [](Handle(Geom2d_Curve) curve) {
        return curve->Value(curve->LastParameter());
      })
    ;

  
  m.def("Sphere", [] (gp_Pnt cc, double r) {
      return BRepPrimAPI_MakeSphere (cc, r).Solid();
    });
  
  m.def("Cylinder", [] (gp_Pnt cpnt, gp_Dir cdir, double r, double h) {
      return BRepPrimAPI_MakeCylinder (gp_Ax2(cpnt, cdir), r, h).Solid();
    }, py::arg("p"), py::arg("d"), py::arg("r"), py::arg("h"));
  m.def("Cylinder", [] (gp_Ax2 ax, double r, double h) {
      return BRepPrimAPI_MakeCylinder (ax, r, h).Solid();
    }, py::arg("axis"), py::arg("r"), py::arg("h"));
  
  m.def("Box", [] (gp_Pnt cp1, gp_Pnt cp2) {
      return BRepPrimAPI_MakeBox (cp1, cp2).Solid();
    });

  m.def("Prism", [] (const TopoDS_Shape & face, gp_Vec vec) {
      return BRepPrimAPI_MakePrism (face, vec).Shape();
    });

  m.def("Pipe", [] (const TopoDS_Wire & spine, const TopoDS_Shape & profile) {
      return BRepOffsetAPI_MakePipe (spine, profile).Shape();
    }, py::arg("spine"), py::arg("profile"));

  // Handle(Geom2d_Ellipse) anEllipse1 = new Geom2d_Ellipse(anAx2d, aMajor, aMinor);
  m.def("Ellipse", [] (const gp_Ax2d & ax, double major, double minor) -> Handle(Geom2d_Curve)
        {
          return new Geom2d_Ellipse(ax, major, minor);
        });
  
  m.def("Segment", [](gp_Pnt2d p1, gp_Pnt2d p2) -> Handle(Geom2d_Curve) { 
      Handle(Geom2d_TrimmedCurve) curve = GCE2d_MakeSegment(p1, p2);
      return curve;
      // return BRepBuilderAPI_MakeEdge(curve).Edge();
      // return GCE2d_MakeSegment(p1, p2);      
    });
  
  m.def("Circle", [](gp_Pnt2d p1, double r) -> Handle(Geom2d_Curve) {
      Handle(Geom2d_Circle) curve = GCE2d_MakeCircle(p1, r);
      return curve;
      // gp_Ax2d ax; ax.SetLocation(p1);
      // return new Geom2d_Circle(ax, r);
    });
  
  m.def("Glue", [] (const std::vector<TopoDS_Shape> shapes) -> TopoDS_Shape
        {
          BOPAlgo_Builder builder;
          for (auto & s : shapes)
            {
              for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next())            
                builder.AddArgument(e.Current());
              if (s.ShapeType() == TopAbs_FACE)
                builder.AddArgument(s);
            }

          builder.Perform();

#ifdef OCC_HAVE_HISTORY          
          Handle(BRepTools_History) history = builder.History ();

          for (auto & s : shapes)
            for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next())
              {
                auto prop = OCCGeometry::global_shape_properties[e.Current().TShape()];
                for (auto mods : history->Modified(e.Current()))
                  OCCGeometry::global_shape_properties[mods.TShape()].Merge(prop);
              }
              /*
              {
                auto name = OCCGeometry::global_shape_names[e.Current().TShape()];
                for (auto mods : history->Modified(e.Current()))
                  OCCGeometry::global_shape_names[mods.TShape()] = name;
              }
              */
#endif // OCC_HAVE_HISTORY
          
          return builder.Shape();
        });

  m.def("Glue", [] (TopoDS_Shape shape) -> TopoDS_Shape
        {
          BOPAlgo_Builder builder;
          
          for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next())
            builder.AddArgument(e.Current());
          
          builder.Perform();
          
          if (builder.HasErrors())
            builder.DumpErrors(cout);
          if (builder.HasWarnings())
            builder.DumpWarnings(cout);

#ifdef OCC_HAVE_HISTORY
          Handle(BRepTools_History) history = builder.History ();

          for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next())
            {
              auto prop = OCCGeometry::global_shape_properties[e.Current().TShape()];
              for (auto mods : history->Modified(e.Current()))
                OCCGeometry::global_shape_properties[mods.TShape()].Merge(prop);
            }
#endif // OCC_HAVE_HISTORY
          
          return builder.Shape();
        });


  // py::class_<Handle(Geom_TrimmedCurve)> (m, "Geom_TrimmedCurve")
  // ;
  
  m.def("Segment", [](gp_Pnt p1, gp_Pnt p2) { 
      Handle(Geom_TrimmedCurve) curve = GC_MakeSegment(p1, p2);
      return BRepBuilderAPI_MakeEdge(curve).Edge();
    });
  m.def("Circle", [](gp_Pnt c, gp_Dir n, double r) {
	Handle(Geom_Circle) curve = GC_MakeCircle (c, n, r);
        return BRepBuilderAPI_MakeEdge(curve).Edge();
    });

  m.def("ArcOfCircle", [](gp_Pnt p1, gp_Pnt p2, gp_Pnt p3) { 
      Handle(Geom_TrimmedCurve) curve = GC_MakeArcOfCircle(p1, p2, p3);
      return BRepBuilderAPI_MakeEdge(curve).Edge();
    }, py::arg("p1"), py::arg("p2"), py::arg("p3"));
  
  m.def("ArcOfCircle", [](gp_Pnt p1, gp_Vec v, gp_Pnt p2) { 
      Handle(Geom_TrimmedCurve) curve = GC_MakeArcOfCircle(p1, v, p2);
      return BRepBuilderAPI_MakeEdge(curve).Edge();
    }, py::arg("p1"), py::arg("v"), py::arg("p2"));


  m.def("Edge", [](Handle(Geom2d_Curve) curve2d, TopoDS_Face face) {
      auto edge = BRepBuilderAPI_MakeEdge(curve2d, BRep_Tool::Surface (face)).Edge();
      BRepLib::BuildCurves3d(edge);
      return edge;
    });
  
  m.def("Wire", [](std::vector<TopoDS_Shape> edges) {
      BRepBuilderAPI_MakeWire builder;
      for (auto s : edges)
        switch (s.ShapeType())
          {
          case TopAbs_EDGE:
            try
              {
              builder.Add(TopoDS::Edge(s)); break;
              }
            catch (Standard_Failure & e)
              {
                e.Print(cout);
                throw NgException("cannot add to wire");
              }
          case TopAbs_WIRE:
            builder.Add(TopoDS::Wire(s)); break;
          default:
            throw Exception("can make wire only from edges and wires");
          }
      try
        {
          return builder.Wire();
        }
      catch (Standard_Failure & e)
        {
          e.Print(cout);
          throw NgException("error in wire builder");
        }
    });

  m.def("Face", [](TopoDS_Wire wire) {
      return BRepBuilderAPI_MakeFace(wire).Face();
    }, py::arg("w"));
  m.def("Face", [](const TopoDS_Face & face, const TopoDS_Wire & wire) {
      // return BRepBuilderAPI_MakeFace(face, wire).Face();
      return BRepBuilderAPI_MakeFace(BRep_Tool::Surface (face), wire).Face();
    }, py::arg("f"), py::arg("w"));
  m.def("Face", [](const TopoDS_Face & face, std::vector<TopoDS_Wire> wires) {
      // return BRepBuilderAPI_MakeFace(face, wire).Face();
      cout << "build from list of wires" << endl;
      auto surf = BRep_Tool::Surface (face);
      BRepBuilderAPI_MakeFace builder(surf, 1e-8);
      for (auto w : wires)
        builder.Add(w);
      return builder.Face();
    }, py::arg("f"), py::arg("w"));
  /*
     not yet working .... ?
  m.def("Face", [](std::vector<TopoDS_Wire> wires) {
      cout << "face from wires" << endl;
      BRepBuilderAPI_MakeFace builder;
      for (auto w : wires)
        {
          cout << "add wire" << endl;
          builder.Add(w);
        }
      return builder.Face();
    }, py::arg("w"));
  */

  m.def("MakeFillet", [](TopoDS_Shape shape, std::vector<TopoDS_Shape> edges, double r) {
      throw Exception("call 'shape.MakeFilled'");
      BRepFilletAPI_MakeFillet mkFillet(shape);
      for (auto e : edges)
        mkFillet.Add (r, TopoDS::Edge(e));
      return mkFillet.Shape();
    });
  
  m.def("MakeThickSolid", [](TopoDS_Shape body, std::vector<TopoDS_Shape> facestoremove,
                             double offset, double tol) {
          throw Exception("call 'shape.MakeThickSolid'");
          TopTools_ListOfShape faces;
          for (auto f : facestoremove)
            faces.Append(f);
          
          BRepOffsetAPI_MakeThickSolid maker;
          maker.MakeThickSolidByJoin(body, faces, offset, tol);
          return maker.Shape();
        });

  m.def("ThruSections", [](std::vector<TopoDS_Shape> wires)
        {
          BRepOffsetAPI_ThruSections aTool(Standard_True);
          for (auto shape : wires)
            aTool.AddWire(TopoDS::Wire(shape));
          aTool.CheckCompatibility(Standard_False);
          return aTool.Shape();
        });



  
  
}

#endif // OCCGEOMETRY
#endif // NG_PYTHON

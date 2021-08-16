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
#include <BRepPrimAPI_MakeRevol.hxx>
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
#include <Geom_Plane.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <GC_MakeSegment.hxx>
#include <GC_MakeCircle.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>

#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepLib.hxx>

#include <Geom2d_Curve.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <GCE2d_MakeCircle.hxx>
#include <GCE2d_MakeArcOfCircle.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <GeomLProp_SLProps.hxx>

#if OCC_VERSION_MAJOR>=7 && OCC_VERSION_MINOR>=4
#define OCC_HAVE_DUMP_JSON
#endif

using namespace netgen;


class ListOfShapes : public std::vector<TopoDS_Shape> { };


void ExtractEdgeData( const TopoDS_Edge & edge, int index, std::vector<double> * p, Box<3> & box )
{
    if (BRep_Tool::Degenerated(edge)) return;

    Handle(Poly_PolygonOnTriangulation) poly;
    Handle(Poly_Triangulation) T;
    TopLoc_Location loc;
    BRep_Tool::PolygonOnTriangulation(edge, poly, T, loc);

    if (poly.IsNull())
      {
        cout << "no edge mesh, do my own sampling" << endl;

        double s0, s1;
        Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);

        constexpr int num = 100;
        for (int i = 0; i < num; i++)
          {
            auto p0 = occ2ng(c->Value (s0 + i*(s1-s0)/num));
            auto p1 = occ2ng(c->Value (s0 + (i+1)*(s1-s0)/num));
            for(auto k : Range(3))
              {
                p[0].push_back(p0[k]);
                p[1].push_back(p1[k]);
              }
            p[0].push_back(index);
            p[1].push_back(index);
            box.Add(p0);
            box.Add(p1);
          }
        return;
      }        

    int nbnodes = poly -> NbNodes();
    for (int j = 1; j < nbnodes; j++)
    {
        auto p0 = occ2ng((T -> Nodes())(poly->Nodes()(j)).Transformed(loc));
        auto p1 = occ2ng((T -> Nodes())(poly->Nodes()(j+1)).Transformed(loc));
        for(auto k : Range(3))
        {
            p[0].push_back(p0[k]);
            p[1].push_back(p1[k]);
        }
        p[0].push_back(index);
        p[1].push_back(index);
        box.Add(p0);
        box.Add(p1);
    }
}

void ExtractFaceData( const TopoDS_Face & face, int index, std::vector<double> * p, Box<3> & box )
{
    TopLoc_Location loc;
    Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);

    bool flip = TopAbs_REVERSED == face.Orientation();

    if (triangulation.IsNull())
      {
        cout << "pls build face triangulation before" << endl;
        return;
      }

    int ntriangles = triangulation -> NbTriangles();
    for (int j = 1; j <= ntriangles; j++)
    {
        Poly_Triangle triangle = (triangulation -> Triangles())(j);
        std::array<Point<3>,3> pts;
        for (int k = 0; k < 3; k++)
            pts[k] = occ2ng( (triangulation -> Nodes())(triangle(k+1)).Transformed(loc) );

        if(flip)
            Swap(pts[1], pts[2]);

        for (int k = 0; k < 3; k++)
        {
            box.Add(pts[k]);
            for (int d = 0; d < 3; d++)
                p[k].push_back( pts[k][d] );
            p[k].push_back( index );
        }
    }
}

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


class WorkPlane : public enable_shared_from_this<WorkPlane>
{
  gp_Ax3 axis;
  gp_Ax2d localpos;
  gp_Pnt2d startpnt;
  TopoDS_Vertex lastvertex, startvertex;
  Handle(Geom_Surface) surf;
  // Geom_Plane surf;

  BRepBuilderAPI_MakeWire wire_builder;
  std::vector<TopoDS_Wire> wires;
  
public:
  
  WorkPlane (const gp_Ax3 & _axis, const gp_Ax2d _localpos = gp_Ax2d())
    : axis(_axis), localpos(_localpos) // , surf(_axis) 
  {
    // surf = GC_MakePlane (gp_Ax1(axis.Location(), axis.Direction()));
    surf = new Geom_Plane(axis);
  }
  

  auto MoveTo (double h, double v)
  {
    startpnt = gp_Pnt2d(h,v);
    localpos.SetLocation(startpnt);
    startvertex.Nullify();
    return shared_from_this();
  }

  auto Direction (double h, double v)
  {
    localpos.SetDirection(gp_Dir2d(h,v));
    return shared_from_this();
  }
  
  auto LineTo (double h, double v, optional<string> name = nullopt)
  {
    gp_Pnt2d old2d = localpos.Location();
    gp_Pnt oldp = axis.Location() . Translated(old2d.X() * axis.XDirection() + old2d.Y() * axis.YDirection());

    // localpos.Translate (gp_Vec2d(h,v));
    localpos.SetLocation (gp_Pnt2d(h,v));
    gp_Pnt2d new2d = localpos.Location();
    gp_Pnt newp = axis.Location() . Translated(new2d.X() * axis.XDirection() + new2d.Y() * axis.YDirection());

    if (new2d.Distance(old2d) < 1e-10) return shared_from_this();    
    bool closing = new2d.Distance(startpnt) < 1e-10;

      
    cout << "lineto, oldp = " << occ2ng(oldp) << endl;
    cout << "lineto, newp = " << occ2ng(newp) << endl;
    gp_Pnt pfromsurf = surf->Value(new2d.X(), new2d.Y());
    cout << "p from plane = " << occ2ng(pfromsurf) << endl;
    
    Handle(Geom_TrimmedCurve) curve = GC_MakeSegment(oldp, newp);

    if (startvertex.IsNull())
      startvertex = lastvertex = BRepBuilderAPI_MakeVertex(oldp);
    auto endv = closing ? startvertex : BRepBuilderAPI_MakeVertex(newp);
    // liefert noch Fehler bei close
    auto edge = BRepBuilderAPI_MakeEdge(curve, lastvertex, endv).Edge();
    lastvertex = endv;

    // auto edge = BRepBuilderAPI_MakeEdge(curve).Edge();
    if (name)
      OCCGeometry::global_shape_properties[edge.TShape()].name = name;
    wire_builder.Add(edge);

    if (closing) Close();
    return shared_from_this();    
  }

  auto Line(double h, double v, optional<string> name = nullopt)
  {
    gp_Pnt2d oldp = localpos.Location();
    oldp.Translate(gp_Vec2d(h,v));
    return LineTo (oldp.X(), oldp.Y(), name);
  }
  
  auto Line(double len, optional<string> name = nullopt)
  {
    gp_Dir2d dir = localpos.Direction();
    cout << "dir = " << dir.X() << ", " << dir.Y() << endl;
    gp_Pnt2d oldp = localpos.Location();
    oldp.Translate(len*dir);
    return LineTo (oldp.X(), oldp.Y(), name);
  }

  auto Rotate (double angle)
  {
    localpos.Rotate(localpos.Location(), angle*M_PI/180);
    return shared_from_this();
  }

  auto ArcTo (double h, double v, const gp_Vec2d t)
  {
    gp_Pnt2d P1 = localpos.Location();

    //check input
    if(P1.X() == h && P1.Y() == v)
        throw Exception("points P1 and P2 must not be congruent");

    localpos.SetLocation (gp_Pnt2d(h,v));
    gp_Pnt2d P2 = localpos.Location();

    cout << "ArcTo:" << endl;
    cout << "P1 = (" << P1.X() <<", " << P1.Y() << ")"<<endl;
    cout << "P2 = (" << P2.X() <<", " << P2.Y() << ")"<<endl;
    cout << "t = (" << t.X() << ", " << t.Y() << ")" << endl;

    //compute circle center point M
    //point midway between p1 and p2
    gp_Pnt2d P12 = gp_Pnt2d((P1.X() + h) / 2, (P1.Y() + v) / 2);
    //vector normal to vector from P1 to P12
    gp_Vec2d p12n = gp_Vec2d( - (P12.Y() - P1.Y()), (P12.X() - P1.X()));
    //M is intersection of p12n and tn (tn ... normalvector to t)
    double k = ((P12.Y()- P1.Y())*p12n.X() + (P1.X() - P12.X())*p12n.Y() )/ (t.X()*p12n.X() + t.Y()*p12n.Y());
    gp_Pnt2d M = gp_Pnt2d(P1.X()-k*t.Y(), P1.Y() + k*t.X());

    cout << "P12 = (" << P12.X() <<", " << P12.Y() << ")"<<endl;
    cout << "p12n = (" << p12n.X() <<", " << p12n.Y() << ")"<<endl;
    cout << "k = " << k <<endl;
    cout << "M = (" << M.X() <<", " << M.Y() << ")"<<endl;

    //radius
    double r = P1.Distance(M);

    //compute point P3 on circle between P1 and P2
    p12n.Normalize();   //docu: reverses direction of p12n ??
    cout << "p12n = (" << p12n.X() <<", " << p12n.Y() << ")"<<endl;

    gp_Pnt2d P3;

    double angletp12n = t.Angle(p12n);
    if(angletp12n > -M_PI/2 && angletp12n < M_PI/2)
        P3 = gp_Pnt2d(M.X() + r * p12n.X() , M.Y() + r * p12n.Y());
    else
        P3 = gp_Pnt2d(M.X() - r * p12n.X() , M.Y() - r * p12n.Y());

    cout << "r = " << r <<endl;
    cout << "angle t,p12n = " << t.Angle(p12n)<<endl;
    cout << "P3 = (" << P3.X() <<", " << P3.Y() << ")"<<endl;
    cout << "dist(M,P3) = " << P3.Distance(M) <<endl;

    //Draw 2d arc of circle from P1 to P2 through P3
    Handle(Geom2d_TrimmedCurve) curve2d = GCE2d_MakeArcOfCircle(P1, P3, P2).Value();

    gp_Pnt P13d = surf->Value(P1.X(), P1.Y());
    gp_Pnt P23d = surf->Value(P2.X(), P2.Y());
    cout << "p13d = " << occ2ng(P13d) << ", p23d = " << occ2ng(P23d) << endl;
    bool closing = P2.Distance(startpnt) < 1e-10;
    if (startvertex.IsNull())
      startvertex = lastvertex = BRepBuilderAPI_MakeVertex(P13d);
    auto endv = closing ? startvertex : BRepBuilderAPI_MakeVertex(P23d);
    // liefert noch Fehler bei close

    cout << "closing = " << closing << endl;
    cout << "startv isnull = " << lastvertex.IsNull() << endl;
    cout << "endv isnull = " << endv.IsNull() << endl;
    //create 3d edge from 2d curve using surf
    auto edge = BRepBuilderAPI_MakeEdge(curve2d, surf, lastvertex, endv).Edge();
    lastvertex = endv;
    cout << "have edge" << endl;
    BRepLib::BuildCurves3d(edge);
    cout << "Have curve3d" << endl;
    wire_builder.Add(edge);

    //compute angle of rotation
    //compute tangent t2 in P2
    gp_Vec2d p2 = gp_Vec2d(P1.X()-P2.X(),P1.Y()-P2.Y());
    gp_Vec2d t2;
    if(t.Angle(p2) >=0)
        t2 = gp_Vec2d((P2.Y()-M.Y()),-(P2.X()-M.X()));
    else
        t2 = gp_Vec2d(-(P2.Y()-M.Y()),(P2.X()-M.X()));
    double angle = -t2.Angle(t);    //angle \in [-pi,pi]
    cout << "angle t2,t = " << angle*180/M_PI << endl;

    //update localpos.Direction()
    Rotate(angle*180/M_PI);
    if (closing)
      {
        cout << "call close from arc" << endl;
        Close();
        cout << "close is back" << endl;
      }
    return shared_from_this();
  }

  auto Arc(double radius, double angle)
  {
    double newAngle = fmod(angle,360)*M_PI/180;

    //check input
    if(newAngle<1e-16 && newAngle>-1e-16)
        throw Exception("angle must not be an integer multiple of 360");

    gp_Dir2d dir = localpos.Direction();
    gp_Dir2d dirn;
    //compute center point of arc
    if(newAngle>=0)
        dirn = gp_Dir2d(-dir.Y(),dir.X());
    else
        dirn = gp_Dir2d(dir.Y(),-dir.X());

    gp_Pnt2d oldp = localpos.Location();

    oldp.Translate(radius*dirn);

    cout << "M = (" << oldp.X() << ", " << oldp.Y() << ")" << endl;

    dirn.Rotate(newAngle-M_PI);
    oldp.Translate(radius*dirn);

    //compute tangent vector in P1
    gp_Vec2d t = gp_Vec2d(dir.X(),dir.Y());

    cout << "t = (" << t.X() << ", " << t.Y() << ")" << endl;

    //add arc
    return ArcTo (oldp.X(), oldp.Y(), t);
  }

  auto Rectangle (double l, double w)
  {
    Line (l);
    Rotate (90);
    Line(w);
    Rotate (90);
    Line (l);
    Rotate (90);
    Line(w);
    Rotate (90);
    // wires.push_back (wire_builder.Wire());
    // wire_builder = BRepBuilderAPI_MakeWire();
    return shared_from_this();            
  }

  auto Circle(double x, double y,  double r)
  {
    
    MoveTo(x+r, y);
    Direction (0, 1);
    Arc(r, 180);
    Arc(r, 180);
    // wires.push_back (wire_builder.Wire());
    // wire_builder = BRepBuilderAPI_MakeWire();
    return shared_from_this();            

    /*

      // could not get it working with MakeCircle 

    cout << "make circle, p = " << p.X() << "/" << p.Y() << ", r = " << r << endl;
    // Handle(Geom2d_Circle) circ_curve = GCE2d_MakeCircle(p, r).Value();
    // Handle(Geom2d_Curve) curve2d = new Geom2d_TrimmedCurve (circ_curve, 0, M_PI);

    gp_Vec2d v(r,0);
    Handle(Geom2d_TrimmedCurve) curve2d = GCE2d_MakeArcOfCircle(p.Translated(v),
                                                                p.Translated(-v),
                                                                p.Translated(v)).Value();
    // Handle(Geom2d_TrimmedCurve) curve2d = GCE2d_MakeCircle(p, r).Value();

    
    auto edge = BRepBuilderAPI_MakeEdge(curve2d, surf).Edge();
    cout << "have edge, is null = " << edge.IsNull() << endl;
    wire_builder.Add(edge);
    wires.push_back (wire_builder.Wire());
    cout << "have wire, is null = " << wires.back().IsNull() << endl;
    wire_builder = BRepBuilderAPI_MakeWire();
    return shared_from_this();    
    */
  }
  
  shared_ptr<WorkPlane> Close ()
  {
    cout << "close called" << endl;

    if (startpnt.Distance(localpos.Location()) > 1e-10)
      {
        cout << "generate closing line" << endl;
        LineTo (startpnt.X(), startpnt.Y());
        return shared_from_this();                    
      }
    
    if (!startvertex.IsNull())
      {
        cout << "I am actually closing" << endl;
        wires.push_back (wire_builder.Wire());
        wire_builder = BRepBuilderAPI_MakeWire();
        startvertex.Nullify();
        cout << "complete" << endl;        
      }
    cout << "close returning" << endl;
    return shared_from_this();            
  }
  
  auto Reverse()
  {
    wires.back().Reverse();
    return shared_from_this();                
  }
  
  auto Offset(double d)
  {
    TopoDS_Wire wire = wires.back();
    wires.pop_back();
    BRepOffsetAPI_MakeOffset builder;
    builder.AddWire(wire);
    cout << "call builder" << endl;
    builder.Perform(d);
    cout << "perform is back" << endl;
    auto shape = builder.Shape();
    cout << "builder is back" << endl;
    cout << "Offset got shape type " << shape.ShapeType() << endl;
    wires.push_back (TopoDS::Wire(shape.Reversed()));
    return shared_from_this();
  }
  
  TopoDS_Wire Last()
  {
    return wires.back();
  }

  TopoDS_Face Face()
  {
    BRepBuilderAPI_MakeFace builder(surf, 1e-8);
    for (auto w : wires)
      builder.Add(w);
    return builder.Face();
  }

  auto Wires()
  {
    ListOfShapes ws;
    for (auto w : wires)
      ws.push_back(w);
    return ws;
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
        auto fused = BRepAlgoAPI_Fuse(shape1, shape2).Shape();
        // return fused;
        
        // make one face when fusing in 2D
        // from https://gitlab.onelab.info/gmsh/gmsh/-/issues/627
        int cntsolid = 0;
        for (TopExp_Explorer e(shape1, TopAbs_SOLID); e.More(); e.Next())
          cntsolid++;
        for (TopExp_Explorer e(shape2, TopAbs_SOLID); e.More(); e.Next())
          cntsolid++;
        if (cntsolid == 0)
          {
            ShapeUpgrade_UnifySameDomain unify(fused, true, true, true);
            unify.Build();
            return unify.Shape();
          }
        else
          return fused;
      })
    .def("__radd__", [] (const TopoDS_Shape & shape, int i) // for sum([shapes])
         { return shape; })
    .def("__mul__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        // return BRepAlgoAPI_Common(shape1, shape2).Shape();
        
        BRepAlgoAPI_Common builder(shape1, shape2);
#ifdef OCC_HAVE_HISTORY
        Handle(BRepTools_History) history = builder.History ();

        
        for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
          for (auto & s : { shape1, shape2 })
            for (TopExp_Explorer e(s, typ); e.More(); e.Next())
              {
                auto prop = OCCGeometry::global_shape_properties[e.Current().TShape()];
                for (auto mods : history->Modified(e.Current()))
                  OCCGeometry::global_shape_properties[mods.TShape()].Merge(prop);
              }

        
        /*
          // work in progress ...
        TopTools_ListOfShape modlist = history->Modified(shape1);
        for (auto s : modlist)
          cout << "modified from list el: " << s.ShapeType() << endl;
        */
        /*
        for (auto & s : { shape1, shape2 })
          for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
            {
              auto & prop = OCCGeometry::global_shape_properties[e.Current().TShape()];
              for (auto smod : history->Modified(e.Current()))            
                OCCGeometry::global_shape_properties[smod.TShape()].Merge(prop);
            }        
        */
#endif // OCC_HAVE_HISTORY
        
        return builder.Shape();
      })
    
    .def("__sub__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        // return BRepAlgoAPI_Cut(shape1, shape2).Shape();
        
        BRepAlgoAPI_Cut builder(shape1, shape2);
#ifdef OCC_HAVE_HISTORY        
        Handle(BRepTools_History) history = builder.History ();

        
        for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
          for (auto & s : { shape1, shape2 })
            for (TopExp_Explorer e(s, typ); e.More(); e.Next())
              {
                auto prop = OCCGeometry::global_shape_properties[e.Current().TShape()];
                for (auto mods : history->Modified(e.Current()))
                  OCCGeometry::global_shape_properties[mods.TShape()].Merge(prop);
              }
        
#ifdef OLD        
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
#endif
        
#endif // OCC_HAVE_HISTORY

        
        return builder.Shape();        
      })

    .def("Reversed", [](const TopoDS_Shape & shape) {
        return CastShape(shape.Reversed()); })

    .def("Extrude", [](const TopoDS_Shape & shape, double h) {
        for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
          {
            Handle(Geom_Surface) surf = BRep_Tool::Surface (TopoDS::Face(e.Current()));
            gp_Vec du, dv;
            gp_Pnt p;
            surf->D1 (0,0,p,du,dv);
            return BRepPrimAPI_MakePrism (shape, h*du^dv).Shape();
          }
        throw Exception("no face found for extrusion");
      })

      .def("Revolve", [](const TopoDS_Shape & shape, const gp_Ax1 &A, const double D) {
        for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
          {
            return BRepPrimAPI_MakeRevol (shape, A, D*M_PI/180).Shape();
          }
        throw Exception("no face found for revolve");
      })
    
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

    .def("Identify", [](const TopoDS_Shape & me, const TopoDS_Shape & you, string name) {
        // only edges supported, by now
        auto me_edge = TopoDS::Edge(me);
        auto you_edge = TopoDS::Edge(you);

        GProp_GProps props;
        BRepGProp::LinearProperties(me, props);
        gp_Pnt cme = props.CentreOfMass();
        BRepGProp::LinearProperties(you, props);
        gp_Pnt cyou = props.CentreOfMass();

        double s0, s1;
        auto curve_me = BRep_Tool::Curve(me_edge, s0, s1);
        auto vme = occ2ng(curve_me->Value(s1))-occ2ng(curve_me->Value(s0));
        auto curve_you = BRep_Tool::Curve(you_edge, s0, s1);
        auto vyou = occ2ng(curve_you->Value(s1))-occ2ng(curve_you->Value(s0));
        
        bool inv = vme*vyou < 0;
        OCCGeometry::identifications[me.TShape()].push_back
                              (OCCIdentification { you, Transformation<3>(occ2ng(cyou) - occ2ng(cme)), inv, name });
        OCCGeometry::identifications[you.TShape()].push_back
                              (OCCIdentification { me, Transformation<3>(occ2ng(cme) - occ2ng(cyou)), inv, name });
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
    .def("_webgui_data", [](const TopoDS_Shape & shape)
         {
           BRepTools::Clean (shape);
           double deflection = 0.01;
           BRepMesh_IncrementalMesh (shape, deflection, true);
           // triangulation = BRep_Tool::Triangulation (face, loc);

           std::vector<double> p[3];
           py::list names, colors;

           int index = 0;

           Box<3> box(Box<3>::EMPTY_BOX);
           for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
           {
               TopoDS_Face face = TopoDS::Face(e.Current());
               // Handle(TopoDS_Face) face = e.Current();
               ExtractFaceData(face, index, p, box);
               auto & props = OCCGeometry::global_shape_properties[face.TShape()];
               if(props.col)
               {
                 auto & c = *props.col;
                 colors.append(py::make_tuple(c[0], c[1], c[2]));
               }
               else
                   colors.append(py::make_tuple(0.0, 1.0, 0.0));
               if(props.name)
               {
                 names.append(*props.name);
               }
               else
                   names.append("");
               index++;
           }

           std::vector<double> edge_p[2];
           py::list edge_names, edge_colors;
           index = 0;
           for (TopExp_Explorer e(shape, TopAbs_EDGE); e.More(); e.Next())
           {
               TopoDS_Edge edge = TopoDS::Edge(e.Current());
               ExtractEdgeData(edge, index, edge_p, box);
               auto & props = OCCGeometry::global_shape_properties[edge.TShape()];
               if(props.col)
               {
                 auto & c = *props.col;
                 edge_colors.append(py::make_tuple(c[0], c[1], c[2]));
               }
               else
                   edge_colors.append(py::make_tuple(0.0, 0.0, 0.0));
               if(props.name)
               {
                 edge_names.append(*props.name);
               }
               else
                   edge_names.append("");
               index++;
           }
           
           
           auto center = box.Center();

           py::list mesh_center;
           mesh_center.append(center[0]);
           mesh_center.append(center[1]);
           mesh_center.append(center[2]);
           py::dict data;
           data["ngsolve_version"] = "Netgen x.x"; // TODO
           data["mesh_dim"] = 3; // TODO
           data["mesh_center"] = mesh_center;
           data["mesh_radius"] = box.Diam()/2;
           data["order2d"] = 1;
           data["order3d"] = 0;
           data["draw_vol"] = false;
           data["draw_surf"] = true;
           data["funcdim"] = 0;
           data["show_wireframe"] = true;
           data["show_mesh"] = true;
           data["Bezier_points"] = py::list{};
           py::list points;
           points.append(p[0]);
           points.append(p[1]);
           points.append(p[2]);
           data["Bezier_trig_points"] = points;
           data["funcmin"] = 0;
           data["funcmax"] = 1;
           data["mesh_regions_2d"] = index;
           data["autoscale"] = false;
           data["colors"] = colors;
           data["names"] = names;

           py::list edges;
           edges.append(edge_p[0]);
           edges.append(edge_p[1]);
           data["edges"] = edges;
           data["edge_names"] = edge_names;
           data["edge_colors"] = edge_colors;
           return data;
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
    .def_property_readonly("start_tangent",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           gp_Pnt p; gp_Vec v;
                           curve->D1(s0, p, v);
                           return v;
                           })
    .def_property_readonly("end_tangent",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           gp_Pnt p; gp_Vec v;
                           curve->D1(s1, p, v);
                           return v;
                           })
    ;
  py::class_<TopoDS_Wire, TopoDS_Shape> (m, "TopoDS_Wire");
  py::class_<TopoDS_Face, TopoDS_Shape> (m, "TopoDS_Face")
    .def(py::init([] (const TopoDS_Shape & shape) {
          return TopoDS::Face(shape);
        }))
    .def_property_readonly("surf", [] (TopoDS_Face face) -> Handle(Geom_Surface)
         {
           Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
           return surf;
         })
    .def("WorkPlane",[] (const TopoDS_Face & face) {
        Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
        gp_Vec du, dv;
        gp_Pnt p;
        surf->D1 (0,0,p,du,dv);
        auto ax = gp_Ax3(p, du^dv, du);
        return make_shared<WorkPlane> (ax);
      })
    ;
  py::class_<TopoDS_Solid, TopoDS_Shape> (m, "TopoDS_Solid");


  
  py::class_<Handle(Geom_Surface)> (m, "Geom_Surface")
    .def("Value", [] (const Handle(Geom_Surface) & surf, double u, double v) {
        return surf->Value(u, v); })
    .def("D1", [] (const Handle(Geom_Surface) & surf, double u, double v) {
        gp_Vec du, dv;
        gp_Pnt p;
        surf->D1 (u,v,p,du,dv);
        return tuple(p,du,dv);
      })
    
    .def("Normal", [] (const Handle(Geom_Surface) & surf, double u, double v) {
        GeomLProp_SLProps lprop(surf,u,v,1,1e-8);
        if (lprop.IsNormalDefined())
          return lprop.Normal();
        throw Exception("normal not defined");
      })
    ;
  
  
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
    .def("__len__", [](const ListOfShapes & self) { return self.size(); })
    .def("__getitem__",[](const ListOfShapes & self, string name)
         {
           ListOfShapes selected;
           for (auto s : self)
             if (auto sname = OCCGeometry::global_shape_properties[s.TShape()].name)
               if (sname == name)
                 selected.push_back(s);
           return selected;
         })
    
    .def("Sorted",[](ListOfShapes self, gp_Vec dir)
         {
           std::map<Handle(TopoDS_TShape), double> sortval;
           for (auto shape : self)
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
               sortval[shape.TShape()] = val;
             }

           std::sort (std::begin(self), std::end(self),
                      [&](TopoDS_Shape a, TopoDS_Shape b)
                      { return sortval[a.TShape()] < sortval[b.TShape()]; });
           return self;
         })
    
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

  m.def("Revolve", [] (const TopoDS_Shape & face,const gp_Ax1 &A, const double D) {
      //comvert angle from deg to rad
      return BRepPrimAPI_MakeRevol (face, A, D*M_PI/180).Shape();
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
              bool has_solid = false;
              for (TopExp_Explorer e(s, TopAbs_SOLID); e.More(); e.Next())
                {
                  builder.AddArgument(e.Current());
                  has_solid = true;
                }
              if (!has_solid)
                for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
                  builder.AddArgument(e.Current());
            }

          builder.Perform();

#ifdef OCC_HAVE_HISTORY          
          Handle(BRepTools_History) history = builder.History ();

          for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
            for (auto & s : shapes)
              for (TopExp_Explorer e(s, typ); e.More(); e.Next())
                {
                  auto prop = OCCGeometry::global_shape_properties[e.Current().TShape()];
                  for (auto mods : history->Modified(e.Current()))
                    OCCGeometry::global_shape_properties[mods.TShape()].Merge(prop);
                }
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


  m.def("BSplineCurve", [](std::vector<gp_Pnt> vpoles, int degree) {
      // not yet working ????
      TColgp_Array1OfPnt poles(0, vpoles.size()-1);
      TColStd_Array1OfReal knots(0, vpoles.size()+degree);
      TColStd_Array1OfInteger mult(0, vpoles.size()+degree);
      int cnt = 0;
      try
        {
          for (int i = 0; i < vpoles.size(); i++)
            {
              poles.SetValue(i, vpoles[i]);
              knots.SetValue(i, i);
              mult.SetValue(i,1);
            }
          for (int i = vpoles.size(); i < vpoles.size()+degree+1; i++)
            {
              knots.SetValue(i, i);
              mult.SetValue(i, 1);
            }

          Handle(Geom_Curve) curve = new Geom_BSplineCurve(poles, knots, mult, degree);
          return BRepBuilderAPI_MakeEdge(curve).Edge();
        }
      catch (Standard_Failure & e)
        {
          stringstream errstr;
          e.Print(errstr);
          throw NgException("cannot create spline: "+errstr.str());
        }
    });
  
  m.def("BezierCurve", [](std::vector<gp_Pnt> vpoles) {
      TColgp_Array1OfPnt poles(0, vpoles.size()-1);
      try
        {
          for (int i = 0; i < vpoles.size(); i++)
            poles.SetValue(i, vpoles[i]);

          Handle(Geom_Curve) curve = new Geom_BezierCurve(poles);
          return BRepBuilderAPI_MakeEdge(curve).Edge();
        }
      catch (Standard_Failure & e)
        {
          stringstream errstr;
          e.Print(errstr);
          throw NgException("cannot create Bezier-spline: "+errstr.str());
        }
    });
  
  m.def("Edge", [](Handle(Geom2d_Curve) curve2d, TopoDS_Face face) {
      auto edge = BRepBuilderAPI_MakeEdge(curve2d, BRep_Tool::Surface (face)).Edge();
      BRepLib::BuildCurves3d(edge);
      return edge;
    });
  
  m.def("Wire", [](std::vector<TopoDS_Shape> edges) {
      BRepBuilderAPI_MakeWire builder;
      try
        {
          for (auto s : edges)
            switch (s.ShapeType())
              {
              case TopAbs_EDGE:
                builder.Add(TopoDS::Edge(s)); break;
              case TopAbs_WIRE:
                builder.Add(TopoDS::Wire(s)); break;
              default:
                throw Exception("can make wire only from edges and wires");
              }
          return builder.Wire();
        }
      catch (Standard_Failure & e)
        {
          stringstream errstr;
          e.Print(errstr);
          throw NgException("error in wire builder: "+errstr.str());
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



  py::class_<WorkPlane, shared_ptr<WorkPlane>> (m, "WorkPlane")
    .def(py::init<gp_Ax3, gp_Ax2d>(), py::arg("axis")=gp_Ax3(), py::arg("pos")=gp_Ax2d())
    .def("MoveTo", &WorkPlane::MoveTo)
    .def("Direction", &WorkPlane::Direction)    
    .def("LineTo", &WorkPlane::LineTo)
    .def("ArcTo", &WorkPlane::ArcTo)
    .def("Arc", &WorkPlane::Arc)
    .def("Rotate", &WorkPlane::Rotate)
    .def("Line", [](WorkPlane&wp,double l, optional<string> name) { return wp.Line(l, name); },
         py::arg("l"), py::arg("name")=nullopt)
    .def("Line", [](WorkPlane&wp,double h,double v, optional<string> name) { return wp.Line(h,v,name); },
         py::arg("dx"), py::arg("dy"), py::arg("name")=nullopt)
    .def("Rectangle", &WorkPlane::Rectangle)
    .def("Circle", &WorkPlane::Circle)
    .def("Offset", &WorkPlane::Offset)
    .def("Reverse", &WorkPlane::Reverse)
    .def("Close", &WorkPlane::Close)
    .def("Last", &WorkPlane::Last)
    .def("Face", &WorkPlane::Face)
    .def("Wires", &WorkPlane::Wires)
    ;


  
  
}

#endif // OCCGEOMETRY
#endif // NG_PYTHON

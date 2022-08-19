#ifdef NG_PYTHON
#ifdef OCCGEOMETRY

#include <regex>

#include <general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include <meshing/python_mesh.hpp>
#include <meshing.hpp>

#include "occgeom.hpp"

#include <BOPAlgo_Builder.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgo_NormalProjection.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepLib.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <GCE2d_MakeArcOfCircle.hxx>
#include <GCE2d_MakeCircle.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeCircle.hxx>
#include <GC_MakeSegment.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2dAPI_PointsToBSpline.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntTools_Context.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax2d.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>

using namespace netgen;

void ExtractEdgeData( const TopoDS_Edge & edge, int index, std::vector<double> * p, Box<3> & box )
{
    if (BRep_Tool::Degenerated(edge)) return;

    Handle(Poly_PolygonOnTriangulation) poly;
    Handle(Poly_Triangulation) T;
    TopLoc_Location loc;
    BRep_Tool::PolygonOnTriangulation(edge, poly, T, loc);

    if (poly.IsNull())
      {
        cout << IM(2) << "no edge mesh, do my own sampling" << endl;

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
        auto p0 = occ2ng((T -> Node(poly->Nodes()(j))).Transformed(loc));
        auto p1 = occ2ng((T -> Node(poly->Nodes()(j+1))).Transformed(loc));
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

void ExtractFaceData( const TopoDS_Face & face, int index, std::vector<double> * p, std::vector<double> * n, Box<3> & box )
{
    TopLoc_Location loc;
    Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);
    Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
    BRepAdaptor_Surface sf(face, Standard_False);
    BRepLProp_SLProps prop(sf, 1, 1e-5);

    bool flip = TopAbs_REVERSED == face.Orientation();

    if (triangulation.IsNull())
      {
        cout << "pls build face triangulation before" << endl;
        return;
      }

    int ntriangles = triangulation -> NbTriangles();
    for (int j = 1; j <= ntriangles; j++)
    {
      Poly_Triangle triangle = triangulation -> Triangle(j);
        std::array<Point<3>,3> pts;
        std::array<Vec<3>,3> normals;
        for (int k = 0; k < 3; k++)
          pts[k] = occ2ng( (triangulation -> Node(triangle(k+1))).Transformed(loc) );

        for (int k = 0; k < 3; k++)
          {
            auto uv = triangulation -> UVNode(triangle(k+1));
            prop.SetParameters (uv.X(), uv.Y());
            if (prop.IsNormalDefined())
              normals[k] = occ2ng (prop.Normal());
            else
              normals[k] = Cross(pts[1]-pts[0], pts[2]-pts[0]);
          }

        if(flip)
        {
            Swap(pts[1], pts[2]);
            Swap(normals[1], normals[2]);
            for (int k = 0; k < 3; k++)
                normals[k] = -normals[k];
        }

        for (int k = 0; k < 3; k++)
        {
            box.Add(pts[k]);
            for (int d = 0; d < 3; d++)
            {
                p[k].push_back( pts[k][d] );
                n[k].push_back( normals[k][d] );
            }
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
    throw Exception("Invalid Shape type");
};


class WorkPlane : public enable_shared_from_this<WorkPlane>
{
  gp_Ax3 axes;
  gp_Ax2d localpos;
  gp_Pnt2d startpnt;
  TopoDS_Vertex lastvertex, startvertex;
  Handle(Geom_Surface) surf;
  // Geom_Plane surf;

  BRepBuilderAPI_MakeWire wire_builder;
  std::vector<TopoDS_Wire> wires;
  
public:
  
  WorkPlane (const gp_Ax3 & _axes, const gp_Ax2d _localpos = gp_Ax2d())
    : axes(_axes), localpos(_localpos) // , surf(_axis) 
  {
    // surf = GC_MakePlane (gp_Ax1(axis.Location(), axis.Direction()));
    surf = new Geom_Plane(axes);
  }


  auto Finish()
  {
    if (!startvertex.IsNull())
      {
        wires.push_back (wire_builder.Wire());
        wire_builder = BRepBuilderAPI_MakeWire();
        startvertex.Nullify();
      }
    return shared_from_this();            
  }

  auto StartPnt() const {
      return startpnt;
  }

  auto CurrentLocation() const
  {
    return localpos.Location();
  }

  auto CurrentDirection() const
  {
      return gp_Vec2d(localpos.Direction());
  }

  auto MoveTo (double h, double v)
  {
    startpnt = gp_Pnt2d(h,v);
    localpos.SetLocation(startpnt);
    startvertex.Nullify();
    return shared_from_this();
  }

  auto Move(double len)
  {
    gp_Dir2d dir = localpos.Direction();
    gp_Pnt2d oldp = localpos.Location();
    auto newp = oldp.Translated(len*dir);
    return MoveTo(newp.X(), newp.Y());
  }
  
  auto Direction (double h, double v)
  {
    localpos.SetDirection(gp_Dir2d(h,v));
    return shared_from_this();
  }
  
  auto LineTo (double h, double v, optional<string> name = nullopt)
  {
    gp_Pnt2d old2d = localpos.Location();
    gp_Pnt oldp = axes.Location() . Translated(old2d.X() * axes.XDirection() + old2d.Y() * axes.YDirection());

    // localpos.Translate (gp_Vec2d(h,v));
    localpos.SetLocation (gp_Pnt2d(h,v));
    gp_Pnt2d new2d = localpos.Location();
    gp_Pnt newp = axes.Location() . Translated(new2d.X() * axes.XDirection() + new2d.Y() * axes.YDirection());

    if (new2d.Distance(old2d) < 1e-10) return shared_from_this();    
    bool closing = new2d.Distance(startpnt) < 1e-10;

      
    cout << IM(6) << "lineto, oldp = " << occ2ng(oldp) << endl;
    cout << IM(6) << "lineto, newp = " << occ2ng(newp) << endl;
    gp_Pnt pfromsurf = surf->Value(new2d.X(), new2d.Y());
    cout << IM(6) << "p from plane = " << occ2ng(pfromsurf) << endl;
    
    Handle(Geom_TrimmedCurve) curve = GC_MakeSegment(oldp, newp);

    if (startvertex.IsNull())
      startvertex = lastvertex = BRepBuilderAPI_MakeVertex(oldp);
    auto endv = closing ? startvertex : BRepBuilderAPI_MakeVertex(newp);
    // liefert noch Fehler bei close
    auto edge = BRepBuilderAPI_MakeEdge(curve, lastvertex, endv).Edge();
    lastvertex = endv;

    // auto edge = BRepBuilderAPI_MakeEdge(curve).Edge();
    if (name)
      OCCGeometry::global_shape_properties[edge].name = name;
    wire_builder.Add(edge);

    if (closing) Finish();
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
    cout << IM(6) << "dir = " << dir.X() << ", " << dir.Y() << endl;
    gp_Pnt2d oldp = localpos.Location();
    oldp.Translate(len*dir);
    return LineTo (oldp.X(), oldp.Y(), name);
  }

  auto Rotate (double angle)
  {
    localpos.Rotate(localpos.Location(), angle*M_PI/180);
    return shared_from_this();
  }

  auto Spline(const std::vector<gp_Pnt2d> &points, bool periodic, double tol, const std::map<int, gp_Vec2d> &tangents,
              bool start_from_localpos)
  {
    gp_Pnt2d P1 = start_from_localpos ? localpos.Location() : points.front();
    gp_Pnt P13d = surf->Value(P1.X(), P1.Y());

    gp_Pnt2d PLast = points.back();
    gp_Pnt PLast3d = surf->Value(PLast.X(), PLast.Y());

    Handle(TColgp_HArray1OfPnt2d) allpoints;
    if (start_from_localpos)
      {
        if (points.front().Distance(P1) <= tol)
          throw Exception("First item of given list of points is too close to current position (distance <= tol).");
        allpoints = new TColgp_HArray1OfPnt2d(1, points.size() + 1);
        allpoints->SetValue(1, P1);
        for (int i = 0; i < points.size(); i++)
          allpoints->SetValue(i + 2, points[i]);
      }
    else
      {
        allpoints = new TColgp_HArray1OfPnt2d(1, points.size());
        for (int i = 0; i < points.size(); i++)
          allpoints->SetValue(i + 1, points[i]);
      }

    Geom2dAPI_Interpolate builder(allpoints, periodic, tol);

    if (tangents.size() > 0)
      {
        const gp_Vec2d dummy_vec = tangents.begin()->second;
        TColgp_Array1OfVec2d tangent_vecs(1, allpoints->Length());
        Handle(TColStd_HArray1OfBoolean) tangent_flags = new TColStd_HArray1OfBoolean(1, allpoints->Length());
        for (int i : Range(allpoints->Length()))
          {
            if (tangents.count(i) > 0)
              {
                tangent_vecs.SetValue(i+1, tangents.at(i));
                tangent_flags->SetValue(i+1, true);
              }
            else
              {
                tangent_vecs.SetValue(i+1, dummy_vec);
                tangent_flags->SetValue(i+1, false);
              }
          }
        builder.Load(tangent_vecs, tangent_flags);
      }


    builder.Perform();
    auto curve2d = builder.Curve();

    const bool closing = periodic || PLast.Distance(startpnt) < 1e-10;
    if (startvertex.IsNull())
        startvertex = lastvertex = BRepBuilderAPI_MakeVertex(P13d).Vertex();
    auto endv = closing ? startvertex : BRepBuilderAPI_MakeVertex(PLast3d).Vertex();

    //create 3d edge from 2d curve using surf
    auto edge = BRepBuilderAPI_MakeEdge(curve2d, surf, lastvertex, endv).Edge();
    lastvertex = endv;
    BRepLib::BuildCurves3d(edge);
    wire_builder.Add(edge);

    // update localpos
    localpos.SetLocation(PLast);

    //compute angle of rotation
    //compute tangent t2 in PLast
    const auto dir = localpos.Direction();
    gp_Vec2d t = gp_Vec2d(dir.X(), dir.Y());
    gp_Vec2d t2 = curve2d->DN(curve2d->LastParameter(), 1);

    double angle = t.Angle(t2);    //angle \in [-pi,pi]
    
    //update localpos.Direction()
    Rotate(angle*180/M_PI);

    if (closing)
        Finish();

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

    cout << IM(6) << "ArcTo:" << endl;
    cout << IM(6) << "P1 = (" << P1.X() <<", " << P1.Y() << ")"<<endl;
    cout << IM(6) << "P2 = (" << P2.X() <<", " << P2.Y() << ")"<<endl;
    cout << IM(6) << "t = (" << t.X() << ", " << t.Y() << ")" << endl;

    //compute circle center point M
    //point midway between p1 and p2
    gp_Pnt2d P12 = gp_Pnt2d((P1.X() + h) / 2, (P1.Y() + v) / 2);
    //vector normal to vector from P1 to P12
    gp_Vec2d p12n = gp_Vec2d( - (P12.Y() - P1.Y()), (P12.X() - P1.X()));
    //M is intersection of p12n and tn (tn ... normalvector to t)
    double k = ((P12.Y()- P1.Y())*p12n.X() + (P1.X() - P12.X())*p12n.Y() )/ (t.X()*p12n.X() + t.Y()*p12n.Y());
    gp_Pnt2d M = gp_Pnt2d(P1.X()-k*t.Y(), P1.Y() + k*t.X());

    cout << IM(6) << "P12 = (" << P12.X() <<", " << P12.Y() << ")"<<endl;
    cout << IM(6) << "p12n = (" << p12n.X() <<", " << p12n.Y() << ")"<<endl;
    cout << IM(6) << "k = " << k <<endl;
    cout << IM(6) << "M = (" << M.X() <<", " << M.Y() << ")"<<endl;

    //radius
    double r = P1.Distance(M);

    //compute point P3 on circle between P1 and P2
    p12n.Normalize();   //docu: reverses direction of p12n ??
    cout << IM(6) << "p12n = (" << p12n.X() <<", " << p12n.Y() << ")"<<endl;

    gp_Pnt2d P3;

    double angletp12n = t.Angle(p12n);
    if(angletp12n > -M_PI/2 && angletp12n < M_PI/2)
        P3 = gp_Pnt2d(M.X() + r * p12n.X() , M.Y() + r * p12n.Y());
    else
        P3 = gp_Pnt2d(M.X() - r * p12n.X() , M.Y() - r * p12n.Y());

    cout << IM(6) << "r = " << r <<endl;
    cout << IM(6) << "angle t,p12n = " << t.Angle(p12n)<<endl;
    cout << IM(6) << "P3 = (" << P3.X() <<", " << P3.Y() << ")"<<endl;
    cout << IM(6) << "dist(M,P3) = " << P3.Distance(M) <<endl;

    //Draw 2d arc of circle from P1 to P2 through P3
    Handle(Geom2d_TrimmedCurve) curve2d = GCE2d_MakeArcOfCircle(P1, P3, P2).Value();

    gp_Pnt P13d = surf->Value(P1.X(), P1.Y());
    gp_Pnt P23d = surf->Value(P2.X(), P2.Y());
    cout << IM(6) << "p13d = " << occ2ng(P13d) << ", p23d = " << occ2ng(P23d) << endl;
    bool closing = P2.Distance(startpnt) < 1e-10;
    if (startvertex.IsNull())
      startvertex = lastvertex = BRepBuilderAPI_MakeVertex(P13d);
    auto endv = closing ? startvertex : BRepBuilderAPI_MakeVertex(P23d);

    //create 3d edge from 2d curve using surf
    auto edge = BRepBuilderAPI_MakeEdge(curve2d, surf, lastvertex, endv).Edge();
    lastvertex = endv;
    BRepLib::BuildCurves3d(edge);
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

    //update localpos.Direction()
    Rotate(angle*180/M_PI);
    if (closing)
      Finish();

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

    cout << IM(6) << "M = (" << oldp.X() << ", " << oldp.Y() << ")" << endl;

    dirn.Rotate(newAngle-M_PI);
    oldp.Translate(radius*dirn);

    //compute tangent vector in P1
    gp_Vec2d t = gp_Vec2d(dir.X(),dir.Y());

    cout << IM(6) << "t = (" << t.X() << ", " << t.Y() << ")" << endl;

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
    return shared_from_this();            
  }

  auto RectangleCentered (double l, double w)
  {
    Move(-l/2);
    Rotate(-90);
    Move(w/2);
    Rotate(90);
    Rectangle(l,w);
    Rotate(-90);
    Move(-w/2);
    Rotate(90);
    Move(l/2);
    return shared_from_this();                
  }

  
  auto Circle(double x, double y,  double r)
  {
    /*
    MoveTo(x+r, y);
    Direction (0, 1);
    Arc(r, 180);
    Arc(r, 180);
    // wires.push_back (wire_builder.Wire());
    // wire_builder = BRepBuilderAPI_MakeWire();
    return shared_from_this();            
    */
    
    gp_Pnt2d p(x,y);
    Handle(Geom2d_Circle) circ_curve = GCE2d_MakeCircle(p, r).Value();
    
    auto edge = BRepBuilderAPI_MakeEdge(circ_curve, surf).Edge();
    BRepLib::BuildCurves3d(edge);

    wire_builder.Add(edge);
    wires.push_back (wire_builder.Wire());
    wire_builder = BRepBuilderAPI_MakeWire();
    return shared_from_this();    
  }

  auto NameVertex (string name)
  {
    if (!lastvertex.IsNull())
      OCCGeometry::global_shape_properties[lastvertex].name = name;
    return shared_from_this();
  }

  auto Circle (double r)
  {
    gp_Pnt2d pos = localpos.Location();
    return Circle (pos.X(), pos.Y(), r);
  }
  
  shared_ptr<WorkPlane> Close ()
  {
    if (startpnt.Distance(localpos.Location()) > 1e-10)
      {
        LineTo (startpnt.X(), startpnt.Y());
        return shared_from_this();                    
      }

    if (!startvertex.IsNull())
      Finish();
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
    builder.Perform(d);
    auto shape = builder.Shape();
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
    wires.clear();
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
           throw Exception ("use 'shape.type' instead of 'shape.ShapeType()'");
         }, "deprecated, use 'shape.type' instead")
    
    .def_property_readonly("type", [](const TopoDS_Shape & shape)
                           { return shape.ShapeType(); }, "returns type of shape, i.e. 'EDGE', 'FACE', ...")    
    
    .def("SubShapes", [] (const TopoDS_Shape & shape, TopAbs_ShapeEnum & type)
         {
           ListOfShapes sub;
           for (TopExp_Explorer e(shape, type); e.More(); e.Next())
             sub.push_back(e.Current());
           return sub;
         }, py::arg("type"), "returns list of sub-shapes of type 'type'")
    
    .def_property_readonly("solids", GetSolids,
            "returns all sub-shapes of type 'SOLID'")
    .def_property_readonly("faces", GetFaces,
            "returns all sub-shapes of type 'FACE'")
    .def_property_readonly("edges", GetEdges,
            "returns all sub-shapes of type 'EDGE'")
    .def_property_readonly("wires", GetWires,
                           "returns all sub-shapes of type 'WIRE'")
    .def_property_readonly("vertices", GetVertices,
            "returns all sub-shapes of type 'VERTEX'")
    .def_property_readonly("bounding_box", [] ( const TopoDS_Shape &shape )
            {
               auto box = GetBoundingBox(shape);
               return py::make_tuple( ng2occ(box.PMin()), ng2occ(box.PMax()) );
            }, "returns bounding box (pmin, pmax)")

    .def("Properties", [] (const TopoDS_Shape & shape)
         {
           GProp_GProps props;
           switch (shape.ShapeType())
             {
             case TopAbs_FACE:
             case TopAbs_SHELL:
               BRepGProp::SurfaceProperties (shape, props); break;
	     case TopAbs_SOLID:
	     case TopAbs_COMPOUND:
	     case TopAbs_COMPSOLID:
               BRepGProp::VolumeProperties (shape, props); break;
             default:
               BRepGProp::LinearProperties(shape, props);
               // throw Exception("Properties implemented only for FACE");
             }
           double mass = props.Mass();
           gp_Pnt center = props.CentreOfMass();
           return tuple( py::cast(mass), py::cast(center) );
         }, "returns tuple of shape properties, currently ('mass', 'center'")
    
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
      }, "returns center of gravity of shape")
    
    .def_property_readonly("mass", [](const TopoDS_Shape & shape) {
           GProp_GProps props;
           switch (shape.ShapeType())
             {
             case TopAbs_FACE:
             case TopAbs_SHELL:
               BRepGProp::SurfaceProperties (shape, props); break;
	     case TopAbs_SOLID:
	     case TopAbs_COMPOUND:
	     case TopAbs_COMPSOLID:
               BRepGProp::VolumeProperties (shape, props); break;
             default:
               BRepGProp::LinearProperties(shape, props);
             }
           return props.Mass();
      }, "returns mass of shape, what is length, face, or volume")

    .def("Move", [](const TopoDS_Shape & shape, const gp_Vec v)
         {
           // which one to choose ? 
           // version 1: Transoformation
           gp_Trsf trafo;
           trafo.SetTranslation(v);
           BRepBuilderAPI_Transform builder(shape, trafo, true);
           PropagateProperties(builder, shape, occ2ng(trafo));
           return CastShape(builder.Shape());
           // version 2: change location
           // ...
         }, py::arg("v"), "copy shape, and translate copy by vector 'v'")


    .def("Rotate", [](const TopoDS_Shape & shape, const gp_Ax1 ax, double ang)
         {
           gp_Trsf trafo;
           trafo.SetRotation(ax, ang*M_PI/180);            
           BRepBuilderAPI_Transform builder(shape, trafo, true);
           PropagateProperties(builder, shape, occ2ng(trafo));
           return builder.Shape();
         }, py::arg("axis"), py::arg("ang"),
         "copy shape, and rotet copy by 'ang' degrees around 'axis'")

    .def("Mirror", [] (const TopoDS_Shape & shape, const gp_Ax3 & ax)
         {
           gp_Trsf trafo;
           trafo.SetMirror(ax.Ax2());
           BRepBuilderAPI_Transform builder(shape, trafo, true);
           PropagateProperties(builder, shape, occ2ng(trafo));
           return builder.Shape();
         }, py::arg("axes"),
         "copy shape, and mirror over plane defined by 'axes'")
    
    .def("Mirror", [] (const TopoDS_Shape & shape, const gp_Ax1 & ax)
         {
           gp_Trsf trafo;
           trafo.SetMirror(ax);
           BRepBuilderAPI_Transform builder(shape, trafo, true);
           PropagateProperties(builder, shape, occ2ng(trafo));
           return builder.Shape();
         }, py::arg("axes"),
         "copy shape, and mirror around axis 'axis'")
    
    .def("Scale", [](const TopoDS_Shape & shape, const gp_Pnt p, double s)
         {
           gp_Trsf trafo;
           trafo.SetScale(p, s);
           BRepBuilderAPI_Transform builder(shape, trafo, true);
           PropagateProperties(builder, shape, occ2ng(trafo));
           return builder.Shape();
         }, py::arg("p"), py::arg("s"),
         "copy shape, and scale copy by factor 's'")

    .def("WriteStep", [](const TopoDS_Shape & shape, string & filename)
            { step_utils::WriteSTEP(shape, filename); }
         , py::arg("filename"), "export shape in STEP - format")

    .def("bc", [](const TopoDS_Shape & shape, const string & name)
         {
           for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
             OCCGeometry::global_shape_properties[e.Current()].name = name;
           return shape;
         }, py::arg("name"), "sets 'name' property for all faces of shape")

    .def("mat", [](const TopoDS_Shape & shape, const string & name)
         {
           for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next())
             OCCGeometry::global_shape_properties[e.Current()].name = name;
           return shape;
         }, py::arg("name"), "sets 'name' property to all solids of shape")
    
    .def_property("name", [](const TopoDS_Shape & self) -> optional<string> {
        if (auto name = OCCGeometry::global_shape_properties[self].name)
          return *name;
        else
          return nullopt;
      }, [](const TopoDS_Shape & self, optional<string> name) {
        OCCGeometry::global_shape_properties[self].name = name;
      }, "'name' of shape")
    
    .def_property("maxh",
                  [](const TopoDS_Shape& self)
                  {
                    return OCCGeometry::global_shape_properties[self].maxh;
                  },
                  [](TopoDS_Shape& self, double val)
                  {
                    for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE, TopAbs_VERTEX })
                      for (TopExp_Explorer e(self, typ); e.More(); e.Next())
                      {
                        auto & maxh = OCCGeometry::global_shape_properties[e.Current()].maxh;
                        maxh = min2(val, maxh);
                      }
                  }, "maximal mesh-size for shape")
    
    .def_property("hpref",
                  [](const TopoDS_Shape& self)
                  {
                    return OCCGeometry::global_shape_properties[self].hpref;
                  },
                  [](TopoDS_Shape& self, double val)
                  {
                    auto & hpref = OCCGeometry::global_shape_properties[self].hpref;
                    hpref = max2(val, hpref);
                  }, "number of refinement levels for geometric refinement")
    
    .def_property("col", [](const TopoDS_Shape & self) {
        auto it = OCCGeometry::global_shape_properties.find(self);
        Vec<4> col(0.2, 0.2, 0.2);
        if (it != OCCGeometry::global_shape_properties.end() && it->second.col)
          col = *it->second.col;
        return std::vector<double> ( { col(0), col(1), col(2) } );
      }, [](const TopoDS_Shape & self, std::vector<double> c) {
        Vec<4> col(c[0], c[1], c[2], 1.0);
        if(c.size() == 4)
          col[3] = c[3];
        OCCGeometry::global_shape_properties[self].col = col;
      }, "color of shape as RGB - tuple")
    .def("UnifySameDomain", [](const TopoDS_Shape& shape,
                               bool edges, bool faces,
                               bool concatBSplines)
    {
      ShapeUpgrade_UnifySameDomain unify(shape, edges, faces,
                                         concatBSplines);
      unify.Build();
      Handle(BRepTools_History) history = unify.History ();
      for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
        for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
          {
            auto prop = OCCGeometry::global_shape_properties[e.Current()];
            for (auto mods : history->Modified(e.Current()))
              OCCGeometry::global_shape_properties[mods].Merge(prop);
          }
      return unify.Shape();
    }, py::arg("unifyEdges")=true, py::arg("unifyFaces")=true,
         py::arg("concatBSplines")=true)
    
    .def_property("location",
                  [](const TopoDS_Shape & shape) { return shape.Location(); },
                  [](TopoDS_Shape & shape, const TopLoc_Location & loc)
                  { shape.Location(loc); }, "Location of shape")
    .def("Located", [](const TopoDS_Shape & shape, const TopLoc_Location & loc)
         { return shape.Located(loc); }, py::arg("loc"), "copy shape and sets location of copy")

    .def("__add__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {

        BRepAlgoAPI_Fuse builder(shape1, shape2);
        PropagateProperties (builder, shape1);
        PropagateProperties (builder, shape2);
        /*
#ifdef OCC_HAVE_HISTORY
        Handle(BRepTools_History) history = builder.History ();
        
        for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
          for (auto & s : { shape1, shape2 })
            for (TopExp_Explorer e(s, typ); e.More(); e.Next())
              {
                auto prop = OCCGeometry::global_shape_properties[e.Current()];
                for (auto mods : history->Modified(e.Current()))
                  OCCGeometry::global_shape_properties[mods].Merge(prop);
              }
#endif        
        */
        auto fused = builder.Shape();        
        
        // make one face when fusing in 2D
        // from https://gitlab.onelab.info/gmsh/gmsh/-/issues/627
        // int cntsolid = 0;
        // for (TopExp_Explorer e(shape1, TopAbs_SOLID); e.More(); e.Next())
        //   cntsolid++;
        // for (TopExp_Explorer e(shape2, TopAbs_SOLID); e.More(); e.Next())
        //   cntsolid++;
        // if (cntsolid == 0)
        //   {
            ShapeUpgrade_UnifySameDomain unify(fused, true, true, true);
            unify.Build();

            // #ifdef OCC_HAVE_HISTORY
            Handle(BRepTools_History) history = unify.History ();
            
            for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
              for (TopExp_Explorer e(fused, typ); e.More(); e.Next())
                {
                  auto prop = OCCGeometry::global_shape_properties[e.Current()];
                  for (auto mods : history->Modified(e.Current()))
                    OCCGeometry::global_shape_properties[mods].Merge(prop);
                }
            // #endif        
            // PropagateProperties (unify, fused);
            
            return unify.Shape();
        //   }
        // else
        //   return fused;
      }, "fuses shapes")
    .def("__radd__", [] (const TopoDS_Shape & shape, int i) // for sum([shapes])
         { return shape; }, "needed for Sum([shapes])")
    .def("__mul__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        
        BRepAlgoAPI_Common builder(shape1, shape2);
        /*
#ifdef OCC_HAVE_HISTORY
        Handle(BRepTools_History) history = builder.History ();

        
        for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
          for (auto & s : { shape1, shape2 })
            for (TopExp_Explorer e(s, typ); e.More(); e.Next())
              {
                auto prop = OCCGeometry::global_shape_properties[e.Current()];
                for (auto mods : history->Modified(e.Current()))
                  OCCGeometry::global_shape_properties[mods].Merge(prop);
              }
#endif // OCC_HAVE_HISTORY
        */
        PropagateProperties (builder, shape1);
        PropagateProperties (builder, shape2);
        
        return builder.Shape();
      }, "common of shapes")
    
    .def("__sub__", [] (const TopoDS_Shape & shape1, const TopoDS_Shape & shape2) {
        
        BRepAlgoAPI_Cut builder(shape1, shape2);
        /*
#ifdef OCC_HAVE_HISTORY        
        Handle(BRepTools_History) history = builder.History ();
        
        for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
          for (auto & s : { shape1, shape2 })
            for (TopExp_Explorer e(s, typ); e.More(); e.Next())
              {
                auto prop = OCCGeometry::global_shape_properties[e.Current()];
                for (auto mods : history->Modified(e.Current()))
                  OCCGeometry::global_shape_properties[mods].Merge(prop);
              }
#endif // OCC_HAVE_HISTORY
        */
        PropagateProperties (builder, shape1);
        PropagateProperties (builder, shape2);
        
        return builder.Shape();        
      }, "cut of shapes")
    .def("__eq__", [] (const TopoDS_Shape& shape1, const TopoDS_Shape& shape2) {
      return shape1.IsSame(shape2);
    })
    .def("__hash__", [] (const TopoDS_Shape& shape) {
      return shape.HashCode(std::numeric_limits<Standard_Integer>::max());
    })

    .def("Reversed", [](const TopoDS_Shape & shape) {
        return CastShape(shape.Reversed()); })

    .def("Extrude", [](const TopoDS_Shape & shape, double h,
                       optional<gp_Vec> dir) {
        for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
          {
            Handle(Geom_Surface) surf = BRep_Tool::Surface (TopoDS::Face(e.Current()));
            gp_Vec edir;
            if(dir.has_value())
              edir = *dir;
            else
              {
                gp_Vec du, dv;
                gp_Pnt p;
                surf->D1 (0,0,p,du,dv);
                edir = du^dv;
              }
            BRepPrimAPI_MakePrism builder(shape, h*edir, true);

            for (auto typ : { TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX })
              for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
                {
                  auto prop = OCCGeometry::global_shape_properties[e.Current()];
                  for (auto mods : builder.Generated(e.Current()))
                    OCCGeometry::global_shape_properties[mods].Merge(prop);
                }
            
            return builder.Shape();
          }
        throw Exception("no face found for extrusion");
    }, py::arg("h"), py::arg("dir")=nullopt, "extrude shape to thickness 'h', shape must contain a plane surface, optionally give an extrusion direction")
    
    .def("Extrude", [] (const TopoDS_Shape & face, gp_Vec vec) {
        return BRepPrimAPI_MakePrism (face, vec).Shape();
      }, py::arg("v"), "extrude shape by vector 'v'")

  .def("Revolve", [](const TopoDS_Shape & shape, const gp_Ax1 &A, const double D) {
      // for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
        {
          // return BRepPrimAPI_MakeRevol (shape, A, D*M_PI/180).Shape();
          BRepPrimAPI_MakeRevol builder(shape, A, D*M_PI/180, true);
            
          for (auto typ : { TopAbs_EDGE, TopAbs_VERTEX })
            for (TopExp_Explorer e(shape, typ); e.More(); e.Next())
              {
                auto prop = OCCGeometry::global_shape_properties[e.Current()];
                for (auto mods : builder.Generated(e.Current()))
                  OCCGeometry::global_shape_properties[mods].Merge(prop);
              }

          return builder.Shape();          
        }
        // throw Exception("no face found for revolve");
    }, py::arg("axis"), py::arg("ang"), "revolve shape around 'axis' by 'ang' degrees")
    
    .def("MakeFillet", [](const TopoDS_Shape & shape, std::vector<TopoDS_Shape> edges, double r) {
        BRepFilletAPI_MakeFillet mkFillet(shape);
        for (auto e : edges)
          mkFillet.Add (r, TopoDS::Edge(e));
        mkFillet.Build();
        PropagateProperties (mkFillet, shape);
        for (auto e : edges)
          for (auto gen : mkFillet.Generated(e))
            OCCGeometry::global_shape_properties[gen].name = "fillet";
        return mkFillet.Shape();
      }, py::arg("edges"), py::arg("r"), "make fillets for edges 'edges' of radius 'r'")
  
    .def("MakeChamfer", [](const TopoDS_Shape & shape, std::vector<TopoDS_Shape> edges, double d) {
#if OCC_VERSION_MAJOR>=7 && OCC_VERSION_MINOR>=4        
        BRepFilletAPI_MakeChamfer mkChamfer(shape);
        for (auto e : edges)
          mkChamfer.Add (d, TopoDS::Edge(e));
        mkChamfer.Build();
        PropagateProperties (mkChamfer, shape);
        for (auto e : edges)
          for (auto gen : mkChamfer.Generated(e))
            OCCGeometry::global_shape_properties[gen].name = "chamfer";
        return mkChamfer.Shape();
#else
        throw Exception("MakeChamfer not available for occ-version < 7.4");
#endif        
      }, py::arg("edges"), py::arg("d"), "make symmetric chamfer for edges 'edges' of distrance 'd'")
  
    .def("MakeThickSolid", [](const TopoDS_Shape & body, std::vector<TopoDS_Shape> facestoremove,
                              double offset, double tol) {
           TopTools_ListOfShape faces;
           for (auto f : facestoremove)
             faces.Append(f);
           
           BRepOffsetAPI_MakeThickSolid maker;
           maker.MakeThickSolidByJoin(body, faces, offset, tol);
           return maker.Shape();
         }, py::arg("facestoremove"), py::arg("offset"), py::arg("tol"), "makes shell-like solid from faces")
    
    .def("MakeTriangulation", [](const TopoDS_Shape & shape)
         {
           BRepTools::Clean (shape);
           double deflection = 0.01;
           BRepMesh_IncrementalMesh (shape, deflection, true);
         })


    .def("Identify", py::overload_cast<const TopoDS_Shape &, const TopoDS_Shape &, string, Identifications::ID_TYPE, std::optional<std::variant<gp_Trsf, gp_GTrsf>>>(&Identify),
            py::arg("other"), py::arg("name"),
            py::arg("type")=Identifications::PERIODIC, py::arg("trafo")=nullopt,
            "Identify shapes for periodic meshing")

    .def("Distance", [](const TopoDS_Shape& self,
                        const TopoDS_Shape& other)
    {
      return BRepExtrema_DistShapeShape(self, other).Value();
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
               Poly_Triangle triangle = triangulation -> Triangle(j);
               std::array<Point<3>,3> pts;
               for (int k = 0; k < 3; k++)
                 pts[k] = occ2ng( (triangulation -> Node(triangle(k+1))).Transformed(loc) );
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
           std::vector<double> n[3];
           py::list names, colors, solid_names;
           std::vector<std::vector<int>> solid_face_map;

           int index = 0;

           Box<3> box(Box<3>::EMPTY_BOX);
           TopTools_IndexedMapOfShape fmap;
           for (TopExp_Explorer e(shape, TopAbs_FACE); e.More(); e.Next())
           {
               TopoDS_Face face = TopoDS::Face(e.Current());
               if(fmap.Contains(face)) continue;
               // Handle(TopoDS_Face) face = e.Current();
               fmap.Add(face);
               ExtractFaceData(face, index, p, n, box);
               auto & props = OCCGeometry::global_shape_properties[face];
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

           for(auto& solid : GetSolids(shape))
             {
               std::vector<int> faces;
               for(auto& face : GetFaces(solid))
                 faces.push_back(fmap.FindIndex(face)-1);
               solid_face_map.push_back(move(faces));
               auto& props = OCCGeometry::global_shape_properties[solid];
               if(props.name)
                 solid_names.append(*props.name);
               else
                 solid_names.append("");
             }

           std::vector<double> edge_p[2];
           py::list edge_names, edge_colors;
           index = 0;
           for (TopExp_Explorer e(shape, TopAbs_EDGE); e.More(); e.Next())
           {
               TopoDS_Edge edge = TopoDS::Edge(e.Current());
               ExtractEdgeData(edge, index, edge_p, box);
               auto & props = OCCGeometry::global_shape_properties[edge];
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
           data["have_normals"] = true;
           data["show_wireframe"] = true;
           data["show_mesh"] = true;
           data["Bezier_points"] = py::list{};
           py::list points;
           points.append(p[0]);
           points.append(p[1]);
           points.append(p[2]);
           points.append(n[0]);
           points.append(n[1]);
           points.append(n[2]);
           data["Bezier_trig_points"] = points;
           data["funcmin"] = 0;
           data["funcmax"] = 1;
           data["mesh_regions_2d"] = index;
           data["autoscale"] = false;
           data["colors"] = colors;
           data["names"] = names;
           data["solid_names"] = solid_names;

           py::list edges;
           edges.append(edge_p[0]);
           edges.append(edge_p[1]);
           data["edges"] = edges;
           data["edge_names"] = edge_names;
           data["edge_colors"] = edge_colors;
           data["solid_face_map"] = solid_face_map;
           return data;
         })
    ;
  
  py::class_<TopoDS_Vertex, TopoDS_Shape> (m, "Vertex")
    .def(py::init([] (const TopoDS_Shape & shape) {
          return TopoDS::Vertex(shape);
        }))
    .def(py::init([] (const gp_Pnt & p) {
          return BRepBuilderAPI_MakeVertex (p).Vertex();
        }))
    .def_property_readonly("p", [] (const TopoDS_Vertex & v) -> gp_Pnt {
        return BRep_Tool::Pnt (v); }, "coordinates of vertex")
    ;
  
  py::class_<TopoDS_Edge, TopoDS_Shape> (m, "Edge")
    .def(py::init([] (const TopoDS_Shape & shape) {
          return TopoDS::Edge(shape);
        }))
    .def(py::init([] (Handle(Geom2d_Curve) curve2d, TopoDS_Face face) {
          auto edge = BRepBuilderAPI_MakeEdge(curve2d, BRep_Tool::Surface (face)).Edge();
          BRepLib::BuildCurves3d(edge);
          return edge;
        }))
    .def("Value", [](const TopoDS_Edge & e, double s) {
        double s0, s1;
        auto curve = BRep_Tool::Curve(e, s0, s1);
        return curve->Value(s);        
      }, py::arg("s"), "evaluate curve for parameters 's'")
    
    .def("Tangent", [](const TopoDS_Edge & e, double s) {
        gp_Pnt p; gp_Vec v;
        double s0, s1;
        auto curve = BRep_Tool::Curve(e, s0, s1);
        curve->D1(s, p, v);
        return v;
      }, py::arg("s"), "tangent vector to curve at parameter 's'")
    
    .def_property_readonly("start",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           return curve->Value(s0);
                           },
                           "start-point of curve")
    .def_property_readonly("end",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           return curve->Value(s1);
                           },
                           "end-point of curve")
    .def_property_readonly("start_tangent",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           gp_Pnt p; gp_Vec v;
                           curve->D1(s0, p, v);
                           return v;
                           },
                           "tangent at start-point")
    .def_property_readonly("end_tangent",
                           [](const TopoDS_Edge & e) {
                           double s0, s1;
                           auto curve = BRep_Tool::Curve(e, s0, s1);
                           gp_Pnt p; gp_Vec v;
                           curve->D1(s1, p, v);
                           return v;
                           },
                           "tangent at end-point")
    .def_property_readonly("parameter_interval",
                           [](const TopoDS_Edge & e) {
                             double s0, s1;
                             auto curve = BRep_Tool::Curve(e, s0, s1);
                             return tuple(s0, s1);
                           },
                           "parameter interval of curve")
    
    
    .def("Split", [](const TopoDS_Edge& self, py::args args)
    {
      ListOfShapes new_edges;
      double s0, s1;
      auto curve = BRep_Tool::Curve(self, s0, s1);
      double tstart, t, dist;
      TopoDS_Vertex vstart, vend;
      vstart = TopExp::FirstVertex(self);
      IntTools_Context context;
      tstart = s0;
      for(auto arg : args)
        {
          if(py::isinstance<py::float_>(arg))
            t = s0 + py::cast<double>(arg) * (s1-s0);
          else
            {
              auto p = py::cast<gp_Pnt>(arg);
              auto result = context.ComputePE(p, 0., self, t, dist);
              if(result != 0)
                throw Exception("Error in finding splitting points on edge!");
            }
          auto p = curve->Value(t);
          vend = BRepBuilderAPI_MakeVertex(p);
          auto newE = TopoDS::Edge(self.EmptyCopied());
          BOPTools_AlgoTools::MakeSplitEdge(self, vstart, tstart, vend, t, newE);
          new_edges.push_back(newE);
          vstart = vend;
          tstart = t;
        }
      auto newE = TopoDS::Edge(self.EmptyCopied());
      t = s1;
      vend = TopExp::LastVertex(self);
      BOPTools_AlgoTools::MakeSplitEdge(self, vstart, tstart, vend, t, newE);
      new_edges.push_back(newE);
      return new_edges;
    }, "Splits edge at given parameters. Parameters can either be floating values in (0,1), then edge parametrization is used. Or it can be points, then the projection of these points are used for splitting the edge.")
    ;
  
  py::class_<TopoDS_Wire, TopoDS_Shape> (m, "Wire")
    .def(py::init([](const TopoDS_Edge & edge) {
          BRepBuilderAPI_MakeWire builder;
          builder.Add(edge); 
          return builder.Wire();
        }))
    .def(py::init([](std::vector<TopoDS_Shape> edges) {
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
        }))
    ;

  py::class_<TopoDS_Face, TopoDS_Shape> (m, "Face")
    .def(py::init([](TopoDS_Wire wire) {
          return BRepBuilderAPI_MakeFace(wire).Face();
        }), py::arg("w"))
    .def(py::init([](const TopoDS_Face & face, const TopoDS_Wire & wire) {
          return BRepBuilderAPI_MakeFace(BRep_Tool::Surface (face), wire).Face();
        }), py::arg("f"), py::arg("w"))
    .def(py::init([](const TopoDS_Face & face, std::vector<TopoDS_Wire> wires) {
          auto surf = BRep_Tool::Surface (face);
          BRepBuilderAPI_MakeFace builder(surf, 1e-8);
          for (auto w : wires)
            builder.Add(w);
          return builder.Face();
        }), py::arg("f"), py::arg("w"))
    .def(py::init([] (const TopoDS_Shape & shape) {
          return TopoDS::Face(shape);
        }))
    .def_property("quad_dominated", [](const TopoDS_Face& self) -> optional<bool>
                  {
                    return OCCGeometry::global_shape_properties[self].quad_dominated;
                  },
                  [](TopoDS_Face& self, optional<bool> quad_dominated)
                  {
                    OCCGeometry::global_shape_properties[self].quad_dominated = quad_dominated;
                  })
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
    .def("ProjectWire", [](const TopoDS_Face& face,
                           const TopoDS_Wire& wire)
    {
      BRepAlgo_NormalProjection builder(face);
      builder.Add(wire);
      builder.Build();
      return builder.Projection();
    })
    ;
  py::class_<TopoDS_Solid, TopoDS_Shape> (m, "Solid")
    .def(py::init([](const TopoDS_Shape& faces)
    {
      BRep_Builder builder;
      TopoDS_Shell shell;
      builder.MakeShell(shell);
      for(auto& face : GetFaces(faces))
        builder.Add(shell, face);
      TopoDS_Solid solid;
      builder.MakeSolid(solid);
      builder.Add(solid, shell);
      return solid;
    }), "Create solid from shell. Shell must consist of topologically closed faces (share vertices and edges).")
    ;
  
  py::class_<TopoDS_Compound, TopoDS_Shape> (m, "Compound")
    .def(py::init([](std::vector<TopoDS_Shape> shapes, bool separate_layers) {
          BRep_Builder builder;
          TopoDS_Compound comp;
          builder.MakeCompound(comp);
          for(auto i : Range(shapes.size()))
          {
            builder.Add(comp, shapes[i]);
            if(separate_layers)
            {
                auto & props = OCCGeometry::global_shape_properties;
                for(auto & s : GetSolids(shapes[i]))
                    props[s].layer = i+1;
                for(auto & s : GetFaces(shapes[i]))
                    props[s].layer = i+1;
                for(auto & s : GetEdges(shapes[i]))
                    props[s].layer = i+1;
                for(auto & s : GetVertices(shapes[i]))
                    props[s].layer = i+1;
            }
          }

          return comp;
        }), py::arg("shapes"), py::arg("separate_layers")=false)
    ;


  
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
  py::implicitly_convertible<TopoDS_Edge, TopoDS_Wire>();
  


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
    .def(py::init<vector<TopoDS_Shape>>())
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
           std::regex pattern(name);
           for (auto s : self)
             if (auto sname = OCCGeometry::global_shape_properties[s].name)
               if (std::regex_match(*sname, pattern))
                 selected.push_back(s);
           return selected;
         }, "returns list of all shapes named 'name'")

    .def("__getitem__",[](const ListOfShapes & self, DirectionalInterval interval)
         {
           ListOfShapes selected;
           for (auto s : self)
             if (interval.Contains(Center(s), GetBoundingBox(s).Diam() * 1e-7))
               selected.push_back(s);
           return selected;
         })
    .def_property_readonly("solids", &ListOfShapes::Solids)
    .def_property_readonly("faces", &ListOfShapes::Faces)
    .def_property_readonly("wires", &ListOfShapes::Wires)
    .def_property_readonly("edges", &ListOfShapes::Edges)
    .def_property_readonly("vertices", &ListOfShapes::Vertices)
    .def(py::self * py::self)

    .def("Sorted",[](ListOfShapes self, gp_Vec dir)
         {
           std::map<TopoDS_Shape, double> sortval;
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
               sortval[shape] = val;
             }

           std::sort (std::begin(self), std::end(self),
                      [&](const TopoDS_Shape& a, const TopoDS_Shape& b)
                      { return sortval[a] < sortval[b]; });
           return self;
         }, py::arg("dir"), "returns list of shapes, where center of gravity is sorted in direction of 'dir'")
    
    .def("Max", [] (ListOfShapes & shapes, gp_Vec dir)
         { return CastShape(shapes.Max(dir)); },
         py::arg("dir"), "returns shape where center of gravity is maximal in the direction 'dir'")
    
    .def("Min", [] (ListOfShapes & shapes, gp_Vec dir) 
         { return CastShape(shapes.Max(-dir)); },
         py::arg("dir"), "returns shape where center of gravity is minimal in the direction 'dir'")

    .def("Nearest", [] (ListOfShapes & shapes, gp_Pnt pnt) 
         { return CastShape(shapes.Nearest(pnt)); },
         py::arg("p"), "returns shape nearest to point 'p'")
    .def("Nearest", [] (ListOfShapes & shapes, gp_Pnt2d pnt) 
         { return CastShape(shapes.Nearest( { pnt.X(), pnt.Y(), 0 })); },
         py::arg("p"), "returns shape nearest to point 'p'")
    
    .def_property("name", [](ListOfShapes& shapes)
    {
      throw Exception("Cannot get property of ListOfShapes, get the property from individual shapes!");
    },
      [](ListOfShapes& shapes, optional<std::string> name)
      {
        for(auto& shape : shapes)
          {
            OCCGeometry::global_shape_properties[shape].name = name;
          }
      }, "set name for all elements of list")
    .def_property("col", [](ListOfShapes& shapes) {
        throw Exception("Cannot get property of ListOfShapes, get the property from individual shapes!");
      }, [](ListOfShapes& shapes, std::vector<double> c) {
        Vec<4> col(c[0], c[1], c[2], 1.0);
        if(c.size() == 4)
          col[3] = c[3];
        for(auto& shape : shapes)
          OCCGeometry::global_shape_properties[shape].col = col;
      }, "set col for all elements of list")
    
    .def_property("maxh", [](ListOfShapes& shapes)
    {
      throw Exception("Cannot get property of ListOfShapes, get the property from individual shapes!");
    },
      [](ListOfShapes& shapes, double maxh)
      {
        for(auto& shape : shapes)
          {
            for(auto& s : GetSolids(shape))
              OCCGeometry::global_shape_properties[s].maxh = maxh;
            for(auto& s : GetFaces(shape))
              OCCGeometry::global_shape_properties[s].maxh = maxh;
            for(auto& s : GetEdges(shape))
              OCCGeometry::global_shape_properties[s].maxh = maxh;
            for(auto& s : GetVertices(shape))
              OCCGeometry::global_shape_properties[s].maxh = maxh;
          }
      }, "set maxh for all elements of list")
    .def_property("hpref", [](ListOfShapes& shapes)
    {
      throw Exception("Cannot get property of ListOfShapes, get the property from individual shapes!");
    },
      [](ListOfShapes& shapes, double hpref)
      {
        for(auto& shape : shapes)
          {
            auto& val = OCCGeometry::global_shape_properties[shape].hpref;
            val = max2(hpref, val);
          }
      }, "set hpref for all elements of list")
    .def_property("quad_dominated", [](ListOfShapes& shapes)
                  {
                    throw Exception("Cannot get property of ListOfShapes, get the property from individual shapes!");
                  },
                  [](ListOfShapes& shapes, optional<bool> quad_dominated)
                  {
                    for(auto& shape : shapes)
                      OCCGeometry::global_shape_properties[shape].quad_dominated = quad_dominated;
                  })
    
    .def("Identify", [](const ListOfShapes& me,
                        const ListOfShapes& other,
                        string name,
                        Identifications::ID_TYPE type,
                        std::variant<gp_Trsf, gp_GTrsf> trafo)
    {
      Identify(me, other, name, type, occ2ng(trafo));
    }, py::arg("other"), py::arg("name"),
         py::arg("type")=Identifications::PERIODIC, py::arg("trafo"),
         "Identify shapes for periodic meshing")

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
    .def("Edge", [](Handle(Geom2d_Curve) curve) {
        // static Geom_Plane surf{gp_Ax3()}; // crashes in nbconvert ???
        static auto surf = new Geom_Plane{gp_Ax3()};          
        auto edge = BRepBuilderAPI_MakeEdge(curve, surf).Edge();
        BRepLib::BuildCurves3d(edge);
        return edge;
      })
    .def("Wire", [](Handle(Geom2d_Curve) curve) {
        // static Geom_Plane surf{gp_Ax3()}; // crashes in nbconvert ???
        static auto surf = new Geom_Plane{gp_Ax3()};          
        auto edge = BRepBuilderAPI_MakeEdge(curve, surf).Edge();
        BRepLib::BuildCurves3d(edge);
        return BRepBuilderAPI_MakeWire(edge).Wire();                
      })
    .def("Face", [](Handle(Geom2d_Curve) curve) {
        // static Geom_Plane surf{gp_Ax3()};  // crashes in nbconvert ???
        static auto surf = new Geom_Plane{gp_Ax3()};
        auto edge = BRepBuilderAPI_MakeEdge(curve, surf).Edge();
        BRepLib::BuildCurves3d(edge);        
        auto wire = BRepBuilderAPI_MakeWire(edge).Wire();        
        return BRepBuilderAPI_MakeFace(wire).Face();
      })
    ;


  py::enum_<GeomAbs_Shape>(m, "ShapeContinuity", "Wrapper for OCC enum GeomAbs_Shape")
    .value("C0", GeomAbs_Shape::GeomAbs_C0)
    .value("C1", GeomAbs_Shape::GeomAbs_C1)
    .value("C2", GeomAbs_Shape::GeomAbs_C2)
    .value("C3", GeomAbs_Shape::GeomAbs_C3)
    .value("CN", GeomAbs_Shape::GeomAbs_CN)
    .value("G1", GeomAbs_Shape::GeomAbs_G1)
    .value("G2", GeomAbs_Shape::GeomAbs_G2);

  py::enum_<Approx_ParametrizationType>(m, "ApproxParamType", "Wrapper for Approx_ParametrizationType")
    .value("Centripetal", Approx_ParametrizationType::Approx_Centripetal)
    .value("ChordLength", Approx_ParametrizationType::Approx_ChordLength)
    .value("IsoParametric", Approx_ParametrizationType::Approx_IsoParametric);


  m.def("HalfSpace", [] (gp_Pnt p, gp_Vec n)
  {
    gp_Pln plane(p, n);
    BRepBuilderAPI_MakeFace bface(plane);
    auto face = bface.Face();
    auto refpnt = p.Translated(-n);
    BRepPrimAPI_MakeHalfSpace builder(face, refpnt);
    return builder.Shape();
  }, py::arg("p"), py::arg("n"), "Create a half space threw point p normal to n");
  m.def("Sphere", [] (gp_Pnt cc, double r) {
      return BRepPrimAPI_MakeSphere (cc, r).Solid();
    }, py::arg("c"), py::arg("r"), "create sphere with center 'c' and radius 'r'");
  
  m.def("Cylinder", [] (gp_Pnt cpnt, gp_Dir cdir, double r, double h,
                        optional<string> bot, optional<string> top, optional<string> mantle) {
    auto builder = BRepPrimAPI_MakeCylinder (gp_Ax2(cpnt, cdir), r, h);
    if(mantle)
      OCCGeometry::global_shape_properties[builder.Face()].name = *mantle;
    auto pyshape = py::cast(builder.Solid());
    gp_Vec v = cdir;
    if(bot)
      pyshape.attr("faces").attr("Min")(v).attr("name") = *bot;
    if(top)
      pyshape.attr("faces").attr("Max")(v).attr("name") = *top;
    return pyshape;
    }, py::arg("p"), py::arg("d"), py::arg("r"), py::arg("h"),
        py::arg("bottom") = nullopt, py::arg("top") = nullopt,
        py::arg("mantle") = nullopt,
    "create cylinder with base point 'p', axis direction 'd', radius 'r', and height 'h'");
  
  m.def("Cylinder", [] (gp_Ax2 ax, double r, double h) {
      return BRepPrimAPI_MakeCylinder (ax, r, h).Solid();
    }, py::arg("axis"), py::arg("r"), py::arg("h"),
    "create cylinder given by axis, radius and height");
  
  m.def("Cone", [] (gp_Ax2 ax, double r1, double r2, double h, double angle) {
     return BRepPrimAPI_MakeCone (ax, r1, r2, h, angle).Solid();
    }, py::arg("axis"), py::arg("r1"), py::arg("r2"), py::arg("h"), py::arg("angle"),
    "create cone given by axis, radius at bottom (z=0) r1, radius at top (z=h) r2, height and angle");

  m.def("Box", [] (gp_Pnt cp1, gp_Pnt cp2) {
      return BRepPrimAPI_MakeBox (cp1, cp2).Solid();
    }, py::arg("p1"), py::arg("p2"),
    "create box with opposite points 'p1' and 'p2'");

  m.def("Prism", [] (const TopoDS_Shape & face, gp_Vec vec) {
      return BRepPrimAPI_MakePrism (face, vec, true).Shape();
    }, py::arg("face"), py::arg("v"),
    "extrude face along the vector 'v'");

  m.def("Revolve", [] (const TopoDS_Shape & face,const gp_Ax1 &A, const double D) {
      //convert angle from deg to rad
      return BRepPrimAPI_MakeRevol (face, A, D*M_PI/180, true).Shape();
    });

  m.def("Pipe", [] (const TopoDS_Wire & spine, const TopoDS_Shape & profile,
                    optional<tuple<gp_Pnt, double>> twist,
                    optional<TopoDS_Wire> auxspine) {
          if (twist)
            {
              auto [pnt, angle] = *twist;

              /*
                cyl = Cylinder((0,0,0), Z, r=1, h=1).faces[0]
                heli = Edge(Segment((0,0), (2*math.pi, 1)), cyl)
                auxspine = Wire( [heli] )
                
                Handle(Geom_Surface) cyl = new Geom_CylindricalSurface (gp_Ax3(pnt, gp_Vec(0,0,1)), 1);
                auto edge = BRepBuilderAPI_MakeEdge(curve2d, cyl).Edge();
                BRepLib::BuildCurves3d(edge);
              */              
              throw Exception("twist not implemented");
            }
          if (auxspine)
            {
              BRepOffsetAPI_MakePipeShell builder(spine);
              builder.SetMode (*auxspine, Standard_True);
              for (TopExp_Explorer e(profile, TopAbs_WIRE); e.More(); e.Next())
                builder.Add (TopoDS::Wire(e.Current()));
              builder.Build();
              builder.MakeSolid();
              return builder.Shape();
            }
          
          return BRepOffsetAPI_MakePipe (spine, profile).Shape();
        }, py::arg("spine"), py::arg("profile"), py::arg("twist")=nullopt, py::arg("auxspine")=nullopt);
  
  m.def("PipeShell", [] (const TopoDS_Wire & spine, const TopoDS_Shape & profile, const TopoDS_Wire & auxspine) {
      try
        {
          BRepOffsetAPI_MakePipeShell builder(spine);
          builder.SetMode (auxspine, Standard_True);
          builder.Add (profile);
          // builder.Build();
          // builder.MakeSolid();
          return builder.Shape();
        }
      catch (Standard_Failure & e)
        {
          stringstream errstr;
          e.Print(errstr);
          throw NgException("cannot create PipeShell: "+errstr.str());
        }
    }, py::arg("spine"), py::arg("profile"), py::arg("auxspine"));


  // Handle(Geom2d_Ellipse) anEllipse1 = new Geom2d_Ellipse(anAx2d, aMajor, aMinor);
  m.def("Ellipse", [] (const gp_Ax2d & ax, double major, double minor) -> Handle(Geom2d_Curve)
        {
          return new Geom2d_Ellipse(ax, major, minor);
        }, py::arg("axes"), py::arg("major"), py::arg("minor"), "create 2d ellipse curve");
  
  m.def("Segment", [](gp_Pnt2d p1, gp_Pnt2d p2) -> Handle(Geom2d_Curve) {
      return Handle(Geom2d_TrimmedCurve)(GCE2d_MakeSegment(p1, p2));   
      /*
      Handle(Geom2d_TrimmedCurve) curve = GCE2d_MakeSegment(p1, p2);
      return curve;
      */
    }, py::arg("p1"), py::arg("p2"), "create 2d line curve");
  
  m.def("Circle", [](gp_Pnt2d p1, double r) -> Handle(Geom2d_Curve) {
      return Handle(Geom2d_Circle)(GCE2d_MakeCircle(p1, r));
      /*
      Handle(Geom2d_Circle) curve = GCE2d_MakeCircle(p1, r);
      return curve;
      */
    }, py::arg("c"), py::arg("r"), "create 2d circle curve");

  m.def("SplineApproximation", [](const std::vector<gp_Pnt2d> &points, Approx_ParametrizationType approx_type, int deg_min,
                                    int deg_max, GeomAbs_Shape continuity, double tol) -> Handle(Geom2d_Curve) {
      TColgp_Array1OfPnt2d hpoints(0, 0);
      hpoints.Resize(0, points.size() - 1, true);
      for (int i = 0; i < points.size(); i++)
          hpoints.SetValue(i, points[i]);

      Geom2dAPI_PointsToBSpline builder(hpoints, approx_type, deg_min, deg_max, continuity, tol);
      return Handle(Geom2d_BSplineCurve)(builder.Curve());
    },
    py::arg("points"),
    py::arg("approx_type") = Approx_ParametrizationType::Approx_ChordLength,
    py::arg("deg_min") = 3,
    py::arg("deg_max") = 8,
    py::arg("continuity") = GeomAbs_Shape::GeomAbs_C2,
    py::arg("tol")=1e-8,
    R"delimiter(
Generate a piecewise continuous spline-curve approximating a list of points in 2d.

Parameters
----------

points : List|Tuple[gp_Pnt2d]
  List (or tuple) of gp_Pnt.

approx_type : ApproxParamType
  Assumption on location of parameters wrt points.

deg_min : int
  Minimum polynomial degree of splines

deg_max : int
  Maximum polynomial degree of splines

continuity : ShapeContinuity
  Continuity requirement on the approximating surface

tol : float
  Tolerance for the distance from individual points to the approximating curve.

)delimiter");

  m.def("SplineInterpolation", [](const std::vector<gp_Pnt2d> &points, bool periodic, double tol, const std::map<int, gp_Vec2d> &tangents) -> Handle(Geom2d_Curve) {
      Handle(TColgp_HArray1OfPnt2d) hpoints = new TColgp_HArray1OfPnt2d(1, points.size());
      for (int i = 0; i < points.size(); i++)
          hpoints->SetValue(i+1, points[i]);
      Geom2dAPI_Interpolate builder(hpoints, periodic, tol);

      if (tangents.size() > 0)
      {
        const gp_Vec2d dummy_vec = tangents.begin()->second;
        TColgp_Array1OfVec2d tangent_vecs(1, points.size());
        Handle(TColStd_HArray1OfBoolean) tangent_flags = new TColStd_HArray1OfBoolean(1, points.size());
        for (int i : Range(points.size()))
        {
          if (tangents.count(i) > 0)
          {
              tangent_vecs.SetValue(i+1, tangents.at(i));
              tangent_flags->SetValue(i+1, true);
          } else{
              tangent_vecs.SetValue(i+1, dummy_vec);
              tangent_flags->SetValue(i+1, false);
          }
        }
        builder.Load(tangent_vecs, tangent_flags);
      }

      builder.Perform();
      return Handle(Geom2d_BSplineCurve)(builder.Curve());
    },
    py::arg("points"),
    py::arg("periodic")=false,
    py::arg("tol")=1e-8,
    py::arg("tangents")=std::map<int, gp_Vec2d>{},
    R"delimiter(
Generate a piecewise continuous spline-curve interpolating a list of points in 2d.

Parameters
----------

points : List|Tuple[gp_Pnt2d]
  List (or tuple) of gp_Pnt2d.

periodic : bool
  Whether the result should be periodic

tol : float
  Tolerance for the distance between points.

tangents : Dict[int, gp_Vec2d]
  Tangent vectors for the points indicated by the key value (0-based).

)delimiter");
  
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
              if (has_solid) continue;
              
              bool has_face = false;
              for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
                {
                  builder.AddArgument(e.Current());
                  has_face = true;
                }
              if (has_face) continue;

              bool has_edge = false;
              for (TopExp_Explorer e(s, TopAbs_EDGE); e.More(); e.Next())
                {
                  builder.AddArgument(e.Current());
                  has_edge = true;
                }
              if (has_edge) continue;

              
              for (TopExp_Explorer e(s, TopAbs_VERTEX); e.More(); e.Next())
                {
                  builder.AddArgument(e.Current());
                }
            }

          builder.Perform();
          
          /*
#ifdef OCC_HAVE_HISTORY          
          Handle(BRepTools_History) history = builder.History ();

          for (auto typ : { TopAbs_SOLID, TopAbs_FACE,  TopAbs_EDGE })
            for (auto & s : shapes)
              for (TopExp_Explorer e(s, typ); e.More(); e.Next())
                {
                  auto prop = OCCGeometry::global_shape_properties[e.Current()];
                  for (auto mods : history->Modified(e.Current()))
                    OCCGeometry::global_shape_properties[mods].Merge(prop);
                }
#endif // OCC_HAVE_HISTORY
          */
          for (auto & s : shapes)
            PropagateProperties (builder, s);          
          return builder.Shape();
        }, py::arg("shapes"), "glue together shapes of list");

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

          /*
#ifdef OCC_HAVE_HISTORY
          Handle(BRepTools_History) history = builder.History ();

          for (TopExp_Explorer e(shape, TopAbs_SOLID); e.More(); e.Next())
            {
              auto prop = OCCGeometry::global_shape_properties[e.Current()];
              for (auto mods : history->Modified(e.Current()))
                OCCGeometry::global_shape_properties[mods].Merge(prop);
            }
#endif // OCC_HAVE_HISTORY
          */
          PropagateProperties (builder, shape);
          
          return builder.Shape();
        }, py::arg("shape"), "glue together shapes from shape, typically a compound");
  m.def("Fuse", [](const vector<TopoDS_Shape>& shapes) -> TopoDS_Shape
  {
    auto s = shapes[0];
    for(auto i : Range(size_t(1), shapes.size()))
      {
        BRepAlgoAPI_Fuse builder(s, shapes[i]);
        PropagateProperties(builder, s);
        PropagateProperties(builder, shapes[i]);
        s = builder.Shape();
      }
    return s;
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
    }, py::arg("p1"), py::arg("p2"), py::arg("p3"),
    "create arc from p1 through p2 to p3");
  
  m.def("ArcOfCircle", [](gp_Pnt p1, gp_Vec v, gp_Pnt p2) { 
      Handle(Geom_TrimmedCurve) curve = GC_MakeArcOfCircle(p1, v, p2);
      return BRepBuilderAPI_MakeEdge(curve).Edge();
    }, py::arg("p1"), py::arg("v"), py::arg("p2"),
    "create arc from p1, with tangent vector v, to point p2");


  m.def("BSplineCurve", [](std::vector<gp_Pnt> vpoles, int degree) {
      // not yet working ????
      TColgp_Array1OfPnt poles(0, vpoles.size()-1);
      TColStd_Array1OfReal knots(0, vpoles.size()+degree);
      TColStd_Array1OfInteger mult(0, vpoles.size()+degree);
      int cnt = 0;

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
    });
  
  m.def("BezierCurve", [](std::vector<gp_Pnt> vpoles) {
      TColgp_Array1OfPnt poles(0, vpoles.size()-1);

      for (int i = 0; i < vpoles.size(); i++)
        poles.SetValue(i, vpoles[i]);
      
      Handle(Geom_Curve) curve = new Geom_BezierCurve(poles);
      return BRepBuilderAPI_MakeEdge(curve).Edge();
    }, py::arg("points"), "create Bezier curve");


  m.def("SplineApproximation", [](const std::vector<gp_Pnt> &points, Approx_ParametrizationType approx_type, int deg_min,
          int deg_max, GeomAbs_Shape continuity, double tol) {
      TColgp_Array1OfPnt hpoints(0, 0);
      hpoints.Resize(0, points.size() - 1, true);
      for (int i = 0; i < points.size(); i++)
        hpoints.SetValue(i, points[i]);

      GeomAPI_PointsToBSpline builder(hpoints, approx_type, deg_min, deg_max, continuity, tol);
      return BRepBuilderAPI_MakeEdge(builder.Curve()).Edge();
    },
    py::arg("points"),
    py::arg("approx_type") = Approx_ParametrizationType::Approx_ChordLength,
    py::arg("deg_min") = 3,
    py::arg("deg_max") = 8,
    py::arg("continuity") = GeomAbs_Shape::GeomAbs_C2,
    py::arg("tol")=1e-8,
    R"delimiter(
Generate a piecewise continuous spline-curve approximating a list of points in 3d.

Parameters
----------

points : List[gp_Pnt] or Tuple[gp_Pnt]
  List (or tuple) of gp_Pnt.

approx_type : ApproxParamType
  Assumption on location of parameters wrt points.

deg_min : int
  Minimum polynomial degree of splines

deg_max : int
  Maxmium polynomial degree of splines

continuity : ShapeContinuity
  Continuity requirement on the approximating surface

tol : float
  Tolerance for the distance from individual points to the approximating curve.

)delimiter");

    m.def("SplineInterpolation", [](const std::vector<gp_Pnt> &points, bool periodic, double tol, const std::map<int, gp_Vec> &tangents) {
        Handle(TColgp_HArray1OfPnt) hpoints = new TColgp_HArray1OfPnt(1, points.size());
        for (int i = 0; i < points.size(); i++)
          hpoints->SetValue(i+1, points[i]);

        GeomAPI_Interpolate builder(hpoints, periodic, tol);

        if (tangents.size() > 0)
        {
          const gp_Vec dummy_vec = tangents.begin()->second;
          TColgp_Array1OfVec tangent_vecs(1, points.size());
          Handle(TColStd_HArray1OfBoolean) tangent_flags = new TColStd_HArray1OfBoolean(1, points.size());
          for (int i : Range(points.size()))
          {
            if (tangents.count(i) > 0)
            {
              tangent_vecs.SetValue(i+1, tangents.at(i));
              tangent_flags->SetValue(i+1, true);
            } else{
              tangent_vecs.SetValue(i+1, dummy_vec);
              tangent_flags->SetValue(i+1, false);
            }
          }
          builder.Load(tangent_vecs, tangent_flags);
        }

        builder.Perform();
        return BRepBuilderAPI_MakeEdge(builder.Curve()).Edge();
      },
      py::arg("points"),
      py::arg("periodic")=false,
      py::arg("tol")=1e-8,
      py::arg("tangents")=std::map<int, gp_Vec>{},
      R"delimiter(
Generate a piecewise continuous spline-curve interpolating a list of points in 3d.

Parameters
----------

points : List|Tuple[gp_Pnt]
  List (or tuple) of gp_Pnt

periodic : bool
  Whether the result should be periodic

tol : float
  Tolerance for the distance between points.

tangents : Dict[int, gp_Vec]
  Tangent vectors for the points indicated by the key value (0-based).

)delimiter");


  m.def("SplineSurfaceApproximation", [](py::array_t<double> pnt_array,
          Approx_ParametrizationType approx_type, int deg_min, int deg_max, GeomAbs_Shape continuity, double tol,
          bool periodic, double degen_tol) {
      if (pnt_array.ndim() != 3)
        throw Exception("`points` array must have dimension 3.");
      if (pnt_array.shape(2) != 3)
        throw Exception("The third dimension must have size 3.");

      auto array = py::extract<py::array_t<double>>(pnt_array)();
      TColgp_Array2OfPnt points(1, pnt_array.shape(0), 1, pnt_array.shape(1));
      auto pnts_unchecked = pnt_array.unchecked<3>();
      for (int i = 0; i < pnt_array.shape(0); ++i)
        for (int j = 0; j < pnt_array.shape(1); ++j)
          points.SetValue(i+1, j+1, gp_Pnt(pnts_unchecked(i, j, 0), pnts_unchecked(i, j, 1), pnts_unchecked(i, j, 2)));

      GeomAPI_PointsToBSplineSurface builder;
#if OCC_VERSION_MAJOR>=7 && OCC_VERSION_MINOR>=4
      builder.Init(points, approx_type, deg_min, deg_max, continuity, tol, periodic);
#else
      if(periodic)
          throw Exception("periodic not supported");
      builder.Init(points, approx_type, deg_min, deg_max, continuity, tol);
#endif
      return BRepBuilderAPI_MakeFace(builder.Surface(), tol).Face();
    },
    py::arg("points"),
    py::arg("approx_type") = Approx_ParametrizationType::Approx_ChordLength,
    py::arg("deg_min") = 3,
    py::arg("deg_max") = 8,
    py::arg("continuity") = GeomAbs_Shape::GeomAbs_C2,
    py::arg("tol") = 1e-3,
    py::arg("periodic") = false,
    py::arg("degen_tol") = 1e-8,
    R"delimiter(
Generate a piecewise continuous spline-surface approximating an array of points.

Parameters
----------

points : np.ndarray
  Array of points coordinates. The first dimension corresponds to the first surface coordinate point
  index, the second dimension to the second surface coordinate point index. The third dimension refers to physical
  coordinates. Such an array can be generated with code like::

      px, py = np.meshgrid(*[np.linspace(0, 1, N)]*2)
      points = np.array([[(px[i,j], py[i,j], px[i,j]*py[i,j]**2) for j in range(N)] for i in range(N)])

approx_type : ApproxParamType
  Assumption on location of parameters wrt points.

deg_min : int
  Minimum polynomial degree of splines

deg_max : int
  Maxmium polynomial degree of splines

continuity : ShapeContinuity
  Continuity requirement on the approximating surface

tol : float
  Tolerance for the distance from individual points to the approximating surface.

periodic : bool
  Whether the result should be periodic in the first surface parameter

degen_tol : double
  Tolerance for resolution of degenerate edges

)delimiter");

    m.def("SplineSurfaceInterpolation", [](
            py::array_t<double> pnt_array, Approx_ParametrizationType approx_type, bool periodic, double degen_tol) {

          if (pnt_array.ndim() != 3)
              throw Exception("`points` array must have dimension 3.");
          if (pnt_array.shape(2) != 3)
              throw Exception("The third dimension must have size 3.");

          auto array = py::extract<py::array_t<double>>(pnt_array)();
          TColgp_Array2OfPnt points(1, pnt_array.shape(0), 1, pnt_array.shape(1));
          auto pnts_unchecked = pnt_array.unchecked<3>();
          for (int i = 0; i < pnt_array.shape(0); ++i)
              for (int j = 0; j < pnt_array.shape(1); ++j)
                  points.SetValue(i+1, j+1, gp_Pnt(pnts_unchecked(i, j, 0), pnts_unchecked(i, j, 1), pnts_unchecked(i, j, 2)));

          GeomAPI_PointsToBSplineSurface builder;
#if OCC_VERSION_MAJOR>=7 && OCC_VERSION_MINOR>=4
          builder.Interpolate(points, approx_type, periodic);
#else
          if(periodic)
              throw Exception("periodic not supported");
          builder.Interpolate(points, approx_type);
#endif
          return BRepBuilderAPI_MakeFace(builder.Surface(), degen_tol).Face();
      },
      py::arg("points"),
      py::arg("approx_type") = Approx_ParametrizationType::Approx_ChordLength,
      py::arg("periodic") = false,
      py::arg("degen_tol") = 1e-8,
      R"delimiter(
Generate a piecewise continuous spline-surface interpolating an array of points.

Parameters
----------

points : np.ndarray
  Array of points coordinates. The first dimension corresponds to the first surface coordinate point
  index, the second dimension to the second surface coordinate point index. The third dimension refers to physical
  coordinates. Such an array can be generated with code like::

      px, py = np.meshgrid(*[np.linspace(0, 1, N)]*2)
      points = np.array([[(px[i,j], py[i,j], px[i,j]*py[i,j]**2) for j in range(N)] for i in range(N)])

approx_type : ApproxParamType
  Assumption on location of parameters wrt points.

periodic : bool
  Whether the result should be periodic in the first surface parameter

degen_tol : double
  Tolerance for resolution of degenerate edges

)delimiter");


  m.def("MakeFillet", [](TopoDS_Shape shape, std::vector<TopoDS_Shape> edges, double r) {
      throw Exception("call 'shape.MakeFilled'");
      BRepFilletAPI_MakeFillet mkFillet(shape);
      for (auto e : edges)
        mkFillet.Add (r, TopoDS::Edge(e));
      return mkFillet.Shape();
    }, "deprecated, use 'shape.MakeFillet'");

  m.def("MakeThickSolid", [](TopoDS_Shape body, std::vector<TopoDS_Shape> facestoremove,
                             double offset, double tol) {
          throw Exception("call 'shape.MakeThickSolid'");
          TopTools_ListOfShape faces;
          for (auto f : facestoremove)
            faces.Append(f);
          
          BRepOffsetAPI_MakeThickSolid maker;
          maker.MakeThickSolidByJoin(body, faces, offset, tol);
          return maker.Shape();
        }, "deprecated, use 'shape.MakeThickSolid'");

  m.def("ThruSections", [](std::vector<TopoDS_Shape> wires, bool solid)
        {
          BRepOffsetAPI_ThruSections aTool(solid); // Standard_True);
          for (auto shape : wires)
            aTool.AddWire(TopoDS::Wire(shape));
          aTool.CheckCompatibility(Standard_False);
          return aTool.Shape();
        }, py::arg("wires"), py::arg("solid")=true,
        "Building a loft. This is a shell or solid passing through a set of sections (wires). "
        "First and last sections may be vertices. See https://dev.opencascade.org/doc/refman/html/class_b_rep_offset_a_p_i___thru_sections.html#details");

  m.def("ConnectEdgesToWires", [](const vector<TopoDS_Shape>& edges,
                                  double tol, bool shared)
  {
    Handle(TopTools_HSequenceOfShape) sedges = new TopTools_HSequenceOfShape;
    Handle(TopTools_HSequenceOfShape) swires = new TopTools_HSequenceOfShape;
    for(auto& e : edges)
      sedges->Append(e);
    ShapeAnalysis_FreeBounds::ConnectEdgesToWires(sedges, tol, shared, swires);
    vector<TopoDS_Wire> wires;
    for(auto& w : *swires)
      wires.push_back(TopoDS::Wire(w));
    return wires;
  }, py::arg("edges"), py::arg("tol")=1e-8, py::arg("shared")=true);

  py::class_<WorkPlane, shared_ptr<WorkPlane>> (m, "WorkPlane")
    .def(py::init<gp_Ax3, gp_Ax2d>(), py::arg("axes")=gp_Ax3(), py::arg("pos")=gp_Ax2d())
    .def_property_readonly("cur_loc", &WorkPlane::CurrentLocation)
    .def_property_readonly("cur_dir", &WorkPlane::CurrentDirection)
    .def_property_readonly("start_pnt", &WorkPlane::StartPnt)
    .def("MoveTo", &WorkPlane::MoveTo, py::arg("h"), py::arg("v"), "moveto (h,v), and start new wire")
    .def("Move", &WorkPlane::Move, py::arg("l"), "move 'l' from current position and direction, start new wire")
    .def("Direction", &WorkPlane::Direction, py::arg("dirh"), py::arg("dirv"), "reset direction to (dirh, dirv)")    
    // .def("LineTo", &WorkPlane::LineTo)
    .def("LineTo", [](WorkPlane&wp, double x, double y, optional<string> name) { return wp.LineTo(x, y, name); },
         py::arg("h"), py::arg("v"), py::arg("name")=nullopt, "draw line to position (h,v)")
    .def("ArcTo", &WorkPlane::ArcTo)
    .def("Arc", &WorkPlane::Arc, py::arg("r"), py::arg("ang"), "draw arc tangential to current pos/dir, of radius 'r' and angle 'ang', draw to the left/right if ang is positive/negative")
    .def("Rotate", &WorkPlane::Rotate, py::arg("ang"), "rotate current direction by 'ang' degrees")
    .def("Line", [](WorkPlane&wp,double l, optional<string> name) { return wp.Line(l, name); },
         py::arg("l"), py::arg("name")=nullopt)
    .def("Line", [](WorkPlane&wp,double h,double v, optional<string> name) { return wp.Line(h,v,name); },
         py::arg("dx"), py::arg("dy"), py::arg("name")=nullopt)
    .def("Spline", &WorkPlane::Spline, py::arg("points"), py::arg("periodic")=false, py::arg("tol")=1e-8,
         py::arg("tangents")=std::map<int, gp_Vec2d>{}, py::arg("start_from_localpos")=true,
         "draw spline (default: starting from current position, which is implicitly added to given list of points), tangents can be specified for each point (0 refers to starting point)")
    .def("Rectangle", &WorkPlane::Rectangle, py::arg("l"), py::arg("w"), "draw rectangle, with current position as corner, use current direction")
    .def("RectangleC", &WorkPlane::RectangleCentered, py::arg("l"), py::arg("w"), "draw rectangle, with current position as center, use current direction")
    .def("Circle", [](WorkPlane&wp, double x, double y, double r) {
        return wp.Circle(x,y,r); }, py::arg("h"), py::arg("v"), py::arg("r"), "draw circle with center (h,v) and radius 'r'")
    .def("Circle", [](WorkPlane&wp, double r) { return wp.Circle(r); }, py::arg("r"), "draw circle with center in current position")
    .def("NameVertex", &WorkPlane::NameVertex, py::arg("name"), "name vertex at current position")
    .def("Offset", &WorkPlane::Offset, py::arg("d"), "replace current wire by offset curve of distance 'd'")
    .def("Reverse", &WorkPlane::Reverse, "revert orientation of current wire")
    .def("Close", &WorkPlane::Close, "draw line to start point of wire, and finish wire")
    .def("Finish", &WorkPlane::Finish, "finish current wire without closing")
    .def("Last", &WorkPlane::Last, "(deprecated) returns current wire")
    .def("Wire", &WorkPlane::Last, "returns current wire")
    .def("Face", &WorkPlane::Face, "generate and return face of all wires, resets list of wires")
    .def("Wires", &WorkPlane::Wires, "returns all wires")
    ;
}

#endif // OCCGEOMETRY
#endif // NG_PYTHON

#ifdef NG_PYTHON

#include "../general/ngpython.hpp"
#include "../core/python_ngcore.hpp"
#include "csg.hpp"
#include "../meshing/python_mesh.hpp"
#include "../general/gzstream.h"


using namespace netgen;
using namespace pybind11::literals;

namespace netgen
{
  extern shared_ptr<NetgenGeometry> ng_geometry;
}



// a shadow solid tree using shared pointers.

class SPSolid
{
  shared_ptr<SPSolid> s1, s2;
  Solid * solid;
  int bc = -1;
  string bcname = "";
  double maxh = -1;
  string material;
  bool owner;
  double red = 0, green = 0, blue = 1;
  bool transp = false;
public:
  enum optyp { TERM, SECTION, UNION, SUB, EXISTING };

  SPSolid (Solid * as) : solid(as), owner(true), op(TERM) { ; }
  SPSolid (Solid * as, int /*dummy*/)
    : solid(as), owner(false), op(EXISTING) { ; }
  ~SPSolid () 
  {
    ; // if (owner) delete solid;
  }  
  SPSolid (optyp aop, shared_ptr<SPSolid> as1, shared_ptr<SPSolid> as2) 
    : s1(as1), s2(as2), owner(true), op(aop) 
  { 
    if (aop == UNION)
      solid = new Solid (Solid::UNION, s1->GetSolid(), s2->GetSolid());
    else if (aop == SECTION)
      solid = new Solid (Solid::SECTION, s1->GetSolid(), s2->GetSolid());
    else if (aop == SUB)
      solid = new Solid (Solid::SUB, s1->GetSolid()); // , s2->GetSolid());
  }

  Solid * GetSolid() { return solid; }
  const Solid * GetSolid() const { return solid; }

  void GiveUpOwner() 
  { 
    owner = false; 
    if (s1) s1 -> GiveUpOwner();
    if (s2) s2 -> GiveUpOwner();
  }

  void AddSurfaces(CSGeometry & geom)
  {
    if (op == TERM)
      geom.AddSurfaces (solid->GetPrimitive());
    if (s1) s1 -> AddSurfaces (geom);
    if (s2) s2 -> AddSurfaces (geom);
  }

  void SetMaterial (string mat)  { material = mat; }

  string GetMaterial ()
  {
    if (!material.empty()) return material;
    if (s1)
      {
        string s1mat = s1->GetMaterial();
        if (!s1mat.empty()) return s1mat;
      }
    if (s2)
      {
        string s2mat = s2->GetMaterial();
        if (!s2mat.empty()) return s2mat;
      }
    return material;
  }

  void SetBC(int abc) 
  {
    if (bc == -1) 
      {
        bc = abc;
        if (s1) s1 -> SetBC(bc);
        if (s2) s2 -> SetBC(bc);
        if (op == TERM)
          {
            Primitive * prim = solid -> GetPrimitive();
            for (int i = 0; i < prim->GetNSurfaces(); i++)
              prim->GetSurface(i).SetBCProperty (abc);
            // cout << "set " << prim->GetNSurfaces() << " surfaces to bc " << bc << endl;
          }
      }
  }

  void SetBCName(string name) 
  {
    if (bcname == "") 
      {
        bcname = name;
        if (s1) s1 -> SetBCName(name);
        if (s2) s2 -> SetBCName(name);
        if (op == TERM)
          {
            Primitive * prim = solid -> GetPrimitive();
            for (int i = 0; i < prim->GetNSurfaces(); i++)
              prim->GetSurface(i).SetBCName (name);
            // cout << "set " << prim->GetNSurfaces() << " surfaces to bc " << bc << endl;
          }
      }
  }



  void SetMaxH(double amaxh) 
  {
    if (maxh == -1) 
      {
        maxh = amaxh;
        if (s1) s1 -> SetMaxH(maxh);
        if (s2) s2 -> SetMaxH(maxh);
        if (op == TERM)
          {
            Primitive * prim = solid -> GetPrimitive();
            for (int i = 0; i < prim->GetNSurfaces(); i++)
              prim->GetSurface(i).SetMaxH (maxh);
          }
      }
  }

  void SetColor(double ared, double agreen, double ablue)
  {
    red = ared;
    green = agreen;
    blue = ablue;
  }

  double GetRed() const { return red; }
  double GetGreen() const { return green; }
  double GetBlue() const { return blue; }

  void SetTransparent() { transp = true; }
  bool IsTransparent() { return transp; }

private:
  optyp op;
};

inline ostream & operator<< (ostream & ost, const SPSolid & sol)
{
  ost << *sol.GetSolid();
  return ost;
}

namespace netgen
{
  extern CSGeometry * ParseCSG (istream & istr, CSGeometry *instance=nullptr);
}



DLL_HEADER void ExportCSG(py::module &m)
{
  py::class_<SplineGeometry<2>, shared_ptr<SplineGeometry<2>>>
    (m, "SplineCurve2d")
    .def(py::init<>())
    .def ("AddPoint", FunctionPointer
          ([] (SplineGeometry<2> & self, double x, double y)
           {
             self.geompoints.Append (GeomPoint<2> (Point<2> (x,y)));
             return self.geompoints.Size()-1;
           }))
    .def ("AddSegment", [] (SplineGeometry<2> & self, int i1, int i2,
                            string bcname, double maxh)
    {
      self.splines.Append (new LineSeg<2> (self.geompoints[i1], self.geompoints[i2], maxh, bcname));
    }, "p1"_a, "p2"_a, "bcname"_a="default", "maxh"_a=1e99)
    .def ("AddSegment", [] (SplineGeometry<2> & self, int i1, int i2,
                            int i3, string bcname, double maxh)
    {
      self.splines.Append (new SplineSeg3<2> (self.geompoints[i1], self.geompoints[i2], self.geompoints[i3], bcname, maxh));
    }, "p1"_a, "p2"_a, "p3"_a, "bcname"_a="default", "maxh"_a=1e99)
    ;

  py::class_<SplineGeometry<3>,shared_ptr<SplineGeometry<3>>> (m,"SplineCurve3d")
    .def(py::init<>())
    .def ("AddPoint", FunctionPointer
          ([] (SplineGeometry<3> & self, double x, double y, double z)
           {
             self.geompoints.Append (GeomPoint<3> (Point<3> (x,y,z)));
             return self.geompoints.Size()-1;
           }))
    .def ("AddSegment", FunctionPointer
          ([] (SplineGeometry<3> & self, int i1, int i2)
           {
             self.splines.Append (new LineSeg<3> (self.geompoints[i1], self.geompoints[i2]));
           }))
    .def ("AddSegment", FunctionPointer
          ([] (SplineGeometry<3> & self, int i1, int i2, int i3)
           {
             self.splines.Append (new SplineSeg3<3> (self.geompoints[i1], self.geompoints[i2], self.geompoints[i3]));
           }))
    ;

  py::class_<SplineSurface, shared_ptr<SplineSurface>> (m, "SplineSurface",
                        "A surface for co dim 2 integrals on the splines")
    .def(py::init([](shared_ptr<SPSolid> base, py::list cuts)
	     {
	       auto primitive = dynamic_cast<OneSurfacePrimitive*> (base->GetSolid()->GetPrimitive());
	       auto acuts = make_shared<NgArray<shared_ptr<OneSurfacePrimitive>>>();
	       for(int i = 0; i<py::len(cuts);i++)
		 {
		   py::extract<shared_ptr<SPSolid>> sps(cuts[i]);
		   if(!sps.check())
		     throw NgException("Cut must be SurfacePrimitive in constructor of SplineSurface!");
		   auto sp = dynamic_cast<OneSurfacePrimitive*>(sps()->GetSolid()->GetPrimitive());
		   if(sp)
		     acuts->Append(shared_ptr<OneSurfacePrimitive>(sp));
		   else
		     throw Exception("Cut must be SurfacePrimitive in constructor of SplineSurface!");
		 }
	       if(!primitive)
		 throw Exception("Base is not a SurfacePrimitive in constructor of SplineSurface!");
	       return make_shared<SplineSurface>(shared_ptr<OneSurfacePrimitive>(primitive),acuts);
	     }),py::arg("base"), py::arg("cuts")=py::list())
    .def("AddPoint", FunctionPointer
	 ([] (SplineSurface & self, double x, double y, double z, bool hpref)
	  {
	    self.AppendPoint(Point<3>(x,y,z),hpref);
	    return self.GetNP()-1;
	  }),
	 py::arg("x"),py::arg("y"),py::arg("z"),py::arg("hpref")=false)
    .def("AddSegment", [] (SplineSurface & self, int i1, int i2, string bcname, double maxh)
	  {
            auto seg = make_shared<LineSeg<3>>(self.GetPoint(i1),self.GetPoint(i2));
	    self.AppendSegment(seg,bcname,maxh);
	  },
	 py::arg("pnt1"),py::arg("pnt2"),py::arg("bcname")="default", py::arg("maxh")=-1.)
    .def("AddSegment", [] (SplineSurface& self, int i1, int i2, int i3, string bcname, double maxh)
         {
           auto seg = make_shared<SplineSeg3<3>>(self.GetPoint(i1), self.GetPoint(i2), self.GetPoint(i3));
           self.AppendSegment(seg, bcname, maxh);
         }, py::arg("pnt1"),py::arg("pnt2"), py::arg("pnt3"),py::arg("bcname")="default", py::arg("maxh")=-1.)
    ;
  
  py::class_<SPSolid, shared_ptr<SPSolid>> (m, "Solid")
    .def ("__str__", &ToString<SPSolid>)
    .def ("__add__", FunctionPointer( [] ( shared_ptr<SPSolid> self, shared_ptr<SPSolid> other) { return make_shared<SPSolid> (SPSolid::UNION, self, other); }))
    .def ("__mul__", FunctionPointer( [] ( shared_ptr<SPSolid> self, shared_ptr<SPSolid> other) { return make_shared<SPSolid> (SPSolid::SECTION, self, other); }))
    .def ("__sub__", FunctionPointer
          ([] ( shared_ptr<SPSolid> self, shared_ptr<SPSolid> other) 
           { return make_shared<SPSolid> (SPSolid::SECTION, self, 
                                          make_shared<SPSolid> (SPSolid::SUB, other, nullptr)); }))

    .def ("bc", FunctionPointer([](shared_ptr<SPSolid> & self, int nr) -> shared_ptr<SPSolid> 
                                { self->SetBC(nr); return self; }))
    .def ("bc", FunctionPointer([](shared_ptr<SPSolid> & self, string name) -> shared_ptr<SPSolid> 
                                { self->SetBCName(name); return self; }))
    .def ("maxh", FunctionPointer([](shared_ptr<SPSolid> & self, double maxh) -> shared_ptr<SPSolid> 
                                { self->SetMaxH(maxh); return self; }))
    .def ("mat", FunctionPointer([](shared_ptr<SPSolid> & self, string mat) -> shared_ptr<SPSolid> 
                                 { self->SetMaterial(mat); return self; }))
    .def ("mat", &SPSolid::GetMaterial)
    .def("col", FunctionPointer([](shared_ptr<SPSolid> & self, py::list rgb) -> shared_ptr<SPSolid>
                                { 
                                  py::extract<double> red(rgb[0]);
                                  py::extract<double> green(rgb[1]);
                                  py::extract<double> blue(rgb[2]);
                                  self->SetColor(red(),green(),blue());
                                  return self; 
                                }))
    .def("transp", FunctionPointer([](shared_ptr<SPSolid> & self)->shared_ptr < SPSolid > { self->SetTransparent(); return self; }))
    ;

  m.def ("Sphere", FunctionPointer([](Point<3> c, double r)
                                     {
                                       Sphere * sp = new Sphere (c, r);
                                       Solid * sol = new Solid (sp);
                                       return make_shared<SPSolid> (sol);
                                     }));
  m.def ("Ellipsoid", FunctionPointer([](Point<3> m, Vec<3> a, Vec<3> b, Vec<3> c)
                                     {
                                       Ellipsoid * ell = new Ellipsoid (m, a, b, c);
                                       Solid * sol = new Solid (ell);
                                       return make_shared<SPSolid> (sol);
                                     }));
  m.def ("Plane", FunctionPointer([](Point<3> p, Vec<3> n)
                                    {
                                      Plane * sp = new Plane (p,n);
                                      Solid * sol = new Solid (sp);
                                      return make_shared<SPSolid> (sol);
                                    }));
  m.def ("Cone", FunctionPointer([](Point<3> a, Point<3> b, double ra, double rb)
                                       {
                                         Cone * cyl = new Cone (a, b, ra, rb);
                                         Solid * sol = new Solid (cyl);
                                         return make_shared<SPSolid> (sol);
                                       }));
  m.def ("Cylinder", FunctionPointer([](Point<3> a, Point<3> b, double r)
                                       {
                                         Cylinder * cyl = new Cylinder (a, b, r);
                                         Solid * sol = new Solid (cyl);
                                         return make_shared<SPSolid> (sol);
                                       }));
  m.def ("OrthoBrick", FunctionPointer([](Point<3> p1, Point<3> p2)
                                         {
                                           OrthoBrick * brick = new OrthoBrick (p1,p2);
                                           Solid * sol = new Solid (brick);
                                           return make_shared<SPSolid> (sol);
                                         }));
  m.def ("Torus", FunctionPointer([](Point<3> c, Vec<3> n, double R, double r)
                                         {
                                           Torus * torus = new Torus (c,n,R,r);
                                           Solid * sol = new Solid (torus);
                                           return make_shared<SPSolid> (sol);
                                         }));
  m.def ("Revolution", [](Point<3> p1, Point<3> p2,
                          shared_ptr<SplineGeometry<2>> spline)
  {
    Revolution * rev = new Revolution (p1, p2, spline);
    Solid * sol = new Solid(rev);
    return make_shared<SPSolid> (sol);
  });
  m.def ("Extrusion", [](shared_ptr<SplineGeometry<3>> path,
                         shared_ptr<SplineGeometry<2>> profile,
                         Vec<3> d)
  {
    Extrusion * extr = new Extrusion (path,profile,d);
    Solid * sol = new Solid(extr);
    return make_shared<SPSolid> (sol);
  }, py::arg("path"), py::arg("profile"), py::arg("d"),
     R"delimiter(A body of extrusion is defined by its profile
(which has to be a closed, clockwiseoriented 2D curve),
 by a path (a 3D curve) and a vector d. It is constructed
 as follows: Take a point p on the path and denote the
 (unit-)tangent of the path in this point by t. If we cut
 the body by the plane given by p and t as normal vector,
 the cut is the profile. The profile is oriented by the
 (local) y-direction `y:=d−(d·t)t` and the (local) x-direction
 `x:=t \times y`.
The following points have to be noticed:
 * If the path is not closed, then also the body is NOT closed.
   In this case e.g. planes or orthobricks have to be used to
   construct a closed body.
 * The path has to be smooth, i.e. the tangents at the end- resp.
   start-point of two consecutive spline or line patches have to
   have the same directions.
)delimiter");
  m.def("EllipticCone", [](const Point<3>& a, const Vec<3>& v, const Vec<3>& w,
                            double h, double r)
        {
          auto ellcone = new EllipticCone(a,v,w,h,r);
          auto sol = new Solid(ellcone);
          return make_shared<SPSolid>(sol);
        }, py::arg("a"), py::arg("vl"), py::arg("vs"), py::arg("h"), py::arg("r"),
        R"raw_string(
An elliptic cone, given by the point 'a' at the base of the cone along the main axis,
the vectors v and w of the long and short axis of the ellipse, respectively,
the height of the cone, h, and ratio of base long axis length to top long axis length, r

Note: The elliptic cone has to be truncated by planes similar to a cone or an elliptic cylinder.
When r =1, the truncated elliptic cone becomes an elliptic cylinder.
When r tends to zero, the truncated elliptic cone tends to a full elliptic cone.
However, when r = 0, the top part becomes a point(tip) and meshing fails!
)raw_string");

  m.def("Polyhedron", [](py::list points, py::list faces)
  {
    auto poly = new Polyhedra();
    for(auto p : points)
      poly->AddPoint(py::cast<Point<3>>(p));
    int fnr = 0;
    for(auto face : faces)
      {
        auto lface = py::cast<py::list>(face);
        if(py::len(lface) == 3)
          poly->AddFace(py::cast<int>(lface[0]),
                        py::cast<int>(lface[1]),
                        py::cast<int>(lface[2]),
                        fnr++);
        else if(py::len(lface) == 4)
          {
            poly->AddFace(py::cast<int>(lface[0]),
                          py::cast<int>(lface[1]),
                          py::cast<int>(lface[2]),
                          fnr);
            poly->AddFace(py::cast<int>(lface[0]),
                          py::cast<int>(lface[2]),
                          py::cast<int>(lface[3]),
                          fnr++);
          }
      }
    return make_shared<SPSolid>(new Solid(poly));
  });
  
  m.def ("Or", FunctionPointer([](shared_ptr<SPSolid> s1, shared_ptr<SPSolid> s2)
                                 {
                                   return make_shared<SPSolid> (SPSolid::UNION, s1, s2);
                                 }));
  m.def ("And", FunctionPointer([](shared_ptr<SPSolid> s1, shared_ptr<SPSolid> s2)
                                  {
                                    return make_shared<SPSolid> (SPSolid::SECTION, s1, s2);
                                  }));


  py::class_<CSGeometry, NetgenGeometry, shared_ptr<CSGeometry>> (m, "CSGeometry")
    .def(py::init<>())
    .def(py::init([](const string& filename)
                  {
                    ifstream ist (filename);
                    auto geo = make_shared<CSGeometry>();
                    ParseCSG(ist, geo.get());
                    geo->FindIdenticSurfaces(1e-8 * geo->MaxSize());
                    return geo;
                  }), py::arg("filename"))
    .def(NGSPickle<CSGeometry>())
    .def("Save", FunctionPointer([] (CSGeometry & self, string filename)
                                 {
                                   cout << "save geometry to file " << filename << endl;
                                   self.Save (filename);
                                 }))
    .def("Add",
         [] (CSGeometry & self, shared_ptr<SPSolid> solid, py::list bcmod, double maxh,
             py::tuple col, bool transparent, int layer)
          {
            solid->AddSurfaces (self);
            solid->GiveUpOwner();
            int tlonr = self.SetTopLevelObject (solid->GetSolid());
            self.GetTopLevelObject(tlonr) -> SetMaterial(solid->GetMaterial());
            self.GetTopLevelObject(tlonr) -> SetRGB(solid->GetRed(),solid->GetGreen(),solid->GetBlue());
            // self.GetTopLevelObject(tlonr)->SetTransparent(solid->IsTransparent());
            self.GetTopLevelObject(tlonr)->SetTransparent(transparent);
            self.GetTopLevelObject(tlonr)->SetMaxH(maxh);
            self.GetTopLevelObject(tlonr)->SetLayer(layer);

            // cout << "rgb = " << py::len(rgb) << endl;
            if (py::len(col)==3)
              self.GetTopLevelObject(tlonr) -> SetRGB(py::cast<double>(col[0]),
                                                      py::cast<double>(col[1]),
                                                      py::cast<double>(col[2]));
            
            // bcmod is list of tuples ( solid, bcnr )
            for (int i = 0; i < py::len(bcmod); i++)
              {
                py::tuple tup = py::extract<py::tuple> (bcmod[i]) ();
                auto mod_solid = py::extract<shared_ptr<SPSolid>> (tup[0]) ();
                int mod_nr = -1;
                string * bcname = nullptr;
                py::object val = tup[1];
                if (py::extract<int>(val).check()) mod_nr = py::extract<int> (val)();
                if (py::extract<string>(val).check()) bcname = new string ( py::extract<string> (val)());

                NgArray<int> si;
                mod_solid -> GetSolid() -> GetSurfaceIndices (si);
                // cout << "change bc on surfaces: " << si << " to " << mod_nr << endl;

                for (int j = 0; j < si.Size(); j++)
                  {
                    CSGeometry::BCModification bcm;
                    bcm.bcname = bcname ? new string (*bcname) : nullptr;
                    bcm.tlonr = tlonr;
                    bcm.si = si[j];
		    bcm.bcnr = mod_nr;
		    self.bcmodifications.Append (bcm);
                  }
                delete bcname;
              }
            return tlonr;
          },
         py::arg("solid"), py::arg("bcmod")=py::list(), py::arg("maxh")=1e99,
         py::arg("col")=py::tuple(), py::arg("transparent")=false, py::arg("layer")=1
         )

    .def("AddSurface", FunctionPointer
         ([] (CSGeometry & self, shared_ptr<SPSolid> surface, shared_ptr<SPSolid> solid)
          {
            solid->AddSurfaces (self);
            solid->GiveUpOwner();
            Surface & surf = surface->GetSolid()->GetPrimitive()->GetSurface();
            int tlonr = self.SetTopLevelObject (solid->GetSolid(), &surf);
            // self.GetTopLevelObject(tlonr) -> SetMaterial(solid->GetMaterial());
            self.GetTopLevelObject(tlonr) -> SetBCProp(surf.GetBCProperty());
            self.GetTopLevelObject(tlonr) -> SetBCName(surf.GetBCName());
            
            self.GetTopLevelObject(tlonr) -> SetRGB(solid->GetRed(),solid->GetGreen(),solid->GetBlue());
            self.GetTopLevelObject(tlonr)->SetTransparent(solid->IsTransparent());
          }),
         py::arg("surface"), py::arg("solid")
         )
    .def("AddSplineSurface", FunctionPointer
	 ([] (CSGeometry & self, shared_ptr<SplineSurface> surf)
	  {
	    auto cuttings = surf->CreateCuttingSurfaces();
	    auto spsol = make_shared<SPSolid>(new Solid(surf.get()));
	    for(auto cut : (*cuttings)){
	      spsol = make_shared<SPSolid>(SPSolid::SECTION,spsol,make_shared<SPSolid>(new Solid(cut.get())));
	    }
	    spsol->AddSurfaces(self);
	    int tlonr = self.SetTopLevelObject(spsol->GetSolid(), surf.get());
	    self.GetTopLevelObject(tlonr) -> SetBCProp(surf->GetBase()->GetBCProperty());
	    self.GetTopLevelObject(tlonr) -> SetBCName(surf->GetBase()->GetBCName());
	    self.GetTopLevelObject(tlonr) -> SetMaxH(surf->GetBase()->GetMaxH());
            NgArray<Point<3>> non_midpoints;
            for(auto spline : surf->GetSplines())
              {
                non_midpoints.Append(spline->GetPoint(0));
              }
	    for(auto p : non_midpoints)
		self.AddUserPoint(p);
            self.AddSplineSurface(surf);
	  }),
	  py::arg("SplineSurface"))
    .def("SingularFace", [] (CSGeometry & self, shared_ptr<SPSolid> sol, shared_ptr<SPSolid> surfaces, double factor)
         {
           int tlonum = -1;
           for (int i = 0; i < self.GetNTopLevelObjects(); i++)
             if (self.GetTopLevelObject(i)->GetSolid() == sol->GetSolid())
               tlonum = i;
           if (tlonum == -1) throw NgException("not a top-level-object");
           if (!surfaces) surfaces = sol;
           auto singface = new SingularFace(tlonum+1, surfaces->GetSolid(), factor);
           self.singfaces.Append(singface);
         }, py::arg("solid"), py::arg("surfaces")=nullptr, py::arg("factor")=0.25)
    .def("SingularEdge", [] (CSGeometry & self, shared_ptr<SPSolid> s1,shared_ptr<SPSolid> s2, double factor)
         {
           auto singedge = new SingularEdge(1, -1, self, s1->GetSolid(), s2->GetSolid(), factor);
           self.singedges.Append (singedge);
         })
    .def("SingularPoint", [] (CSGeometry & self, shared_ptr<SPSolid> s1,shared_ptr<SPSolid> s2,
                             shared_ptr<SPSolid> s3, double factor)
         {
           auto singpoint = new SingularPoint(1, s1->GetSolid(), s2->GetSolid(), s3->GetSolid(), factor);
           self.singpoints.Append (singpoint);
         })
    .def("CloseSurfaces", FunctionPointer
         ([] (CSGeometry & self, shared_ptr<SPSolid> s1, shared_ptr<SPSolid> s2, py::list aslices )
          {
            NgArray<int> si1, si2;
            s1->GetSolid()->GetSurfaceIndices (si1);
            s2->GetSolid()->GetSurfaceIndices (si2);
            Flags flags;

            try
            {
                int n = py::len(aslices);
                Array<double> slices(n);
                for(int i=0; i<n; i++)
                {
                    slices[i]= py::extract<double>(aslices[i])();
                }
                flags.SetFlag("slices", slices);
            }
            catch( py::error_already_set const & ) {
                cout << "caught python error:" << endl;
                PyErr_Print();
            }

            const TopLevelObject * domain = nullptr;
            self.AddIdentification
              (new CloseSurfaceIdentification
               (self.GetNIdentifications()+1, self,
                self.GetSurface (si1[0]), self.GetSurface (si2[0]),
                domain,
                flags));
          }),
         py::arg("solid1"), py::arg("solid2"), py::arg("slices")
         )
    .def("CloseSurfaces", FunctionPointer
         ([] (CSGeometry & self, shared_ptr<SPSolid> s1, shared_ptr<SPSolid> s2,
              int reflevels, shared_ptr<SPSolid> domain_solid)
          {
            NgArray<int> si1, si2;
            s1->GetSolid()->GetSurfaceIndices (si1);
            s2->GetSolid()->GetSurfaceIndices (si2);
            cout << "surface ids1 = " << si1 << endl;
            cout << "surface ids2 = " << si2 << endl;

            Flags flags;
            const TopLevelObject * domain = nullptr;
            if (domain_solid)
              domain = self.GetTopLevelObject(domain_solid->GetSolid());
              
            self.AddIdentification 
              (new CloseSurfaceIdentification 
               (self.GetNIdentifications()+1, self, 
                self.GetSurface (si1[0]), self.GetSurface (si2[0]),
                domain,
                flags));
          }),
         py::arg("solid1"), py::arg("solid2"), py::arg("reflevels")=2, py::arg("domain")=nullptr
         )
    
    .def("PeriodicSurfaces", FunctionPointer
         ([] (CSGeometry & self, shared_ptr<SPSolid> s1, shared_ptr<SPSolid> s2,
              Transformation<3> trafo)
          {
            NgArray<int> si1, si2;
            s1->GetSolid()->GetSurfaceIndices (si1);
            s2->GetSolid()->GetSurfaceIndices (si2);
            cout << "identify surfaces " << si1[0] << " and " << si2[0] << endl;
            self.AddIdentification 
              (new PeriodicIdentification 
               (self.GetNIdentifications()+1, self, 
                self.GetSurface (si1[0]), self.GetSurface (si2[0]),
                trafo));
          }),
         py::arg("solid1"), py::arg("solid2"),
         py::arg("trafo")=Transformation<3>(Vec<3>(0,0,0))
         )
    .def("NameEdge", [] (CSGeometry & self, shared_ptr<SPSolid> s1, shared_ptr<SPSolid> s2, string name)
         {
           Array<Surface*> surfs1, surfs2;
           s1->GetSolid()->ForEachSurface( [&surfs1] (Surface * s, bool inv) { surfs1.Append(s); });
           s2->GetSolid()->ForEachSurface( [&surfs2] (Surface * s, bool inv) { surfs2.Append(s); });
           for (auto s1 : surfs1)
             for (auto s2 : surfs2)
               self.named_edges[tuple(s1,s2)] = name;
         })
         
    .def("AddPoint", [] (CSGeometry & self, Point<3> p, variant<int,string> index) -> CSGeometry&
         {
           if (auto pint = std::get_if<int> (&index))
             self.AddUserPoint(CSGeometry::UserPoint(p, *pint));
           if (auto pstr = std::get_if<string> (&index))
             self.AddUserPoint(CSGeometry::UserPoint(p, *pstr));
           return self;
         })
    
    .def("GetTransparent", FunctionPointer
         ([] (CSGeometry & self, int tlonr)
          {
            return self.GetTopLevelObject(tlonr)->GetTransparent();
          }),
         py::arg("tlonr")
         )
    .def("SetTransparent", FunctionPointer
         ([] (CSGeometry & self, int tlonr, bool transparent)
          {
            self.GetTopLevelObject(tlonr)->SetTransparent(transparent);
          }),
         py::arg("tlonr"), py::arg("transparent")
         )

    .def("GetVisible", FunctionPointer
         ([] (CSGeometry & self, int tlonr)
          {
            return self.GetTopLevelObject(tlonr)->GetVisible();
          }),
         py::arg("tlonr")
         )
    .def("SetVisible", FunctionPointer
         ([] (CSGeometry & self, int tlonr, bool visible)
          {
            self.GetTopLevelObject(tlonr)->SetVisible(visible);
          }),
         py::arg("tlonr"), py::arg("visible")
         )
    .def("SetBoundingBox", FunctionPointer
         ([] (CSGeometry & self, Point<3> pmin, Point<3> pmax)
          {
            self.SetBoundingBox(Box<3> (pmin, pmax));
          }),
         py::arg("pmin"), py::arg("pmax")
         )
    .def("Draw", FunctionPointer
         ([] (shared_ptr<CSGeometry> self)
          {
             self->FindIdenticSurfaces(1e-8 * self->MaxSize());
             self->CalcTriangleApproximation(0.01, 20);
             ng_geometry = self;
          })
         )
    .def("GetSolids", [](CSGeometry& self)
                      {
                        py::list lst;
                        for(auto i : Range(self.GetSolids().Size()))
                          lst.append(make_shared<SPSolid>(self.GetSolids()[i], 1234));
                        return lst;
                      })
    .def_property_readonly ("ntlo", &CSGeometry::GetNTopLevelObjects)
    .def("_visualizationData", [](shared_ptr<CSGeometry> csg_geo)
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
           for (int i = 0; i < csg_geo->GetNSurf(); i++)
             {
               auto surf = csg_geo->GetSurface(i);
               surfnames.push_back(surf->GetBCName());
             }
           csg_geo->FindIdenticSurfaces(1e-8 * csg_geo->MaxSize());
           csg_geo->CalcTriangleApproximation(0.01,20);
           auto nto = csg_geo->GetNTopLevelObjects();
           size_t np = 0;
           size_t ntrig = 0;
           for (int i = 0; i < nto; i++){
             np += csg_geo->GetTriApprox(i)->GetNP();
             ntrig += csg_geo->GetTriApprox(i)->GetNT();
           }
           vertices.reserve(np*3);
           trigs.reserve(ntrig*4);
           normals.reserve(np*3);
           int offset_points = 0;
           for (int i = 0; i < nto; i++)
             {
               auto triapprox = csg_geo->GetTriApprox(i);
               for (int j = 0; j < triapprox->GetNP(); j++)
                 for(int k = 0; k < 3; k++) {
                   float val = triapprox->GetPoint(j)[k];
                   vertices.push_back(val);
                   min[k] = min2(min[k], val);
                   max[k] = max2(max[k],val);
                   normals.push_back(triapprox->GetNormal(j)[k]);
                 }
               for (int j = 0; j < triapprox->GetNT(); j++)
                 {
                   for(int k = 0; k < 3; k++)
                     trigs.push_back(triapprox->GetTriangle(j)[k]+offset_points);
                   trigs.push_back(triapprox->GetTriangle(j).SurfaceIndex());
                 }
               offset_points += triapprox->GetNP();
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
  .def("GenerateMesh", [](shared_ptr<CSGeometry> geo,
                          MeshingParameters* pars, py::kwargs kwargs)
           {
             MeshingParameters mp;
             if(pars) mp = *pars;
             {
               py::gil_scoped_acquire aq;
               CreateMPfromKwargs(mp, kwargs);
             }
             auto mesh = make_shared<Mesh>();
             SetGlobalMesh (mesh);
             mesh->SetGeometry(geo);
	     ng_geometry = geo;
             geo->FindIdenticSurfaces(1e-8 * geo->MaxSize());
             auto result = geo->GenerateMesh (mesh, mp);
             if(result != 0)
               throw Exception("Meshing failed!");
             return mesh;
           }, py::arg("mp") = nullptr,
       meshingparameter_description.c_str(),
    py::call_guard<py::gil_scoped_release>())
    ;

  m.def("Save", FunctionPointer 
          ([](const Mesh & self, const string & filename, const CSGeometry & geom)
           {
             ostream * outfile;
             if (filename.substr (filename.length()-3, 3) == ".gz")
               outfile = new ogzstream (filename.c_str());
             else
               outfile = new ofstream (filename.c_str());
             
             self.Save (*outfile);
             *outfile << endl << endl << "endmesh" << endl << endl;
             geom.SaveToMeshFile (*outfile);
             delete outfile;
           }),py::call_guard<py::gil_scoped_release>())
    ;



  m.def("ZRefinement", FunctionPointer
          ([](Mesh & mesh, CSGeometry & geom)
          {
            ZRefinementOptions opt;
            opt.minref = 5;
            ZRefinement (mesh, &geom, opt);
          }),py::call_guard<py::gil_scoped_release>())
    ;
}

PYBIND11_MODULE(libcsg, m) {
  ExportCSG(m);
}
#endif


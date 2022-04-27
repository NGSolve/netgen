#ifdef NG_PYTHON

#include "../general/ngpython.hpp"
#include "../core/python_ngcore.hpp"
#include "../meshing/python_mesh.hpp"

#include "../include/meshing.hpp"
#include "../include/geometry2d.hpp"
#include "csg2d.hpp"

using namespace netgen;
using namespace pybind11::literals;

namespace netgen
{
  extern std::shared_ptr<NetgenGeometry> ng_geometry;
}


NGCORE_API_EXPORT void ExportGeom2d(py::module &m) 
{
  py::class_<SplineSegExt, shared_ptr<SplineSegExt>>
    (m, "Spline", "Spline of a SplineGeometry object")
    .def_property("leftdom", [] (SplineSegExt& self) { return self.leftdom; },
                  [](SplineSegExt& self, int dom) { self.leftdom = dom; })
    .def_property("rightdom", [] (SplineSegExt& self) { return self.rightdom; },
                  [](SplineSegExt& self, int dom) { self.rightdom = dom; })
    .def_property_readonly("bc", [] (SplineSegExt& self) { return self.bc; })
    .def("GetNormal", [](SplineSegExt& self, double t)
                      {
                        auto tang = self.GetTangent(t).Normalize();
                        return Vec<2>(tang[1], -tang[0]);
                      })
    .def("StartPoint", [](SplineSegExt& self) { return Point<2>(self.StartPI()); })
    .def("EndPoint", [](SplineSegExt& self) { return Point<2>(self.EndPI()); })
    ;
    
  py::class_<SplineGeometry2d, NetgenGeometry, shared_ptr<SplineGeometry2d>>
    (m, "SplineGeometry",
     "a 2d boundary representation geometry model by lines and splines",
     py::multiple_inheritance())
    .def(py::init<>())
    .def(py::init([](const string& filename)
                  {
                    auto geo = make_shared<SplineGeometry2d>();
                    geo->Load(filename.c_str());
                    ng_geometry = geo;
                    return geo;
                  }))
    .def(NGSPickle<SplineGeometry2d>())
    .def("Load",&SplineGeometry2d::Load)
    .def("SetDomainLayer", &SplineGeometry2d::SetDomainLayer)
    .def("AppendPoint", FunctionPointer
         ([](SplineGeometry2d &self, double px, double py, double maxh, double hpref, string name)
          {
            Point<2> p;
            p(0) = px;
            p(1) = py;
            GeomPoint<2> gp(p);
            gp.hmax = maxh;
            gp.hpref = hpref;
            gp.name = name;
            self.geompoints.Append(gp);
            return self.geompoints.Size()-1;
	  }),
         py::arg("x"), py::arg("y"), py::arg("maxh") = 1e99, py::arg("hpref")=0, py::arg("name")="")
    .def("Append", FunctionPointer([](SplineGeometry2d &self, py::list segment, int leftdomain, int rightdomain,
                                      optional<variant<int, string>> bc, optional<int> copy, double maxh,
                                      double hpref, double hprefleft, double hprefright)
	  {
            SplineSegExt * seg;
            if(py::isinstance<py::str>(segment[0]))
              {
                auto segtype = py::cast<std::string>(segment[0]);
            
                if (segtype == "line")
                  {
                    LineSeg<2> * l = new LineSeg<2>(self.GetPoint(py::cast<int>(segment[1])),
                                                    self.GetPoint(py::cast<int>(segment[2])));
                    seg = new SplineSegExt(*l);
                  }
                else if (segtype == "spline3")
                  {
                    SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(py::cast<int>(segment[1])),
                                                             self.GetPoint(py::cast<int>(segment[2])),
                                                             self.GetPoint(py::cast<int>(segment[3])));
                    seg = new SplineSegExt(*seg3);
                  }
                else
                  throw Exception("Appended segment is not a line or a spline3");
              }
            else
              {
                if(py::len(segment) == 2)
                  {
                    auto l = new LineSeg<2>(self.GetPoint(py::cast<int>(segment[0])),
                                            self.GetPoint(py::cast<int>(segment[1])));
                    seg = new SplineSegExt(*l);
                  }
                else if(py::len(segment) == 3)
                  {
                    SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(py::cast<int>(segment[0])),
                                                             self.GetPoint(py::cast<int>(segment[1])),
                                                             self.GetPoint(py::cast<int>(segment[2])));
                    seg = new SplineSegExt(*seg3);
                  }
                else
                  throw Exception("Appended segment must either have 2 or 3 points");
              }
            seg->leftdom = leftdomain;
            seg->rightdom = rightdomain;
            seg->hmax = maxh;
            seg->hpref_left = max(hpref, hprefleft);
            seg->hpref_right = max(hpref,hprefright);
            seg->reffak = 1;
            seg->copyfrom = -1;
            if (copy.has_value())
              seg->copyfrom = *copy+1;

            if (bc.has_value())
              {
                if(auto intptr = get_if<int>(&*bc); intptr)
                  seg->bc = *intptr;
                else
                  {
                    auto bcname = get_if<string>(&*bc);
                    seg->bc = self.GetNSplines() + 1;
                    self.SetBCName(seg->bc, *bcname);
                  }
              }
            else
              seg->bc = self.GetNSplines()+1;
            self.AppendSegment(seg);
            return self.GetNSplines()-1;
	  }), py::arg("point_indices"), py::arg("leftdomain") = 1, py::arg("rightdomain") = py::int_(0),
         py::arg("bc")=nullopt, py::arg("copy")=nullopt, py::arg("maxh")=1e99,
         py::arg("hpref")=0,py::arg("hprefleft")=0,py::arg("hprefright")=0)

    
    .def("AppendSegment", FunctionPointer([](SplineGeometry2d &self, py::list point_indices, int leftdomain, int rightdomain)
                                          {
		  int npts = py::len(point_indices);
		  SplineSegExt * seg;
		  //int a = py::extract<int>(point_indices[0]);
		  if (npts == 2)
		  {
			  LineSeg<2> * l = new LineSeg<2>(self.GetPoint(py::extract<int>(point_indices[0])()), self.GetPoint(py::extract<int>(point_indices[1])()));
			  seg = new SplineSegExt(*l);
			  
		  }
		  else if (npts == 3)
		  {
			  SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(py::extract<int>(point_indices[0])()), self.GetPoint(py::extract<int>(point_indices[1])()), self.GetPoint(py::extract<int>(point_indices[2])()));
			  seg = new SplineSegExt(*seg3);

		  }
                  else
                    throw Exception("Can only append segments with 2 or 3 points!");
		  seg->leftdom = leftdomain;
		  seg->rightdom = rightdomain;
		  seg->hmax = 1e99;
		  seg->reffak = 1;
		  seg->copyfrom = -1;
		  self.AppendSegment(seg);
                  }), py::arg("point_indices"), py::arg("leftdomain") = 1, py::arg("rightdomain") = py::int_(0))


    .def("AddCurve", 
         [] (SplineGeometry2d & self, py::object func,
                         int leftdomain, int rightdomain, py::object bc, double maxh)
         {
           int n = 1000;
           NgArray<Point<2>> points;
           for (int i = 0; i <= n; i++)
             {
               double t = double(i)/n;
               py::tuple xy = func(t);
               double x = py::cast<double>(xy[0]);
               double y = py::cast<double>(xy[1]);
               points.Append (Point<2>(x,y));
             }
           auto spline = new DiscretePointsSeg<2> (points);
           SplineSegExt * spex = new SplineSegExt (*spline);
           
           spex -> leftdom = leftdomain;
           spex -> rightdom = rightdomain;
           spex->hmax = maxh;
           spex->reffak = 1;
           spex->copyfrom = -1;
           
           if (py::extract<int>(bc).check())
             spex->bc = py::extract<int>(bc)();
           else if (py::extract<string>(bc).check())
             {
               string bcname = py::extract<string>(bc)();
               spex->bc = self.GetNSplines()+1;
               self.SetBCName(spex->bc, bcname);
             }
           else
             spex->bc = self.GetNSplines()+1;

           
           self.AppendSegment (spex);
         }, py::arg("func"), py::arg("leftdomain") = 1, py::arg("rightdomain") = py::int_(0),
         py::arg("bc")=NGDummyArgument(), py::arg("maxh")=1e99,
         "Curve is given as parametrization on the interval [0,1]")
    
    .def("SetMaterial", &SplineGeometry2d::SetMaterial)
    .def("SetDomainMaxH", &SplineGeometry2d::SetDomainMaxh)

    .def("GetBCName", [](SplineGeometry2d& self, size_t index) { return self.GetBCName(index); })

    .def("GetNDomains", [](SplineGeometry2d& self) { return self.GetNDomains(); })

    .def("GetNSplines", [](SplineGeometry2d& self) { return self.splines.Size(); })
    .def("GetSpline", [](SplineGeometry2d& self, size_t index)
                      { return shared_ptr<SplineSegExt>(&self.GetSpline(index), NOOP_Deleter); },
         py::return_value_policy::reference_internal)
    .def("GetNPoints", [](SplineGeometry2d& self) { return self.GetNP(); })
    .def("GetPoint", [](SplineGeometry2d& self, size_t index) { return Point<2>(self.GetPoint(index)); })

	.def("PlotData", FunctionPointer([](SplineGeometry2d &self)
	  {
		  Box<2> box(self.GetBoundingBox());
		  double xdist = box.PMax()(0) - box.PMin()(0);
		  double ydist = box.PMax()(1) - box.PMin()(1);
		  py::tuple xlim = py::make_tuple(box.PMin()(0) - 0.1*xdist, box.PMax()(0) + 0.1*xdist);
		  py::tuple ylim = py::make_tuple(box.PMin()(1) - 0.1*ydist, box.PMax()(1) + 0.1*ydist);

		  py::list xpoints, ypoints;

		  for (int i = 0; i < self.splines.Size(); i++)
		  {
			  py::list xp, yp;
			  if (self.splines[i]->GetType().compare("line")==0)
			  {
				  GeomPoint<2> p1 = self.splines[i]->StartPI();
				  GeomPoint<2> p2 = self.splines[i]->EndPI();
				  xp.append(py::cast(p1(0)));
				  xp.append(py::cast(p2(0)));
				  yp.append(py::cast(p1(1)));
				  yp.append(py::cast(p2(1)));
			  }
			  else if (self.splines[i]->GetType().compare("spline3")==0)
			  {
				  double len = self.splines[i]->Length();
				  int n = floor(len/(0.05*min(xdist,ydist)));
				  
				  for (int j = 0; j <= n; j++)
				  {
					  GeomPoint<2> point = self.splines[i]->GetPoint(j*1./n);
					  xp.append(py::cast(point(0)));
					  yp.append(py::cast(point(1)));
				  }
			  }
			  else
			  {
				  cout << "spline is neither line nor spline3" << endl;
			  }
			  xpoints.append(xp);
			  ypoints.append(yp);
				  
		  }
		  return py::tuple(py::make_tuple(xlim, ylim, xpoints, ypoints));

	  }))
    .def("_visualizationData", [](SplineGeometry2d &self)
         {
           Box<2> box(self.GetBoundingBox());
           double xdist = box.PMax()(0) - box.PMin()(0);
           double ydist = box.PMax()(1) - box.PMin()(1);
           py::dict data;
           py::dict segment_data;
           auto min_val = py::make_tuple(box.PMin()(0), box.PMin()(1),0);
           auto max_val = py::make_tuple(box.PMax()(1),box.PMax()(1),0);
           py::list vertices;
           py::list domains;
           py::list segment_points;
           py::list segment_normals;
           py::list leftdom;
           py::list rightdom;
           int max_bcnr = 0;
           for(int i = 0; i < self.splines.Size(); i++)
             {
               std::vector<netgen::GeomPoint<2>> lst;
               if (self.splines[i]->GetType().compare("line") == 0)
                 lst = { self.splines[i]->StartPI(), self.splines[i]->EndPI() };
               else if(self.splines[i]->GetType().compare("spline3") == 0)
                 {
                   double len = self.splines[i]->Length();
                   int n = floor(len/(0.05*min(xdist,ydist)));
                   n = max(3, n);
                   lst.push_back(self.splines[i]->StartPI());
                   for (int j = 1; j < n; j++){
                     lst.push_back(self.splines[i]->GetPoint(j*1./n));
                     lst.push_back(self.splines[i]->GetPoint(j*1./n));
                   }
                   lst.push_back(self.splines[i]->EndPI());
                   }
               else
                 {
                   throw NgException("Spline is neither line nor spline3");
                 }
               for (auto point : lst)
                 {
                   for(auto val : {point(0), point(1), 0.})
                      vertices.append(val);
                    int bcnr = self.GetSpline(i).bc;
                    max_bcnr = max2(max_bcnr, bcnr);
                    domains.append(bcnr);
                    domains.append(self.GetSpline(i).leftdom);
                    domains.append(self.GetSpline(i).rightdom);
                 }

               // segment data
               auto pnt = self.splines[i]->GetPoint(0.5);
               segment_points.append(py::make_tuple(pnt(0),pnt(1)));
               auto normal = self.GetSpline(i).GetTangent(0.5);
               std::swap(normal(0),normal(1));
               normal(1) *= -1;
               normal *= 1./sqrt(normal(0) * normal(0) + normal(1)*normal(1));
               segment_normals.append(py::make_tuple(normal(0),normal(1)));
               leftdom.append(self.GetSpline(i).leftdom);
               rightdom.append(self.GetSpline(i).rightdom);
             }
           py::list bcnames;
           for (int i = 1; i<max_bcnr + 1; i++)
             bcnames.append(self.GetBCName(i));
           segment_data["midpoints"] = segment_points;
           segment_data["normals"] = segment_normals;
           segment_data["leftdom"] = leftdom;
           segment_data["rightdom"] = rightdom;
           data["segment_data"] = segment_data;
           data["vertices"] = vertices;
           data["domains"] = domains;
           data["min"] = min_val;
           data["max"] = max_val;
           data["bcnames"] = bcnames;
           return data;
         })
	.def("PointData", FunctionPointer([](SplineGeometry2d &self)
	  {
		  py::list xpoints, ypoints, pointindex;
		  
		  for (int i = 0; i < self.geompoints.Size(); i++)
		  {
			  pointindex.append(py::cast(i));
			  xpoints.append(py::cast(self.geompoints[i][0]));
			  ypoints.append(py::cast(self.geompoints[i][1]));
		  }
		  return py::tuple(py::make_tuple(xpoints, ypoints, pointindex));
		  
	  }))
	.def("SegmentData", FunctionPointer([](SplineGeometry2d &self)
	  {
		  py::list leftpoints, rightpoints, leftdom, rightdom;

		  for (int i = 0; i < self.splines.Size(); i++)
		  {
			  GeomPoint<2> point = self.splines[i]->GetPoint(0.5);
			  Vec<2> normal = self.GetSpline(i).GetTangent(0.5);
			  double temp = normal(0);
			  normal(0) = normal(1);
			  normal(1) = -temp;

			  leftdom.append(py::cast(self.GetSpline(i).leftdom));
			  rightdom.append(py::cast(self.GetSpline(i).rightdom));

			  rightpoints.append(py::make_tuple(point(0), point(1), normal(0)<0, normal(1)<0));
			  leftpoints.append(py::make_tuple(point(0), point(1), normal(0)<0, normal(1)<0));
		  }
		  return py::tuple(py::make_tuple(leftpoints, rightpoints, leftdom, rightdom));

	  }))
	.def("Print", FunctionPointer([](SplineGeometry2d &self)
	  {
		  for (int i = 0; i < self.geompoints.Size(); i++)
		  {
			  cout << i << " : " << self.geompoints[i][0] << " , " << self.geompoints[i][1] << endl;
		  }
		  //Box<2> box(self.GetBoundingBox());
		  //cout << box.PMin() << endl;
		  //cout << box.PMax() << endl;
		  cout << self.splines.Size() << endl;
		  for (int i = 0; i < self.splines.Size(); i++)
		  {
			  cout << self.splines[i]->GetType() << endl;
			  //cout << i << " : " << self.splines[i]->GetPoint(0.1) << " , " << self.splines[i]->GetPoint(0.5) << endl;
		  }
	  }))
    .def("Draw", FunctionPointer
         ([] (shared_ptr<SplineGeometry2d> self)
          {
             ng_geometry = self;
             py::module::import("netgen").attr("Redraw")();
          })
         )
    
    .def("GenerateMesh", [](shared_ptr<SplineGeometry2d> self,
                            optional<MeshingParameters> pars, py::kwargs kwargs)
		{
                  MeshingParameters mp;
                  if(pars) mp = *pars;
                  {
                    py::gil_scoped_acquire aq;
                    CreateMPfromKwargs(mp, kwargs);
                  }
		  auto mesh = make_shared<Mesh>();
                  mesh->SetGeometry(self);
                  SetGlobalMesh (mesh);
                  ng_geometry = self;
		  auto result = self->GenerateMesh(mesh, mp);
                  if(result != 0)
                    throw Exception("Meshing failed!");
		  return mesh;
                }, py::arg("mp") = nullopt,
      py::call_guard<py::gil_scoped_release>(),
      meshingparameter_description.c_str())
    .def("_SetDomainTensorMeshing", &SplineGeometry2d::SetDomainTensorMeshing)
    ;
  
  py::class_<Solid2d>(m, "Solid2d")
    .def(py::init<>())
    .def(py::init<Array<std::variant<Point<2>, EdgeInfo, PointInfo>>, std::string, std::string>(), py::arg("points"), py::arg("mat")=MAT_DEFAULT, py::arg("bc")=BC_DEFAULT)

    .def(py::self+py::self)
    .def(py::self-py::self)
    .def(py::self*py::self)
    .def(py::self+=py::self)
    .def(py::self-=py::self)
    .def(py::self*=py::self)

    .def("Mat", &Solid2d::Mat)
    .def("BC", &Solid2d::BC)
    .def("Maxh", &Solid2d::Maxh)
    .def("Layer", &Solid2d::Layer)

    .def("Copy", [](Solid2d & self) -> Solid2d { return self; })
    .def("Move", &Solid2d::Move)
    .def("Scale", static_cast<Solid2d& (Solid2d::*)(double)>(&Solid2d::Scale))
    .def("Scale", static_cast<Solid2d& (Solid2d::*)(Vec<2>)>(&Solid2d::Scale))
    .def("Rotate", &Solid2d::RotateDeg, py::arg("angle"), py::arg("center")=Point<2>{0,0})
    ;
  

  m.def("Rectangle", [](Point<2> p0, Point<2> p1, string mat, string bc, optional<string> bottom, optional<string> right, optional<string> top, optional<string> left) -> Solid2d
		  {
                      using P = Point<2>;
                      return { {
                              p0,    EdgeInfo{bottom ? *bottom : bc},
                              P{p1[0],p0[1]}, EdgeInfo {right  ? *right  : bc},
                              p1,             EdgeInfo {top    ? *top    : bc},
                              P{p0[0],p1[1]}, EdgeInfo {left   ? *left   : bc},
                             }, mat};
                  },
		  "pmin"_a, "pmax"_a, "mat"_a=MAT_DEFAULT, "bc"_a=BC_DEFAULT,
                  "bottom"_a=nullopt, "right"_a=nullopt, "top"_a=nullopt, "left"_a=nullopt
       );
  m.def("Circle", Circle, py::arg("center"), py::arg("radius"), py::arg("mat")=MAT_DEFAULT, py::arg("bc")=BC_DEFAULT);

  py::class_<CSG2d>(m, "CSG2d")
    .def(py::init<>())
    .def("GenerateSplineGeometry", &CSG2d::GenerateSplineGeometry)
    .def("Add", &CSG2d::Add)
    .def("GenerateMesh", [](CSG2d & self, optional<MeshingParameters> pars, py::kwargs kwargs)
		{
                  MeshingParameters mp;
                  if(pars) mp = *pars;
                  {
                    py::gil_scoped_acquire aq;
                    CreateMPfromKwargs(mp, kwargs);
                  }
		  auto mesh = make_shared<Mesh>();
                  auto geo = self.GenerateSplineGeometry();
                  mesh->SetGeometry(geo);
                  SetGlobalMesh (mesh);
                  ng_geometry = geo;
		  auto result = geo->GenerateMesh(mesh, mp);
                  if(result != 0)
                    throw Exception("Meshing failed!");
		  return mesh;
                }, py::arg("mp") = nullopt,
      py::call_guard<py::gil_scoped_release>(),
      meshingparameter_description.c_str())
    ;

  py::class_<EdgeInfo>(m, "EdgeInfo")
    .def(py::init<>())
    .def(py::init<const Point<2>&>(), py::arg("control_point"))
    .def(py::init<double>(), py::arg("maxh"))
    .def(py::init<string>(), py::arg("bc"))
    .def(py::init<optional<Point<2>>, double, string>(), py::arg("control_point")=nullopt, py::arg("maxh")=MAXH_DEFAULT, py::arg("bc")=BC_DEFAULT)
    ;
  py::class_<PointInfo>(m, "PointInfo")
    .def(py::init<>())
    .def(py::init<double>(), "maxh"_a)
    .def(py::init<string>(), "name"_a)
    .def(py::init<double, string>(), "maxh"_a, "name"_a)
    ;
}

PYBIND11_MODULE(libgeom2d, m) {
  ExportGeom2d(m);
}

#endif


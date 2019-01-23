#ifdef NG_PYTHON

#include <../general/ngpython.hpp>

#include <meshing.hpp>
#include <geometry2d.hpp>

using namespace netgen;

namespace netgen
{
  extern std::shared_ptr<NetgenGeometry> ng_geometry;
}


DLL_HEADER void ExportGeom2d(py::module &m) 
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
                                      py::object bc, py::object copy, double maxh, double hpref)
	  {
            py::extract<std::string> segtype(segment[0]);
            
            SplineSegExt * seg;
            if (segtype().compare("line") == 0)
              {
                py::extract<int> point_index1(segment[1]);
                py::extract<int> point_index2(segment[2]);
                //point_index1.check()
                
                LineSeg<2> * l = new LineSeg<2>(self.GetPoint(point_index1()), self.GetPoint(point_index2()));
                seg = new SplineSegExt(*l);
              }
            else if (segtype().compare("spline3") == 0)
              {
                py::extract<int> point_index1(segment[1]);
                py::extract<int> point_index2(segment[2]);
                py::extract<int> point_index3(segment[3]);
                
                SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(point_index1()), self.GetPoint(point_index2()), self.GetPoint(point_index3()));
                seg = new SplineSegExt(*seg3);
              }
            else
              {
                cout << "Appended segment is not a line or a spline3" << endl;
              }
            seg->leftdom = leftdomain;
            seg->rightdom = rightdomain;
            seg->hmax = maxh;
            seg->hpref_left = hpref;
            seg->hpref_right = hpref;
            seg->reffak = 1;
            seg->copyfrom = -1;
            if (py::extract<int>(copy).check())
              seg->copyfrom = py::extract<int>(copy)()+1;
              
            if (py::extract<int>(bc).check())
              seg->bc = py::extract<int>(bc)();
            else if (py::extract<string>(bc).check())
              {
                string bcname = py::extract<string>(bc)();
                seg->bc = self.GetNSplines()+1;
                self.SetBCName(seg->bc, bcname);
              }
            else
              seg->bc = self.GetNSplines()+1;
            self.AppendSegment(seg);
            return self.GetNSplines()-1;
	  }), py::arg("point_indices"), py::arg("leftdomain") = 1, py::arg("rightdomain") = py::int_(0),
               py::arg("bc")=NGDummyArgument(), py::arg("copy")=NGDummyArgument(), py::arg("maxh")=1e99, py::arg("hpref")=0)

    
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
		  seg->leftdom = leftdomain;
		  seg->rightdom = rightdomain;
		  seg->hmax = 1e99;
		  seg->reffak = 1;
		  seg->copyfrom = -1;
		  self.AppendSegment(seg);
                  }), py::arg("point_indices"), py::arg("leftdomain") = 1, py::arg("rightdomain") = py::int_(0))

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
	  .def("GenerateMesh", [](shared_ptr<SplineGeometry2d> self, MeshingParameters & mparam)
		{
		  shared_ptr<Mesh> mesh = make_shared<Mesh> ();
                  mesh->SetGeometry(self);
                  SetGlobalMesh (mesh);
                  ng_geometry = self;
		  self->GenerateMesh(mesh, mparam);
		  return mesh;
                },py::call_guard<py::gil_scoped_release>())
	  
	  ;
  
}

PYBIND11_MODULE(libgeom2d, m) {
  ExportGeom2d(m);
}

#endif


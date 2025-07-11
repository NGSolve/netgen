#ifdef NG_PYTHON

#include <regex>

#include <../general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include "python_mesh.hpp"

#include <mystdlib.h>
#include "meshing.hpp"
#include "boundarylayer.hpp"
// #include <csg.hpp>
// #include <geometry2d.hpp>
#include <../interface/rw_medit.hpp>
#include <../interface/writeuser.hpp>
#include <../include/nginterface.h>
#include <../general/gzstream.h>


class ClearSolutionClass
{
public:
  ClearSolutionClass() { } 
  ~ClearSolutionClass() { Ng_ClearSolutionData(); }
};



using namespace netgen;

extern const char *ngscript[];

namespace netgen
{
  extern bool netgen_executable_started;
  extern shared_ptr<NetgenGeometry> ng_geometry;
  extern void Optimize2d (Mesh & mesh, MeshingParameters & mp, int faceindex=0);
#ifdef NG_CGNS
  extern tuple<shared_ptr<Mesh>, vector<string>, vector<Array<double>>, vector<int>> ReadCGNSFile(const filesystem::path & filename, int base);
  extern void WriteCGNSFile(shared_ptr<Mesh> mesh, const filesystem::path & filename, vector<string> fields, vector<Array<double>> values, vector<int> locations);
#endif // NG_CGNS
}


void TranslateException (const NgException & ex)
{
  string err = string("Netgen exception: ")+ex.What();
  PyErr_SetString(PyExc_RuntimeError, err.c_str());
}

static Transformation<3> global_trafo(Vec<3> (0,0,0));





DLL_HEADER void ExportNetgenMeshing(py::module &m) 
{
  py::register_exception<NgException>(m, "NgException");
  m.attr("_netgen_executable_started") = py::cast(netgen::netgen_executable_started);
  string script;
  const char ** hcp = ngscript;
  while (*hcp)
      script += *hcp++;

  m.attr("_ngscript") = py::cast(script);

  m.def("_GetStatus", []()
        {
          MyStr s; double percent;
          GetStatus(s, percent);
          return py::make_tuple(s.c_str(), percent);
        });
  m.def("_PushStatus", [](string s) { PushStatus(MyStr(s)); });
  m.def("_SetThreadPercentage", [](double percent) { SetThreadPercent(percent); });

  py::enum_<Identifications::ID_TYPE>(m,"IdentificationType")
    .value("UNDEFINED", Identifications::UNDEFINED)
    .value("PERIODIC", Identifications::PERIODIC)
    .value("CLOSESURFACES", Identifications::CLOSESURFACES)
    .value("CLOSEEDGES", Identifications::CLOSEEDGES)
    ;

  py::implicitly_convertible<int, Identifications::ID_TYPE>();

  
  py::class_<NGDummyArgument>(m, "NGDummyArgument")
    .def("__bool__", []( NGDummyArgument &self ) { return false; } )
    ;

  py::class_<LocalH, shared_ptr<LocalH>>(m, "LocalH");
  
  py::class_<Point<2>> (m, "Point2d")
    .def(py::init<double,double>())
    .def(py::init( [] (std::pair<double,double> xy)
            {
                return Point<2>{xy.first, xy.second};
            }))
    .def ("__str__", &ToString<Point<2>>)
    .def(py::self-py::self)
    .def(py::self+Vec<2>())
    .def(py::self-Vec<2>())
    .def("__getitem__", [](Point<2>& self, int index) { return self[index]; })
    .def("__len__", [](Point<2>& /*unused*/) { return 2; })
    ;

  py::implicitly_convertible<py::tuple, Point<2>>();

  py::class_<Point<3>> (m, "Point3d")
    .def(py::init<double,double,double>())
    .def(py::init([](py::tuple p)
    {
      return Point<3> { p[0].cast<double>(), p[1].cast<double>(),
        p[2].cast<double>() };
    }))
    .def ("__str__", &ToString<Point<3>>)
    .def(py::self-py::self)
    .def(py::self+Vec<3>())
    .def(py::self-Vec<3>())
    .def("__getitem__", [](Point<3>& self, int index) { return self[index]; })
    .def("__len__", [](Point<3>& /*unused*/) { return 3; })
    ;

  py::implicitly_convertible<py::tuple, Point<3>>();

  m.def("Pnt", [](double x, double y, double z)
               { return global_trafo(Point<3>(x,y,z)); });
  m.def("Pnt", [](double x, double y) { return Point<2>(x,y); });
  m.def("Pnt", [](py::array_t<double> np_array)
               {
                 int dim = np_array.size();
                 if(!(dim == 2 || dim == 3))
                   throw Exception("Invalid dimension of input array!");
                 if(dim == 2)
                   return py::cast(Point<2>(np_array.at(0),
                                            np_array.at(1)));
                 return py::cast(global_trafo(Point<3>(np_array.at(0),
                                                       np_array.at(1),
                                                       np_array.at(2))));
               });

  py::class_<Vec<2>> (m, "Vec2d")
    .def(py::init<double,double>())
    .def(py::init( [] (std::pair<double,double> xy)
            {
                return Vec<2>{xy.first, xy.second};
            }))
    .def ("__str__", &ToString<Vec<3>>)
    .def(py::self==py::self)
    .def(py::self+py::self)
    .def(py::self-py::self)
    .def(-py::self)
    .def(double()*py::self)
    .def(py::self*double())
    .def("Norm", &Vec<2>::Length)
    .def("__getitem__", [](Vec<2>& vec, int index) { return vec[index]; })
    .def("__len__", [](Vec<2>& /*unused*/) { return 2; })
    ;

  py::implicitly_convertible<py::tuple, Vec<2>>();

  py::class_<Vec<3>> (m, "Vec3d")
    .def(py::init<double,double,double>())
    .def(py::init([](py::tuple v)
    {
      return Vec<3> { v[0].cast<double>(), v[1].cast<double>(),
        v[2].cast<double>() };
    }))
    .def ("__str__", &ToString<Vec<3>>)
    .def(py::self==py::self)
    .def(py::self+py::self)
    .def(py::self-py::self)
    .def(-py::self)
    .def(double()*py::self)
    .def(py::self*double())
    .def("Norm", &Vec<3>::Length)
    .def("__getitem__", [](Vec<3>& vec, int index) { return vec[index]; })
    .def("__len__", [](Vec<3>& /*unused*/) { return 3; })
    ;

  py::implicitly_convertible<py::tuple, Vec<3>>();

  py::class_<Mat<3,3>>(m, "Mat33")
    .def(py::init([](py::tuple m)
    {
      if(m.size() != 9)
        throw std::length_error("Invalid dimension of input array!");
      Mat<3,3> mat;
      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
          mat(i,j) = m[i*3+j].cast<double>();
      return mat;
    }))
    .def("__getitem__", [](Mat<3,3>& mat, py::tuple index)
    {
      if(index.size() != 2)
        throw std::length_error("Invalid dimension of input array!");
      return mat(index[0].cast<int>(), index[1].cast<int>());
    })
    .def("__setitem__", [](Mat<3,3>& mat, py::tuple index, double val)
    {
      if(index.size() != 2)
        throw std::length_error("Invalid dimension of input array!");
      mat(index[0].cast<int>(), index[1].cast<int>()) = val;
    })
    .def("__str__", &ToString<Mat<3,3>>)
    ;

  py::implicitly_convertible<py::tuple, Mat<3,3>>();

  m.def ("Vec", FunctionPointer
           ([] (double x, double y, double z) { return global_trafo(Vec<3>(x,y,z)); }));
  m.def("Vec", [](py::array_t<double> np_array)
               {
                 int dim = np_array.size();
                 if(!(dim == 2 || dim == 3))
                   throw Exception("Invalid dimension of input array!");
                 if(dim == 2)
                   return py::cast(Vec<2>(np_array.at(0),
                                          np_array.at(1)));
                 return py::cast(global_trafo(Vec<3>(np_array.at(0),
                                                     np_array.at(1),
                                                     np_array.at(2))));
               });

  m.def ("Vec", FunctionPointer
           ([] (double x, double y) { return Vec<2>(x,y); }));

  py::class_<Transformation<3>> (m, "Trafo")
    .def(py::init<Vec<3>>(), "a translation")
    .def(py::init<Point<3>,Vec<3>,double>(), "a rotation given by point on axes, direction of axes, angle")
    .def("__mul__", [](Transformation<3> a, Transformation<3> b)->Transformation<3>
         { Transformation<3> res; res.Combine(a,b); return res; })
    .def("__call__", [] (Transformation<3> trafo, Point<3> p) { return trafo(p); })
    .def_property("mat", &Transformation<3>::GetMatrix,
                  [](Transformation<3>& self, const Mat<3,3>& mat)
                  {
                    self.GetMatrix() = mat;
                  })
    ;

  m.def ("GetTransformation", [] () { return global_trafo; });
  m.def ("SetTransformation", [] (Transformation<3> trafo) { global_trafo = trafo; });
  m.def ("SetTransformation", 
         [](int dir, double angle)
         {
           if (dir > 0)
             global_trafo.SetAxisRotation (dir, angle*M_PI/180);
           else
             global_trafo = Transformation<3> (Vec<3>(0,0,0));
         },
         py::arg("dir")=int(0), py::arg("angle")=int(0));
  m.def ("SetTransformation", 
         [](Point<3> p0, Vec<3> ex, Vec<3> ey, Vec<3> ez)
            {
              Point<3> pnts[4];
              pnts[0] = p0;
              pnts[1] = p0 + ex;
              pnts[2] = p0 + ey;
              pnts[3] = p0 + ez;
              global_trafo = Transformation<3> (pnts);
            },
         py::arg("p0"), py::arg("ex"), py::arg("ey"), py::arg("ez"));


  
  py::class_<PointIndex>(m, "PointId")
    .def(py::init<int>())
    .def("__repr__", &ToString<PointIndex>)
    .def("__str__", &ToString<PointIndex>)
    .def_property_readonly("nr", &PointIndex::operator int)
    .def("__eq__" , FunctionPointer( [](PointIndex &self, PointIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](PointIndex &self ) { return static_cast<int>(self); }) )
    ;

  py::class_<ElementIndex>(m, "ElementId3D")
    .def(py::init<int>())
    .def("__repr__", &ToString<ElementIndex>)
    .def("__str__", &ToString<ElementIndex>)
    .def_property_readonly("nr", &ElementIndex::operator int)
    .def("__eq__" , FunctionPointer( [](ElementIndex &self, ElementIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](ElementIndex &self ) { return static_cast<int>(self); }) )
    ;


  py::class_<SurfaceElementIndex>(m, "ElementId2D")
    .def(py::init<int>())
    .def("__repr__", &ToString<SurfaceElementIndex>)
    .def("__str__", &ToString<SurfaceElementIndex>)
    .def_property_readonly("nr", &SurfaceElementIndex::operator int)
    .def("__eq__" , FunctionPointer( [](SurfaceElementIndex &self, SurfaceElementIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](SurfaceElementIndex &self ) { return static_cast<int>(self); }) )
    ;

  py::class_<SegmentIndex>(m, "ElementId1D")
    .def(py::init<int>())
    .def("__repr__", &ToString<SegmentIndex>)
    .def("__str__", &ToString<SegmentIndex>)
    .def_property_readonly("nr", &SegmentIndex::operator int)
    .def("__eq__" , FunctionPointer( [](SegmentIndex &self, SegmentIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](SegmentIndex &self ) { return static_cast<int>(self); }) )
    ;



  /*  
  py::class_<Point<3>> ("Point")
    .def(py::init<double,double,double>())
    ;
  */

  py::class_<MeshPoint /* ,py::bases<Point<3>> */ >(m, "MeshPoint")
    .def(py::init<Point<3>>())
    .def("__str__", &ToString<MeshPoint>)
    .def("__repr__", &ToString<MeshPoint>)
    .def_property_readonly("p", [](const MeshPoint & self)
                           {
                             py::list l;
                             l.append ( py::cast(self[0]) );
                             l.append ( py::cast(self[1]) );
                             l.append ( py::cast(self[2]) );
                             return py::tuple(l);
                           })
    .def("__getitem__", [](const MeshPoint & self, int index) {
	  if(index<0 || index>2)
              throw py::index_error();
	  return self[index];
	})
    .def("__setitem__", [](MeshPoint & self, int index, double val) {
	  if(index<0 || index>2)
              throw py::index_error();
	  self(index) = val;
    })
    .def_property("singular",
                  [](const MeshPoint & pnt) { return pnt.Singularity(); },
                  [](MeshPoint & pnt, double sing) { pnt.Singularity(sing); })
    ;

  py::class_<Element>(m, "Element3D")
    .def(py::init([](int index, std::vector<PointIndex> vertices)
                  {
                    int np = vertices.size();
                    ELEMENT_TYPE et;
                    switch (np)
                      {
                      case 4: et = TET; break;
                      case 5: et = PYRAMID; break;
                      case 6: et = PRISM; break;
                      case 8: et = HEX; break;
                      case 10: et = TET10; break;
                      case 13: et = PYRAMID13; break;
                      case 15: et = PRISM15; break;
                      case 20: et = HEX20; break;
                      default:
                        throw Exception ("no Element3D with " + ToString(np) +
                                         " points");
                      }

                    auto newel = new Element(et);
                    for(int i=0; i<np; i++)
                      (*newel)[i] = vertices[i];
                    newel->SetIndex(index);
                    return newel;
                  }),
          py::arg("index")=1,py::arg("vertices"),
         "create volume element"
         )
    .def("__repr__", &ToString<Element>)
    .def_property("index", &Element::GetIndex, &Element::SetIndex)
    .def_property("curved", &Element::IsCurved, &Element::SetCurved)
    .def_property("refine", &Element::TestRefinementFlag, &Element::SetRefinementFlag)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Element & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNV(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    .def_property_readonly("points", 
                  FunctionPointer ([](const Element & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNP(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    ;

  if(ngcore_have_numpy)
  {
    auto data_layout = Element::GetDataLayout();

    py::detail::npy_format_descriptor<Element>::register_dtype({
        py::detail::field_descriptor {
          "nodes", data_layout["pnum"],
          ELEMENT_MAXPOINTS * sizeof(PointIndex),
          py::format_descriptor<int[ELEMENT_MAXPOINTS]>::format(),
          py::detail::npy_format_descriptor<int[ELEMENT_MAXPOINTS]>::dtype() },
        py::detail::field_descriptor {
          "index", data_layout["index"], sizeof(int),
          py::format_descriptor<int>::format(),
          py::detail::npy_format_descriptor<int>::dtype() },
        py::detail::field_descriptor {
          "np", data_layout["np"], sizeof(int8_t),
          py::format_descriptor<signed char>::format(),
            pybind11::dtype("int8") },
        py::detail::field_descriptor {
          "refine", data_layout["refine"], sizeof(bool),
          py::format_descriptor<bool>::format(),
          py::detail::npy_format_descriptor<bool>::dtype() },            
        py::detail::field_descriptor {
          "curved", data_layout["curved"], sizeof(bool),
          py::format_descriptor<bool>::format(),
          py::detail::npy_format_descriptor<bool>::dtype()}            
      });
  }

  py::class_<Element2d>(m, "Element2D")
    .def(py::init ([](int index, std::vector<PointIndex> vertices)
                   {
                     Element2d * newel = nullptr;
                     if (vertices.size() == 3)
                       {
                         newel = new Element2d(TRIG);
                         for (int i = 0; i < 3; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else if (vertices.size() == 4)
                       {
                         newel = new Element2d(QUAD);
                         for (int i = 0; i < 4; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else if (vertices.size() == 6)
                       {
                         newel = new Element2d(TRIG6);
                         for(int i = 0; i<6; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else if (vertices.size() == 8)
                       {
                         newel = new Element2d(QUAD8);
                         for(int i = 0; i<8; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else 
                       throw NgException("Inconsistent number of vertices in Element2D");
                     return newel;
                   }),
                   py::arg("index")=1,py::arg("vertices"),
         "create surface element"
         )
    .def_property("index", &Element2d::GetIndex, &Element2d::SetIndex)
    .def_property("curved", &Element2d::IsCurved, &Element2d::SetCurved)
    .def_property("refine", &Element2d::TestRefinementFlag, &Element2d::SetRefinementFlag)
    .def_property_readonly("geominfo", [](const Element2d& self) -> py::list
    {
      py::list li;
      for (const auto &pgi : self.GeomInfo())
        li.append(py::make_tuple(pgi.trignum, pgi.u, pgi.v));
      return li;
    })
    .def_property_readonly("vertices",
                  FunctionPointer([](const Element2d & self) -> py::list
                                  {
                                    py::list li;
                                    for (int i = 0; i < self.GetNV(); i++)
                                      li.append(py::cast(self[i]));
                                    return li;
                                  }))
    .def_property_readonly("points", 
                  FunctionPointer ([](const Element2d & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNP(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))    
    ;

  if(ngcore_have_numpy)
  {
    auto data_layout = Element2d::GetDataLayout();
    py::detail::npy_format_descriptor<Element2d>::register_dtype({
        py::detail::field_descriptor {
          "nodes", data_layout["pnum"],
          ELEMENT2D_MAXPOINTS * sizeof(PointIndex),
          py::format_descriptor<int[ELEMENT2D_MAXPOINTS]>::format(),
          py::detail::npy_format_descriptor<int[ELEMENT2D_MAXPOINTS]>::dtype() },
        py::detail::field_descriptor {
          "index", data_layout["index"], sizeof(int),
          py::format_descriptor<int>::format(),
          py::detail::npy_format_descriptor<int>::dtype() },
        py::detail::field_descriptor {
          "np", data_layout["np"], sizeof(int8_t),
          py::format_descriptor<signed char>::format(),
        pybind11::dtype("int8") },
        py::detail::field_descriptor {
          "refine", data_layout["refine"], sizeof(bool),
          py::format_descriptor<bool>::format(),
          py::detail::npy_format_descriptor<bool>::dtype() },
        py::detail::field_descriptor {
          "curved", data_layout["curved"], sizeof(bool),
          py::format_descriptor<bool>::format(),
          py::detail::npy_format_descriptor<bool>::dtype() }
      });
  }

  py::class_<Segment>(m, "Element1D")
    .def(py::init([](py::list vertices, py::list surfaces, int index, int edgenr,
                     py::list trignums)
                  {
                    Segment * newel = new Segment();
                    for (int i = 0; i < 2; i++)
                      (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                    newel -> si = index;
                    newel -> epgeominfo[0].edgenr = edgenr;
                    newel -> epgeominfo[1].edgenr = edgenr;
                    newel -> edgenr = index;
                    for(auto i : Range(len(trignums)))
                      newel->geominfo[i].trignum = py::cast<int>(trignums[i]);
                    if (len(surfaces))
                      {
                        newel->surfnr1 = py::extract<int>(surfaces[0])();
                        newel->surfnr2 = py::extract<int>(surfaces[1])();
                      }
                    return newel;
                  }),
          py::arg("vertices"),
           py::arg("surfaces")=py::list(),
           py::arg("index")=1,
           py::arg("edgenr")=1,
           py::arg("trignums")=py::list(), // for stl
         "create segment element"
         )
    .def("__repr__", &ToString<Segment>)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < 2; i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    .def_property_readonly("points", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNP(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    .def_property_readonly("surfaces", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     li.append (py::cast(self.surfnr1));
                                     li.append (py::cast(self.surfnr2));
                                     return li;
                                   }))
    .def_property("index",
                  [](const Segment &self) -> size_t
                  {
                    return self.si;
                  },
                  [](Segment& self, int index)
                  {
                    self.si = index;
                  })
    .def_property("edgenr",
                  [](const Segment & self) -> size_t
                  {
                    return self.edgenr;
                  },
                  [](Segment& self, int edgenr)
                  {
                    self.edgenr = edgenr;
                  })
    .def_property("singular",
                  [](const Segment & seg) { return seg.singedge_left; },
                  [](Segment & seg, double sing) { seg.singedge_left = sing; seg.singedge_right=sing; })
    ;

  if(ngcore_have_numpy)
  {
    py::detail::npy_format_descriptor<Segment>::register_dtype({
        py::detail::field_descriptor {
          "nodes", offsetof(Segment, pnums),
          3 * sizeof(PointIndex),
          py::format_descriptor<int[3]>::format(),
          py::detail::npy_format_descriptor<int[3]>::dtype() },
        py::detail::field_descriptor {
          "index", offsetof(Segment, si), sizeof(int),
          py::format_descriptor<int>::format(),
          py::detail::npy_format_descriptor<int>::dtype() },
        py::detail::field_descriptor {
          "edgenr", offsetof(Segment, edgenr), sizeof(int),
          py::format_descriptor<int>::format(),
          py::detail::npy_format_descriptor<int>::dtype() },
      });
  }

  py::class_<Element0d>(m, "Element0D")
    .def(py::init([](PointIndex vertex, int index)
                  {
                    Element0d * instance = new Element0d;
                    instance->pnum = vertex;
                    instance->index = index;
                    return instance;
                  }),
         py::arg("vertex"),
         py::arg("index")=1,
         "create point element"
         )
    .def("__repr__", &ToString<Element0d>)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Element0d & self) -> py::list
                                   {
                                     py::list li;
                                     li.append (py::cast(self.pnum));
                                     return li;
                                   }))
    ;
  
  
  


  py::class_<FaceDescriptor>(m, "FaceDescriptor")
    .def(py::init<const FaceDescriptor&>())
    .def(py::init([](int surfnr, int domin, int domout, int bc)
                  {
                    FaceDescriptor * instance = new FaceDescriptor();
                    instance->SetSurfNr(surfnr);
                    instance->SetDomainIn(domin);
                    instance->SetDomainOut(domout);
                    instance->SetBCProperty(bc);
                             return instance;
                  }),
         py::arg("surfnr")=1, 
         py::arg("domin")=1,
         py::arg("domout")=py::int_(0),
         py::arg("bc")=py::int_(0),
         "create facedescriptor")
    .def("__str__", &ToString<FaceDescriptor>)
    .def("__repr__", &ToString<FaceDescriptor>)
    .def_property("surfnr", &FaceDescriptor::SurfNr, &FaceDescriptor::SetSurfNr)
    .def_property("domin", &FaceDescriptor::DomainIn, &FaceDescriptor::SetDomainIn)
    .def_property("domout", &FaceDescriptor::DomainOut, &FaceDescriptor::SetDomainOut)
    .def_property("domin_singular", &FaceDescriptor::DomainInSingular, &FaceDescriptor::SetDomainInSingular)
    .def_property("domout_singular", &FaceDescriptor::DomainOutSingular, &FaceDescriptor::SetDomainOutSingular)

    .def_property("bc", &FaceDescriptor::BCProperty, &FaceDescriptor::SetBCProperty)
    .def_property("bcname",
                  [](FaceDescriptor & self) -> string { return self.GetBCName(); },
                  [](FaceDescriptor & self, string name) { self.SetBCName(new string(name)); } // memleak
                  )
    .def_property("color",
                  [](const FaceDescriptor& self)
                  {
                    auto sc = self.SurfColour();
                    return py::make_tuple(sc[0], sc[1], sc[2], sc[3]);
                  },
                  [](FaceDescriptor& self, py::tuple col)
                  {
                    Vec<4> sc = 1;
                    sc[0] = py::cast<double>(col[0]);
                    sc[1] = py::cast<double>(col[1]);
                    sc[2] = py::cast<double>(col[2]);
                    if(py::len(col) > 3)
                      sc[3] = py::cast<double>(col[3]);
                    self.SetSurfColour(sc);
                  }
                  )
    .def_property("transparency",
                  [](const FaceDescriptor& self)
                  {
                    return self.SurfColour()[3];
                  },
                  [](FaceDescriptor& self, double val)
                  {
                    auto sc = self.SurfColour();
                    sc[3] = val;
                    self.SetSurfColour(sc);
                  })
    ;

  py::implicitly_convertible< int, SurfaceElementIndex>();
  PYBIND11_NUMPY_DTYPE(SurfaceElementIndex, i);
  ExportArray<SurfaceElementIndex, SurfaceElementIndex>(m);
  
  py::implicitly_convertible< int, ElementIndex>();
  PYBIND11_NUMPY_DTYPE(ElementIndex, i);
  ExportArray<ElementIndex, ElementIndex>(m);
  
  ExportArray<Element,ElementIndex>(m);
  ExportArray<Element2d,SurfaceElementIndex>(m);
  ExportArray<Segment,SegmentIndex>(m);
  ExportArray<Element0d>(m);
  ExportArray<MeshPoint,PointIndex>(m);
  ExportArray<FaceDescriptor>(m);

  string export_docu = "Export mesh to other file format. Supported formats are:\n";
  Array<string> export_formats;
  for(auto & e : UserFormatRegister::entries)
    if(e.write) {
      string s = '\t'+e.format+"\t("+e.extensions[0];
      for(auto & ext : e.extensions.Range(1, e.extensions.Size()))
        s += ", "+ext;
      s += ")\n";
      export_formats.Append(s);
    }
  QuickSort(export_formats);
  for(const auto & s : export_formats)
    export_docu += s;

  py::implicitly_convertible< int, PointIndex>();

  py::class_<NetgenGeometry, shared_ptr<NetgenGeometry>> (m, "NetgenGeometry", py::dynamic_attr())
    .def("RestrictH", &NetgenGeometry::RestrictH)
             ;
  
  py::class_<Mesh,shared_ptr<Mesh>>(m, "Mesh")
    // .def(py::init<>("create empty mesh"))

    .def(py::init( [] (int dim, NgMPI_Comm comm)
                   {
                     auto mesh = make_shared<Mesh>();
		     mesh->SetCommunicator(comm);
                     mesh -> SetDimension(dim);
                     SetGlobalMesh(mesh);  // for visualization
                     mesh -> SetGeometry (nullptr);
                     return mesh;
                   } ),
         py::arg("dim")=3, py::arg("comm")=NgMPI_Comm{}
         )
    .def(NGSPickle<Mesh>())
    .def_property_readonly("comm", [](const Mesh & amesh) -> NgMPI_Comm
			   { return amesh.GetCommunicator(); },
                           "MPI-communicator the Mesh lives in")
    /*
    .def("__init__",
         [](Mesh *instance, int dim)
                           {
                             new (instance) Mesh();
                             instance->SetDimension(dim);
                           },
           py::arg("dim")=3
          )
    */
    
    .def_property_readonly("_timestamp", &Mesh::GetTimeStamp)
    .def_property_readonly("ne", [](Mesh& m) { return m.GetNE(); })
    .def_property_readonly("bounding_box", [](Mesh& m) {
          Point3d pmin, pmax;
          m.GetBox(pmin, pmax);
          return py::make_tuple( Point<3>(pmin),Point<3>(pmax));
    })
    .def("Partition", [](shared_ptr<Mesh> self, int numproc) {
        self->ParallelMetis(numproc);
      }, py::arg("numproc"))
    
    .def("Distribute", [](shared_ptr<Mesh> self, NgMPI_Comm comm) {
	self->SetCommunicator(comm);
	if(comm.Size()==1) return self;
	// if(MyMPI_GetNTasks(comm)==2) throw NgException("Sorry, cannot handle communicators with NP=2!");
	// cout << " rank " << MyMPI_GetId(comm) << " of " << MyMPI_GetNTasks(comm) << " called Distribute " << endl;
	if(comm.Rank()==0) self->Distribute();
	else self->SendRecvMesh();
	return self;
      }, py::arg("comm"))
    .def_static("Receive", [](NgMPI_Comm comm) -> shared_ptr<Mesh> {
        auto mesh = make_shared<Mesh>();
        mesh->SetCommunicator(comm);
        mesh->SendRecvMesh();
        return mesh;
      }, py::arg("comm"))
    .def("Load",  FunctionPointer 
	 ([](shared_ptr<Mesh> self, const string & filename)
	  {

	    auto comm = self->GetCommunicator();
	    int id = comm.Rank();
	    int ntasks = comm.Size();
	    auto & mesh = self;

	    {
	      ifstream infile(filename.c_str());
	      if(!infile.good())
		throw NgException(string("Error opening file ") + filename);
	    }

	    if ( filename.find(".vol") == string::npos )
	      {
		if(ntasks>1)
		  throw NgException("Not sure what to do with this?? Does this work with MPI??");
		mesh->SetCommunicator(comm);
		ReadFile(*mesh,filename.c_str());
		//mesh->SetGlobalH (mparam.maxh);
		//mesh->CalcLocalH();
		return;
	      }

	    istream * infile = nullptr;
	    Array<char> buf; // for distributing geometry!
	    int strs;

	    if( id == 0) {
	      if (filename.length() > 8 && filename.substr (filename.length()-8, 8) == ".vol.bin")
                mesh -> Load(filename);
              else if (filename.substr (filename.length()-3, 3) == ".gz")
		infile = new igzstream (filename.c_str());
	      else
		infile = new ifstream (filename.c_str());

              if(infile)
                {
                  mesh -> Load(*infile);
                  // make string from rest of file (for geometry info!)
                  // (this might be empty, in which case we take the global ng_geometry)
                  stringstream geom_part;
                  geom_part << infile->rdbuf();
                  string geom_part_string = geom_part.str();
                  strs = geom_part_string.size();
                  // buf = new char[strs];
                  buf.SetSize(strs);
                  memcpy(buf.Data(), geom_part_string.c_str(), strs*sizeof(char));
                  delete infile;
                }


	      if (ntasks > 1)
		{

		  char * weightsfilename = new char [filename.size()+1];
		  strcpy (weightsfilename, filename.c_str());            
		  weightsfilename[strlen (weightsfilename)-3] = 'w';
		  weightsfilename[strlen (weightsfilename)-2] = 'e';
		  weightsfilename[strlen (weightsfilename)-1] = 'i';

		  ifstream weightsfile(weightsfilename);      
		  delete [] weightsfilename;  
	  
		  if (!(weightsfile.good()))
		    {
		      // cout << "regular distribute" << endl;
		      mesh -> Distribute();
		    }
		  else
		    {
		      char str[20];   
		      bool endfile = false;
		      int n, dummy;
	      
		      NgArray<int> segment_weights;
		      NgArray<int> surface_weights;
		      NgArray<int> volume_weights;
	      
		      while (weightsfile.good() && !endfile)
			{
			  weightsfile >> str;
		  
			  if (strcmp (str, "edgeweights") == 0)
			    {
			      weightsfile >> n;
			      segment_weights.SetSize(n);
			      for (int i = 0; i < n; i++)
				weightsfile >> dummy >> segment_weights[i];
			    }
		  
			  if (strcmp (str, "surfaceweights") == 0)
			    {
			      weightsfile >> n;
			      surface_weights.SetSize(n);
			      for (int i=0; i<n; i++)
				weightsfile >> dummy >> surface_weights[i];
			    }
		  
			  if (strcmp (str, "volumeweights") == 0)
			    {
			      weightsfile >> n;
			      volume_weights.SetSize(n);
			      for (int i=0; i<n; i++)
				weightsfile >> dummy >> volume_weights[i];
			    }
		  
			  if (strcmp (str, "endfile") == 0)
			    endfile = true;  
			}     
	      
		      mesh -> Distribute(volume_weights, surface_weights, segment_weights);
		    }
		} // ntasks>1 end
	    } // id==0 end
	    else {
	      mesh->SendRecvMesh();
	    }

	    if(ntasks>1) {
              // #ifdef PARALLEL
	      /** Scatter the geometry-string (no dummy-implementation in mpi_interface) **/
              /*
	      int strs = buf.Size();
	      MyMPI_Bcast(strs, comm);
	      if(strs>0)
		MyMPI_Bcast(buf, comm);
              */
              comm.Bcast(buf);
              // #endif
	    }

	    shared_ptr<NetgenGeometry> geo;
	    if(buf.Size()) { // if we had geom-info in the file, take it
	      istringstream geom_infile(string((const char*)buf.Data(), buf.Size()));
	      geo = GeometryRegister().LoadFromMeshFile(geom_infile);
	    }
	    if(geo!=nullptr) mesh->SetGeometry(geo);
	    else if(ng_geometry!=nullptr) mesh->SetGeometry(ng_geometry);
	  }),py::call_guard<py::gil_scoped_release>())
    .def("Save", static_cast<void(Mesh::*)(const filesystem::path & name)const>(&Mesh::Save),py::call_guard<py::gil_scoped_release>())
    .def("Export",
         [] (Mesh & self, string filename, string format)
          {
            if (WriteUserFormat (format, self, filename))
                throw Exception ("Nothing known about format"+format);
          },
         py::arg("filename"), py::arg("format"), export_docu.c_str(),
         py::call_guard<py::gil_scoped_release>())
    
    .def_property("dim", &Mesh::GetDimension, &Mesh::SetDimension)

    .def("Elements3D", 
         static_cast<Array<Element,ElementIndex>&(Mesh::*)()> (&Mesh::VolumeElements),
         py::return_value_policy::reference)

    .def("Elements2D", 
         static_cast<Array<Element2d,SurfaceElementIndex>&(Mesh::*)()> (&Mesh::SurfaceElements),
         py::return_value_policy::reference)

    .def("Elements1D", 
         static_cast<Array<Segment, SegmentIndex>&(Mesh::*)()> (&Mesh::LineSegments),
         py::return_value_policy::reference)

    .def("Elements0D", FunctionPointer([] (Mesh & self) -> Array<Element0d>&
                                       {
                                         return self.pointelements;
                                       } ),
         py::return_value_policy::reference)

    .def("Points", 
         static_cast<Mesh::T_POINTS&(Mesh::*)()> (&Mesh::Points),
         py::return_value_policy::reference)

    .def("Coordinates", [](Mesh & self) {
        return py::array
          (
           py::memoryview::from_buffer
           (&self.Points()[PointIndex::BASE](0), sizeof(double),
            py::format_descriptor<double>::value,
            { self.Points().Size(), size_t(self.GetDimension())  }, 
            { sizeof(self.Points()[PointIndex::BASE]), sizeof(double) } )
           );
      })
    .def_property_readonly("parentelements", [](Mesh & self) {
      // return FlatArray<int>(self.mlparentelement.Size(), &self.mlparentelement[0]);
      return FlatArray(self.mlparentelement);
    }, py::keep_alive<0,1>())
    .def_property_readonly("parentsurfaceelements", [](Mesh & self) {
      // return FlatArray<int>(self.mlparentsurfaceelement.Size(),
      // &self.mlparentsurfaceelement[0]);
      return FlatArray(self.mlparentsurfaceelement);      
    }, py::keep_alive<0,1>())
    .def_property_readonly("macromesh", [](Mesh & self) {
      auto coarsemesh = make_shared<Mesh>();
      *coarsemesh = *self.coarsemesh;
      return coarsemesh;
    }, "mesh before hp-refinement")
    .def("MacroElementNr", [](Mesh & self, int elnr, optional<int> dim) {
      // cout << "hpels = " << self.hpelements->Size() << endl;
      // return self[ElementIndex(elnr)].GetHpElnr();
      if (!dim) dim = self.GetDimension();
      switch (*dim)
        {
        case 2:
          return (*self.hpelements)[self[SurfaceElementIndex(elnr)].GetHpElnr()].coarse_elnr;
        case 3:
          return (*self.hpelements)[self[ElementIndex(elnr)].GetHpElnr()].coarse_elnr;
        }
      throw Exception ("MacroElementNr not implemented for dim");
    }, py::arg("elnr"), py::arg("dim")=nullopt, "number of macro element of element number elnr")
    .def("FaceDescriptor", static_cast<FaceDescriptor&(Mesh::*)(int)> (&Mesh::GetFaceDescriptor),
         py::return_value_policy::reference)
    .def("GetNFaceDescriptors", &Mesh::GetNFD)
    .def("RestrictLocalH", [](Mesh& self, const Point<3>& pnt, double maxh,
                              int layer)
    {
      self.RestrictLocalH(pnt, maxh, layer);
    }, py::arg("p"), py::arg("h"), py::arg("layer")=1)
    .def("FaceDescriptors", 
         // static_cast<Array<Element>&(Mesh::*)()> (&Mesh::FaceDescriptors),
         &Mesh::FaceDescriptors,         
         py::return_value_policy::reference)
    
    
    .def("GetNDomains", &Mesh::GetNDomains)

    .def("GetVolumeNeighboursOfSurfaceElement", [](Mesh & self, size_t sel)
                                                {
                                                  int elnr1, elnr2;
                                                  self.GetTopology().GetSurface2VolumeElement(sel+1, elnr1, elnr2);
                                                  return py::make_tuple(elnr1, elnr2);
                                                }, "Returns element nrs of volume element connected to surface element, -1 if no volume element")

    .def("GetNCD2Names", &Mesh::GetNCD2Names)
    

    .def("__getitem__", [](const Mesh & self, PointIndex id) { return self[id]; })
    .def("__getitem__", [](const Mesh & self, ElementIndex id) { return self[id]; })
    .def("__getitem__", [](const Mesh & self, SurfaceElementIndex id) { return self[id]; })
    .def("__getitem__", [](const Mesh & self, SegmentIndex id) { return self[id]; })

    .def("__setitem__", [](Mesh & self, PointIndex id, const MeshPoint & mp) { return self[id] = mp; })
    
    .def ("Add", [](Mesh & self, MeshPoint p)
          {
            return self.AddPoint (Point3d(p));
          })
          
    .def ("Add", [](Mesh & self, const Element & el)
          {
            return self.AddVolumeElement (el);
          })
          
    .def ("Add", [](Mesh & self, const Element2d & el)
          {
            return self.AddSurfaceElement (el);
          })

    .def ("Add", [](Mesh & self, const Segment & el, bool project_geominfo)
          {
            if (project_geominfo)
              {
                auto &p1 = self[el[0]];
                auto &p2 = self[el[1]];
                auto geo = self.GetGeometry();
                geo->ProjectPointEdge
                  (0,0,p1,
                   const_cast<EdgePointGeomInfo*>(&el.epgeominfo[0]));
                geo->ProjectPointEdge
                  (0,0,p2,
                   const_cast<EdgePointGeomInfo*>(&el.epgeominfo[1]));
              }
            return self.AddSegment (el);
          }, py::arg("el"), py::arg("project_geominfo")=false)
          
    .def ("Add", [](Mesh & self, const Element0d & el)
          {
            return self.pointelements.Append (el);
          })

    .def ("Add", [](Mesh & self, const FaceDescriptor & fd)
          {
            return self.AddFaceDescriptor (fd);
          })

    .def ("AddSingularity", [](Mesh & self, PointIndex pi, double factor)
         {
	   self[pi].Singularity(factor);
         })

    .def ("AddPoints", [](Mesh & self, py::buffer b1)
          {
            static Timer timer("Mesh::AddPoints");
            static Timer timercast("Mesh::AddPoints - casting");            
            RegionTimer reg(timer);

            timercast.Start();
            // casting from here: https://github.com/pybind/pybind11/issues/1908
            auto b = b1.cast<py::array_t<double_t, py::array::c_style | py::array::forcecast>>();
            timercast.Stop();
            
            py::buffer_info info = b.request();
            // cout << "data format = " << info.format << endl;
            if (info.ndim != 2)
              throw std::runtime_error("AddPoints needs buffer of dimension 2");
            // if (info.format != py::format_descriptor<double>::format())
            // throw std::runtime_error("AddPoints needs buffer of type double");
            if (info.strides[0] != sizeof(double)*info.shape[1])
              throw std::runtime_error("AddPoints needs packed array");              
            double * ptr = static_cast<double*> (info.ptr);
            
            self.Points().SetAllocSize(self.Points().Size()+info.shape[0]);
            if (info.shape[1]==2)
              for ([[maybe_unused]] auto i : Range(info.shape[0]))
                {
                  self.AddPoint (Point<3>(ptr[0], ptr[1], 0));
                  ptr += 2;
                }
            if (info.shape[1]==3)
              for ([[maybe_unused]] auto i : Range(info.shape[0]))
                {
                  self.AddPoint (Point<3>(ptr[0], ptr[1], ptr[2]));
                  ptr += 3;
                }
          })
    .def ("AddElements", [](Mesh & self, int dim, int index, py::buffer b1, int base,
                            bool project_geometry)
          {
            static Timer timer("Mesh::AddElements");
            static Timer timercast("Mesh::AddElements casting");
            RegionTimer reg(timer);

            timercast.Start();
            auto b = b1.cast<py::array_t<int, py::array::c_style | py::array::forcecast>>();
            timercast.Stop();
            
            py::buffer_info info = b.request();
            if (info.ndim != 2)
              throw std::runtime_error("AddElements needs buffer of dimension 2");
            // if (info.format != py::format_descriptor<int>::format())
            // throw std::runtime_error("AddPoints needs buffer of type int");

            int * ptr = static_cast<int*> (info.ptr);
            if (dim == 1)
              {
                // ELEMENT_TYPE type;
                int np = info.shape[1];
                self.LineSegments().SetAllocSize(self.LineSegments().Size()+info.shape[0]);                
                for ([[maybe_unused]] auto i : Range(info.shape[0]))
                  {
                    Segment el;
                    for (int j = 0; j < np; j++)
                      el[j] = ptr[j]+PointIndex::BASE-base;
                    el.si = index;
                    self.AddSegment(el);
                    ptr += info.strides[0]/sizeof(int);
                  }
              }
            if (dim == 2)
              {
                ELEMENT_TYPE type;
                int np = info.shape[1];
                switch (np)
                  {
                  case 3: type = TRIG; break;
                  case 4: type = QUAD; break;
                  case 6: type = TRIG6; break;
                  case 8: type = QUAD8; break;
                  default:
                    throw Exception("unsupported 2D element with "+ToString(np)+" points");
                  }
                self.SurfaceElements().SetAllocSize(self.SurfaceElements().Size()+info.shape[0]);                
                for ([[maybe_unused]] auto i : Range(info.shape[0]))
                  {
                    Element2d el(type);
                    for (int j = 0; j < np; j++)
                      el[j] = ptr[j]+PointIndex::BASE-base;
                    el.SetIndex(index);
                    if(project_geometry)
                      {
                        // find some point in the mid of trig/quad for
                        // quick + stable uv-projection of all points
                        auto startp = Center(self[el[0]], self[el[1]], self[el[2]]);
                        PointGeomInfo gi = self.GetGeometry()->ProjectPoint(index,
                                                                            startp);
                        for(auto i : Range(np))
                          {
                            el.GeomInfo()[i] = gi;
                            self.GetGeometry()->ProjectPointGI(index,
                                                               self[el[i]],
                                                               el.GeomInfo()[i]);
                          }
                      }
                    self.AddSurfaceElement (el);
                    ptr += info.strides[0]/sizeof(int);
                  }
              }
            if (dim == 3)
              {
                ELEMENT_TYPE type;
                int np = info.shape[1];
                switch (np)
                  {
                  case 4: type = TET; break;
                    /* // have to check ordering of points
                       case 10: type = TET10; break;
                       case 8: type = HEX; break;
                       case 6: type = PRISM; break;
                    */
                  default:
                    throw Exception("unsupported 3D element with "+ToString(np)+" points");
                  }
                self.VolumeElements().SetAllocSize(self.VolumeElements().Size()+info.shape[0]);
                for ([[maybe_unused]] auto i : Range(info.shape[0]))
                  {
                    Element el(type);
                    for (int j = 0; j < np;j ++)
                      el[j] = ptr[j]+PointIndex::BASE-base;
                    el.SetIndex(index);
                    self.AddVolumeElement (el);
                    ptr += info.strides[0]/sizeof(int);
                  }
              }
            
          }, py::arg("dim"), py::arg("index"), py::arg("data"), py::arg("base")=0,
          py::arg("project_geometry")=false)
    
    .def ("DeleteSurfaceElement",
          [](Mesh & self, SurfaceElementIndex i)
          {
            return self.Delete(i);
          })
          
    .def ("Compress", [](Mesh & self)
          {
            return self.Compress ();
          } ,py::call_guard<py::gil_scoped_release>())
          
    .def ("AddRegion", [] (Mesh & self, string name, int dim) -> int
         {
           auto & regionnames = self.GetRegionNamesCD(self.GetDimension()-dim);
           regionnames.Append (new string(name));
           int idx = regionnames.Size();
           if (dim == 2)
             {
               FaceDescriptor fd;
               fd.SetBCName(regionnames.Last());
               fd.SetBCProperty(idx);
               self.AddFaceDescriptor(fd);
             }
           return idx;
         }, py::arg("name"), py::arg("dim"))

    .def ("GetRegionNames", [] (Mesh & self, optional<int> optdim, optional<int> optcodim)
          {
            int codim;
            if (optdim)
              codim = self.GetDimension() - *optdim;
            else if (optcodim)
              codim = *optcodim;
            else
              throw Exception("either 'dim' or 'codim' must be specified");
            
            Array<string*> & codimnames = self.GetRegionNamesCD (codim);
            
            std::vector<string> names;
            for (auto name : codimnames)
              {
                if (name)
                  names.push_back(*name);
                else
                  names.push_back("");
              }
            return names;               
          }, py::arg("dim")=nullopt, py::arg("codim")=nullopt)
    
    .def ("SetBCName", &Mesh::SetBCName)
    .def ("GetBCName", FunctionPointer([](Mesh & self, int bc)->string 
                                       { return self.GetBCName(bc); }))
    .def ("SetMaterial", &Mesh::SetMaterial)
    .def ("GetMaterial", FunctionPointer([](Mesh & self, int domnr)
                                         { return string(self.GetMaterial(domnr)); }))

    .def ("GetCD2Name", &Mesh::GetCD2Name)
    .def ("SetCD2Name", &Mesh::SetCD2Name)

    .def ("GetCD3Name", &Mesh::GetCD3Name)
    .def ("SetCD3Name", &Mesh::SetCD3Name)
    .def ("SplitFacesByAdjacentDomains", &Mesh::SplitFacesByAdjacentDomains)
    .def ("GetSubMesh", &Mesh::GetSubMesh, py::arg("domains")="", py::arg("faces")="")
    .def("GetIdentifications", [](Mesh & self) -> py::list
         {
           py::list points;
           for(const auto& pair : self.GetIdentifications().GetIdentifiedPoints())
             {
               // py::tuple pnts = py::make_tuple(pair.first.I1(), pair.first.I2());
               
               auto [pi1, pi2] = get<0> (pair.first);
               py::tuple pnts = py::make_tuple(pi1, pi2);
               points.append(pnts);
             }
           return points;
         })
    .def ("AddPointIdentification", [](Mesh & self, py::object pindex1, py::object pindex2, int identnr, Identifications::ID_TYPE type)
                           {
			     if(py::extract<PointIndex>(pindex1).check() && py::extract<PointIndex>(pindex2).check())
			       {
				 self.GetIdentifications().Add (py::extract<PointIndex>(pindex1)(), py::extract<PointIndex>(pindex2)(), identnr);
				 self.GetIdentifications().SetType(identnr, type); // type = 2 ... periodic
			       }
                           },
          //py::default_call_policies(),
          py::arg("pid1"),
           py::arg("pid2"),
           py::arg("identnr"),
           py::arg("type")=Identifications::PERIODIC)
    .def("IdentifyPeriodicBoundaries", &Mesh::IdentifyPeriodicBoundaries,
         py::arg("identification_name"), py::arg("face1"), py::arg("mapping"),
py::arg("point_tolerance") = -1.)
    .def("GetCurveOrder", [] (Mesh & self)
          {
            return self.GetCurvedElements().GetOrder();
          })
    .def("GetNrIdentifications", [](Mesh& self)
                                 {
                                   return self.GetIdentifications().GetMaxNr();
                                 })
    .def ("CalcLocalH", &Mesh::CalcLocalH)
    .def ("SetMaxHDomain", [] (Mesh& self, py::list maxhlist)
          {
            NgArray<double> maxh;
            for(auto el : maxhlist)
              maxh.Append(py::cast<double>(el));
            self.SetMaxHDomain(maxh);
          })
    .def ("GenerateVolumeMesh", 
          [](Mesh & self, MeshingParameters* pars,
             py::kwargs kwargs)
           {
             MeshingParameters mp;
             if(pars) mp = *pars;
             CreateMPfromKwargs(mp, kwargs);
             py::gil_scoped_release gil_release;
             MeshVolume (mp, self);
             OptimizeVolume (mp, self);
           }, py::arg("mp")=nullptr,
          meshingparameter_description.c_str())

    .def ("OptimizeVolumeMesh", [](Mesh & self, MeshingParameters* pars)
          {
            MeshingParameters mp;
            if(pars) mp = *pars;
            else mp.optsteps3d = 5;
            OptimizeVolume (mp, self);
          }, py::arg("mp"), py::call_guard<py::gil_scoped_release>())
    .def("SetLocalH",[](Mesh& self, shared_ptr<LocalH> localh, int layer)
         {
           self.SetLocalH(localh, layer);
         }, py::arg("localh"), py::arg("layer")=1)
    .def("GetLocalH", &Mesh::GetLocalH)
    .def ("OptimizeMesh2d", [](Mesh & self, MeshingParameters* pars, int faceindex)
          {
            self.CalcLocalH(0.5);
            MeshingParameters mp;
            if(pars) mp = *pars;
            else mp.optsteps2d = 5;
            if(!self.GetGeometry())
              throw Exception("Cannot optimize surface mesh without geometry!");
            Optimize2d (self, mp, faceindex);
          }, py::arg("mp")=nullptr, py::arg("faceindex")=0, py::call_guard<py::gil_scoped_release>())
    
    .def ("Refine", FunctionPointer
          ([](Mesh & self, bool adaptive)
           {
             if (!adaptive)
               {
                 self.GetGeometry()->GetRefinement().Refine(self);
                 self.UpdateTopology();
               }
             else
               {
                 BisectionOptions biopt;
                 biopt.usemarkedelements = 1;
                 biopt.refine_p = 0;
                 biopt.refine_hp = 0;
                 /*
                   biopt.onlyonce = onlyonce;
                   if (reftype == NG_REFINE_P)
                   biopt.refine_p = 1;
                   if (reftype == NG_REFINE_HP)
                   biopt.refine_hp = 1;
                 */
                 self.GetGeometry()->GetRefinement().Bisect (self, biopt);
                 self.UpdateTopology();
                 self.GetCurvedElements().SetIsHighOrder (false);
               }
           }), py::arg("adaptive")=false, py::call_guard<py::gil_scoped_release>())
    
    .def("ZRefine", &Mesh::ZRefine)
    .def("Split2Tets", &Mesh::Split2Tets)
    .def ("SplitAlfeld", FunctionPointer
          ([](Mesh & self)
           {
            NgLock meshlock (self.MajorMutex(), true);
            Refinement & ref = const_cast<Refinement&> (self.GetGeometry()->GetRefinement());
            ::netgen::HPRefinement (self, &ref, SPLIT_ALFELD, 1, 0.5, true, true);
           }
           ), py::call_guard<py::gil_scoped_release>())
    .def ("SplitPowellSabin", FunctionPointer
          ([](Mesh & self)
           {
            NgLock meshlock (self.MajorMutex(), true);
            Refinement & ref = const_cast<Refinement&> (self.GetGeometry()->GetRefinement());
            ::netgen::HPRefinement (self, &ref, SPLIT_POWELL, 1, 0.5, true, true);
           }
           ), py::call_guard<py::gil_scoped_release>())
    .def ("SecondOrder", [](Mesh & self)
          {
            self.GetGeometry()->GetRefinement().MakeSecondOrder(self);
          })
    
    .def ("Curve", [](Mesh & self, int order)
          {
            self.BuildCurvedElements(order);
          })
    .def ("CalcElementMapping", [](Mesh & self, py::buffer refpts1, py::buffer physpts1)
          {
            auto refpts = refpts1.cast<py::array_t<double_t, py::array::c_style | py::array::forcecast>>();
            auto physpts = physpts1.cast<py::array_t<double_t, py::array::c_style | py::array::forcecast>>();            
            
            py::buffer_info ref_info = refpts.request();
            py::buffer_info phys_info = physpts.request();            
            double * ref_ptr = static_cast<double*> (ref_info.ptr);
            double * phys_ptr = static_cast<double*> (phys_info.ptr);
            
            if (ref_info.ndim != 2)
              throw std::runtime_error("Reference points need buffer of dimension 2");
            if (phys_info.ndim != 3)
              throw std::runtime_error("Physical points need buffer of dimension 3");

            /*
            cout << "ref_info.shape = " << FlatArray(2, &ref_info.shape[0]) << endl;
            cout << "ref_info.stride = " << FlatArray(2, &ref_info.strides[0]) << endl;
            cout << "phys_info.shape = " << FlatArray(3, &phys_info.shape[0]) << endl;
            cout << "phys_info.stride = " << FlatArray(3, &phys_info.strides[0]) << endl;
            */
            
            size_t npts = ref_info.shape[0];
            size_t dim = ref_info.shape[1];
            // size_t nel = phys_info.shape[0];
            size_t dim_phys = phys_info.shape[2];            

            size_t stride_refpts = ref_info.strides[0]/sizeof(double);
            size_t stride_physels = phys_info.strides[0]/sizeof(double);
            size_t stride_physpts = phys_info.strides[1]/sizeof(double);
            
            auto & curved = self.GetCurvedElements();

            if (dim == 2)  // mapping of 2D elements
              {
                for (SurfaceElementIndex i = 0; i < self.GetNSE(); i++)
                  for (size_t j = 0; j < npts; j++)
                    {
                      Point<2> xref;
                      Point<3> xphys;
                      for (size_t k = 0; k < 2; k++)
                        xref(k) = ref_ptr[j*stride_refpts+k];
                      curved.CalcSurfaceTransformation(xref, i, xphys);
                      for (size_t k = 0; k < dim_phys; k++)
                        phys_ptr[i*stride_physels+j*stride_physpts+k] = xphys(k);
                    }
              }
            
            if (dim == 3)  // mapping of 3D elements
              {
                for (ElementIndex i = 0; i < self.GetNE(); i++)
                  for (size_t j = 0; j < npts; j++)
                    {
                      Point<3> xref;
                      Point<3> xphys;
                      for (size_t k = 0; k < 3; k++)
                        xref(k) = ref_ptr[j*stride_refpts+k];
                      curved.CalcElementTransformation(xref, i, xphys);
                      for (size_t k = 0; k < 3; k++)
                        phys_ptr[i*stride_physels+j*stride_physpts+k] = xphys(k);
                    }
              }
          })
    
    .def ("GetGeometry", [](Mesh & self) { return self.GetGeometry(); })
    .def ("SetGeometry", [](Mesh & self, shared_ptr<NetgenGeometry> geo)
           {
             self.SetGeometry(geo);
           })

    /*
    .def ("SetGeometry", FunctionPointer
          ([](Mesh & self, shared_ptr<CSGeometry> geo)
           {
             self.SetGeometry(geo);
           }))
    */
    
    .def ("BuildSearchTree", &Mesh::BuildElementSearchTree,py::call_guard<py::gil_scoped_release>(),
          py::arg("dim")=3)

    .def ("BoundaryLayer2", GenerateBoundaryLayer2, py::arg("domain"), py::arg("thicknesses"), py::arg("make_new_domain")=true, py::arg("boundaries")=Array<int>{})
    .def ("BoundaryLayer", [](Mesh & self, variant<string, int, std::vector<int>> boundary,
                              variant<double, std::vector<double>> thickness,
                              optional<variant<string, map<string, string>>> material,
                              variant<string, int, std::vector<int>> domain, bool outside,
                              optional<variant<string, std::vector<int>>> project_boundaries,
                              bool grow_edges, bool limit_growth_vectors,
                              bool sides_keep_surfaceindex,
                              bool disable_curving)
           {
             throw Exception("Call syntax has changed! Pass a list of BoundaryLayerParameters to the GenerateMesh call instead: \ngeo.GenerateMesh(..., boundary_layers=[BoundaryLayerParameters(...), BoundaryLayerParameters(...), ...])");
             BoundaryLayerParameters blp;
             blp.boundary = boundary;
             blp.thickness = thickness;
             blp.new_material = material;
             blp.domain = domain;
             blp.outside = outside;
             blp.project_boundaries = project_boundaries;
             blp.grow_edges = grow_edges;
             blp.limit_growth_vectors = limit_growth_vectors;
             blp.sides_keep_surfaceindex = sides_keep_surfaceindex;
             blp.disable_curving = disable_curving;
             GenerateBoundaryLayer (self, blp);
             self.UpdateTopology();
           }, py::arg("boundary"), py::arg("thickness"), py::arg("material")=nullopt,
          py::arg("domains") = ".*", py::arg("outside") = false,
          py::arg("project_boundaries")=nullopt, py::arg("grow_edges")=true, py::arg("limit_growth_vectors") = false, py::arg("sides_keep_surfaceindex")=false,
          py::arg("disable_curving")=true, "Add boundary layer to mesh. see help(BoundaryLayerParameters) for details.")

    .def_static ("EnableTableClass", [] (string name, bool set)
          {
            MeshTopology::EnableTableStatic(name, set);
          },
          py::arg("name"), py::arg("set")=true)
    .def ("EnableTable", [] (Mesh & self, string name, bool set)
          {
            const_cast<MeshTopology&>(self.GetTopology()).EnableTable(name, set);
          },
          py::arg("name"), py::arg("set")=true)
    
    .def ("Scale", [](Mesh & self, double factor)
          {
            for(auto & pnt : self.Points())
	      pnt.Scale(factor);
          })
    .def ("Copy", [](Mesh & self)
          {
            auto m2 = make_shared<Mesh> ();
            *m2 = self;
            return m2;
          })
    .def ("CalcMinMaxAngle", [](Mesh & self, double badel_limit)
          {
            double values[4];
            self.CalcMinMaxAngle (badel_limit, values);
            py::dict res;
            res["trig"] = py::make_tuple( values[0], values[1] );
            res["tet"] = py::make_tuple( values[2], values[3] );
            return res;
          }, py::arg("badelement_limit")=175.0)
    .def ("Update", [](Mesh & self)
          {
            self.SetNextTimeStamp();
          })
    .def ("CalcTotalBadness", &Mesh::CalcTotalBad)
    .def ("GetQualityHistogram", &Mesh::GetQualityHistogram)
    .def("Mirror", &Mesh::Mirror)
    .def("_getVertices", [](Mesh & self)
          {
            // std::vector<float> verts(3*self.GetNV());
            Array<float> verts(3*self.GetNV());
            ParallelForRange( self.GetNV(), [&](auto myrange) {
                const auto & points = self.Points();
                for(auto i : myrange)
                {
                    auto p = points[PointIndex::BASE+i];
                    auto * v = &verts[3*i];
                    for(auto k : Range(3))
                        v[k] = p[k];
                } });
            return verts;
          })
    .def("_getSegments", [](Mesh & self)
          {
            // std::vector<int> output;
            // output.resize(2*self.GetNSeg());
            Array<int> output(2*self.GetNSeg());
            ParallelForRange( self.GetNSeg(), [&](auto myrange) {
                const auto & segs = self.LineSegments();
                for(auto i : myrange)
                {
                    const auto & seg = segs[i];
                    for(auto k : Range(2))
                      output[2*i+k] = seg[k]-IndexBASE<PointIndex>();
                } });
            return output;
          })
    .def("_getWireframe", [](Mesh & self)
          {
            const auto & topo = self.GetTopology();
            size_t n = topo.GetNEdges();
            /*
            std::vector<int> output;
            output.resize(2*n);
            */
            Array<int> output(2*n);
            ParallelForRange( n, [&](auto myrange) {
                for(auto i : myrange)
                {
                  // PointIndex p0,p1;
                  // topo.GetEdgeVertices(i+1, p0, p1);
                  auto [p0,p1] = topo.GetEdgeVertices(i);
                    output[2*i] = p0-IndexBASE<PointIndex>();
                    output[2*i+1] = p1-IndexBASE<PointIndex>();
                } });
            return output;
          })
    .def("_get2dElementsAsTriangles", [](Mesh & self)
          {
            /*
            std::vector<int> trigs;
            trigs.resize(3*self.GetNSE());
            */
            Array<int> trigs(3*self.GetNSE());
            ParallelForRange( self.GetNSE(), [&](auto myrange) {
                const auto & surfels = self.SurfaceElements();
                for(auto i : myrange)
                {
                    const auto & sel = surfels[i];
                    auto * trig = &trigs[3*i];
                    for(auto k : Range(3))
                        trig[k] = sel[k]-IndexBASE<PointIndex>();
                        // todo: quads (store the second trig in thread-local extra array, merge them at the end (mutex)
                } });
            return trigs;
          })
    .def("_get3dElementsAsTets", [](Mesh & self)
        {
          // std::vector<int> tets;
          // tets.resize(4*self.GetNE());

            Array<int> tets(4*self.GetNE());
            ParallelForRange( self.GetNE(), [&](auto myrange) {
                const auto & els = self.VolumeElements();
                for(auto i : myrange)
                {
                    const auto & el = els[i];
                    auto * trig = &tets[4*i];
                    for(auto k : Range(4))
                        trig[k] = el[k]-IndexBASE<PointIndex>();
                        // todo: prisms etc (store the extra tets in thread-local extra array, merge them at the end (mutex)
                } });
            return tets;
        })
    ;

  string import_docu = "Import mesh from other file format. Leaving format parameter empty guesses based on file extension.\nSupported formats are:\n";
  UserFormatRegister::IterateFormats([&](auto & e) {
      string s = '\t'+e.format+"\t("+e.extensions[0];
      for(auto & ext : e.extensions.Range(1, e.extensions.Size()))
        s += ", "+ext;
      s += ")\n";
      import_docu += s;
  }, true);

  m.def("ReadMedit", [](const string& filename) {
          map<tuple<int, int>, int> index_map;
          auto mesh = make_shared<Mesh>();
          ReadMeditFormat(*mesh, filename, index_map);
          return py::make_tuple(mesh, index_map);
  });
  m.def("WriteMedit", [](const Mesh& mesh, const string& filename) {
          map<tuple<int,int>, int> index_map;
          WriteMeditFormat(mesh, filename, index_map);
          return index_map;
  });
  m.def("ImportMesh", [](const string& filename, const string & format)
                      {
                        auto mesh = make_shared<Mesh>();
                        ReadUserFormat(*mesh, filename, format);
                        return mesh;
                      }, py::arg("filename"), py::arg("format")="", import_docu.c_str());
  py::enum_<MESHING_STEP>(m,"MeshingStep")
    .value("ANALYSE", MESHCONST_ANALYSE)
    .value("MESHEDGES", MESHCONST_MESHEDGES)
    .value("MESHSURFACE", MESHCONST_OPTSURFACE)
    .value("MESHVOLUME", MESHCONST_OPTVOLUME)
    ;
         
  typedef MeshingParameters MP;
  auto mp = py::class_<MP> (m, "MeshingParameters")
    .def(py::init<>())
            .def(py::init([](MeshingParameters* other, py::kwargs kwargs)
                  {
                    MeshingParameters mp;
                    if(other) mp = *other;
                    CreateMPfromKwargs(mp, kwargs, false);
                    return mp;
                  }), py::arg("mp")=nullptr, meshingparameter_description.c_str())
    .def("__str__", &ToString<MP>)
    .def("RestrictH", [](MP & mp, double x, double y, double z, double h, int layer)
          {
            mp.meshsize_points.Append ( MeshingParameters::MeshSizePoint(Point<3> (x,y,z), h, layer));
          }, py::arg("x"), py::arg("y"), py::arg("z"), py::arg("h"), py::arg("layer")=1
         )
    .def("RestrictH", [](MP & mp, const Point<3>& p, double h, int layer)
    {
      mp.meshsize_points.Append ({p, h, layer});
    }, py::arg("p"), py::arg("h"), py::arg("layer")=1)
    .def("RestrictHLine", [](MP& mp, const Point<3>& p1, const Point<3>& p2,
                             double maxh, int layer)
    {
      int steps = int(Dist(p1, p2) / maxh) + 2;
      auto v = p2 - p1;
      for (int i = 0; i <= steps; i++)
        {
          mp.meshsize_points.Append({p1 + double(i)/steps * v, maxh, layer});
        }
    }, py::arg("p1"), py::arg("p2"), py::arg("maxh"), py::arg("layer")=1)
    ;

  m.def("SetTestoutFile", FunctionPointer ([] (const string & filename)
                                             {
                                               delete testout;
                                               testout = new ofstream (filename);
                                             }));

  m.def("SetMessageImportance", FunctionPointer ([] (int importance)
                                                   {
                                                     int old = printmessage_importance;
                                                     printmessage_importance = importance;
                                                     return old;
                                                   }));

  py::class_<DebugParameters> (m, "_DebugParameters")
      .def_readwrite("debugoutput", &DebugParameters::debugoutput)
      .def_readwrite("slowchecks", &DebugParameters::slowchecks)
      .def_readwrite("haltsuccess", &DebugParameters::haltsuccess)
      .def_readwrite("haltnosuccess", &DebugParameters::haltnosuccess)
      .def_readwrite("haltlargequalclass", &DebugParameters::haltlargequalclass)
      .def_readwrite("haltsegment", &DebugParameters::haltsegment)
      .def_readwrite("haltnode", &DebugParameters::haltnode)
      .def_readwrite("haltsegmentp1", &DebugParameters::haltsegmentp1)
      .def_readwrite("haltsegmentp2", &DebugParameters::haltsegmentp2)
      .def_readwrite("haltexistingline", &DebugParameters::haltexistingline)
      .def_readwrite("haltoverlap", &DebugParameters::haltoverlap)
      .def_readwrite("haltface", &DebugParameters::haltface)
      .def_readwrite("haltfacenr", &DebugParameters::haltfacenr)
      .def_readwrite("write_mesh_on_error", &DebugParameters::write_mesh_on_error)
      ;

  m.attr("debugparam") = py::cast(&debugparam);

  py::class_<BoundaryLayerParameters>(m, "BoundaryLayerParameters")
    .def(py::init([]( 
        std::variant<string, int, std::vector<int>> boundary,
        std::variant<double, std::vector<double>> thickness,
        std::optional<std::variant<string, std::map<string, string>>> new_material,
        std::variant<string, int, std::vector<int>> domain,
        bool outside,
        std::optional<std::variant<string, std::vector<int>>> project_boundaries,
        bool grow_edges,
        bool limit_growth_vectors,
        std::optional<bool> sides_keep_surfaceindex,
        bool disable_curving)
        {
          BoundaryLayerParameters blp;
          blp.boundary = boundary;
          blp.thickness = thickness;
          blp.new_material = new_material;
          blp.domain = domain;
          blp.outside = outside;
          blp.project_boundaries = project_boundaries;
          blp.grow_edges = grow_edges;
          blp.limit_growth_vectors = limit_growth_vectors;
          blp.sides_keep_surfaceindex = sides_keep_surfaceindex;
          blp.disable_curving = disable_curving;
          return blp;
        }),
           py::arg("boundary"), py::arg("thickness"), py::arg("new_material")=nullopt,
           py::arg("domain") = ".*", py::arg("outside") = false,
           py::arg("project_boundaries")=nullopt, py::arg("grow_edges")=true,
           py::arg("limit_growth_vectors") = false, py::arg("sides_keep_surfaceindex")=nullopt,
           py::arg("disable_curving")=true,
           R"delimiter(
Add boundary layer to mesh.

Parameters
----------

boundary : string or int
  Boundary name or number.

thickness : float or List[float]
  Thickness of boundary layer(s).

material : str or List[str]
  Material name of boundary layer(s).

domain : str or int
  Regexp for domain boundarylayer is going into.

outside : bool = False
  If true add the layer on the outside

grow_edges : bool = False
  Grow boundary layer over edges.

project_boundaries : Optional[str] = None
  Project boundarylayer to these boundaries if they meet them. Set
  to boundaries that meet boundarylayer at a non-orthogonal edge and
  layer-ending should be projected to that boundary.

)delimiter")
    .def(py::init([]( const py::dict & d ) {
      try {
        // Call other constructor with named arguments by unpacking the dictionary
        py::object cls = py::type::of<BoundaryLayerParameters>();
        return cls(**d).cast<BoundaryLayerParameters>();
      }
      catch (py::error_already_set & e) {
        cerr << "Error creating BoundaryLayerParameters from dict:" << endl;
        cerr << e.what() << endl;
        throw;
      }
    }))
    .def_readwrite("boundary", &BoundaryLayerParameters::boundary)
    .def_readwrite("thickness", &BoundaryLayerParameters::thickness)
    .def_readwrite("new_material", &BoundaryLayerParameters::new_material)
    .def_readwrite("domain", &BoundaryLayerParameters::domain)
    .def_readwrite("outside", &BoundaryLayerParameters::outside)
    .def_readwrite("project_boundaries", &BoundaryLayerParameters::project_boundaries)
    .def_readwrite("grow_edges", &BoundaryLayerParameters::grow_edges)
    .def_readwrite("limit_growth_vectors", &BoundaryLayerParameters::limit_growth_vectors)
    .def_readwrite("sides_keep_surfaceindex", &BoundaryLayerParameters::sides_keep_surfaceindex)
    .def_readwrite("disable_curving", &BoundaryLayerParameters::disable_curving)
    ;
  py::implicitly_convertible<py::dict, BoundaryLayerParameters>();

#ifdef NG_CGNS
  m.def("ReadCGNSFile", &ReadCGNSFile, py::arg("filename"), py::arg("base")=1, "Read mesh and solution vectors from CGNS file");
  m.def("WriteCGNSFile", &WriteCGNSFile, py::arg("mesh"), py::arg("filename"), py::arg("names"), py::arg("values"), py::arg("locations"),
      R"(Write mesh and solution vectors to CGNS file, possible values for locations:
      Vertex     = 0
      EdgeCenter = 1
      FaceCenter = 2
      CellCenter = 3
      )");
#endif // NG_CGNS

    py::class_<SurfaceGeometry, NetgenGeometry, shared_ptr<SurfaceGeometry>> (m, "SurfaceGeometry")
    .def(py::init<>())
    .def(py::init([](py::object pyfunc)
                  {
                    std::function<Vec<3> (Point<2>)> func = [pyfunc](Point<2> p)
                                                    {
                                                      py::gil_scoped_acquire aq;
                                                      py::tuple pyres = py::extract<py::tuple>(pyfunc(p[0],p[1],0.0)) ();
                                                      return Vec<3>(py::extract<double>(pyres[0])(),py::extract<double>(pyres[1])(),py::extract<double>(pyres[2])());
                                                    };
                    auto geo = make_shared<SurfaceGeometry>(func);
                    return geo;
                  }), py::arg("mapping"))
    .def(NGSPickle<SurfaceGeometry>())
    .def("GenerateMesh", [](shared_ptr<SurfaceGeometry> geo,
                            bool quads, int nx, int ny, bool flip_triangles, py::list py_bbbpts, py::list py_bbbnames, py::list py_hppnts, py::dict/*list*/ py_hpbnd, py::dict py_layers)
           {
             if (py::len(py_bbbpts) != py::len(py_bbbnames))
               throw Exception("In SurfaceGeometry::GenerateMesh bbbpts and bbbnames do not have same lengths.");
             Array<Point<3>> bbbpts(py::len(py_bbbpts));
             Array<string> bbbname(py::len(py_bbbpts));
             Array<Point<3>> hppnts(py::len(py_hppnts));
             Array<float> hppntsfac(py::len(py_hppnts));
             Array<string> hpbnd(py::len(py_hpbnd));
             Array<float> hpbndfac(py::len(py_hpbnd));
             for(int i = 0; i<py::len(py_bbbpts);i++)
		 {
                   py::tuple pnt = py::extract<py::tuple>(py_bbbpts[i])();
                   bbbpts[i] = Point<3>(py::extract<double>(pnt[0])(),py::extract<double>(pnt[1])(),py::extract<double>(pnt[2])());
                   bbbname[i] = py::extract<string>(py_bbbnames[i])();
                 }
             for(int i = 0; i<py::len(py_hppnts);i++)
		 {
                   py::tuple pnt = py::extract<py::tuple>(py_hppnts[i])();
                   hppnts[i] = Point<3>(py::extract<double>(pnt[0])(),py::extract<double>(pnt[1])(),py::extract<double>(pnt[2])());
                   hppntsfac[i] = py::extract<double>(pnt[3])();
                 }

	     int ii=0;
             for(auto val : py_hpbnd)
               {
                 hpbnd[ii] = py::cast<string>(val.first);
		 hpbndfac[ii] = py::cast<float>(val.second);
		 ii++;
	       }

             
             Array<double> layer_thickness[4];
             bool layer_quad = false;

             for(auto val : py_layers)
               {
		 int index = -1;
                 if (py::cast<string>(val.first) == "left") index = 0;
                 else if (py::cast<string>(val.first) == "top") index = 3;
                 else if (py::cast<string>(val.first) == "right") index = 2;
                 else if (py::cast<string>(val.first) == "bottom") index = 1;
		 else if (py::cast<string>(val.first) == "quads") layer_quad = py::cast<bool>(val.second);
		 else throw Exception("Unknown parameter " + string(py::cast<string>(val.first)));
		 if (index < 0) continue;

		 auto list = py::cast<py::list>(val.second);
		 layer_thickness[index] = Array<double>(py::len(list));
		 for (size_t i = 0; i < py::len(list); i++)
		   layer_thickness[index][i] = py::cast<double>(list[i]);
               }
                   
             auto mesh = make_shared<Mesh>();
             SetGlobalMesh (mesh);
             mesh->SetGeometry(geo);
	     ng_geometry = geo;
             auto result = geo->GenerateStructuredMesh (mesh, quads, nx, ny, flip_triangles, bbbpts, bbbname, hppnts, hppntsfac, hpbnd, hpbndfac, layer_thickness, layer_quad);
             if(result != 0)
               throw Exception("SurfaceGeometry: Meshing failed!");
             return mesh;
           }, py::arg("quads")=true, py::arg("nx")=10, py::arg("ny")=10, py::arg("flip_triangles")=false, py::arg("bbbpts")=py::list(), py::arg("bbbnames")=py::list(), py::arg("hppnts")=py::list(), py::arg("hpbnd")=py::dict(), py::arg("boundarylayer")=py::dict());/*, R"raw_string(
      Generate a structured 2D surface mesh

    Parameters:
    
    quads : bool
      If True, a quadrilateral mesh is generated. If False, the quads are split to triangles.

    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    flip_triangles : bool
      If set to True together with quads=False the quads are cut the other way round

    bbbpts : list
      List of points which should be handled as BBBND and are named with bbbnames. The mesh must be constructed in such a way that the bbbpts coincide with generated points.

    bbbnames : list
      List of bbbnd names as strings. Size must coincide with size of bbbpts.

    hppnts : list
      If not None it expects a list of the form [ (px1,py1,pz1, hpref1), (px2,py2,pz2, hpref2), ... ] where px,py,pz are the point coordinates which have to be resolved in the mesh and hpref the refinement factor.

    hpbnd : dict
      If not None it expects a dictionary of the form {"boundaryname" : hpref } where boundaryname in [left, right, top, bottom] and hpref the refinement factor.

    boundarylayer : dict
      If not None it expects a dictionary of the form { "boundaryname" : [t1,...,tn], "quads" : False } where ti denote the thickness of layer i. The number of layers are included in nx/ny. After the layers are placed the remaining number of cells are used to divide the remaining grid uniformly. If quads are set to True quadrilaterals are used inside the boundarylayer. If set False the value of "quads" of the function call is used.
      )raw_string");*/
    ;

    py::class_<ClearSolutionClass> (m, "ClearSolutionClass")
      .def(py::init<>())
      ;
    m.def("SetParallelPickling", [](bool par) { parallel_pickling = par; });
    m.def ("_Redraw",
        ([](bool blocking, double fr)
          {
            static auto last_time = std::chrono::system_clock::now()-std::chrono::seconds(10);
            auto now = std::chrono::system_clock::now();
            double elapsed = std::chrono::duration<double>(now-last_time).count();
            if (blocking || elapsed * fr > 1)
              {
                Ng_Redraw(blocking);
                last_time = std::chrono::system_clock::now();
                return true;
              }
            return false;
          }),
        py::arg("blocking")=false, py::arg("fr") = 25, R"raw_string(
  Redraw all

  Parameters:

  blocking : bool
    input blocking

  fr : double
    input framerate

  )raw_string");
}

PYBIND11_MODULE(libmesh, m) {
  ExportNetgenMeshing(m);
}
#endif





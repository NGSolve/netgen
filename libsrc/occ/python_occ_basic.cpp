#ifdef NG_PYTHON
#ifdef OCCGEOMETRY

#include <general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include <meshing/python_mesh.hpp>
#include <meshing.hpp>

#include "occgeom.hpp"

#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax2d.hxx>
#include <gp_Ax3.hxx>
#include <gp_Trsf.hxx>
#include <gp_GTrsf.hxx>

using namespace netgen;

DLL_HEADER void ExportNgOCCBasic(py::module &m) 
{
  py::class_<gp_Pnt>(m, "gp_Pnt", "3d OCC point")
    .def(py::init([] (py::tuple pnt)
                  {
                    if (py::len(pnt) != 3)
                      throw Exception("need 3-tuple to create gp_Pnt");
                    
                    return gp_Pnt(py::cast<double>(pnt[0]),
                                  py::cast<double>(pnt[1]),
                                  py::cast<double>(pnt[2]));
                  }))
    .def(py::init([] (double x, double y, double z) {
          return gp_Pnt(x, y, z);
        }), py::arg("x"), py::arg("y"), py::arg("z"))
    .def_property("x", [](gp_Pnt&p) { return p.X(); }, [](gp_Pnt&p,double x) { p.SetX(x); })
    .def_property("y", [](gp_Pnt&p) { return p.Y(); }, [](gp_Pnt&p,double y) { p.SetY(y); })
    .def_property("z", [](gp_Pnt&p) { return p.Z(); }, [](gp_Pnt&p,double z) { p.SetZ(z); })
    .def("__str__", [] (const gp_Pnt & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")";
        return str.str();
      })
    .def("__repr__", [] (const gp_Pnt & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")";
        return str.str();
      })

    .def("__sub__", [](gp_Pnt p1, gp_Pnt p2) { return gp_Vec(p2, p1); }) 
    .def("__add__", [](gp_Pnt p, gp_Vec v) { return p.Translated(v); }) // gp_Pnt(p.X()+v.X(), p.Y()+v.Y(), p.Z()+v.Z()); })
    .def("__sub__", [](gp_Pnt p, gp_Vec v) { return p.Translated(-v); }) // gp_Pnt(p.X()-v.X(), p.Y()-v.Y(), p.Z()-v.Z()); })
    .def("__getitem__", [](const gp_Pnt& p, int index)
    {
      if(index == 0)
        return p.X();
      if(index == 1)
        return p.Y();
      if(index == 2)
        return p.Z();
      throw std::out_of_range("Point index must be in range [0,3)!");
    })
    ;
  
  py::class_<gp_Vec>(m, "gp_Vec", "3d OCC vector")
    .def(py::init([] (py::tuple vec)
                  {
                    return gp_Vec(py::cast<double>(vec[0]),
                                  py::cast<double>(vec[1]),
                                  py::cast<double>(vec[2]));
                  }))
    .def(py::init([] (double x, double y, double z) {
          return gp_Vec(x, y, z);
        }), py::arg("x"), py::arg("y"), py::arg("z"))
    .def_property("x", [](gp_Vec&p) { return p.X(); }, [](gp_Vec&p,double x) { p.SetX(x); })
    .def_property("y", [](gp_Vec&p) { return p.Y(); }, [](gp_Vec&p,double y) { p.SetY(y); })
    .def_property("z", [](gp_Vec&p) { return p.Z(); }, [](gp_Vec&p,double z) { p.SetZ(z); })
    .def("Norm", [](const gp_Vec& v)
    { return v.Magnitude(); })
    .def("__str__", [] (const gp_Vec & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")";
        return str.str();
      })
    .def("__repr__", [] (const gp_Vec & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")";
        return str.str();
      })
    .def("__add__", [](gp_Vec v1, gp_Vec v2) { return v1+v2; }) 
    .def("__sub__", [](gp_Vec v1, gp_Vec v2) { return v1-v2; }) 
    .def("__rmul__", [](gp_Vec v, double s) { return s*v; }) 
    .def("__mul__", [](gp_Vec v1, gp_Vec v2) { return v1*v2; })
    .def("__neg__", [](gp_Vec v) { return -v; }) 
    .def("__xor__", [](gp_Vec v1, gp_Vec v2) { return v1^v2; })
    
    .def("__lt__", [](gp_Vec v, double val)
         {
           cout << IM(6) << "vec, lt v - " << netgen::occ2ng(v) << ", val = " << val << endl;
           return DirectionalInterval(v) < val;
         })
    .def("__gt__", [](gp_Vec v, double val)
         {
           cout << IM(6) << "vec, gt v - " << netgen::occ2ng(v) << ", val = " << val << endl;           
           return DirectionalInterval(v) > val;
         })
    .def("__le__", [](gp_Vec v, double val)
        {
          return DirectionalInterval(v) <= val;
        })
    .def("__ge__", [](gp_Vec v, double val)
        {
          return DirectionalInterval(v) >= val;
        })
    ;

  py::class_<gp_Dir>(m, "gp_Dir", "3d OCC direction")
    .def(py::init([] (py::tuple dir)
                  {
                    return gp_Dir(py::cast<double>(dir[0]),
                                  py::cast<double>(dir[1]),
                                  py::cast<double>(dir[2]));
                  }))
    .def(py::init([] (double x, double y, double z) {
          return gp_Dir(x, y, z);
        }), py::arg("x"), py::arg("y"), py::arg("z"))
    .def(py::init<gp_Vec>())
    .def("__str__", [] (const gp_Dir & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")";
        return str.str();
      })
    ;
  
  py::class_<gp_Ax1>(m, "Axis", "an OCC axis in 3d") 
    .def(py::init([](gp_Pnt p, gp_Dir d) {
          return gp_Ax1(p,d);
        }), py::arg("p"), py::arg("d"))
    ;
  py::class_<gp_Ax2>(m, "gp_Ax2")
    .def(py::init([](gp_Pnt p, gp_Dir d) {
          return gp_Ax2(p,d);
        }))
    .def(py::init([](const gp_Ax3 & ax3) {
          return gp_Ax2(ax3.Ax2());
        }))
    ;

  py::class_<gp_Ax3>(m, "Axes", "an OCC coordinate system in 3d")
    .def(py::init([](gp_Pnt p, gp_Dir N, gp_Dir Vx) {
          return gp_Ax3(p,N, Vx);
        }), py::arg("p")=gp_Pnt(0,0,0), py::arg("n")=gp_Vec(0,0,1), py::arg("h")=gp_Vec(1,0,0))
    .def(py::init<gp_Ax2>())
    .def_property("p", [](gp_Ax3 & ax) { return ax.Location(); }, [](gp_Ax3&ax, gp_Pnt p) { ax.SetLocation(p); })
    ;


  py::class_<gp_Pnt2d>(m, "gp_Pnt2d", "2d OCC point")
    .def(py::init([] (py::tuple pnt)
                  {
                    if (py::len(pnt) != 2)
                      throw Exception("need 2-tuple to create gp_Pnt2d");
                    return gp_Pnt2d(py::cast<double>(pnt[0]),
                                    py::cast<double>(pnt[1]));
                  }))
    .def(py::init([] (double x, double y) {
          return gp_Pnt2d(x, y);
        }), py::arg("x"), py::arg("y"))
    .def_property("x", [](gp_Pnt2d&p) { return p.X(); }, [](gp_Pnt2d&p,double x) { p.SetX(x); })
    .def_property("y", [](gp_Pnt2d&p) { return p.Y(); }, [](gp_Pnt2d&p,double y) { p.SetY(y); })
    .def("__str__", [] (const gp_Pnt2d & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ")";
        return str.str();
      })
    .def("__repr__", [] (const gp_Pnt2d & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ")";
        return str.str();
      })

    .def("__sub__", [](gp_Pnt2d p1, gp_Pnt2d p2) { return gp_Vec2d(p1.X()-p2.X(), p1.Y()-p2.Y()); })
    .def("__add__", [](gp_Pnt2d p, gp_Vec2d v) { return p.Translated(v); })
    .def("__sub__", [](gp_Pnt2d p, gp_Vec2d v) { return p.Translated(-v); })
    ;
  
  py::class_<gp_Vec2d>(m, "gp_Vec2d", "2d OCC vector")
    .def(py::init([] (py::tuple vec)
                  {
                    if (py::len(vec) != 2)
                      throw Exception("need 2-tuple to create gp_Vec2d");                    
                    return gp_Vec2d(py::cast<double>(vec[0]),
                                    py::cast<double>(vec[1]));
                  }))
    .def(py::init([] (double x, double y) {
          return gp_Vec2d(x, y);
        }), py::arg("x"), py::arg("y"))
    .def_property("x", [](gp_Vec2d&p) { return p.X(); }, [](gp_Vec2d&p,double x) { p.SetX(x); })
    .def_property("y", [](gp_Vec2d&p) { return p.Y(); }, [](gp_Vec2d&p,double y) { p.SetY(y); })
    .def("__str__", [] (const gp_Vec & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ")";
        return str.str();
      })
    .def("__repr__", [] (const gp_Vec & p) {
        stringstream str;
        str << "(" << p.X() << ", " << p.Y() << ")";
        return str.str();
      })
    .def("__add__", [](gp_Vec2d v1, gp_Vec2d v2) { return v1+v2; }) 
    .def("__sub__", [](gp_Vec2d v1, gp_Vec2d v2) { return v1-v2; })
    .def("__rmul__", [](gp_Vec2d v, double s) { return s*v; })
    .def("__neg__", [](gp_Vec2d v) { return -v; }) 
    .def("__xor__", [](gp_Vec2d v1, gp_Vec2d v2) { return v1^v2; }) 
    ;

  py::class_<gp_Dir2d>(m, "gp_Dir2d", "2d OCC direction")
    .def(py::init([] (py::tuple dir)
                  {
                    if (py::len(dir) != 2)
                      throw Exception("need 2-tuple to create gp_Dir2d");                    
                    return gp_Dir2d(py::cast<double>(dir[0]),
                                    py::cast<double>(dir[1]));
                  }))
    .def(py::init([] (double x, double y) {
          return gp_Dir2d(x, y);
        }), py::arg("x"), py::arg("y"))
    ;

  m.def("Pnt", [](double x, double y) { return gp_Pnt2d(x,y); },
        py::arg("x"), py::arg("y"), "create 2d OCC point");
  m.def("Pnt", [](double x, double y, double z) { return gp_Pnt(x,y,z); },
        py::arg("x"), py::arg("y"), py::arg("z"), "create 3d OCC point");
  m.def("Pnt", [](std::vector<double> p)
        {
          if (p.size() == 2)
            return py::cast(gp_Pnt2d(p[0], p[1]));
          if (p.size() == 3)
            return py::cast(gp_Pnt(p[0], p[1], p[2]));
          throw Exception("OCC-Points only in 2D or 3D");
        }, py::arg("p"), "create 2d or 3d OCC point");

  m.def("Vec", [](double x, double y) { return gp_Vec2d(x,y); },
        py::arg("x"), py::arg("y"), "create 2d OCC point");
  m.def("Vec", [](double x, double y, double z) { return gp_Vec(x,y,z); },
        py::arg("x"), py::arg("y"), py::arg("z"), "create 3d OCC point");        
  m.def("Vec", [](std::vector<double> p)
        {
          if (p.size() == 2)
            return py::cast(gp_Vec2d(p[0], p[1]));
          if (p.size() == 3)
            return py::cast(gp_Vec(p[0], p[1], p[2]));
          throw Exception("OCC-Vecs only in 2D or 3D");
        }, py::arg("v"), "create 2d or 3d OCC vector");          

  m.def("Dir", [](double x, double y) { return gp_Dir2d(x,y); },
        py::arg("x"), py::arg("y"), "create 2d OCC direction");
  m.def("Dir", [](double x, double y, double z) { return gp_Dir(x,y,z); },
        py::arg("x"), py::arg("y"), py::arg("z"), "create 3d OCC direction");        
  m.def("Dir", [](std::vector<double> p)
        {
          if (p.size() == 2)
            return py::cast(gp_Dir2d(p[0], p[1]));
          if (p.size() == 3)
            return py::cast(gp_Dir(p[0], p[1], p[2]));
          throw Exception("OCC-Dirs only in 2D or 3D");
        }, py::arg("d"), "create 2d or 3d OCC direction");                    



  
  py::class_<gp_Ax2d>(m, "gp_Ax2d", "2d OCC coordinate system")
    .def(py::init([](gp_Pnt2d p, gp_Dir2d d) {
          return gp_Ax2d(p,d);
        }), py::arg("p")=gp_Pnt2d(0,0), py::arg("d")=gp_Dir2d(1,0))
    ;

  py::class_<gp_GTrsf>(m, "gp_GTrsf")
    .def(py::init([](const std::vector<double>& mat,
                     const std::vector<double>& vec)
         {
           if(mat.size() != 9)
             throw Exception("Need 9 matrix values for construction of gp_GTrsf");
           if(vec.size() != 3)
             throw Exception("Need 3 vector values for construction of gp_GTrsf");
           gp_GTrsf trafo;
           trafo.SetVectorialPart({ mat[0], mat[1], mat[2],
               mat[3], mat[4], mat[5],
               mat[6], mat[7], mat[8] });
           trafo.SetTranslationPart( { vec[0], vec[1], vec[2] });
           return trafo;
         }), py::arg("mat"), py::arg("vec") = std::vector<double>{ 0., 0., 0. })
    .def("__call__", [] (gp_GTrsf & trafo, const TopoDS_Shape & shape) {
        BRepBuilderAPI_GTransform builder(shape, trafo, true);
        PropagateProperties(builder, shape, occ2ng(trafo));
        return builder.Shape();
      })
    ;
  
  py::class_<gp_Trsf>(m, "gp_Trsf")
    .def(py::init<>())    
    .def("SetMirror", [] (gp_Trsf & trafo, const gp_Ax1 & ax) { trafo.SetMirror(ax); return trafo; })
    .def("Inverted", &gp_Trsf::Inverted)
    .def_static("Translation", [] (const gp_Vec & v) { gp_Trsf trafo; trafo.SetTranslation(v); return trafo; })
    .def_static("Scale", [] (const gp_Pnt & p, double s) { gp_Trsf trafo; trafo.SetScale(p,s); return trafo; })    
    .def_static("Mirror", [] (const gp_Ax1 & ax) { gp_Trsf trafo; trafo.SetMirror(ax); return trafo; })
    .def_static("Rotation", [] (const gp_Ax1 & ax, double ang) { gp_Trsf trafo; trafo.SetRotation(ax, ang*M_PI/180); return trafo; })
    .def_static("Rotation", [] (const gp_Pnt & p, const gp_Dir & d, double ang)
                { gp_Trsf trafo; trafo.SetRotation(gp_Ax1(p,d), ang*M_PI/180); return trafo; })    
    .def_static("Transformation", [] (const gp_Ax3 & ax) { gp_Trsf trafo; trafo.SetTransformation(ax); return trafo; })
    .def_static("Transformation", [] (const gp_Ax3 & from, const gp_Ax3 to)
                { gp_Trsf trafo; trafo.SetTransformation(from, to); return trafo; })
    .def(py::self * py::self)
    .def("__call__", [] (gp_Trsf & trafo, const TopoDS_Shape & shape) {
        BRepBuilderAPI_Transform builder(shape, trafo, true);
        PropagateProperties(builder, shape, occ2ng(trafo));
        return builder.Shape();
      })
    .def("__str__", [](gp_Trsf & trafo)
    {
      stringstream str;
      gp_XYZ xyz = trafo.TranslationPart();
      str << xyz.X() << ", " << xyz.Y() << ", " << xyz.Z();
      return str.str();
    })
    ;

  py::class_<TopLoc_Location>(m, "TopLoc_Location")
    .def(py::init<gp_Trsf>())
    .def("Transformation", [](const TopLoc_Location & loc) { return loc.Transformation(); })
    ;


  py::class_<DirectionalInterval> (m, "DirectionalInterval")
    .def("__str__", [](DirectionalInterval self)
         {
           stringstream str;
           str << "(" << self.minval << ", " << self.maxval << ")";
           return str.str();
         })
         
    .def("__lt__", [](DirectionalInterval i, double val)
         {
           cout << "directionalinterval, lt, imin/max = " << i.minval << " / " << i.maxval << endl;
           return i < val;
         })
    .def("__gt__", [](DirectionalInterval i, double val)
         {
           cout << "directionalinterval, gt, imin/max = " << i.minval << " / " << i.maxval << endl;
           return i > val;
         })
    .def("__and__", [](DirectionalInterval self, DirectionalInterval other) {
        cout << "and of intervals" << endl;
        return self.Intersect(other);
      })
    ;    
  
  py::implicitly_convertible<py::tuple, gp_Pnt>();
  py::implicitly_convertible<py::tuple, gp_Vec>();
  py::implicitly_convertible<py::tuple, gp_Dir>();
  py::implicitly_convertible<gp_Vec, gp_Dir>();  
  py::implicitly_convertible<py::tuple, gp_Pnt2d>();  
  py::implicitly_convertible<py::tuple, gp_Vec2d>();  
  py::implicitly_convertible<py::tuple, gp_Dir2d>();


  py::implicitly_convertible<gp_Ax3, gp_Ax2>();

  m.attr("X") = py::cast(gp_Vec(1,0,0));
  m.attr("Y") = py::cast(gp_Vec(0,1,0));
  m.attr("Z") = py::cast(gp_Vec(0,0,1));
}


#endif // OCCGEOMETRY
#endif // NG_PYTHON
 

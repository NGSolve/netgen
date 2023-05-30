// a wrapper to load netgen-dll into python

#include <iostream>
#include <../general/ngpython.hpp>
#include <core/ngcore_api.hpp>

void NGCORE_API_IMPORT ExportNetgenMeshing(py::module &m);
void NGCORE_API_IMPORT ExportCSG(py::module &m);
void NGCORE_API_IMPORT ExportGeom2d(py::module &m);
void NGCORE_API_IMPORT ExportSTL(py::module &m);
#ifdef OCCGEOMETRY
void NGCORE_API_IMPORT ExportNgOCC(py::module &m);
#endif // OCCGEOMETRY

PYBIND11_MODULE(libngpy, ngpy)
{
  py::module::import("pyngcore");
    py::module meshing = ngpy.def_submodule("_meshing", "pybind meshing module");
    ExportNetgenMeshing(meshing);
    py::module csg = ngpy.def_submodule("_csg", "pybind csg module");
    ExportCSG(csg);
    py::module geom2d = ngpy.def_submodule("_geom2d", "pybind geom2d module");
    ExportGeom2d(geom2d);
    py::module stl = ngpy.def_submodule("_stl", "pybind stl module");
    ExportSTL(stl);
#ifdef OCCGEOMETRY
    py::module NgOCC = ngpy.def_submodule("_NgOCC", "pybind NgOCC module");
    ExportNgOCC(NgOCC);
#endif // OCCGEOMETRY
}

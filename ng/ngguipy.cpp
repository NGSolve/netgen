#include <iostream>
#include <../general/ngpython.hpp>
#include <core/ngcore_api.hpp>

void NGCORE_API_IMPORT ExportMeshVis(py::module &m);
void NGCORE_API_IMPORT ExportCSGVis(py::module &m);
void NGCORE_API_IMPORT ExportSTLVis(py::module &m);
namespace netgen
{
  std::vector<unsigned char> NGCORE_API_IMPORT Snapshot( int w, int h );
}

PYBIND11_MODULE(libngguipy, ngpy)
{
    py::module::import("pyngcore");
    py::module meshvis = ngpy.def_submodule("meshvis", "pybind meshvis module");
    ExportMeshVis(meshvis);
    py::module csgvis = ngpy.def_submodule("csgvis", "pybind csgvis module");
    ExportCSGVis(csgvis);
    py::module stlvis = ngpy.def_submodule("stlvis", "pybind stlvis module");
    ExportSTLVis(stlvis);
    ngpy.def("Snapshot", netgen::Snapshot);
}

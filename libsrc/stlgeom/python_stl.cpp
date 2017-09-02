
#ifdef NG_PYTHON

#include <../general/ngpython.hpp>
#include <stlgeom.hpp>

#ifdef WIN32
   #define DLL_HEADER   __declspec(dllexport)
#endif

using namespace netgen;
namespace netgen
{
  //extern shared_ptr<Mesh> mesh;
  extern shared_ptr<NetgenGeometry> ng_geometry;
}


DLL_HEADER void ExportSTL(py::module & m)
{
  py::class_<STLGeometry,shared_ptr<STLGeometry>> (m,"STLGeometry")
    .def(py::init<>())
    ;
  m.def("LoadSTLGeometry", FunctionPointer([] (const string & filename)
					   {
					     ifstream ist(filename);
					     return shared_ptr<STLGeometry>(STLGeometry::Load(ist));
					   }));
  m.def("GenerateMesh", FunctionPointer([] (shared_ptr<STLGeometry> geo, MeshingParameters &param)
					{
					  auto mesh = make_shared<Mesh>();
					  SetGlobalMesh(mesh);
					  mesh->SetGeometry(geo);
					  ng_geometry = geo;
					  try
					    {
					      geo->GenerateMesh(mesh,param);
					    }
					  catch (NgException ex)
					    {
					      cout << "Caught NgException: " << ex.What() << endl;
					    }
					  return mesh;
					}))
    ;
}

PYBIND11_MODULE(libstl, m) {
  ExportSTL(m);
}

#endif

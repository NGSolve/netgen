#ifdef NG_PYTHON
#ifdef OCCGEOMETRY

#include <../general/ngpython.hpp>

#include <meshing.hpp>
#include <occgeom.hpp>

using namespace netgen;

namespace netgen
{
  extern std::shared_ptr<NetgenGeometry> ng_geometry;
}


DLL_HEADER void ExportNgOCC(py::module &m) 
{
  py::class_<OCCGeometry, shared_ptr<OCCGeometry>> (m, "OCCGeometry", R"raw_string(Use LoadOCCGeometry to load the geometry from a *.step file.)raw_string")
    .def(py::init<>())
    .def("Heal",[](OCCGeometry & self, double tolerance, bool fixsmalledges, bool fixspotstripfaces, bool sewfaces, bool makesolids, bool splitpartitions)
         {
           self.tolerance = tolerance;
           self.fixsmalledges = fixsmalledges;
           self.fixspotstripfaces = fixspotstripfaces;
           self.sewfaces = sewfaces;
           self.makesolids = makesolids;
           self.splitpartitions = splitpartitions;

           self.HealGeometry();
           self.BuildFMap();
         },py::arg("tolerance")=1e-3, py::arg("fixsmalledges")=true, py::arg("fixspotstripfaces")=true, py::arg("sewfaces")=true, py::arg("makesolids")=true, py::arg("splitpartitions")=false,R"raw_string(Heal the OCCGeometry.)raw_string")
      ;
    m.def("LoadOCCGeometry",FunctionPointer([] (const string & filename)
           {
             cout << "load OCC geometry";
             ifstream ist(filename);
             OCCGeometry * instance = new OCCGeometry();
             instance = LoadOCC_STEP(filename.c_str());
             return shared_ptr<OCCGeometry>(instance, NOOP_Deleter);
           }));
  m.def("GenerateMesh", FunctionPointer([] (shared_ptr<OCCGeometry> geo, MeshingParameters &param)
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

PYBIND11_MODULE(libNgOCC, m) {
  ExportNgOCC(m);
}

#endif // OCCGEOMETRY
#endif // NG_PYTHON

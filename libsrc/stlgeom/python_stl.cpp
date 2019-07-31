
#ifdef NG_PYTHON

#include <../general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include <stlgeom.hpp>
#include "../meshing/python_mesh.hpp"

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
  py::class_<STLGeometry,shared_ptr<STLGeometry>, NetgenGeometry> (m,"STLGeometry")
    .def(py::init<>())
    .def(py::init<>([](const string& filename)
                    {
                      ifstream ist(filename);
                      return shared_ptr<STLGeometry>(STLGeometry::Load(ist));
                    }), py::arg("filename"),
      py::call_guard<py::gil_scoped_release>())
    .def(NGSPickle<STLGeometry>())
    .def("_visualizationData", [](shared_ptr<STLGeometry> stl_geo)
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

           surfnames.push_back("stl");
           vertices.reserve(stl_geo->GetNT()*3*3);
           trigs.reserve(stl_geo->GetNT()*4);
           normals.reserve(stl_geo->GetNT()*3*3);
           size_t ii = 0;
           for(int i = 0; i < stl_geo->GetNT(); i++)
             {
               auto& trig = stl_geo->GetTriangle(i+1);
               for(int k = 0; k < 3; k++)
                 {
                   trigs.push_back(ii++);
                   auto& pnt = stl_geo->GetPoint(trig[k]);
                   for (int l = 0; l < 3; l++)
                     {
                       float val = pnt[l];
                       vertices.push_back(val);
                       min[l] = min2(min[l], val);
                       max[l] = max2(max[l], val);
                       normals.push_back(trig.Normal()[l]);
                     }
                 }
               trigs.push_back(0);
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
    .def("GenerateMesh", [] (shared_ptr<STLGeometry> geo,
                             MeshingParameters* pars, py::kwargs kwargs)
                         {
                           MeshingParameters mp;
                           if(pars) mp = *pars;
                           {
                             py::gil_scoped_acquire aq;
                             CreateMPfromKwargs(mp, kwargs);
                           }
                           auto mesh = make_shared<Mesh>();
                           SetGlobalMesh(mesh);
                           mesh->SetGeometry(geo);
                           ng_geometry = geo;
                           geo->GenerateMesh(mesh,mp);
                           return mesh;
                         }, py::arg("mp") = nullptr,
      py::call_guard<py::gil_scoped_release>(),
      meshingparameter_description.c_str())
    ;
  m.def("LoadSTLGeometry", [] (const string & filename)
                           {
                             cout << "WARNING: LoadSTLGeometry is deprecated, use the STLGeometry(filename) constructor instead!" << endl;
                             ifstream ist(filename);
                             return shared_ptr<STLGeometry>(STLGeometry::Load(ist));
                           },py::call_guard<py::gil_scoped_release>());
}

PYBIND11_MODULE(libstl, m) {
  ExportSTL(m);
}

#endif

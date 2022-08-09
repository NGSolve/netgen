
#ifdef NG_PYTHON

#include <../general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include <stlgeom.hpp>
#include "../meshing/python_mesh.hpp"

using namespace netgen;
namespace netgen
{
  extern shared_ptr<Mesh> mesh;
  extern shared_ptr<NetgenGeometry> ng_geometry;
}

static string stlparameter_description = R"delimiter(
STL Specific Meshing Parameters
-------------------------------

yangle: float = 30.
  Angle for edge detection

contyangle: float = 20.
  Edges continue if angle > contyangle

edgecornerangle: float = 60.
  Angle of geometry edge at which the mesher should set a point.

closeedgefac: Optional[float] = 1.
  Factor for meshing close edges, if None it is disabled.

minedgelen: Optional[float] = 0.001
  Minimum edge length to be used for dividing edges to mesh points. If
  None this is disabled.
)delimiter";

void CreateSTLParametersFromKwargs(STLParameters& stlparam, py::dict kwargs)
{
  if(kwargs.contains("yangle"))
    stlparam.yangle = py::cast<double>(kwargs.attr("pop")("yangle"));
  if(kwargs.contains("contyangle"))
    stlparam.contyangle = py::cast<double>(kwargs.attr("pop")("contyangle"));
  if(kwargs.contains("edgecornerangle"))
    stlparam.edgecornerangle = py::cast<double>(kwargs.attr("pop")("edgecornerangle"));
  if(kwargs.contains("chartangle"))
    stlparam.chartangle = py::cast<double>(kwargs.attr("pop")("chartangle"));
  if(kwargs.contains("outerchartangle"))
    stlparam.outerchartangle = py::cast<double>(kwargs.attr("pop")("outerchartangle"));
  if(kwargs.contains("usesearchtree"))
    stlparam.usesearchtree = py::cast<int>(kwargs.attr("pop")("usesearchtree"));
  if(kwargs.contains("atlasfac"))
  {
    auto val = kwargs.attr("pop")("resthatlasfac");
    if(val.is_none())
      stlparam.resthatlasenable = false;
    else
    {
      stlparam.resthatlasenable = true;
      stlparam.resthatlasfac = py::cast<double>(val);
    }
  }
  if(kwargs.contains("atlasminh"))
    stlparam.atlasminh = py::cast<double>(kwargs.attr("pop")("atlasminh"));
  if(kwargs.contains("surfcurvfac"))
  {
    auto val = kwargs.attr("pop")("surfcurvfac");
    if(val.is_none())
      stlparam.resthsurfcurvenable = false;
    else
    {
      stlparam.resthsurfcurvenable = true;
      stlparam.resthsurfcurvfac = py::cast<double>(val);
    }
  }
  if(kwargs.contains("chartdistfac"))
  {
    auto val = kwargs.attr("pop")("chartdistfac");
    if(val.is_none())
      stlparam.resthchartdistenable = false;
    else
    {
      stlparam.resthchartdistenable = true;
      stlparam.resthchartdistfac = py::cast<double>(val);
    }
  }
  if(kwargs.contains("edgeanglefac"))
  {
    auto val = kwargs.attr("pop")("edgeanglefac");
    if(val.is_none())
      stlparam.resthedgeangleenable = false;
    else
    {
      stlparam.resthedgeangleenable = true;
      stlparam.resthedgeanglefac = py::cast<double>(val);
    }
  }
  if(kwargs.contains("surfmeshcurvfac"))
  {
    auto val = kwargs.attr("pop")("surfmeshcurvfac");
    if(val.is_none())
      stlparam.resthsurfmeshcurvenable = false;
    else
    {
      stlparam.resthsurfmeshcurvenable = true;
      stlparam.resthsurfmeshcurvfac = py::cast<double>(val);
    }
  }
  if(kwargs.contains("linelengthfac"))
  {
    auto val = kwargs.attr("pop")("linelengthfac");
    if(val.is_none())
      stlparam.resthlinelengthenable = false;
    else
    {
      stlparam.resthlinelengthenable = true;
      stlparam.resthlinelengthfac = py::cast<double>(val);
    }
  }
  if(kwargs.contains("recalc_h_opt"))
    stlparam.recalc_h_opt = py::cast<bool>(kwargs.attr("pop")("recalc_h_opt"));
}


NGCORE_API_EXPORT void ExportSTL(py::module & m)
{
  py::class_<STLGeometry,shared_ptr<STLGeometry>, NetgenGeometry> (m,"STLGeometry")
    .def(py::init<>())
    .def(py::init<>([](const string& filename, bool surface)
                    {
                      ifstream ist(filename);
                      return shared_ptr<STLGeometry>(STLGeometry::Load(ist,
                                                                       surface));
                    }), py::arg("filename"), py::arg("surface")=false,
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
                           STLParameters stlparam;
                           { py::gil_scoped_acquire aq;
                             if(pars)
                             {
                               auto mp_flags = pars->geometrySpecificParameters;
                               auto mp_kwargs = CreateDictFromFlags(mp_flags);
                               CreateSTLParametersFromKwargs(stlparam, mp_kwargs);
                               mp = *pars;
                             }
                             CreateSTLParametersFromKwargs(stlparam, kwargs);
                             CreateMPfromKwargs(mp, kwargs); // this will throw if any kwargs are not passed
                           }
                           auto mesh = make_shared<Mesh>();
                           mesh->SetGeometry(geo);
                           ng_geometry = geo;
                           SetGlobalMesh(mesh);
                           auto result = STLMeshingDummy(geo.get(), mesh, mp, stlparam);
                           if(result != 0)
                             {
                               netgen::mesh = mesh;
                               throw Exception("Meshing failed!");
                             }

                           return mesh;
                         }, py::arg("mp") = nullptr,
      py::call_guard<py::gil_scoped_release>(),
         (meshingparameter_description + stlparameter_description).c_str())
    .def("Draw", FunctionPointer
         ([] (shared_ptr<STLGeometry> self)
          {
             ng_geometry = self;
          })
         )
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

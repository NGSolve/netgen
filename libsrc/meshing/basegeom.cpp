#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

  DLL_HEADER GeometryRegisterArray geometryregister; 
  //DLL_HEADER NgArray<GeometryRegister*> geometryregister;

  GeometryRegister :: ~GeometryRegister()
  { ; }

  void NetgenGeometry :: Analyse(Mesh& mesh,
                                 const MeshingParameters& mparam)
  {
    static Timer t1("SetLocalMeshsize"); RegionTimer regt(t1);
    mesh.SetGlobalH(mparam.maxh);
    mesh.SetMinimalH(mparam.minh);

    mesh.SetLocalH(bounding_box.PMin(), bounding_box.PMax(),
                   mparam.grading);

    if(mparam.uselocalh)
      RestrictLocalMeshsize(mesh, mparam);
    mesh.LoadLocalMeshSize(mparam.meshsizefilename);
  }

  void NetgenGeometry :: OptimizeSurface(Mesh& mesh, const MeshingParameters& mparam)
  {
    const auto savetask = multithread.task;
    multithread.task = "Optimizing surface";

    static Timer timer_opt2d("Optimization 2D");
    RegionTimer reg(timer_opt2d);
    auto meshopt = MeshOptimize2d(mesh);
    for(auto i : Range(mparam.optsteps2d))
      {
        PrintMessage(2, "Optimization step ", i);
        for(auto optstep : mparam.optimize2d)
          {
            switch(optstep)
              {
              case 's':
                meshopt.EdgeSwapping(0);
                break;
              case 'S':
                meshopt.EdgeSwapping(1);
                break;
              case 'm':
                meshopt.ImproveMesh(mparam);
                break;
              case 'c':
                meshopt.CombineImprove();
                break;
              }
          }
      }
    mesh.CalcSurfacesOfNode();
    mesh.Compress();
    multithread.task = savetask;
  }
  
  shared_ptr<NetgenGeometry> GeometryRegisterArray :: LoadFromMeshFile (istream & ist) const
  {
    for (int i = 0; i < Size(); i++)
      {
        NetgenGeometry * hgeom = (*this)[i]->LoadFromMeshFile (ist);
        if (hgeom)
          return shared_ptr<NetgenGeometry>(hgeom);
      }
    return nullptr;
  }



  
  int NetgenGeometry :: GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam)
  {
    multithread.percent = 0;

    if(mparam.perfstepsstart <= MESHCONST_ANALYSE)
      {
        if(!mesh)
          mesh = make_shared<Mesh>();
        mesh->geomtype = GetGeomType();
        Analyse(*mesh, mparam);
      }

    if(multithread.terminate || mparam.perfstepsend <= MESHCONST_ANALYSE)
      return 0;

    if(mparam.perfstepsstart <= MESHCONST_MESHEDGES)
      FindEdges(*mesh, mparam);

    if(multithread.terminate || mparam.perfstepsend <= MESHCONST_MESHEDGES)
      return 0;

    if (mparam.perfstepsstart <= MESHCONST_MESHSURFACE)
      {
        MeshSurface(*mesh, mparam);
        mesh->CalcSurfacesOfNode();
      }

    if (multithread.terminate || mparam.perfstepsend <= MESHCONST_MESHSURFACE)
      return 0;

    if (mparam.perfstepsstart <= MESHCONST_OPTSURFACE)
      OptimizeSurface(*mesh, mparam);

    if (multithread.terminate || mparam.perfstepsend <= MESHCONST_OPTSURFACE)
      return 0;


    if(mparam.perfstepsstart <= MESHCONST_MESHVOLUME)
      {
        multithread.task = "Volume meshing";

        MESHING3_RESULT res = MeshVolume (mparam, *mesh);

        if (res != MESHING3_OK) return 1;
        if (multithread.terminate) return 0;

        RemoveIllegalElements (*mesh);
        if (multithread.terminate) return 0;

        MeshQuality3d (*mesh);
      }

    if (multithread.terminate || mparam.perfstepsend <= MESHCONST_MESHVOLUME)
      return 0;


    if (mparam.perfstepsstart <= MESHCONST_OPTVOLUME)
      {
	multithread.task = "Volume optimization";

	OptimizeVolume (mparam, *mesh);
	if (multithread.terminate) return 0;
      }
    FinalizeMesh(*mesh);
    return 0;
  }

  void NetgenGeometry :: Save (string filename) const
  {
    throw NgException("Cannot save geometry - no geometry available");
  }

  static RegisterClassForArchive<NetgenGeometry> regnggeo;
}

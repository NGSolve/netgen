#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

  DLL_HEADER void Optimize2d (Mesh & mesh, MeshingParameters & mp)
  {
    static Timer timer("optimize2d"); RegionTimer reg(timer);

    mesh.CalcSurfacesOfNode();

    bool secondorder = mesh.GetNP() > mesh.GetNV();


    if (secondorder)
      {
      for (SurfaceElementIndex ei = 0; ei < mesh.GetNSE(); ei++)
        mesh[ei].SetType(TRIG);
      }
    mesh.Compress();

    const char * optstr = mp.optimize2d.c_str();
    int optsteps = mp.optsteps2d;

    for (int i = 1; i <= optsteps; i++)
      for (size_t j = 1; j <= strlen(optstr); j++)
	{
	  if (multithread.terminate) break;
	  switch (optstr[j-1])
	    {
	    case 's': 
	      {  // topological swap
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.EdgeSwapping (mesh, 0);
		break;
	      }
	    case 'S': 
	      {  // metric swap
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.EdgeSwapping (mesh, 1);
		break;
	      }
	    case 'm': 
	      {
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.ImproveMesh(mesh, mp);
		break;
	      }
	    case 'c': 
	      {
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.CombineImprove(mesh);
		break;
	      }
	    default:
	      cerr << "Optimization code " << optstr[j-1] << " not defined" << endl;
	    }  
	}
    if (secondorder)
      {
        mesh.GetGeometry()->GetRefinement().MakeSecondOrder(mesh);
      }
  }

}

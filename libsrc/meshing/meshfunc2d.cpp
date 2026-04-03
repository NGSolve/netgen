#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

  DLL_HEADER void Optimize2d (Mesh & mesh, MeshingParameters & mp, int faceindex)
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

    bool optimize_swap_separate_faces = false;
    if(!mp.quad)
    {
      bool mixed = false;
      ParallelFor( Range(mesh.GetNSE()), [&] (auto i) NETGEN_LAMBDA_INLINE
          {
            if (mesh[SurfaceElementIndex(i)].GetNP() != 3)
                mixed = true;
          });
      if(mixed)
        optimize_swap_separate_faces = true;
    }

    if(faceindex)
      optimize_swap_separate_faces = false;

    const char * optstr = mp.optimize2d.c_str();
    int optsteps = mp.optsteps2d;

    // reset topology
    mesh.GetTopology() = MeshTopology(mesh);
    for (int i = 1; i <= optsteps; i++)
      for (size_t j = 1; j <= strlen(optstr); j++)
	{
	  if (multithread.terminate) break;
          MeshOptimize2d meshopt(mesh);
          meshopt.SetMetricWeight (mp.elsizeweight);
          meshopt.SetFaceIndex(faceindex);
	  switch (optstr[j-1])
	    {
	    case 's': 
	      {  // topological swap

                if(optimize_swap_separate_faces)
                {
                  for(auto i : Range(1, mesh.GetNFD()+1))
                  {
                    meshopt.SetFaceIndex(i);
                    meshopt.EdgeSwapping (0);
                  }
                }
                else
                {
                  meshopt.EdgeSwapping (0);
                }
		break;
	      }
	    case 'S': 
	      {  // metric swap
                if(optimize_swap_separate_faces)
                {
                  for(auto i : Range(1, mesh.GetNFD()+1))
                  {
                    meshopt.SetFaceIndex(i);
                    meshopt.EdgeSwapping (1);
                  }
                }
                else
                {
                  meshopt.EdgeSwapping (1);
                }
		break;
	      }
	    case 'm': 
	      {
		meshopt.ImproveMesh(mp);
		break;
	      }
	    case 'c': 
	      {
		meshopt.CombineImprove();
		break;
	      }
	    default:
	      cerr << "Optimization code " << optstr[j-1] << " not defined" << endl;
	    }  
	}
    mesh.Compress(); // better: compress in individual steps, if necessary
    if (secondorder)
      {
        mesh.GetGeometry()->GetRefinement().MakeSecondOrder(mesh);
      }
  }

}

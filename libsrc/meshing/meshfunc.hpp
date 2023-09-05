#ifndef FILE_MESHFUNC
#define FILE_MESHFUNC

/**************************************************************************/
/* File:   meshfunc.hpp                                                   */
/* Author: Johannes Gerstmayr, Joachim Schoeberl                          */
/* Date:   26. Jan. 98                                                    */
/**************************************************************************/

#include <mydefs.hpp>
#include "meshing3.hpp"
#include "meshtype.hpp"

namespace netgen
{
/*
  Functions for mesh-generations strategies
 */

class Mesh;
// class CSGeometry;

/// Build tet-mesh
DLL_HEADER MESHING3_RESULT MeshVolume (const MeshingParameters & mp, Mesh& mesh3d);

/// Build mixed-element mesh
// MESHING3_RESULT MeshMixedVolume (MeshingParameters & mp, Mesh& mesh3d);

/// Optimize tet-mesh
DLL_HEADER MESHING3_RESULT OptimizeVolume (const MeshingParameters & mp, Mesh& mesh3d);
//			       const CSGeometry * geometry = NULL);

DLL_HEADER void RemoveIllegalElements (Mesh & mesh3d);


enum MESHING_STEP { 
  MESHCONST_ANALYSE = 1,
  MESHCONST_MESHEDGES = 2,
  MESHCONST_MESHSURFACE = 3,
  MESHCONST_OPTSURFACE = 4,
  MESHCONST_MESHVOLUME = 5,
  MESHCONST_OPTVOLUME = 6
};
} // namespace netgen

#endif

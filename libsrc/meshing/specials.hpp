#ifndef FILE_SPECIALS
#define FILE_SPECIALS

/*

  Very special implementations ..
  
 */

namespace netgen {
///
DLL_HEADER extern void CutOffAndCombine (Mesh & mesh, const Mesh & othermesh);

DLL_HEADER extern void HelmholtzMesh (Mesh & mesh);

} // namespace netgen
#endif

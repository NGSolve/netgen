#ifndef NETGEN_CLUSTERS_HPP
#define NETGEN_CLUSTERS_HPP

/**************************************************************************/
/* File:   clusers.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   28. Apr. 01                                                    */
/**************************************************************************/

/*
  Anisotropic clusters

  nodes, edges, faces, elements
*/

#include "meshclass.hpp"

namespace netgen
{

class AnisotropicClusters
{
  const Mesh & mesh;

  int nv, ned, nfa, ne;

  // connected nodes, nodes = vertices, edges, faces, elements
  NgArray<int> cluster_reps;

public:
  AnisotropicClusters (const Mesh & amesh);
  ~AnisotropicClusters();

  void Update();

  int GetVertexRepresentant (int vnr) const
  { return cluster_reps[vnr-1]; }
  int GetEdgeRepresentant (int ednr) const
  { return cluster_reps[nv+ednr-1]; }
  int GetFaceRepresentant (int fnr) const
  { return cluster_reps[nv+ned+fnr-1]; }
  int GetElementRepresentant (int enr) const
  { return cluster_reps[nv+ned+nfa+enr-1]; }
};
} // namespace netgen
#endif // NETGEN_CLUSTERS_HPP

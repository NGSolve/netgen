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
  { return cluster_reps.Get(vnr); }
  int GetEdgeRepresentant (int ednr) const
  { return cluster_reps.Get(nv+ednr); }
  int GetFaceRepresentant (int fnr) const
  { return cluster_reps.Get(nv+ned+fnr); }
  int GetElementRepresentant (int enr) const
  { return cluster_reps.Get(nv+ned+nfa+enr); }
};
} // namespace netgen
#endif // NETGEN_CLUSTERS_HPP

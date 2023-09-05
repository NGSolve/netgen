#ifndef NETGEN_IMPROVE2_HPP
#define NETGEN_IMPROVE2_HPP

#include "meshtype.hpp"

namespace netgen
{

inline void AppendEdges( const Element2d & elem, PointIndex pi, Array<std::tuple<PointIndex,PointIndex>> & edges )
{
  for (int j = 0; j < 3; j++)
  {
      PointIndex pi0 = elem[j];
      PointIndex pi1 = elem[(j+1)%3];
      if (pi1 < pi0) Swap(pi0, pi1);
      if(pi0==pi)
          edges.Append(std::make_tuple(pi0, pi1));
  }
}

inline void AppendEdges( const Element & elem, PointIndex pi, Array<std::tuple<PointIndex,PointIndex>> & edges )
{
  static constexpr int tetedges[6][2] =
  { { 0, 1 }, { 0, 2 }, { 0, 3 },
      { 1, 2 }, { 1, 3 }, { 2, 3 } };

  if(elem.Flags().fixed)
      return;
  for (int j = 0; j < 6; j++)
  {
      PointIndex pi0 = elem[tetedges[j][0]];
      PointIndex pi1 = elem[tetedges[j][1]];
      if (pi1 < pi0) Swap(pi0, pi1);
      if(pi0==pi)
          edges.Append(std::make_tuple(pi0, pi1));
  }
}

template<typename TINDEX>
void BuildEdgeList( const Mesh & mesh, const Table<TINDEX, PointIndex> & elementsonnode, Array<std::tuple<PointIndex, PointIndex>> & edges )
{
  static_assert(is_same_v<TINDEX, ElementIndex>||is_same_v<TINDEX,SurfaceElementIndex>, "Invalid type for TINDEX");
  static Timer tbuild_edges("Build edges"); RegionTimer reg(tbuild_edges);

  int ntasks = 4*ngcore::TaskManager::GetMaxThreads();
  Array<Array<std::tuple<PointIndex,PointIndex>>> task_edges(ntasks);

  ParallelFor(IntRange(ntasks), [&] (int ti)
    {
      auto myrange = mesh.Points().Range().Split(ti, ntasks);
      ArrayMem<std::tuple<PointIndex,PointIndex>, 100> local_edges;
      for (auto pi : myrange)
      {
        local_edges.SetSize(0);

        for(auto ei : elementsonnode[pi])
        {
            const auto & elem = mesh[ei];
            if (elem.IsDeleted()) continue;

            AppendEdges(elem, pi, local_edges);
        }
        QuickSort(local_edges);

        auto edge_prev = std::make_tuple<PointIndex, PointIndex>(-1,-1);

        for(auto edge : local_edges)
            if(edge != edge_prev)
            {
                task_edges[ti].Append(edge);
                edge_prev = edge;
            }
      }
    }, ntasks);

  int num_edges = 0;
  for (auto & edg : task_edges)
      num_edges += edg.Size();
  edges.SetAllocSize(num_edges);
  for (auto & edg : task_edges)
      edges.Append(edg);
}


class Neighbour
{
  int nr[3];
  int orient[3];

public:
  Neighbour () { ; }

  void SetNr (int side, int anr) { nr[side] = anr; }
  int GetNr (int side) { return nr[side]; }

  void SetOrientation (int side, int aorient) { orient[side] = aorient; }
  int GetOrientation (int side) { return orient[side]; }
};

///
class MeshOptimize2d
{
  int faceindex = 0;
  int improveedges = 0;
  double metricweight = 0.;
  int writestatus = 1;
  Mesh& mesh;
  const NetgenGeometry& geo;
public:
  ///
  MeshOptimize2d(Mesh& amesh) : mesh(amesh), geo(*mesh.GetGeometry())
  {}
  virtual ~MeshOptimize2d() { ; }
  ///
  DLL_HEADER void ImproveMesh (const MeshingParameters & mp);
  DLL_HEADER void ImproveMeshJacobian (const MeshingParameters & mp);
  DLL_HEADER void ImproveVolumeMesh ();
  DLL_HEADER void ProjectBoundaryPoints(NgArray<int> & surfaceindex, 
			     const NgArray<Point<3>* > & from, NgArray<Point<3>* > & dest);

  DLL_HEADER bool EdgeSwapping (const int usemetric, Array<Neighbour> &neighbors, Array<bool> &swapped,
                                const SurfaceElementIndex t1, const int edge, const int t, Array<int,PointIndex> &pdef, const bool check_only=false);
  DLL_HEADER void EdgeSwapping (int usemetric);
  DLL_HEADER void CombineImprove ();
  DLL_HEADER void SplitImprove ();

  DLL_HEADER void GenericImprove ();


  void SetFaceIndex (int fi) { faceindex = fi; }
  void SetImproveEdges (int ie) { improveedges = ie; }
  void SetMetricWeight (double mw) { metricweight = mw; }
  void SetWriteStatus (int ws) { writestatus = ws; }


  /// liefert zu einem 3d-Punkt die geominfo (Dreieck) und liefert 1, wenn erfolgreich, 
  /// 0, wenn nicht (Punkt ausserhalb von chart)
  ///

  void CheckMeshApproximation (Mesh & mesh);


  ///
  friend class Opti2SurfaceMinFunction;
  ///
  friend class Opti2EdgeMinFunction;
  ///
  friend double Opti2FunctionValueGrad (const Vector & x, Vector & grad);
  ///
  friend double Opti2EdgeFunctionValueGrad (const Vector & x, Vector & grad);



};


extern void CalcTriangleBadness (double x2, double x3, double y3, 
				 double metricweight,
				 double h, double & badness, 
				 double & g1x, double & g1y);




extern double CalcTriangleBadness (const Point<3> & p1, 
				   const Point<3> & p2, 
				   const Point<3> & p3,
				   double metricweight,
				   double h);

extern double CalcTriangleBadness (const Point<3> & p1, 
				   const Point<3> & p2, 
				   const Point<3> & p3,
				   const Vec<3> & n,
				   double metricweight,
				   double h);
} // namespace netgen
#endif // NETGEN_IMPROVE2_HPP

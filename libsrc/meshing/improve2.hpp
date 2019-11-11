#ifndef FILE_IMPROVE2
#define FILE_IMPROVE2

template<typename TINDEX>
void BuildEdgeList( const Mesh & mesh, const Table<TINDEX, PointIndex> & elementsonnode, Array<std::tuple<PointIndex, PointIndex>> & edges )
{
  static Timer tbuild_edges("Build edges"); RegionTimer reg(tbuild_edges);

  static constexpr int tetedges[6][2] =
    { { 0, 1 }, { 0, 2 }, { 0, 3 },
        { 1, 2 }, { 1, 3 }, { 2, 3 } };

  int ntasks = 2*ngcore::TaskManager::GetMaxThreads();
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

            for (int j = 0; j < 6; j++)
            {
                PointIndex pi0 = elem[tetedges[j][0]];
                PointIndex pi1 = elem[tetedges[j][1]];
                if (pi1 < pi0) Swap(pi0, pi1);
                if(pi0==pi)
                    local_edges.Append(std::make_tuple(pi0, pi1));
            }
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
  void ImproveMesh (const MeshingParameters & mp);
  void ImproveMeshJacobian (const MeshingParameters & mp);
  void ImproveVolumeMesh ();
  void ProjectBoundaryPoints(NgArray<int> & surfaceindex, 
			     const NgArray<Point<3>* > & from, NgArray<Point<3>* > & dest);

  bool EdgeSwapping (const int usemetric, Array<Neighbour> &neighbors, Array<bool> &swapped,
    const SurfaceElementIndex t1, const int edge, const int t, Array<int,PointIndex> &pdef, const bool check_only=false);
  void EdgeSwapping (int usemetric);
  void CombineImprove ();
  void SplitImprove ();

  void GenericImprove ();


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

#endif



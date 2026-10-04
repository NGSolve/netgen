#ifndef FILE_PARALLELTOP
#define FILE_PARALLELTOP

namespace netgen
{

  /*
    Which other ranks share a vertex, edge or face of the local mesh, and the
    global vertex numbers. All update functions are collective over the mesh
    communicator.
  */
  class ParallelMeshTopology
  {
    const Mesh & mesh;

    /// row per local vertex/edge/face: the other ranks sharing it, sorted
    DynamicTable<int> loc2distvert;
    DynamicTable<int> loc2distedge, loc2distface;

    /// global vertex numbers (1-based)
    Array<int> glob_vert;

  public:
    ParallelMeshTopology (const Mesh & amesh);

    /// vertices created by refinement (mesh.mlbetweennodes) are shared where both parents are
    void IdentifyNewVertices ();
    /// numbers new vertices globally and reorders the points by global number
    void EnumeratePointsGlobally ();
    /// shared edges and faces from the shared vertices, needs the mesh topology
    void UpdateEdgesAndFaces ();

    void AddDistantProc     (PointIndex pi, int proc) { loc2distvert.AddUnique (pi-IndexBASE<PointIndex>(), proc); }
    void AddDistantEdgeProc (int edge, int proc) { loc2distedge.AddUnique (edge, proc); }
    void AddDistantFaceProc (int face, int proc) { loc2distface.AddUnique (face, proc); }

    FlatArray<int> GetDistantProcs     (PointIndex pi) const { return loc2distvert[pi-IndexBASE<PointIndex>()]; }
    FlatArray<int> GetDistantEdgeProcs (EdgeIndex locnum) const { return loc2distedge[locnum.Nr0()]; }
    FlatArray<int> GetDistantFaceProcs (FaceIndex locnum) const { return loc2distface[locnum.Nr0()]; }

    auto & L2G (PointIndex pi) { return glob_vert[pi-IndexBASE<PointIndex>()]; }
    auto L2G (PointIndex pi) const { return glob_vert[pi-IndexBASE<PointIndex>()]; }
    int GetGlobalPNum (PointIndex locnum) const { return glob_vert[locnum-IndexBASE<PointIndex>()]; }

    /// number of local vertices: resizes the distant-proc table, keeps existing rows
    void SetNV (int anv);
    /// number of local vertices: resets the global numbers
    void SetNV_Loc2Glob (int anv);
  };

}

#endif

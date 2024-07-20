#ifndef FILE_PARALLELTOP
#define FILE_PARALLELTOP

namespace netgen
{


  class ParallelMeshTopology
  {
    const Mesh & mesh;

    /**
       mapping from local to distant vertex number
       each row of the table corresponds to one vertex
       each row contains a list of pairs (procnr, dist_vnum)
    */
    
    DynamicTable<int> loc2distvert;
    DynamicTable<int> loc2distedge, loc2distface;

    Array<int> glob_vert;

    // will get rid of them
    NgArray<int> glob_edge, glob_face;
    NgArray<int> glob_el, glob_surfel, glob_segm;

    bool is_updated;

  public:

    ParallelMeshTopology (const Mesh & amesh);
    ~ParallelMeshTopology ();

    void Reset ();
    void Print() const;

    
    void UpdateCoarseGrid();
    // [[deprecated("should not need it anymore")]]                    
    // void UpdateCoarseGridGlobal();
    void IdentifyVerticesAfterRefinement();
    void EnumeratePointsGlobally ();
        
    void AddDistantProc    (PointIndex pi, int proc) { loc2distvert.AddUnique (pi-PointIndex::BASE, proc); }
    void AddDistantFaceProc (int edge, int proc) { loc2distface.AddUnique (edge, proc); }
    void AddDistantEdgeProc (int face, int proc) { loc2distedge.AddUnique (face, proc); }
    
    FlatArray<int> GetDistantProcs (PointIndex pi) const { return loc2distvert[pi-PointIndex::BASE]; }
    FlatArray<int> GetDistantFaceProcs (int locnum) const { return loc2distface[locnum]; }
    FlatArray<int> GetDistantEdgeProcs (int locnum) const { return loc2distedge[locnum]; }


    
    auto & L2G (PointIndex pi) { return glob_vert[pi-PointIndex::BASE]; } 
    auto L2G (PointIndex pi) const { return glob_vert[pi-PointIndex::BASE]; } 


    /// set number of local vertices, reset sizes of loc2dist_vert, isexchangevert...
    void SetNV (int anv);
    void SetNV_Loc2Glob (int anv);
    void SetNE (int ane);
    void SetNSE (int anse);
    void SetNSegm (int anseg);

    [[deprecated("Use AddDistantFaceProc instead!")]]                
    void SetDistantFaceNum (int dest, int locnum) { loc2distface.AddUnique (locnum-1, dest); }
    [[deprecated("Use AddDistantProc instead!")]]                
    void SetDistantPNum    (int dest, int locnum) { loc2distvert.AddUnique (locnum-1, dest); }
    [[deprecated("Use AddDistantEdgeProc instead!")]]                
    void SetDistantEdgeNum (int dest, int locnum) { loc2distedge.AddUnique (locnum-1, dest); }

    [[deprecated("Use GetDistantFaceProcx instead!")]]                    
    FlatArray<int> GetDistantFaceNums (int locnum) const { return loc2distface[locnum]; }
    [[deprecated("Use GetDistantEdgeProcx instead!")]]
    FlatArray<int> GetDistantEdgeNums (int locnum) const { return loc2distedge[locnum]; }


    
    [[deprecated("Use L2G(pi) instead!")]]                
    void SetLoc2Glob_Vert   (int locnum, int globnum) { glob_vert[locnum-1] = globnum; }
    // [[deprecated("Try to avoid global enumration!")]]                
    void SetLoc2Glob_Edge   (int locnum, int globnum) { glob_edge[locnum-1] = globnum; }
    // [[deprecated("Try to avoid global enumration!")]]                
    void SetLoc2Glob_Face   (int locnum, int globnum) { glob_face[locnum-1] = globnum; }
    // [[deprecated("Try to avoid global enumration!")]]                
    void SetLoc2Glob_VolEl  (int locnum, int globnum) { glob_el[locnum-1] = globnum; }
    // [[deprecated("Try to avoid global enumration!")]]                
    void SetLoc2Glob_SurfEl (int locnum, int globnum) { glob_surfel[locnum-1] = globnum; }
    // [[deprecated("Try to avoid global enumration!")]]                
    void SetLoc2Glob_Segm   (int locnum, int globnum) { glob_segm[locnum-1] = globnum; }

    // [[deprecated("Try to avoid global enumration!")]]                    
    int GetGlobalPNum    (PointIndex locnum) const { return glob_vert[locnum-PointIndex::BASE]; }
    [[deprecated("Try to avoid global enumration!")]]                
    int GetGlobalEdgeNum (int locnum) const { return glob_edge[locnum-1]; }
    [[deprecated("Try to avoid global enumration!")]]                
    int GetGlobalFaceNum (int locnum) const { return glob_face[locnum-1]; }
    [[deprecated("Try to avoid global enumration!")]]                
    int GetGlobalElNum   (int locnum) const { return glob_el[locnum-1]; }
    [[deprecated("Try to avoid global enumration!")]]                
    int GetGlobalSElNum  (int locnum) const { return glob_surfel[locnum-1]; }

    

    [[deprecated("Use GetDistantPNums(locnum).Size() instead!")]]            
    int GetNDistantPNums (int locpnum) const { return loc2distvert[locpnum-1].Size(); }

    [[deprecated("Use GetDistantFaceNums(locnum).Size() instead!")]]                
    int GetNDistantFaceNums (int locfacenum) const { return loc2distface[locfacenum-1].Size(); }

    [[deprecated("Use GetDistantEdgeNums(locnum).Size() instead!")]]                    
    int GetNDistantEdgeNums ( int locedgenum) const { return loc2distedge[locedgenum-1].Size(); }

    [[deprecated("Use GetDistantPNums(locnum) -> FlatArray instead!")]]                
    void GetDistantPNums (int locpnum, int * distpnums ) const
    {
      for (int i = 0; i < loc2distvert[locpnum-1].Size(); i++ )
	distpnums[i] = loc2distvert[locpnum-1][i];
    } 

    [[deprecated("Use GetDistantFaceNums(locnum) -> FlatArray instead!")]]                    
    void GetDistantFaceNums (int locfacenum, int * distfacenums ) const
    {
      for ( int i = 0; i < loc2distface[locfacenum-1].Size(); i++ )
	distfacenums[i] = loc2distface[locfacenum-1][i];
    } 

    [[deprecated("Use GetDistantFaceNums(locnum) -> FlatArray instead!")]]                        
    void GetDistantFaceNums (int locfacenum, NgArray<int> & distfacenums ) const
    {
      // distfacenums = loc2distface[locfacenum-1];
      auto loc = loc2distface[locfacenum-1];
      distfacenums.SetSize (loc.Size());
      for (int i = 0; i < loc.Size(); i++)
        distfacenums[i] = loc[i];
    }

    [[deprecated("Use GetDistantEdgeNums(locnum) -> FlatArray instead!")]]                            
    void GetDistantEdgeNums (int locedgenum, int * distedgenums ) const
    {
      for (int i = 0; i < loc2distedge[locedgenum-1].Size(); i++ )
	distedgenums[i] = loc2distedge[locedgenum-1][i];
    } 

    [[deprecated("Use GetDistantEdgeNums(locnum) -> FlatArray instead!")]]                                
    void GetDistantEdgeNums (int locedgenum, NgArray<int> & distedgenums ) const
    {
      // distedgenums = loc2distedge[locedgenum-1];
      auto loc = loc2distedge[locedgenum-1];
      distedgenums.SetSize (loc.Size());
      for (int i = 0; i < loc.Size(); i++)
        distedgenums[i] = loc[i];
    } 

    [[deprecated("Use GetDistantProcs(..)!")]]                    
    FlatArray<int> GetDistantPNums (int locnum) const { return loc2distvert[locnum]; }


    
    [[deprecated("Use GetDistantProcs(..).Contains instead!")]]                
    bool IsExchangeVert (int dest, int vnum) const
    {
      return loc2distvert[vnum-1].Contains (dest);
    }
  };

}




#endif

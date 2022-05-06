// #ifdef PARALLEL


#include <meshing.hpp>
#include "paralleltop.hpp"


namespace netgen
{

  ParallelMeshTopology :: ParallelMeshTopology (const Mesh & amesh)
    : mesh(amesh)
  {
    is_updated = false;
  }
 

  ParallelMeshTopology :: ~ParallelMeshTopology ()
  {
    ;
  }


  void ParallelMeshTopology :: Reset ()
  {
    *testout << "ParallelMeshTopology::Reset" << endl;

    if ( mesh.GetCommunicator().Size() == 1 ) return;

    size_t ned = mesh.GetTopology().GetNEdges();
    size_t nfa = mesh.GetTopology().GetNFaces();

    if (glob_edge.Size() != ned)
      {
	glob_edge.SetSize(ned);
	glob_face.SetSize(nfa);
	glob_edge = -1;
	glob_face = -1;

	loc2distedge.ChangeSize (ned);
	loc2distface.ChangeSize (nfa);
      }

    if (glob_vert.Size() != mesh.GetNV())
      {
	SetNV(mesh.GetNV());
	SetNE(mesh.GetNE());
      }
  }


  void ParallelMeshTopology ::  Print() const
  {
    ;
  }


  void ParallelMeshTopology :: EnumeratePointsGlobally ()
  {
    auto comm = mesh.GetCommunicator();
    auto rank = comm.Rank();

    size_t oldnv = glob_vert.Size();
    size_t nv = loc2distvert.Size();
    *testout << "enumerate globally, loc2distvert.size = " << loc2distvert.Size()
             << ", glob_vert.size = " << glob_vert.Size() << endl;

    
    if (rank == 0)
      nv = 0;

    // IntRange newvr(oldnv, nv); // new vertex range
    auto new_pir = Range(PointIndex(oldnv+PointIndex::BASE),
                         PointIndex(nv+PointIndex::BASE));
    
    glob_vert.SetSize (nv);
    for (auto pi : new_pir)
      L2G(pi) = -1;

    int num_master_points = 0;

    for (auto pi : new_pir)
      {
        auto dps = GetDistantProcs(pi);
        // check sorted:
        for (int j = 0; j+1 < dps.Size(); j++)
          if (dps[j+1] < dps[j]) cout << "wrong sort" << endl;
        
        if (dps.Size() == 0 || dps[0] > comm.Rank())
          L2G(pi) = num_master_points++;
      }
    
    // *testout << "nummaster = " << num_master_points << endl;

    Array<int> first_master_point(comm.Size());
    comm.AllGather (num_master_points, first_master_point);
    auto max_oldv = comm.AllReduce (Max (glob_vert.Range(0, oldnv)), MPI_MAX);
    if (comm.AllReduce (oldnv, MPI_SUM) == 0)
      max_oldv = PointIndex::BASE-1;
    
    size_t num_glob_points = max_oldv+1;
    for (int i = 0; i < comm.Size(); i++)
      {
        int cur = first_master_point[i];
        first_master_point[i] = num_glob_points;
        num_glob_points += cur;
      }
    
    for (auto pi : new_pir)
      if (L2G(pi) != -1)
        L2G(pi) += first_master_point[comm.Rank()];
    
    // ScatterDofData (global_nums); 
    
    Array<int> nsend(comm.Size()), nrecv(comm.Size());
    nsend = 0;
    nrecv = 0;

    /** Count send/recv size **/
    for (auto pi : new_pir)
      if (auto dps = GetDistantProcs(pi); dps.Size())
        {
          if (rank < dps[0])
            for (auto p : dps)
              nsend[p]++;
          else
            nrecv[dps[0]]++;
        }
    
    Table<PointIndex> send_data(nsend);
    Table<PointIndex> recv_data(nrecv);
    
    /** Fill send_data **/
    nsend = 0;
    for (auto pi : new_pir)
      if (auto dps = GetDistantProcs(pi); dps.Size())
        if (rank < dps[0])
          for (auto p : dps)
            send_data[p][nsend[p]++] = L2G(pi);

    Array<MPI_Request> requests;
    for (int i = 0; i < comm.Size(); i++)
      {
        if (nsend[i])
          requests.Append (comm.ISend (send_data[i], i, 200));
        if (nrecv[i])
          requests.Append (comm.IRecv (recv_data[i], i, 200));
      }
    
    MyMPI_WaitAll (requests);
    
    Array<int> cnt(comm.Size());
    cnt = 0;

    for (auto pi : new_pir)
      if (auto dps = GetDistantProcs(pi); dps.Size())
        if (int master = dps[0]; master < comm.Rank())
          L2G(pi) = recv_data[master][cnt[master]++];
    
    // reorder following global ordering:
    Array<int> index0(glob_vert.Size());
    for (int pi : Range(index0))
      index0[pi] = pi;
    QuickSortI (glob_vert, index0);

    if (rank != 0)
      {
        Array<PointIndex, PointIndex> inv_index(index0.Size());
        for (int i = 0; i < index0.Size(); i++)
          inv_index[index0[i]+PointIndex::BASE] = i+PointIndex::BASE;
        
        for (auto & el : mesh.VolumeElements())
          for (PointIndex & pi : el.PNums())
            pi = inv_index[pi];
        for (auto & el : mesh.SurfaceElements())
          for (PointIndex & pi : el.PNums())
            pi = inv_index[pi];
        for (auto & el : mesh.LineSegments())
          for (PointIndex & pi : el.PNums())
            pi = inv_index[pi];
        
        // auto hpoints (mesh.Points());
        Array<MeshPoint, PointIndex> hpoints { mesh.Points() };
        for (PointIndex pi : Range(mesh.Points()))
          mesh.Points()[inv_index[pi]] = hpoints[pi];

        if (mesh.mlbetweennodes.Size() == mesh.Points().Size())
          {
            NgArray<PointIndices<2>,PointIndex::BASE> hml { mesh.mlbetweennodes };
            for (PointIndex pi : Range(mesh.Points()))
              mesh.mlbetweennodes[inv_index[pi]] = hml[pi];
          }

        // *testout << "index0 = " << endl << index0 << endl;
        // *testout << "loc2distvertold = " << endl;
        // for (auto i : Range(index0))
        // *testout << "l " << i << " globi "<< glob_vert[i]  << " dist = " << loc2distvert[i] << endl;

        DynamicTable<int> oldtable = std::move(loc2distvert);        
        loc2distvert = DynamicTable<int> (oldtable.Size());
        for (size_t i = 0; i < oldtable.Size(); i++)
          for (auto val : oldtable[index0[i]])
            loc2distvert.Add (i, val);

        Array<int> hglob_vert(glob_vert);
        for (int i = 0; i < index0.Size(); i++)
          glob_vert[i] = hglob_vert[index0[i]];

        // *testout << "loc2distvertnew = " << endl;
        // for (auto i : Range(index0))
        // *testout << "l " << i << " globi "<< glob_vert[i]  << " dist = " << loc2distvert[i] << endl;
      }

    /*
    for (size_t i = 0; i+1 < glob_vert.Size(); i++)
      if (glob_vert[i] > glob_vert[i+1])
        cout << "wrong ordering of globvert" << endl;
    */
    if (glob_vert.Size() > 1)
      for (auto i : Range(glob_vert).Modify(0,-1))
        if (glob_vert[i] > glob_vert[i+1])
          cout << "wrong ordering of globvert" << endl;
  }
  
  /*
  void ParallelMeshTopology :: SetDistantFaceNum (int dest, int locnum)
  {
    for ( int i = 0; i < loc2distface[locnum-1].Size(); i+=1 )
      if ( loc2distface[locnum-1][i] == dest )
	return;
    loc2distface.Add(locnum-1, dest);
  }

  void ParallelMeshTopology :: SetDistantPNum (int dest, int locnum)
  {
    for ( int i = 0;  i < loc2distvert[locnum-1].Size(); i+=1 )
      if ( loc2distvert[locnum-1][i] == dest )
	return;
    loc2distvert.Add (locnum-1, dest);  
  }


  void ParallelMeshTopology :: SetDistantEdgeNum (int dest, int locnum)
  {
    for ( int i = 0; i < loc2distedge[locnum-1].Size(); i+=1 )
      if ( loc2distedge[locnum-1][i] == dest )
	return;
    loc2distedge.Add (locnum-1, dest);
  }
  */
  
  void ParallelMeshTopology :: SetNV_Loc2Glob (int anv)
  {
    glob_vert.SetSize(anv);
    glob_vert = -1;
  }
  
  void ParallelMeshTopology :: SetNV (int anv)
  {
    // glob_vert.SetSize(anv);
    // glob_vert = -1;
    // loc2distvert.ChangeSize (anv);

    DynamicTable<int> oldtable(loc2distvert.Size());    
    for (size_t i = 0; i < loc2distvert.Size(); i++)
      for (auto val : loc2distvert[i])
        oldtable.Add (i, val);
    loc2distvert = DynamicTable<int> (anv);
    for (size_t i = 0; i < min(size_t(anv), oldtable.Size()); i++)
      for (auto val : oldtable[i])
        loc2distvert.Add (i, val);
  }

  void ParallelMeshTopology :: SetNE ( int ane )
  {
    glob_el.SetSize (ane);
    glob_el = -1;
  }

  void ParallelMeshTopology :: SetNSE ( int anse )
  {
    glob_surfel.SetSize(anse);
    glob_surfel = -1;
  }

  void ParallelMeshTopology :: SetNSegm ( int anseg )
  {
    glob_segm.SetSize (anseg);
    glob_segm = -1;
  }






  void ParallelMeshTopology :: UpdateCoarseGridGlobal ()
  {
    // cout << "updatecoarsegridglobal called" << endl;
    if (id == 0)
      PrintMessage ( 3, "UPDATE GLOBAL COARSEGRID STARTS" );      

    int timer = NgProfiler::CreateTimer ("UpdateCoarseGridGlobal");
    NgProfiler::RegionTimer reg(timer);

    *testout << "ParallelMeshTopology :: UpdateCoarseGridGlobal" << endl;

    const MeshTopology & topology = mesh.GetTopology();
    auto comm = mesh.GetCommunicator();

    if ( id == 0 )
      {
	NgArray<NgArray<int>*> sendarrays(ntasks);
	for (int dest = 1; dest < ntasks; dest++)
	  sendarrays[dest] = new NgArray<int>;

	NgArray<int> edges, faces;
	for (int el = 1; el <= mesh.GetNE(); el++)
	  {
	    topology.GetElementFaces (el, faces);
	    topology.GetElementEdges (el, edges);
	    const Element & volel = mesh.VolumeElement (el);

	    // NgArray<int> & sendarray = *sendarrays[volel.GetPartition()];
            NgArray<int> & sendarray = *sendarrays[mesh.vol_partition[el-1]];

	    for ( int i = 0; i < edges.Size(); i++ )
	      sendarray.Append (edges[i]);
	    for ( int i = 0; i < faces.Size(); i++ )
	      sendarray.Append (faces[i]);
	  }

	for (int el = 1; el <= mesh.GetNSE(); el++)
	  {
	    topology.GetSurfaceElementEdges (el, edges);
	    const Element2d & surfel = mesh.SurfaceElement (el);
	    // NgArray<int> & sendarray = *sendarrays[surfel.GetPartition()];
            NgArray<int> & sendarray = *sendarrays[mesh.surf_partition[el-1]];

	    for ( int i = 0; i < edges.Size(); i++ )
	      sendarray.Append (edges[i]);
	    sendarray.Append (topology.GetSurfaceElementFace (el));
	  }

	Array<MPI_Request> sendrequests;
	for (int dest = 1; dest < ntasks; dest++)
	  // sendrequests.Append (MyMPI_ISend (*sendarrays[dest], dest, MPI_TAG_MESH+10, comm));
          sendrequests.Append (comm.ISend (FlatArray<int>(*sendarrays[dest]), dest, MPI_TAG_MESH+10));
	MyMPI_WaitAll (sendrequests);

	for (int dest = 1; dest < ntasks; dest++)
	  delete sendarrays[dest];
      }

    else

      {
	// NgArray<int> recvarray;
	// MyMPI_Recv (recvarray, 0, MPI_TAG_MESH+10, comm);
	Array<int> recvarray;
	comm.Recv (recvarray, 0, MPI_TAG_MESH+10); // MyMPI_Recv (recvarray, 0, MPI_TAG_MESH+10, comm);

	int ii = 0;

	NgArray<int> faces, edges;

	for (int volel = 1; volel <= mesh.GetNE(); volel++)
	  {
	    topology.GetElementEdges ( volel, edges);
	    for ( int i = 0; i  < edges.Size(); i++)
	      SetLoc2Glob_Edge ( edges[i], recvarray[ii++]);

	    topology.GetElementFaces( volel, faces);
	    for ( int i = 0; i  < faces.Size(); i++)
	      SetLoc2Glob_Face ( faces[i], recvarray[ii++]);
	  }

	for (int surfel = 1; surfel <= mesh.GetNSE(); surfel++)
	  {
	    topology.GetSurfaceElementEdges (surfel, edges);
	    for (int i = 0; i  < edges.Size(); i++)
	      SetLoc2Glob_Edge (edges[i], recvarray[ii++]);
	    int face = topology.GetSurfaceElementFace (surfel);
	    SetLoc2Glob_Face ( face, recvarray[ii++]);
	  }
      }
    
    is_updated = true;
  }

  
  void ParallelMeshTopology :: IdentifyVerticesAfterRefinement()
  {
    static Timer t("ParallelTopology::UpdateCoarseGrid"); RegionTimer r(t);

    NgMPI_Comm comm = mesh.GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    if (ntasks == 1) return;
    
    Reset();
    static int timer = NgProfiler::CreateTimer ("UpdateCoarseGrid");
    NgProfiler::RegionTimer reg(timer);


    (*testout) << "UPDATE COARSE GRID PARALLEL TOPOLOGY " << endl;
    if (id == 0)
      PrintMessage (1, "update parallel topology");
    
    
    const MeshTopology & topology = mesh.GetTopology();

    Array<int> cnt_send(ntasks);

    int maxsize = comm.AllReduce (mesh.mlbetweennodes.Size(), MPI_MAX);
    // update new vertices after mesh-refinement
    if (maxsize > 0)
      {
        int newnv = mesh.mlbetweennodes.Size();
        
        loc2distvert.ChangeSize(mesh.mlbetweennodes.Size());

	bool changed = true;
	while (changed)
	  {
	    changed = false;

	    // build exchange vertices
	    cnt_send = 0;
	    for (PointIndex pi : mesh.Points().Range())
	      for (int dist : GetDistantProcs(pi))
		cnt_send[dist]++;
            // TABLE<int> dest2vert(cnt_send);    
	    DynamicTable<int> dest2vert(cnt_send);    
	    for (PointIndex pi : mesh.Points().Range())
	      for (int dist : GetDistantProcs(pi))
		dest2vert.Add (dist, pi);
            
	    for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
              if (auto [v1,v2] = mesh.mlbetweennodes[pi]; v1.IsValid())              
                {
                  auto procs1 = GetDistantProcs(v1);
                  auto procs2 = GetDistantProcs(v2);
                  for (int p : procs1)
                    if (procs2.Contains(p))
                      cnt_send[p]++;
                }

	    // TABLE<int> dest2pair(cnt_send);
            DynamicTable<int> dest2pair(cnt_send);            
            
            for (PointIndex pi : mesh.mlbetweennodes.Range())
              if (auto [v1,v2] = mesh.mlbetweennodes[pi]; v1.IsValid())
                {
                  auto procs1 = GetDistantProcs(v1);
                  auto procs2 = GetDistantProcs(v2);
                  for (int p : procs1)
                    if (procs2.Contains(p))
                      dest2pair.Add (p, pi);
                }

	    cnt_send = 0;
	    for (PointIndex pi : mesh.mlbetweennodes.Range())
              if (auto [v1,v2] = mesh.mlbetweennodes[pi]; v1.IsValid())
                {
                  auto procs1 = GetDistantProcs(v1);
                  auto procs2 = GetDistantProcs(v2);
                  
                  for (int p : procs1)
                    if (procs2.Contains(p))
                      cnt_send[p]+=2;
                }
	    
	    // TABLE<int> send_verts(cnt_send);
            DynamicTable<int> send_verts(cnt_send);

	    NgArray<int, PointIndex::BASE> loc2exchange(mesh.GetNV());

	    for (int dest = 0; dest < ntasks; dest++)
	      if (dest != id)
		{
		  loc2exchange = -1;
		  int cnt = 0;
		  for (PointIndex pi : dest2vert[dest])
		    loc2exchange[pi] = cnt++;
		  
		  for (PointIndex pi : dest2pair[dest])
                    if (auto [v1,v2] = mesh.mlbetweennodes[pi]; v1.IsValid())                    
                      {
                        auto procs1 = GetDistantProcs(v1);
                        auto procs2 = GetDistantProcs(v2);
                        
                        if (procs1.Contains(dest) && procs2.Contains(dest))
                          {
                            send_verts.Add (dest, loc2exchange[v1]);
                            send_verts.Add (dest, loc2exchange[v2]);
                          }
                      }
		}

	    DynamicTable<int> recv_verts(ntasks);
	    // MyMPI_ExchangeTable (send_verts, recv_verts, MPI_TAG_MESH+9, comm);
            comm.ExchangeTable (send_verts, recv_verts, MPI_TAG_MESH+9);

	    for (int dest = 0; dest < ntasks; dest++)
	      if (dest != id)
		{
		  loc2exchange = -1;
		  int cnt = 0;

		  for (PointIndex pi : dest2vert[dest])
		    loc2exchange[pi] = cnt++;
		  
		  FlatArray<int> recvarray = recv_verts[dest];
		  for (int ii = 0; ii < recvarray.Size(); ii+=2)
		    for (PointIndex pi : dest2pair[dest])
		      {
			PointIndex v1 = mesh.mlbetweennodes[pi][0];
			PointIndex v2 = mesh.mlbetweennodes[pi][1];
			if (v1.IsValid())
			  {
			    INDEX_2 re(recvarray[ii], recvarray[ii+1]);
			    INDEX_2 es(loc2exchange[v1], loc2exchange[v2]);
			    // if (es == re && !IsExchangeVert(dest, pi))
                            if (es == re && !GetDistantProcs(pi).Contains(dest))
			      {
				// SetDistantPNum(dest, pi);
                                AddDistantProc (pi, dest);
				changed = true;
			      }
			  }
		      }
		}

            changed = comm.AllReduce (changed, MPI_LOR);
	  }
      }

    NgArray<int> sendarray, recvarray;
    // cout << "UpdateCoarseGrid - edges" << endl;

    // static int timerv = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex vertices");
    static int timere = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex edges");
    static int timerf = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex faces");

    
    NgProfiler::StartTimer (timere);

    // build exchange vertices
    cnt_send = 0;
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantProcs(pi))
	cnt_send[dist]++;
    // TABLE<int> dest2vert(cnt_send);
    DynamicTable<int> dest2vert(cnt_send);    
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantProcs(pi))
	dest2vert.Add (dist, pi);
    
    // MPI_Group_free(&MPI_LocalGroup);
    // MPI_Comm_free(&MPI_LocalComm);
  }


  void ParallelMeshTopology :: UpdateCoarseGrid ()
  {
    static Timer t("ParallelTopology::UpdateCoarseGrid"); RegionTimer r(t);
    // cout << "UpdateCoarseGrid" << endl;
    // if (is_updated) return;

    NgMPI_Comm comm = mesh.GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    if (ntasks == 1) return;
    
    Reset();
    static int timer = NgProfiler::CreateTimer ("UpdateCoarseGrid");
    NgProfiler::RegionTimer reg(timer);


    (*testout) << "UPDATE COARSE GRID PARALLEL TOPOLOGY " << endl;
    if (id == 0)
      PrintMessage (1, "update parallel topology");


    // UpdateCoarseGridGlobal();


    /*
    // MPI_Barrier (MPI_COMM_WORLD);

    MPI_Group MPI_GROUP_comm;
    MPI_Group MPI_LocalGroup;
    MPI_Comm MPI_LocalComm;

    int process_ranks[] = { 0 };
    MPI_Comm_group (comm, &MPI_GROUP_comm);
    MPI_Group_excl (MPI_GROUP_comm, 1, process_ranks, &MPI_LocalGroup);
    MPI_Comm_create (comm, MPI_LocalGroup, &MPI_LocalComm);

    if (id == 0)
      {
        // SetNV(0);
        // EnumeratePointsGlobally();
        return;
      }
    */
    
    const MeshTopology & topology = mesh.GetTopology();

    Array<int> cnt_send(ntasks);

    // NgArray<int> sendarray, recvarray;
    // cout << "UpdateCoarseGrid - edges" << endl;

    // static int timerv = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex vertices");
    static int timere = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex edges");
    static int timerf = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex faces");


    NgProfiler::StartTimer (timere);


    int nfa = topology . GetNFaces();
    int ned = topology . GetNEdges();
    
    // build exchange vertices
    cnt_send = 0;
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantProcs(pi))
	cnt_send[dist]++;
    // TABLE<int> dest2vert(cnt_send);
    DynamicTable<int> dest2vert(cnt_send);    
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantProcs(pi))
	dest2vert.Add (dist, pi);

    // exchange edges
    cnt_send = 0;
    int v1, v2;
    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
        /*
	for (int dest = 1; dest < ntasks; dest++)
	  // if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
          if (GetDistantProcs(v1).Contains(dest) && GetDistantProcs(v2).Contains(dest))
	    cnt_send[dest-1]+=1;
        */
        for (auto p : GetDistantProcs(v1))
          if (GetDistantProcs(v2).Contains(p))
	    cnt_send[p]+=1;
      }
    
    // TABLE<int> dest2edge(cnt_send);
    DynamicTable<int> dest2edge(cnt_send);
    for (int & v : cnt_send) v *= 2;
    // TABLE<int> send_edges(cnt_send);
    DynamicTable<int> send_edges(cnt_send);

    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
	for (int dest = 0; dest < ntasks; dest++)
	  // if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
          if (GetDistantProcs(v1).Contains(dest) && GetDistantProcs(v2).Contains(dest))
	    dest2edge.Add (dest, edge);
      }


    NgArray<int, PointIndex::BASE> loc2exchange(mesh.GetNV());
    for (int dest = 0; dest < ntasks; dest++)
      {
        loc2exchange = -1;
        int cnt = 0;
        for (PointIndex pi : dest2vert[dest])
	  loc2exchange[pi] = cnt++;

	for (int edge : dest2edge[dest])
          {
            topology.GetEdgeVertices (edge, v1, v2);
            // if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
            if (GetDistantProcs(v1).Contains(dest) && GetDistantProcs(v2).Contains(dest))            
              {
                send_edges.Add (dest, loc2exchange[v1]);
                send_edges.Add (dest, loc2exchange[v2]);
              }
          }
      }

    // cout << "UpdateCoarseGrid - edges mpi-exchange" << endl;
    // TABLE<int> recv_edges(ntasks);
    DynamicTable<int> recv_edges(ntasks);
    // MyMPI_ExchangeTable (send_edges, recv_edges, MPI_TAG_MESH+9, comm);
    comm.ExchangeTable (send_edges, recv_edges, MPI_TAG_MESH+9);
    // cout << "UpdateCoarseGrid - edges mpi-exchange done" << endl;

    for (int dest = 0; dest < ntasks; dest++)
      {
	auto ex2loc = dest2vert[dest];
	if (ex2loc.Size() == 0) continue;

	INDEX_2_CLOSED_HASHTABLE<int> vert2edge(2*dest2edge[dest].Size()+10); 
	for (int edge : dest2edge[dest])
	  {
	    topology.GetEdgeVertices (edge, v1, v2);
	    vert2edge.Set(INDEX_2(v1,v2), edge);
	  }

	FlatArray<int> recvarray = recv_edges[dest];
        for (int ii = 0; ii < recvarray.Size(); ii+=2)
	  {
	    INDEX_2 re(ex2loc[recvarray[ii]], 
		       ex2loc[recvarray[ii+1]]);
	    if (vert2edge.Used(re))
	      // SetDistantEdgeNum(dest, vert2edge.Get(re));
              AddDistantEdgeProc(vert2edge.Get(re)-1, dest);
	  }
      }



    NgProfiler::StopTimer (timere);

    // cout << "UpdateCoarseGrid - faces" << endl;
    if (mesh.GetDimension() == 3)
      {
	NgProfiler::StartTimer (timerf);
	NgArray<int> verts;

	// exchange faces
	cnt_send = 0;
	for (int face = 1; face <= nfa; face++)
	  {
	    topology.GetFaceVertices (face, verts);
	    for (int dest = 0; dest < ntasks; dest++)
	      if (dest != id)
                /*
		if (IsExchangeVert (dest, verts[0]) && 
		    IsExchangeVert (dest, verts[1]) &&
		    IsExchangeVert (dest, verts[2]))
                */
                if (GetDistantProcs (verts[0]).Contains(dest) &&
                    GetDistantProcs (verts[1]).Contains(dest) &&
                    GetDistantProcs (verts[2]).Contains(dest))
		  cnt_send[dest]++;
	  }
	
	// TABLE<int> dest2face(cnt_send);
        DynamicTable<int> dest2face(cnt_send);
	for (int face = 1; face <= nfa; face++)
	  {
	    topology.GetFaceVertices (face, verts);
	    for (int dest = 0; dest < ntasks; dest++)
	      if (dest != id)
                /*
		if (IsExchangeVert (dest, verts[0]) && 
		    IsExchangeVert (dest, verts[1]) &&
		    IsExchangeVert (dest, verts[2]))
                */
                if (GetDistantProcs (verts[0]).Contains(dest) && 
                    GetDistantProcs (verts[1]).Contains(dest) &&
                    GetDistantProcs (verts[2]).Contains(dest))
		  dest2face.Add(dest, face);
	  }

	for (int & c : cnt_send) c*=3;
	// TABLE<int> send_faces(cnt_send);
        DynamicTable<int> send_faces(cnt_send);
	NgArray<int, PointIndex::BASE> loc2exchange(mesh.GetNV());
	for (int dest = 0; dest < ntasks; dest++)
	  if (dest != id)
	    {
	      if (dest2vert[dest].Size() == 0) continue;

	      loc2exchange = -1;
	      int cnt = 0;
	      for (PointIndex pi : dest2vert[dest])
		loc2exchange[pi] = cnt++;
	      
	      for (int face : dest2face[dest])
		{
		  topology.GetFaceVertices (face, verts);
                  /*
		  if (IsExchangeVert (dest, verts[0]) && 
		      IsExchangeVert (dest, verts[1]) &&
		      IsExchangeVert (dest, verts[2]))
                  */
                  if (GetDistantProcs (verts[0]).Contains(dest) &&
                      GetDistantProcs (verts[1]).Contains(dest) &&
                      GetDistantProcs (verts[2]).Contains(dest))
		    {
		      send_faces.Add (dest, loc2exchange[verts[0]]);
		      send_faces.Add (dest, loc2exchange[verts[1]]);
		      send_faces.Add (dest, loc2exchange[verts[2]]);
		    }
		}
	    }
	
	// cout << "UpdateCoarseGrid - faces mpi-exchange" << endl;
	// TABLE<int> recv_faces(ntasks);
	DynamicTable<int> recv_faces(ntasks);
	// MyMPI_ExchangeTable (send_faces, recv_faces, MPI_TAG_MESH+9, comm);
        comm.ExchangeTable (send_faces, recv_faces, MPI_TAG_MESH+9);
	// cout << "UpdateCoarseGrid - faces mpi-exchange done" << endl;
	
	for (int dest = 0; dest < ntasks; dest++)
	  {
	    auto ex2loc = dest2vert[dest];
	    if (ex2loc.Size() == 0) continue;
	    
	    INDEX_3_CLOSED_HASHTABLE<int> vert2face(2*dest2face[dest].Size()+10); 
	    for (int face : dest2face[dest])
	      {
		topology.GetFaceVertices (face, verts);
		vert2face.Set(INDEX_3(verts[0], verts[1], verts[2]), face);
	      }
	    
	    FlatArray<int> recvarray = recv_faces[dest];
	    for (int ii = 0; ii < recvarray.Size(); ii+=3)
	      {
		INDEX_3 re(ex2loc[recvarray[ii]], 
			   ex2loc[recvarray[ii+1]],
			   ex2loc[recvarray[ii+2]]);
		if (vert2face.Used(re))
		  AddDistantFaceProc(vert2face.Get(re)-1, dest);
	      }
	  }
	
	NgProfiler::StopTimer (timerf);
      }
    // cout << "UpdateCoarseGrid - done" << endl;
    // EnumeratePointsGlobally();
    is_updated = true;

    // MPI_Group_free(&MPI_LocalGroup);
    // MPI_Comm_free(&MPI_LocalComm);
  }
}


// #endif

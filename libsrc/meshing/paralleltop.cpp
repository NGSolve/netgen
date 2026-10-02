#include <meshing.hpp>
#include "paralleltop.hpp"


namespace netgen
{

  ParallelMeshTopology :: ParallelMeshTopology (const Mesh & amesh)
    : mesh(amesh)
  { ; }


  void ParallelMeshTopology :: SetNV_Loc2Glob (int anv)
  {
    glob_vert.SetSize(anv);
    glob_vert = -1;
  }

  void ParallelMeshTopology :: SetNV (int anv)
  {
    DynamicTable<int> oldtable(loc2distvert.Size());
    for (size_t i = 0; i < loc2distvert.Size(); i++)
      for (auto val : loc2distvert[i])
        oldtable.Add (i, val);
    loc2distvert = DynamicTable<int> (anv);
    for (size_t i = 0; i < min(size_t(anv), oldtable.Size()); i++)
      for (auto val : oldtable[i])
        loc2distvert.Add (i, val);
  }


  void ParallelMeshTopology :: EnumeratePointsGlobally ()
  {
    static Timer t("ParallelTopology::EnumeratePointsGlobally"); RegionTimer r(t);

    auto comm = mesh.GetCommunicator();
    auto rank = comm.Rank();

    size_t oldnv = glob_vert.Size();
    size_t nv = loc2distvert.Size();
    *testout << "enumerate globally, loc2distvert.size = " << loc2distvert.Size()
             << ", glob_vert.size = " << glob_vert.Size() << endl;


    // IntRange newvr(oldnv, nv); // new vertex range
    auto new_pir = Range(PointIndex::FromNr0(oldnv), PointIndex::FromNr0(nv));
    
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
    auto max_oldv = comm.AllReduce (Max (glob_vert.Range(0, oldnv)), NG_MPI_MAX);
    if (comm.AllReduce (oldnv, NG_MPI_SUM) == 0)
      max_oldv = long(PointIndex::BASE)-1;
    
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
    
    Table<int> send_data(nsend);   // global point numbers
    Table<int> recv_data(nrecv);
    
    /** Fill send_data **/
    nsend = 0;
    for (auto pi : new_pir)
      if (auto dps = GetDistantProcs(pi); dps.Size())
        if (rank < dps[0])
          for (auto p : dps)
            send_data[p][nsend[p]++] = L2G(pi);

    NgMPI_Requests requests;
    for (int i = 0; i < comm.Size(); i++)
      {
        if (nsend[i])
          requests += comm.ISend (send_data[i], i, 200);
        if (nrecv[i])
          requests += comm.IRecv (recv_data[i], i, 200);
      }
    
    requests.WaitAll();
    
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

    {
        Array<PointIndex, PointIndex> inv_index(index0.Size());
        for (int i = 0; i < index0.Size(); i++)
          inv_index[PointIndex::FromNr0(index0[i])] = PointIndex::FromNr0(i);
        
        for (auto el : mesh.VolumeElements())
          for (PointIndex & pi : el.PNums())
            pi = inv_index[pi];
        for (auto el : mesh.SurfaceElements())
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
            Array<PointIndices<2>,PointIndex> hml (mesh.mlbetweennodes);
            for (PointIndex pi : Range(mesh.Points()))
              mesh.mlbetweennodes[inv_index[pi]] = hml[pi];
          }


        DynamicTable<int> oldtable = std::move(loc2distvert);        
        loc2distvert = DynamicTable<int> (oldtable.Size());
        for (size_t i = 0; i < oldtable.Size(); i++)
          for (auto val : oldtable[index0[i]])
            loc2distvert.Add (i, val);

        Array<int> hglob_vert(glob_vert);
        for (int i = 0; i < index0.Size(); i++)
          glob_vert[i] = hglob_vert[index0[i]];

    }

    if (glob_vert.Size() > 1)
      for (auto i : Range(glob_vert).Modify(0,-1))
        if (glob_vert[i] > glob_vert[i+1])
          cout << "wrong ordering of globvert" << endl;
  }


  void ParallelMeshTopology :: IdentifyNewVertices ()
  {
    static Timer t("ParallelTopology::IdentifyNewVertices"); RegionTimer r(t);

    NgMPI_Comm comm = mesh.GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();
    if (ntasks == 1) return;

    if (loc2distvert.Size() != mesh.GetNV())
      SetNV (mesh.GetNV());

    Array<int> cnt_send(ntasks);

    int maxsize = comm.AllReduce (mesh.mlbetweennodes.Size(), NG_MPI_MAX);
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
            DynamicTable<PointIndex> dest2vert(cnt_send);    
            for (PointIndex pi : mesh.Points().Range())
              for (int dist : GetDistantProcs(pi))
                dest2vert.Add (dist, pi);
            
            for (PointIndex pi = IndexBASE<PointIndex>(); pi < newnv+IndexBASE<PointIndex>(); pi++)
              if (auto [v1,v2] = mesh.mlbetweennodes[pi]; v1.IsValid())              
                {
                  auto procs1 = GetDistantProcs(v1);
                  auto procs2 = GetDistantProcs(v2);
                  for (int p : procs1)
                    if (procs2.Contains(p))
                      cnt_send[p]++;
                }

            DynamicTable<PointIndex> dest2pair(cnt_send);            
            
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
            
            DynamicTable<int> send_verts(cnt_send);

            Array<int, PointIndex> loc2exchange(mesh.GetNV());

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
            comm.ExchangeTable (send_verts, recv_verts, NG_MPI_TAG_MESH+9);

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
                            IVec<2> re(recvarray[ii], recvarray[ii+1]);
                            IVec<2> es(loc2exchange[v1], loc2exchange[v2]);
                            if (es == re && !GetDistantProcs(pi).Contains(dest))
                              {
                                AddDistantProc (pi, dest);
                                changed = true;
                              }
                          }
                      }
                }

            changed = comm.AllReduce (changed, NG_MPI_LOR);
          }
      }
  }


  void ParallelMeshTopology :: UpdateEdgesAndFaces ()
  {
    static Timer t("ParallelTopology::UpdateEdgesAndFaces"); RegionTimer r(t);

    NgMPI_Comm comm = mesh.GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();
    if (ntasks == 1) return;

    if (id == 0)
      PrintMessage (3, "update parallel topology");

    const MeshTopology & topology = mesh.GetTopology();
    Array<int> cnt_send(ntasks);

    static Timer timere("UpdateEdgesAndFaces - edges");
    static Timer timerf("UpdateEdgesAndFaces - faces");
    timere.Start();

    // recomputed from scratch, the topology may have renumbered edges and faces
    loc2distedge = DynamicTable<int> (topology.GetNEdges());
    loc2distface = DynamicTable<int> (topology.GetNFaces());

    int nfa = topology . GetNFaces();
    int ned = topology . GetNEdges();
    
    // build exchange vertices
    cnt_send = 0;
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantProcs(pi))
        cnt_send[dist]++;
    DynamicTable<PointIndex> dest2vert(cnt_send);    
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantProcs(pi))
        dest2vert.Add (dist, pi);

    // exchange edges
    cnt_send = 0;
    for (int edge = 1; edge <= ned; edge++)
      {
        auto [v1,v2] = topology.GetEdgeVertices(EdgeIndex::FromNr1(edge));
        /*
        for (int dest = 1; dest < ntasks; dest++)
          if (GetDistantProcs(v1).Contains(dest) && GetDistantProcs(v2).Contains(dest))
            cnt_send[dest-1]+=1;
        */
        for (auto p : GetDistantProcs(v1))
          if (GetDistantProcs(v2).Contains(p))
            cnt_send[p]+=1;
      }
    
    DynamicTable<int> dest2edge(cnt_send);
    for (int & v : cnt_send) v *= 2;
    DynamicTable<int> send_edges(cnt_send);

    for (int edge = 1; edge <= ned; edge++)
      {
        auto [v1,v2] = topology.GetEdgeVertices(EdgeIndex::FromNr1(edge));        
        for (int dest = 0; dest < ntasks; dest++)
          if (GetDistantProcs(v1).Contains(dest) && GetDistantProcs(v2).Contains(dest))
            dest2edge.Add (dest, edge);
      }


    Array<int, PointIndex> loc2exchange(mesh.GetNV());
    for (int dest = 0; dest < ntasks; dest++)
      {
        loc2exchange = -1;
        int cnt = 0;
        for (PointIndex pi : dest2vert[dest])
          loc2exchange[pi] = cnt++;

        for (int edge : dest2edge[dest])
          {
            auto [v1,v2] = topology.GetEdgeVertices(EdgeIndex::FromNr1(edge));            
            if (GetDistantProcs(v1).Contains(dest) && GetDistantProcs(v2).Contains(dest))            
              {
                send_edges.Add (dest, loc2exchange[v1]);
                send_edges.Add (dest, loc2exchange[v2]);
              }
          }
      }

    DynamicTable<int> recv_edges(ntasks);
    comm.ExchangeTable (send_edges, recv_edges, NG_MPI_TAG_MESH+9);

    for (int dest = 0; dest < ntasks; dest++)
      {
        auto ex2loc = dest2vert[dest];
        if (ex2loc.Size() == 0) continue;

        ClosedHashTable<PointIndices<2>, int> vert2edge(4*dest2edge[dest].Size()+16);
        for (int edge : dest2edge[dest])
          {
            auto [v1,v2] = topology.GetEdgeVertices(EdgeIndex::FromNr1(edge));            
            vert2edge.Set(PointIndices<2>(v1,v2), edge);
          }

        FlatArray<int> recvarray = recv_edges[dest];
        for (int ii = 0; ii < recvarray.Size(); ii+=2)
          {
            PointIndices<2> re(ex2loc[recvarray[ii]], 
                               ex2loc[recvarray[ii+1]]);
            if (vert2edge.Used(re))
              AddDistantEdgeProc (vert2edge.Get(re)-1, dest);
          }
      }



    timere.Stop();

    if (mesh.GetDimension() == 3)
      {
        timerf.Start();

        // exchange faces
        cnt_send = 0;
        for (int face = 0; face < nfa; face++)
          {
            auto verts = topology.GetFaceVertices (FaceIndex::FromNr0(face));
            for (int dest = 0; dest < ntasks; dest++)
              if (dest != id)
                if (GetDistantProcs (verts[0]).Contains(dest) &&
                    GetDistantProcs (verts[1]).Contains(dest) &&
                    GetDistantProcs (verts[2]).Contains(dest))
                  cnt_send[dest]++;
          }
        
        DynamicTable<int> dest2face(cnt_send);
        for (int face = 1; face <= nfa; face++)
          {
            auto verts = topology.GetFaceVertices (FaceIndex::FromNr1(face));
            for (int dest = 0; dest < ntasks; dest++)
              if (dest != id)
                if (GetDistantProcs (verts[0]).Contains(dest) && 
                    GetDistantProcs (verts[1]).Contains(dest) &&
                    GetDistantProcs (verts[2]).Contains(dest))
                  dest2face.Add(dest, face);
          }

        for (int & c : cnt_send) c*=3;
        DynamicTable<int> send_faces(cnt_send);
        Array<int, PointIndex> loc2exchange(mesh.GetNV());
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
                  auto verts = topology.GetFaceVertices (FaceIndex::FromNr1(face));
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
        
        DynamicTable<int> recv_faces(ntasks);
        comm.ExchangeTable (send_faces, recv_faces, NG_MPI_TAG_MESH+9);
        
        for (int dest = 0; dest < ntasks; dest++)
          {
            auto ex2loc = dest2vert[dest];
            if (ex2loc.Size() == 0) continue;
            
            ClosedHashTable<PointIndices<3>, int> vert2face(4*dest2face[dest].Size()+16);
            for (int face : dest2face[dest])
              {
                auto verts = topology.GetFaceVertices (FaceIndex::FromNr1(face));
                vert2face.Set(PointIndices<3>(verts[0], verts[1], verts[2]), face);
              }
            
            FlatArray<int> recvarray = recv_faces[dest];
            for (int ii = 0; ii < recvarray.Size(); ii+=3)
              {
                PointIndices<3> re(ex2loc[recvarray[ii]], 
                                   ex2loc[recvarray[ii+1]],
                                   ex2loc[recvarray[ii+2]]);
                if (vert2face.Used(re))
                  AddDistantFaceProc(vert2face.Get(re)-1, dest);
              }
          }
        
        timerf.Stop();
      }
  }
}


#include <meshing.hpp>
#include "paralleltop.hpp"

// #define METIS4


#ifdef METIS
namespace metis {
  extern "C" {

#include <metis.h>

#if METIS_VER_MAJOR >= 5
#define METIS5
    typedef idx_t idxtype;   
#else
#define METIS4
    typedef idxtype idx_t;  
#endif
  } 
}

using namespace metis;
#endif

namespace netgen
{

  void Mesh :: SendRecvMesh ()
  {
    int id = GetCommunicator().Rank();
    int np = GetCommunicator().Size();

    if (np == 1) {
      throw NgException("SendRecvMesh called, but only one rank in communicator!!");
    }
    
    if (id == 0)
      PrintMessage (1, "Send/Receive mesh");

    if (id == 0)
      SendMesh ();
    else
      ReceiveParallelMesh();
  }


  void Mesh :: UpdateParallelTopology ()
  {
    if (GetCommunicator().Size() == 1) return;
    paralleltop->IdentifyNewVertices();
    paralleltop->EnumeratePointsGlobally();
  }


  /*
    The master walks the mesh once and writes, for every destination rank, one
    binary MemoryOutArchive with the part of the mesh that rank gets:

      dim
      global vertex numbers (0-based)      -> defines the local vertex numbering
      points                               -> Array<MeshPoint>
      identifications                      -> int table, local vertex numbers
      distant procs of vertices            -> pairs (local vertex, proc)
      volume elements                      (StridedElementArray::DoArchiveCurrent)
      regions of dimension 2, 1, 3, 0      -> descriptors including names
      surface elements, their geometry info
      segments, their geometry info
      point elements

    Element point numbers are mapped to the destination's local numbers before
    archiving, so the receiver reads straight into its arrays with the same
    DoArchive methods used for files and pickling.
  */
  void Mesh :: SendMesh () const   
  {
    static Timer tsend("SendMesh"); RegionTimer reg(tsend);
    static Timer tbuildvertex("SendMesh::BuildVertex");
    static Timer tbuildvertexa("SendMesh::BuildVertex a");
    static Timer tbuildvertexb("SendMesh::BuildVertex b");
    static Timer tarchive("SendMesh::Archive");
    
    NgMPI_Comm comm = GetCommunicator();
    int ntasks = comm.Size();
    auto & self = const_cast<Mesh&>(*this);

    int dim = GetDimension();

    // periodic identifications need the vertex2element tables
    bool has_periodic = false;
    {
      auto & idents = GetIdentifications();
      for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
        if (idents.GetType(idnr) == Identifications::PERIODIC)
          { has_periodic = true; break; }
    }

    // If the topology is not already updated, we do not need to
    // build edges/faces.
    auto & top = const_cast<MeshTopology&>(GetTopology());
    if(top.NeedsUpdate()) {
      top.SetBuildVertex2Element(has_periodic);
      top.SetBuildEdges(false);
      top.SetBuildFaces(false);
      top.Update();
    }
    
    PrintMessage ( 3, "Building element tables");
    
    Array<int> num_els_on_proc(ntasks);
    num_els_on_proc = 0;
    for (ElementIndex ei : VolumeElements().Range())
      num_els_on_proc[vol_partition[ei]]++;

    Table<ElementIndex> els_of_proc (num_els_on_proc);
    num_els_on_proc = 0;
    for (ElementIndex ei : VolumeElements().Range())
      {
        auto nr = vol_partition[ei];
        els_of_proc[nr][num_els_on_proc[nr]++] = ei;
      }
    
    Array<int> num_sels_on_proc(ntasks);
    num_sels_on_proc = 0;
    for (SurfaceElementIndex ei : SurfaceElements().Range())
      num_sels_on_proc[surf_partition[ei]]++;

    Table<SurfaceElementIndex> sels_of_proc (num_sels_on_proc);
    num_sels_on_proc = 0;
    for (SurfaceElementIndex ei : SurfaceElements().Range())
      {
        auto nr = surf_partition[ei];
        sels_of_proc[nr][num_sels_on_proc[nr]++] = ei;
      }
    

    Array<int> num_segs_on_proc(ntasks);
    num_segs_on_proc = 0;
    for (SegmentIndex ei : LineSegments().Range())
      // num_segs_on_proc[(*this)[ei].GetPartition()]++;
      num_segs_on_proc[seg_partition[ei]]++;

    DynamicTable<SegmentIndex> segs_of_proc (num_segs_on_proc);
    for (SegmentIndex ei : LineSegments().Range())
      segs_of_proc.Add (seg_partition[ei], ei);


    /**
            ----- STRATEGY FOR PERIODIC MESHES -----

       Whenever two vertices are identified by periodicity, any proc 
       that gets one of the vertices actually gets both of them.
       This has to be transitive, that is, if
       a <-> b and  b <-> c,
       then any proc that has vertex a also has vertices b and c!

       Surfaceelements and Segments that are identified by
       periodicity are treated the same way.
       
       We need to duplicate these so we have containers to
       hold the edges/facets. Afaik, a mesh cannot have nodes 
       that are not part of some sort of element.

     **/

    /** First, we build tables for vertex identification. **/
    Array<PointIndices<2>> per_pairs;
    Array<PointIndices<2>> pp2;
    auto & idents = GetIdentifications();
    for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
      {
        if(idents.GetType(idnr)!=Identifications::PERIODIC) continue;
        idents.GetPairs(idnr, pp2);
        per_pairs += pp2;
      }
    Array<int, PointIndex> npvs(GetNV());
    npvs = 0;
    for (auto [p1, p2] : per_pairs) {
      npvs[p1]++;
      npvs[p2]++;
    }

    /** for each vertex, gives us all identified vertices **/
    DynamicTable<PointIndex, PointIndex> per_verts(GetNV());
    for (auto [p1, p2] : per_pairs) {
      per_verts.Add(p1, p2);
      per_verts.Add(p2, p1);
    }
    for (PointIndex k : Range(PointIndex::FromNr0(0), PointIndex::FromNr0(GetNV()))) {
      BubbleSort(per_verts[k]);
    }

    /** The same table as per_verts, but TRANSITIVE!! **/
    auto iterate_per_verts_trans = [&](auto f){
      Array<PointIndex> allvs;
      // for (int k = PointIndex::BASE; k < GetNV()+PointIndex::BASE; k++)
      for (PointIndex k = IndexBASE<PointIndex>();
           k < GetNV()+IndexBASE<PointIndex>(); k++)      
        {
          allvs.SetSize(0);
          allvs.Append(per_verts[k]);
          bool changed = true;
          while(changed) {
            changed = false;
            for (int j = 0; j<allvs.Size(); j++)
              {
                auto pervs2 = per_verts[allvs[j]];
                for (int l = 0; l < pervs2.Size(); l++)
                  {
                    auto addv = pervs2[l];
                    if (allvs.Contains(addv) || addv==k) continue;
                    changed = true;
                    allvs.Append(addv);
                  }
              }
          }
          f(k, allvs);
        }
    };
    iterate_per_verts_trans([&](auto k, auto & allvs) {
        npvs[k] = allvs.Size();
      });
    DynamicTable<PointIndex, PointIndex> per_verts_trans(GetNV());
    iterate_per_verts_trans([&](auto k, auto & allvs) {
        for (int j = 0; j<allvs.Size(); j++)
          per_verts_trans.Add(k, allvs[j]);
      });
    for (PointIndex k : Range(PointIndex::FromNr0(0), PointIndex::FromNr0(GetNV()))) {
      BubbleSort(per_verts_trans[k]);
    }

    /** Now we build the vertex-data to send to the workers. **/
    tbuildvertex.Start();
    Array<int, PointIndex> vert_flag (GetNV());
    Array<int, PointIndex> num_procs_on_vert (GetNV());
    Array<int> num_verts_on_proc (ntasks);
    num_verts_on_proc = 0;
    num_procs_on_vert = 0;
    
    auto iterate_vertices = [&](auto f) {
      vert_flag = -1;
      for (int dest = 0; dest < ntasks; dest++)
        {
          for (auto ei : els_of_proc[dest])
            for (auto pnum : (*this)[ei].PNums())
              f(pnum, dest);

          for (auto ei : sels_of_proc[dest])
            for (auto pnum : (*this)[ei].PNums())
              f(pnum, dest);

          /*
          FlatArray<SegmentIndex> segs = segs_of_proc[dest];
          for (int hi = 0; hi < segs.Size(); hi++)
            {
              const Segment & el = (*this) [segs[hi]];
              for (int i = 0; i < 2; i++)
                f(el[i], dest);
            }
          */
          for (auto segi : segs_of_proc[dest])
            for (auto pnum : (*this)[segi].PNums())
              f(pnum, dest);
        }
    };
    /** count vertices per proc and procs per vertex **/
    tbuildvertexa.Start();    
    iterate_vertices([&](auto vertex, auto dest){
        auto countit = [&] (auto vertex, auto dest) {
          if (vert_flag[vertex] < dest)
            {
              vert_flag[vertex] = dest;
              num_verts_on_proc[dest]++;
              num_procs_on_vert[vertex]++;
            }
        };
        countit(vertex, dest);
        for (auto v : per_verts_trans[vertex])
          countit(v, dest);
      });
    tbuildvertexa.Stop();    

    
    tbuildvertexb.Start();    
    
    DynamicTable<int> verts_of_proc (num_verts_on_proc);   // 0-based offsets into points
    DynamicTable<int, PointIndex> procs_of_vert (GetNV());
    DynamicTable<PointIndex, PointIndex> loc_num_of_vert (GetNV());
    /** Write vertex/proc mappingfs to tables **/
    iterate_vertices([&](auto vertex, auto dest) {
        auto addit = [&] (auto vertex, auto dest) {
          if (vert_flag[vertex] < dest)
            {
              vert_flag[vertex] = dest;
              procs_of_vert.Add (vertex, dest);
            }
        };
        addit(vertex, dest);
        for (auto v : per_verts_trans[vertex])
          addit(v, dest);
      });
    tbuildvertexb.Stop();        
    /** 
        local vertex numbers on distant procs 
        (I think this was only used for debugging??) 
    **/
    // for (int vert = 1; vert <= GetNP(); vert++ )
    for (PointIndex vert : Points().Range())
      {
        FlatArray<int> procs = procs_of_vert[vert];
        for (int j = 0; j < procs.Size(); j++)
          {
            int dest = procs[j];
            // !! we also use this as offsets for MPI-type, if this is changed, also change ReceiveParallelMesh
            verts_of_proc.Add (dest, vert - IndexBASE<T_POINTS::index_type>());
            loc_num_of_vert.Add (vert, verts_of_proc[dest].Size() -1+IndexBASE<T_POINTS::index_type>());
          }
      }
    tbuildvertex.Stop();    

    // local number of a global vertex on dest, INVALID if the vertex is not sent there
    auto loc_num = [&] (PointIndex vert, int dest) -> PointIndex
    {
      auto procs = procs_of_vert[vert];
      for (int j = 0; j < procs.Size(); j++)
        if (procs[j] == dest) return loc_num_of_vert[vert][j];
      return PointIndex::INVALID;
    };

    /** periodic identifications, per destination, in local numbers:
        maxidentnr, type of each ident, nr of pairs of each ident, pairs **/
    PrintMessage ( 3, "Building identifications");
    int maxidentnr = idents.GetMaxNr();
    Array<int> ppd_sizes(ntasks);
    ppd_sizes = 1 + 2*maxidentnr;
    for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
      {
        if(idents.GetType(idnr)!=Identifications::PERIODIC) continue;
        idents.GetPairs(idnr, pp2);
        for (auto pair : pp2)
          for (auto p : procs_of_vert[pair[0]])   // both are on the same procs
            ppd_sizes[p] += 2;
      }
    DynamicTable<int> pp_data(ppd_sizes);
    for (int dest = 0; dest < ntasks; dest++)
      {
        pp_data.Add(dest, maxidentnr);
        for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
          pp_data.Add(dest, idents.GetType(idnr));
        for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
          pp_data.Add(dest, 0);
      }
    for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
      {
        if(idents.GetType(idnr)!=Identifications::PERIODIC) continue;
        idents.GetPairs(idnr, pp2);
        for (auto pair : pp2)
          for (auto p : procs_of_vert[pair[0]])
            {
              PointIndex l0 = loc_num(pair[0], p), l1 = loc_num(pair[1], p);
              if (!l0.IsValid() || !l1.IsValid()) continue;
              pp_data[p][maxidentnr + idnr]++;
              pp_data.Add(p, l0.Nr0());
              pp_data.Add(p, l1.Nr0());
            }
      }

    /** distant procs of the vertices: pairs (local vertex, proc) **/
    PrintMessage ( 3, "Building distant procs");
    Array<int> num_distpnums(ntasks);
    num_distpnums = 0;
    for (PointIndex vert : Points().Range())
      {
        FlatArray<int> procs = procs_of_vert[vert];
        for (auto p : procs)
          num_distpnums[p] += 2 * (procs.Size()-1);
      }
    DynamicTable<int> distpnums (num_distpnums);
    for (PointIndex vert : Points().Range())
      {
        FlatArray<int> procs = procs_of_vert[vert];
        for (int j = 0; j < procs.Size(); j++)
          for (int k = 0; k < procs.Size(); k++)
            if (j != k)
              {
                distpnums.Add (procs[j], loc_num_of_vert[vert][j].Nr0());
                distpnums.Add (procs[j], procs[k]);
              }
      }

    /** surface elements, segments and point elements per destination (periodic copies included) **/
    PrintMessage ( 3, "Building surface element and segment tables");
    // build sel-identification
    size_t nse = GetNSE();
    Array<SurfaceElementIndex, SurfaceElementIndex> ided_sel(nse);
    ided_sel = SurfaceElementIndex::INVALID;
    [[maybe_unused]] bool has_ided_sels = false;
    if(GetNE() && has_periodic) //we can only have identified surf-els if we have vol-els (right?)
      {
        Array<SurfaceElementIndex> os1, os2;
        for (SurfaceElementIndex sei : SurfaceElements().Range())
          {
            if(ided_sel[sei].IsValid()) continue;
            const Element2dRef & sel = (*this)[sei];
            auto points = sel.PNums();
            auto ided1 = per_verts[points[0]];
            os1.SetSize(0);
            for (int j = 0; j < ided1.Size(); j++)
              os1.Append(GetTopology().GetVertexSurfaceElements(ided1[j]));
            for (int j = 1; j < points.Size(); j++)
              {
                os2.SetSize(0);
                auto p2 = points[j];
                auto ided2 = per_verts[p2];
                for (int l = 0; l < ided2.Size(); l++)
                  os2.Append(GetTopology().GetVertexSurfaceElements(ided2[l]));
                for (int m = 0; m<os1.Size(); m++) {
                  if(!os2.Contains(os1[m])) {
                    os1.DeleteElement(m);
                    m--;
                  }
                }
              }
            if(!os1.Size()) continue;
            if(os1.Size()>1) {
              throw NgException("SurfaceElement identified with more than one other??");
            }
            // const Element2dRef & sel2 = (*this)[sei];
            // auto points2 = sel2.PNums();
            has_ided_sels = true;
            ided_sel[sei] = os1[0];
            ided_sel[os1[0]] = sei;
          }
      }
    auto iterate_sels = [&](auto f) {
      for (SurfaceElementIndex sei : SurfaceElements().Range())
        {
          const Element2dRef & sel = (*this)[sei];
          // int dest = (*this)[sei].GetPartition();
          int dest = surf_partition[sei];
          f(sei, sel, dest);
          if(ided_sel[sei].IsValid())
            {
              // int dest2 = (*this)[ided_sel[sei]].GetPartition();
              int dest2 = surf_partition[ided_sel[sei]];
              f(sei, sel, dest2);
            }
        }      
    };
    DynamicTable<SurfaceElementIndex> sels_to_send(ntasks);
    iterate_sels([&](SurfaceElementIndex sei, const Element2dRef & sel, int dest)
                 { sels_to_send.Add (dest, sei); });

    auto iterate_segs1 = [&](auto f) {
      Array<SegmentIndex> osegs1, osegs2, osegs_both;
      Array<int> type1, type2;
      for (SegmentIndex segi : LineSegments().Range())
        {
          const Segment & seg = (*this)[segi];
          int segnp = seg.GetNP();
          PointIndex pi1 = seg[0];
          auto ided1 = per_verts[pi1];
          PointIndex pi2 = seg[1];
          auto ided2 = per_verts[pi2];
          if (!(ided1.Size() && ided2.Size())) continue;
          osegs1.SetSize(0);
          type1.SetSize(0);
          for (int l = 0; l<ided1.Size(); l++)
            {
              auto ospart = GetTopology().GetVertexSegments(ided1[l]);
              for(int j=0; j<ospart.Size(); j++)
                {
                  if(osegs1.Contains(ospart[j]))
                    throw NgException("Periodic Mesh did something weird.");
                  osegs1.Append(ospart[j]);
                  type1.Append(idents.GetSymmetric(pi1, ided1[l]));
                }
            }
          osegs2.SetSize(0);
          type2.SetSize(0);
          for (int l = 0; l<ided2.Size(); l++)
            {
              auto ospart = GetTopology().GetVertexSegments(ided2[l]);
              for(int j=0; j<ospart.Size(); j++)
                {
                  if(osegs2.Contains(ospart[j]))
                    throw NgException("Periodic Mesh did something weird.");
                  osegs2.Append(ospart[j]);
                  type2.Append(idents.GetSymmetric(pi2, ided2[l]));
                }
            }
          osegs_both.SetSize(0);
          for (int l = 0; l<osegs1.Size(); l++) {
            auto pos = osegs2.Pos(osegs1[l]);
            if (pos == -1) continue;
            if (type1[l] != type2[pos]) continue;
            osegs_both.Append(osegs1[l]);
          }
          for(int l = 0; l<osegs_both.Size(); l++) {
            int segnp2 = (*this)[osegs_both[l]].GetNP();
            if(segnp!=segnp2)
              throw NgException("Tried to identify non-curved and curved Segment!");
          }
          for(int l = 0; l<osegs_both.Size(); l++) {
            f(segi, osegs_both[l]);
          }
        }
    };
    Array<int, SegmentIndex> per_seg_size(GetNSeg());
    per_seg_size = 0;
    iterate_segs1([&](SegmentIndex segi1, SegmentIndex segi2)
                  { per_seg_size[segi1]++; });
    DynamicTable<SegmentIndex, SegmentIndex> per_seg(per_seg_size);
    iterate_segs1([&](SegmentIndex segi1, SegmentIndex segi2)
                  { per_seg.Add(segi1, segi2); });
    // make per_seg transitive
    auto iterate_per_seg_trans = [&](auto f){
      Array<SegmentIndex> allsegs;
      for (SegmentIndex segi : LineSegments().Range())
        {
          allsegs.SetSize(0);
          allsegs.Append(per_seg[segi]);
          bool changed = true;
          while (changed)
            {
              changed = false;
              for (int j = 0; j<allsegs.Size(); j++)
                {
                  auto persegs2 = per_seg[allsegs[j]];
                  for (int l = 0; l<persegs2.Size(); l++)
                    {
                      auto addseg = persegs2[l];
                      if (allsegs.Contains(addseg) || addseg==segi) continue;
                      allsegs.Append(addseg);
                      changed = true;
                    }
                }
            }
          f(segi, allsegs);
        }
    };
    iterate_per_seg_trans([&](SegmentIndex segi, Array<SegmentIndex> & segs){
        for (int j = 0; j < segs.Size(); j++)
          per_seg_size[segi] = segs.Size();
      });
    DynamicTable<SegmentIndex, SegmentIndex> per_seg_trans(per_seg_size);
    iterate_per_seg_trans([&](SegmentIndex segi, Array<SegmentIndex> & segs){
        for (int j = 0; j < segs.Size(); j++)
          per_seg_trans.Add(segi, segs[j]);
      });
    Array<int> dests;
    auto iterate_segs2 = [&](auto f)
      {
        for (SegmentIndex segi : LineSegments().Range())
          {
            const Segment & seg = (*this)[segi];
            dests.SetSize(0);
            // dests.Append(seg.GetPartition());
            dests.Append(seg_partition[segi]);
            for (int l = 0; l < per_seg_trans[segi].Size(); l++)
              {
                // int dest2 = (*this)[per_seg_trans[segi][l]].GetPartition();
                int dest2 = seg_partition[per_seg_trans[segi][l]];
                if(!dests.Contains(dest2))
                  dests.Append(dest2);
              }
            for (int l = 0; l < dests.Size(); l++)
              f(segi, seg, dests[l]);
          }
      };
    DynamicTable<SegmentIndex> segs_to_send(ntasks);
    iterate_segs2([&](auto segi, const auto & seg, int dest)
                  {
                    for (auto pi : seg.PNums())
                      if (!loc_num(pi, dest).IsValid()) return;
                    segs_to_send.Add (dest, segi);
                  });

    DynamicTable<int> pels_to_send(ntasks);
    for (auto k : Range(pointelements))
      for (auto dest : procs_of_vert[pointelements[k].pnum])
        pels_to_send.Add (dest, k);

    /** serialize and send **/
    PrintMessage ( 3, "Sending mesh");
    tarchive.Start();
    Array<unique_ptr<MemoryOutArchive>> archives(ntasks);   // outlive the requests
    NgMPI_Requests sendrequests;

    auto archive_array = [] (Archive & ar, auto data)   // compatible with Array<T>::DoArchive
    {
      size_t n = data.Size();
      ar & n;
      ar.Do (data.Data(), n);
    };

    for (int dest = 0; dest < ntasks; dest++)
      {
        archives[dest] = make_unique<MemoryOutArchive>();
        auto & ar = *archives[dest];

        ar & dim;

        // vertices
        FlatArray<int> verts = verts_of_proc[dest];
        archive_array (ar, verts);
        size_t nv = verts.Size();
        ar & nv;
        for (int v : verts)
          ar & self.points[PointIndex::FromNr0(v)];

        archive_array (ar, pp_data[dest]);
        archive_array (ar, distpnums[dest]);

        // volume elements
        {
          auto els = els_of_proc[dest];
          size_t n = els.Size(), w = volelements.Width();   // as StridedElementArray::DoArchiveCurrent
          ar & n & w;
          for (auto ei : els)
            {
              Element el ((*this)[ei]);
              for (auto & pi : el.PNums()) pi = loc_num (pi, dest);
              el.DoArchive (ar);
            }
        }

        ar & self.Regions<2>() & self.Regions<1>() & self.Regions<3>() & self.Regions<0>();

        // surface elements
        {
          auto sels = sels_to_send[dest];
          size_t n = sels.Size(), w = surfelements.Width();
          ar & n & w;
          for (auto sei : sels)
            {
              Element2d el ((*this)[sei]);
              for (auto & pi : el.PNums()) pi = loc_num (pi, dest);
              el.DoArchive (ar);
            }
          for (auto sei : sels)
            Element2d((*this)[sei]).DoArchiveGeomInfo (ar);
        }

        // segments
        {
          auto segs = segs_to_send[dest];
          size_t n = segs.Size();
          ar & n;
          for (auto segi : segs)
            {
              Segment seg = (*this)[segi];
              for (auto & pi : seg.PNums()) pi = loc_num (pi, dest);
              seg.DoArchive (ar);
            }
          for (auto segi : segs)
            Segment((*this)[segi]).DoArchiveGeomInfo (ar);
        }

        // point elements
        {
          auto pels = pels_to_send[dest];
          size_t n = pels.Size();
          ar & n;
          for (auto k : pels)
            {
              Element0d el = pointelements[k];
              el.pnum = loc_num (el.pnum, dest);
              ar & el;
            }
        }

        auto & data = ar.Data();
        if (data.size() > size_t(std::numeric_limits<int>::max()))
          throw NgException("SendMesh: mesh part for rank " + ToString(dest) + " exceeds the MPI message size limit");
        if (dest != comm.Rank())   // own part is unpacked below
          sendrequests += comm.ISend (FlatArray<std::byte>(data.size(), data.data()), dest, NG_MPI_TAG_MESH);
      }
    tarchive.Stop();

    PrintMessage ( 3, "now wait ...");
    sendrequests.WaitAll();

    PrintMessage( 3, "Clean up local memory");

    self.points = T_POINTS(0);
    self.surfelements = T_SURFELEMENTS();
    self.volelements = T_VOLELEMENTS();
    self.segments = Array<Segment>(0);
    self.hp_surfinfo.SetSize(0);
    self.hp_volinfo.SetSize(0);
    self.hp_seginfo.SetSize(0);
    self.pointelements = Array<Element0d>(0);
    self.lockedpoints = Array<PointIndex>(0);
    /*
    auto cleanup_ptr = [](auto & ptr) {
      if (ptr != nullptr) {
        delete ptr;
        ptr = nullptr;
      }
    };
    cleanup_ptr(self.boundaryedges);
    cleanup_ptr(self.segmentht);
    cleanup_ptr(self.surfelementht);
    */
    self.boundaryedges = nullptr;
    self.segmentht = nullptr;
    self.surfelementht = nullptr;
    
    self.openelements = Array<Element2d>(0);
    self.opensegments = Array<Segment>(0);
    self.numvertices = -1;   // all points are vertices until ComputeNVertices
    self.mlbetweennodes = Array<PointIndices<2>,PointIndex> (0);
    self.mlparentelement = Array<ElementIndex, ElementIndex>(0);
    self.mlparentsurfaceelement = Array<SurfaceElementIndex, SurfaceElementIndex>(0);
    self.curvedelems = make_unique<CurvedElements> (self);
    self.clusters = make_unique<AnisotropicClusters> (self);
    self.ident = make_unique<Identifications> (self);
    self.topology = MeshTopology(*this);
    self.vol_partition.SetSize(0);
    self.surf_partition.SetSize(0);
    self.seg_partition.SetSize(0);

    // the own part: empty unless root_participates
    auto & own = archives[comm.Rank()]->Data();
    self.UnpackMeshPart (FlatArray<std::byte>(own.size(), own.data()));

    PrintMessage( 3, "send mesh complete");
  }




  // workers receive the mesh from the master
  void Mesh :: ReceiveParallelMesh ( )
  {
    static Timer timer("ReceiveParallelMesh"); RegionTimer reg(timer);

    Array<std::byte> buffer;
    GetCommunicator().Recv (buffer, 0, NG_MPI_TAG_MESH);
    UnpackMeshPart (buffer);
  }


  // collective: every rank unpacks its part, the global enumeration needs all of them
  void Mesh :: UnpackMeshPart (FlatArray<std::byte> data)
  {
    static Timer timer("UnpackMeshPart"); RegionTimer reg(timer);
    static Timer timer_unpack("Unpack mesh");

    NgMPI_Comm comm = GetCommunicator();
    int id = comm.Rank();

    timer_unpack.Start();
    MemoryInArchive ar(data.Data(), data.Size());

    int dim;
    ar & dim;
    SetDimension(dim);

    // vertices
    Array<int> verts;
    ar & verts;
    int numvert = verts.Size();
    paralleltop -> SetNV (numvert);
    paralleltop -> SetNV_Loc2Glob (numvert);
    for (int vert = 0; vert < numvert; vert++)
      paralleltop->L2G (PointIndex::FromNr0(vert)) = PointIndex::FromNr0(verts[vert]).Nr1();
    ar & points;

    // identifications
    Array<int> pp_data;
    ar & pp_data;
    int maxidentnr = pp_data[0];
    auto & idents = GetIdentifications();
    for (int idnr = 1; idnr < maxidentnr+1; idnr++)
      idents.SetType(idnr, (Identifications::ID_TYPE)pp_data[idnr]);
    int offset = 2*maxidentnr+1;
    for (int idnr = 1; idnr < maxidentnr+1; idnr++)
      {
        int npairs = pp_data[maxidentnr+idnr];
        FlatArray<int> pairdata(2*npairs, &pp_data[offset]);
        offset += 2*npairs;
        for (int k = 0; k < npairs; k++)
          idents.Add (PointIndex::FromNr0(pairdata[2*k]), PointIndex::FromNr0(pairdata[2*k+1]), idnr);
      }

    // distant procs
    Array<int> dist_pnums; 
    ar & dist_pnums;
    for (int hi = 0; hi < dist_pnums.Size(); hi += 2)
      paralleltop -> AddDistantProc (PointIndex::FromNr0(dist_pnums[hi]), dist_pnums[hi+1]);
    *testout << "got " << numvert << " vertices" << endl;

    volelements.DoArchiveCurrent (ar);

    ar & Regions<2>() & Regions<1>() & Regions<3>() & Regions<0>();

    surfelements.DoArchiveCurrent (ar);
    for (auto el : surfelements)
      el.DoArchiveGeomInfo (ar);

    ar & segments;
    for (auto & seg : segments)
      seg.DoArchiveGeomInfo (ar);

    ar & pointelements;
    timer_unpack.Stop();

    RebuildSurfaceElementLists();
    RebuildFDIndices();

    UpdateParallelTopology();

    static Timer timerloc("Update local mesh");
    static Timer timerloc2("CalcSurfacesOfNode");

    RegionTimer regloc(timerloc);
    stringstream str;
    str << "p" << id << ": got " << GetNE() << " elements and " 
         << GetNSE() << " surface elements";
    PrintMessage(2, str.str());
    // cout << str.str() << endl;
    // PrintMessage (2, "Got ", GetNE(), " elements and ", GetNSE(), " surface elements");
    // PrintMessage (2, "Got ", GetNSE(), " surface elements");

    timerloc2.Start();

    CalcSurfacesOfNode ();

    timerloc2.Stop();

    UpdateTopology();   // includes the shared edges and faces
    SetNextMajorTimeStamp();
  }




  /*
    Every rank packs its part with global point numbers (points, elements with
    geometry info, identification pairs) into a MemoryOutArchive; the root
    unpacks all parts into one mesh. Partition arrays and curving are not
    gathered: the result is a plain serial mesh, curve it again if needed.
  */
  shared_ptr<Mesh> Mesh :: GatherToRoot (int root) const
  {
    static Timer t("Mesh::GatherToRoot"); RegionTimer r(t);
    NgMPI_Comm comm = GetCommunicator();

    if (comm.Size() == 1)
      {
        auto m = make_shared<Mesh>();
        *m = *this;
        return m;
      }

    auto & self = const_cast<Mesh&>(*this);
    auto & partop = GetParallelTopology();

    Array<PointIndex, PointIndex> globnum(points.Size());
    PointIndex maxglob = PointIndex::INVALID;
    for (auto pi : Range(points))
      {
        globnum[pi] = PointIndex::FromNr1(partop.GetGlobalPNum(pi));
        maxglob = max(globnum[pi], maxglob);
      }
    maxglob = comm.AllReduce (maxglob, NG_MPI_MAX);
    int numglob = maxglob+1-IndexBASE<PointIndex>();

    auto pack = [&] (Archive & ar)
    {
      auto renum = [&] (auto pnums) { for (auto & pi : pnums) if (pi.IsValid()) pi = globnum[pi]; };
      ar & globnum;
      ar & self.points;
      T_VOLELEMENTS el3d (volelements);
      for (auto el : el3d) renum (el.PNums());
      el3d.DoArchiveCurrent (ar);
      T_SURFELEMENTS el2d (surfelements);
      for (auto el : el2d) renum (el.PNums());
      el2d.DoArchiveCurrent (ar);
      for (auto el : el2d) el.DoArchiveGeomInfo (ar);
      Array<Segment, SegmentIndex> el1d (segments);
      for (auto & seg : el1d) renum (seg.PNums());
      ar & el1d;
      for (auto & seg : el1d) seg.DoArchiveGeomInfo (ar);
      Array<Element0d> el0d (pointelements);
      for (auto & el : el0d) if (el.pnum.IsValid()) el.pnum = globnum[el.pnum];
      ar & el0d;
      // identification pairs: idnr, p1, p2 (1-based global)
      Array<int> idpairs;
      Array<PointIndices<2>> pairs;
      for (int idnr = 1; idnr <= ident->GetMaxNr(); idnr++)
        {
          ident->GetPairs (idnr, pairs);
          for (auto [p1, p2] : pairs)
            { idpairs += idnr; idpairs += globnum[p1].Nr1(); idpairs += globnum[p2].Nr1(); }
        }
      ar & idpairs;
    };

    if (comm.Rank() != root)
      {
        MemoryOutArchive out;
        pack (out);
        auto & data = out.Data();
        if (data.size() > size_t(std::numeric_limits<int>::max()))
          throw NgException("GatherToRoot: mesh part exceeds the MPI message size limit");
        comm.Send (FlatArray<std::byte>(data.size(), data.data()), root, NG_MPI_TAG_MESH+1);
        return nullptr;
      }

    auto result = make_shared<Mesh>();
    Mesh & m = *result;
    m.SetDimension (GetDimension());
    for (int idnr = 1; idnr <= ident->GetMaxNr(); idnr++)
      m.ident->SetType (idnr, ident->GetType(idnr));
    Array<MeshPoint, PointIndex> globpoints(numglob);

    auto unpack = [&] (Archive & ar)
    {
      Array<PointIndex, PointIndex> gn;
      Array<MeshPoint, PointIndex> pts;
      ar & gn & pts;
      for (auto i : Range(gn)) globpoints[gn[i]] = pts[i];
      T_VOLELEMENTS el3d;
      el3d.DoArchiveCurrent (ar);
      for (auto el : el3d) m.volelements.Append (el);
      T_SURFELEMENTS el2d;
      el2d.DoArchiveCurrent (ar);
      for (auto el : el2d) el.DoArchiveGeomInfo (ar);
      for (auto el : el2d) m.surfelements.Append (el);
      Array<Segment, SegmentIndex> el1d;
      ar & el1d;
      for (auto & seg : el1d) seg.DoArchiveGeomInfo (ar);
      for (auto & seg : el1d) m.segments.Append (seg);
      Array<Element0d> el0d;
      ar & el0d;
      for (auto & el : el0d) m.pointelements.Append (el);
      Array<int> idpairs;
      ar & idpairs;
      for (int k = 0; k < idpairs.Size(); k += 3)
        m.ident->Add (PointIndex::FromNr1(idpairs[k+1]), PointIndex::FromNr1(idpairs[k+2]), idpairs[k]);
    };

    {  // own part
      MemoryOutArchive out;
      pack (out);
      auto & data = out.Data();
      MemoryInArchive in (data.data(), data.size());
      unpack (in);
    }
    for (int src = 0; src < comm.Size(); src++)
      if (src != root)
        {
          Array<std::byte> buffer;
          comm.Recv (buffer, src, NG_MPI_TAG_MESH+1);
          MemoryInArchive in (buffer.Data(), buffer.Size());
          unpack (in);
        }

    m.points = std::move(globpoints);
    m.numvertices = numglob;
    m.Regions<3>() = Regions<3>();
    m.Regions<2>() = Regions<2>();
    m.Regions<1>() = Regions<1>();
    m.Regions<0>() = Regions<0>();
    m.SetGeometry (GetGeometry());

    m.RebuildSurfaceElementLists();
    m.RebuildFDIndices();
    m.CalcSurfacesOfNode();
    m.topology.Update();
    m.clusters->Update();
    m.SetNextMajorTimeStamp();
    return result;
  }


  // distribute the mesh to the worker processors
  // call it only for the master !
  void Mesh :: Distribute (bool root_participates)
  {
    NgMPI_Comm comm = GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    if (id != 0 || ntasks == 1 ) return;

    if (vol_partition.Size() < GetNE() || surf_partition.Size() < GetNSE() ||
        seg_partition.Size() < GetNSeg())
      ParallelMetis (comm.Size(), root_participates);

    /*
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
    for (SurfaceElementIndex ei : SurfaceElements().Range())
      *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
      */
    
    // MyMPI_SendCmd ("mesh");
    SendRecvMesh (); 
  }
  

#ifdef METIS5
  void Mesh :: ParallelMetis (int nproc, bool root_participates)
  {
    PrintMessage (3, "call metis 5 ...");

    static Timer timer("Mesh::Partition");
    RegionTimer reg(timer);

    idx_t ne = GetNE() + GetNSE() + GetNSeg();
    idx_t nn = GetNP();

    Array<idx_t> eptr, eind;
    for (int i = 0; i < GetNE(); i++)
      {
        eptr.Append (eind.Size());
        auto el = (*this)[ElementIndex::FromNr1(i+1)];
        for (int j = 0; j < el.GetNP(); j++)
          eind.Append (el[j].Nr0());
      }
    for (int i = 0; i < GetNSE(); i++)
      {
        eptr.Append (eind.Size());
        const Element2dRef & el = (*this)[SurfaceElementIndex::FromNr1(i+1)];
        for (int j = 0; j < el.GetNP(); j++)
          eind.Append (el[j].Nr0());
      }
    for (int i = 0; i < GetNSeg(); i++)
      {
        eptr.Append (eind.Size());
        const Segment & el = (*this)[SegmentIndex::FromNr1(i+1)];
        eind.Append (el[0].Nr0());
        eind.Append (el[1].Nr0());
      }
    eptr.Append (eind.Size());
    Array<idx_t> epart(ne), npart(nn);

    // partition numbers are destination ranks
    int first_rank = root_participates ? 0 : 1;
    idxtype nparts = nproc - first_rank;

    vol_partition.SetSize(GetNE());
    surf_partition.SetSize(GetNSE());
    seg_partition.SetSize(GetNSeg());
    if (nparts == 1)
      {
        for (int i = 0; i < GetNE(); i++)
          vol_partition[ElementIndex::FromNr0(i)]= first_rank;
        for (int i = 0; i < GetNSE(); i++)
          surf_partition[SurfaceElementIndex::FromNr0(i)] = first_rank;
        for (int i = 0; i < GetNSeg(); i++)
          seg_partition[SegmentIndex::FromNr0(i)] = first_rank;
      }

    else
      
      {

        idxtype edgecut;
        
        idxtype ncommon = GetDimension();
        PrintMessage (3, "metis start");

        static Timer tm("metis library");
        tm.Start();
        METIS_PartMeshDual (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &ncommon, &nparts,
                            NULL, NULL,
                            &edgecut, &epart[0], &npart[0]);
        tm.Stop();

        PrintMessage (3, "metis complete");
        
        for (int i = 0; i < GetNE(); i++)
          vol_partition[ElementIndex::FromNr0(i)]= epart[i] + first_rank;
        for (int i = 0; i < GetNSE(); i++)
          surf_partition[SurfaceElementIndex::FromNr0(i)] = epart[i+GetNE()] + first_rank;
        for (int i = 0; i < GetNSeg(); i++)
          seg_partition[SegmentIndex::FromNr0(i)] = epart[i+GetNE()+GetNSE()] + first_rank;
      }
    
        
    // surface elements attached to volume elements
    Array<bool, PointIndex> boundarypoints (GetNP());
    boundarypoints = false;

    if(GetDimension() == 3)
      for (auto sel : SurfaceElements())
        {
          const Element2dRef & el = sel;
          for (int j = 0; j < el.GetNP(); j++)
            boundarypoints[el[j]] = true;
        }
    else
      for (auto & seg : LineSegments())
        {
          for (int j = 0; j < 2; j++)
            boundarypoints[seg[j]] = true;
        }

    
    // Build Pnt2Element table, boundary points only
    Array<int, PointIndex> cnt(GetNP());
    cnt = 0;

    auto loop_els_2d = [&](auto f) {
      for (SurfaceElementIndex sei : SurfaceElements().Range())
        {
          const Element2dRef & el = (*this)[sei];
          for (int j = 0; j < el.GetNP(); j++) {
            f(el[j], sei);
          }
        }
    };
    auto loop_els_3d = [&](auto f) {
      for (ElementIndex ei : VolumeElements().Range())
        {
          auto el = (*this)[ei];
          for (int j = 0; j < el.GetNP(); j++)
            f(el[j], ei);
        }
    };
    auto loop_els = [&](auto f)
      {
        if (GetDimension() == 3 ) 
          loop_els_3d(f);
        else
          loop_els_2d(f);
      };

    
    loop_els([&](auto vertex, auto index)
        {
          if(boundarypoints[vertex])
            cnt[vertex]++;
        });
    DynamicTable<int, PointIndex> pnt2el(GetNP());
    loop_els([&](auto vertex, auto index)
        {
          if(boundarypoints[vertex])
            pnt2el.Add(vertex, index.Nr0());
        });


    if (GetDimension() == 3)
      {
        for (SurfaceElementIndex sei : SurfaceElements().Range())
          {
            Element2dRef sel = (*this)[sei];
            PointIndex pi1 = sel[0];
            // FlatArray<ElementIndex> els = pnt2el[pi1];
            FlatArray<int> els = pnt2el[pi1];
            
            // sel.SetPartition (-1);
            surf_partition[sei] = -1;
            
            for (int j = 0; j < els.Size(); j++)
              {
                auto el = (*this)[ElementIndex::FromNr0(els[j])];
                
                bool hasall = true;
                
                for (int k = 0; k < sel.GetNP(); k++)
                  {
                    bool haspi = false;
                    for (int l = 0; l < el.GetNP(); l++)
                      if (sel[k] == el[l])
                        haspi = true;

                    if (!haspi) hasall = false;
                  }
                
                if (hasall)
                  {
                    // sel.SetPartition (el.GetPartition());
                    surf_partition[sei] = vol_partition[ElementIndex::FromNr0(els[j])];
                    break;
                  }
              }
            // if (sel.GetPartition() == -1)
            if (surf_partition[sei] == -1)
              cerr << "no volume element found" << endl;
          }


        for (SegmentIndex si : LineSegments().Range())
          {
            Segment & sel = (*this)[si];
            PointIndex pi1 = sel[0];
            FlatArray<int> els = pnt2el[pi1];
            
            // sel.SetPartition (-1);
            seg_partition[si] = -1;
            
            for (int j = 0; j < els.Size(); j++)
              {
                auto el = (*this)[ElementIndex::FromNr0(els[j])];
                
                bool haspi[9] = { false };  // max surfnp
                
                for (int k = 0; k < 2; k++)
                  for (int l = 0; l < el.GetNP(); l++)
                    if (sel[k] == el[l])
                      haspi[k] = true;
                
                bool hasall = true;
                for (int k = 0; k < sel.GetNP(); k++)
                  if (!haspi[k]) hasall = false;
                
                if (hasall)
                  {
                    // sel.SetPartition (el.GetPartition());
                    seg_partition[si] = vol_partition[ElementIndex::FromNr0(els[j])];
                    break;
                  }
              }
            // if (sel.GetPartition() == -1)
            if (seg_partition[si] == -1)
              cerr << "no volume element found" << endl;
          }
      }
    else
      {
        for (SegmentIndex segi : LineSegments().Range())
          {
            Segment & seg = (*this)[segi];
            // seg.SetPartition(-1);
            seg_partition[segi] = -1;
            PointIndex pi1 = seg[0];

            FlatArray<int> sels = pnt2el[pi1];
            for (int j = 0; j < sels.Size(); j++)
              {
                SurfaceElementIndex sei = SurfaceElementIndex::FromNr0(sels[j]);
                Element2dRef se = (*this)[sei];
                bool found = false;
                for (int l = 0; l < se.GetNP(); l++ && !found)
                  found |= (se[l]==seg[1]);
                if(found) {
                  // seg.SetPartition(se.GetPartition());
                  seg_partition[segi] = surf_partition[sei];
                  break;
                }
              }
            
            // if (seg.GetPartition() == -1) {
            if (seg_partition[segi] == -1) {
              cout << endl << "segi: " << segi << endl;
              cout << "points: " << seg[0] << " " << seg[1] << endl;
              cout << "surfels: " << endl << sels << endl;
              throw NgException("no surface element found");
            }
          }
        
      }
  }

#endif





//========================== weights =================================================================



  // distribute the mesh to the worker processors
  // call it only for the master !
  void Mesh :: Distribute (Array<int> & volume_weights , Array<int>  & surface_weights, Array<int>  & segment_weights,
                           bool root_participates)
  {
    NgMPI_Comm comm = GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    if (id != 0 || ntasks == 1 ) return;

    ParallelMetis (volume_weights, surface_weights, segment_weights, root_participates);

    /*
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
    for (SurfaceElementIndex ei : SurfaceElements().Range())
      *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
      */
    
    // MyMPI_SendCmd ("mesh");
    SendRecvMesh (); 
  }
  

#ifdef METIS5
  void Mesh :: ParallelMetis (Array<int> & volume_weights , Array<int> & surface_weights, Array<int> & segment_weights,
                              bool root_participates)
  {
    PrintMessage (3, "call metis 5 with weights ...");
    
    // cout << "segment_weights " << segment_weights << endl;
    // cout << "surface_weights " << surface_weights << endl;
    // cout << "volume_weights " << volume_weights << endl;

    Timer timer("Mesh::Partition");
    RegionTimer reg(timer);

    idx_t ne = GetNE() + GetNSE() + GetNSeg();
    idx_t nn = GetNP();
    
    Array<idx_t> eptr, eind , nwgt;
    for (int i = 0; i < GetNE(); i++)
      {
        eptr.Append (eind.Size());
        
        auto el = (*this)[ElementIndex::FromNr1(i+1)];
        
        int ind = el.GetIndex().Nr1();        
        if (volume_weights.Size()<ind)
            nwgt.Append(0);
        else
            nwgt.Append (volume_weights[ind -1]);
        
        for (int j = 0; j < el.GetNP(); j++)
          eind.Append (el[j].Nr0());
      }
    for (int i = 0; i < GetNSE(); i++)
      {
        eptr.Append (eind.Size());
        const Element2dRef & el = (*this)[SurfaceElementIndex::FromNr1(i+1)];
        
        
        int ind = GetFaceDescriptor(el).BCProperty();
        if (surface_weights.Size()<ind)
            nwgt.Append(0);
        else
            nwgt.Append (surface_weights[ind -1]);

        
        for (int j = 0; j < el.GetNP(); j++)
          eind.Append (el[j].Nr0());
      }
    for (int i = 0; i < GetNSeg(); i++)
      {
        eptr.Append (eind.Size());
        
        const Segment & el = (*this)[SegmentIndex::FromNr1(i+1)];       
        
        int ind = el.GetIndex().Nr1();
        if (segment_weights.Size()<ind)
            nwgt.Append(0);
        else
            nwgt.Append (segment_weights[ind -1]);
        
        eind.Append (el[0].Nr0());
        eind.Append (el[1].Nr0());
      }
      
    eptr.Append (eind.Size());
    Array<idx_t> epart(ne), npart(nn);

    int first_rank = root_participates ? 0 : 1;
    idxtype nparts = GetCommunicator().Size() - first_rank;
    vol_partition.SetSize(GetNE());
    surf_partition.SetSize(GetNSE());
    seg_partition.SetSize(GetNSeg());
    
    if (nparts == 1)
      {
        for (int i = 0; i < GetNE(); i++)
          // VolumeElement(i+1).SetPartition(1);
          vol_partition[ElementIndex::FromNr0(i)] = first_rank;
        for (int i = 0; i < GetNSE(); i++)
          // SurfaceElement(i+1).SetPartition(1);
          surf_partition[SurfaceElementIndex::FromNr0(i)] = first_rank;
        for (int i = 0; i < GetNSeg(); i++)
          // LineSegment(i+1).SetPartition(1);
          seg_partition[SegmentIndex::FromNr0(i)] = first_rank;
        return;
      }

    
    idxtype edgecut;


    idxtype ncommon = 3;
    METIS_PartMeshDual (&ne, &nn, &eptr[0], &eind[0], &nwgt[0], NULL, &ncommon, &nparts,
                        NULL, NULL,
                        &edgecut, &epart[0], &npart[0]);
    /*
    METIS_PartMeshNodal (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &nparts,
                         NULL, NULL,
                         &edgecut, &epart[0], &npart[0]);
    */
    PrintMessage (3, "metis complete");
    // cout << "done" << endl;

    for (int i = 0; i < GetNE(); i++)
      // VolumeElement(i+1).SetPartition(epart[i] + 1);
      vol_partition[ElementIndex::FromNr0(i)] = epart[i] + first_rank;
    for (int i = 0; i < GetNSE(); i++)
      // SurfaceElement(i+1).SetPartition(epart[i+GetNE()] + 1);
      surf_partition[SurfaceElementIndex::FromNr0(i)] = epart[i+GetNE()] + first_rank;
    for (int i = 0; i < GetNSeg(); i++)
      // LineSegment(i+1).SetPartition(epart[i+GetNE()+GetNSE()] + 1);
      seg_partition[SegmentIndex::FromNr0(i)] = epart[i+GetNE()+GetNSE()] + first_rank;
  }
#endif

#ifndef METIS5
  void Mesh :: ParallelMetis (int /* nproc */)
  {
    throw NgException("Mesh::ParallelMetis: Netgen was built without METIS");
  }
  void Mesh :: ParallelMetis (Array<int> &, Array<int> &, Array<int> &)
  {
    throw NgException("Mesh::ParallelMetis: Netgen was built without METIS");
  }
#endif
 



//===========================================================================================









#ifdef METIS4
  void Mesh :: ParallelMetis ( )  
  {
    static Timer timer("Mesh::Partition");
    RegionTimer reg(timer);

    PrintMessage (3, "Metis called");
      
    if (GetDimension() == 2) 
      {
        PartDualHybridMesh2D ( ); // neloc );
        return;
      }


    idx_t ne = GetNE();
    idx_t nn = GetNP();

    if (ntasks <= 2 || ne <= 1)
      {
        if (ntasks == 1) return;
        
        for (int i=1; i<=ne; i++)
          VolumeElement(i).SetPartition(1);

        for (int i=1; i<=GetNSE(); i++)
          SurfaceElement(i).SetPartition(1);

        return;
      }


    bool uniform_els = true;

    ELEMENT_TYPE elementtype = TET; 
    for (int el = 1; el <= GetNE(); el++)
      if (VolumeElement(el).GetType() != elementtype)
        {
          uniform_els = false;
          break;
        }


    if (!uniform_els)
      {
        PartHybridMesh ();  
      }
    else
      {
        
        // uniform (TET) mesh,  JS
        int npe = VolumeElement(1).GetNP();
        Array<idxtype> elmnts(ne*npe);
        
        int etype;
        if (elementtype == TET)
          etype = 2;
        else if (elementtype == HEX)
          etype = 3;
        
    
        for (int i=1; i<=ne; i++)
          for (int j = 0; j < npe; j++)
            elmnts[(i-1)*npe+(j)] = VolumeElement(i)[j]-1;
        
        int numflag = 0;
        int nparts = ntasks-1;
        int ncommon = 3;
        int edgecut;
        Array<idxtype> epart(ne), npart(nn);
        
        //     if ( ntasks == 1 ) 
        //       {
        //      (*this) = *mastermesh;
        //      nparts = 4;        
        //      metis :: METIS_PartMeshDual (&ne, &nn, elmnts, &etype, &numflag, &nparts,
        //                                   &edgecut, epart, npart);
        //      cout << "done" << endl;
        
        //      cout << "edge-cut: " << edgecut << ", balance: " << metis :: ComputeElementBalance(ne, nparts, epart) << endl;
        
        //      for (int i=1; i<=ne; i++)
        //        {
        //          mastermesh->VolumeElement(i).SetPartition(epart[i-1]);
        //        }
        
        //      return;
        //       }
        
        
        static Timer timermetis("Metis itself");
        timermetis.Start();
        
#ifdef METIS4
        cout << "call metis(4)_PartMeshDual ... " << flush;
        METIS_PartMeshDual (&ne, &nn, &elmnts[0], &etype, &numflag, &nparts,
                            &edgecut, &epart[0], &npart[0]);
#else
        cout << "call metis(5)_PartMeshDual ... " << endl;
        // idx_t options[METIS_NOPTIONS];
        
        Array<idx_t> eptr(ne+1);
        for (int j = 0; j < ne+1; j++)
          eptr[j] = 4*j;
        
        METIS_PartMeshDual (&ne, &nn, &eptr[0], &elmnts[0], NULL, NULL, &ncommon, &nparts,
                            NULL, NULL,
                            &edgecut, &epart[0], &npart[0]);
#endif
        
        timermetis.Stop();
        
        cout << "complete" << endl;
#ifdef METIS4
        cout << "edge-cut: " << edgecut << ", balance: " 
             << ComputeElementBalance(ne, nparts, &epart[0]) << endl;
#endif
        
        // partition numbering by metis : 0 ...  ntasks - 1
        // we want:                       1 ...  ntasks
        for (int i=1; i<=ne; i++)
          VolumeElement(i).SetPartition(epart[i-1] + 1);
      }
    

    for (SurfaceElementIndex sei : SurfaceElements().Range())
      {
        ElementIndex ei1, ei2;
        GetTopology().GetSurface2VolumeElement (sei, ei1, ei2);
        Element2dRef sel = (*this)[sei];

        for (int j = 0; j < 2; j++)
          {
            ElementIndex ei = (j == 0) ? ei1 : ei2;
            if ( ei.IsValid() && ei.Nr0() < GetNE() )
              {
                sel.SetPartition ((*this)[ei].GetPartition());
                break;
              }
          }     
      }
    
  }
#endif


  void Mesh :: PartHybridMesh () 
  {
    throw Exception("PartHybridMesh not supported");    
#ifdef METISxxx
    int ne = GetNE();
    
    int nn = GetNP();
    int nedges = topology.GetNEdges();

    idxtype *xadj, * adjacency;
    // idxtype *v_weights = NULL, *e_weights = NULL;

    int weightflag = 0;
    int numflag = 0;
    int nparts = ntasks - 1;

    int options[5];
    options[0] = 0;
    int edgecut;
    idxtype * part;

    xadj = new idxtype[nn+1];
    part = new idxtype[nn];

    Array<int> cnt(nn+1);
    cnt = 0;

    for ( int edge = 0; edge < nedges; edge++ )
      {
        // int v1, v2;
        // topology.GetEdgeVertices ( edge, v1, v2);
        auto [v1,v2] = topology.GetEdgeVertices(edge);
        cnt[v1-1] ++;
        cnt[v2-1] ++;
      }

    xadj[0] = 0;
    for ( int n = 1; n <= nn; n++ )
      {
        xadj[n] = idxtype(xadj[n-1] + cnt[n-1]); 
      }

    adjacency = new idxtype[xadj[nn]];
    cnt = 0;

    for ( int edge = 0; edge < nedges; edge++ )
      {
        // int v1, v2;
        // topology.GetEdgeVertices ( edge, v1, v2);
        auto [v1,v2] = topology.GetEdgeVertices(edge);        
        adjacency[ xadj[v1-1] + cnt[v1-1] ] = v2-1;
        adjacency[ xadj[v2-1] + cnt[v2-1] ] = v1-1;
        cnt[v1-1]++;
        cnt[v2-1]++;
      }

    for ( int vert = 0; vert < nn; vert++ )
      {
        FlatArray<idxtype> array ( cnt[vert], &adjacency[ xadj[vert] ] );
        BubbleSort(array);
      }

#ifdef METIS4
    METIS_PartGraphKway ( &nn, xadj, adjacency, v_weights, e_weights, &weightflag, 
                          &numflag, &nparts, options, &edgecut, part );
#else
    cout << "currently not supported (metis5), A" << endl;
#endif

    Array<int> nodesinpart(ntasks);
    vol_partition.SetSize(ne);
    for ( int el = 1; el <= ne; el++ )
      {
        Element & volel = VolumeElement(el);
        nodesinpart = 0;

        
        int el_np = volel.GetNP();
        int partition = 0; 
        for ( int i = 0; i < el_np; i++ )
          nodesinpart[ part[volel[i]-1]+1 ] ++;

        for ( int i = 1; i < ntasks; i++ )
          if ( nodesinpart[i] > nodesinpart[partition] ) 
            partition = i;

        // volel.SetPartition(partition);
        vol_partition[el-1] = partition;
      }

    delete [] xadj;
    delete [] part;
    delete [] adjacency;
#else
    cout << "parthybridmesh not available" << endl;
#endif
  }


  void Mesh :: PartDualHybridMesh ( ) // Array<int> & neloc ) 
  {
    throw Exception("PartDualHybridMesh not supported");
#ifdef OLD      
#ifdef METIS
    int ne = GetNE();
    
    // int nn = GetNP();
    // int nedges = topology->GetNEdges();
    int nfaces = topology.GetNFaces();

    idxtype  *xadj, * adjacency, *v_weights = NULL, *e_weights = NULL;

    int weightflag = 0;
    // int numflag = 0;
    int nparts = ntasks - 1;

    int options[5];
    options[0] = 0;
    int edgecut;
    idxtype * part;

    Array<int> facevolels1(nfaces), facevolels2(nfaces);
    facevolels1 = -1;
    facevolels2 = -1;

    // Array<int, 0> elfaces;
    xadj = new idxtype[ne+1];
    part = new idxtype[ne];

    Array<int> cnt(ne+1);
    cnt = 0;

    for ( int el=1; el <= ne; el++ )
      {
        Element volel = VolumeElement(el);
        // topology.GetElementFaces(el, elfaces);
        auto elfaces = topology.GetFaces (ElementIndex(el-1));
        for ( int i = 0; i < elfaces.Size(); i++ )
          {
            if ( facevolels1[elfaces[i]] == -1 )
              facevolels1[elfaces[i]] = el;
            else
              {
                facevolels2[elfaces[i]] = el;
                cnt[facevolels1[elfaces[i]]-1]++;
                cnt[facevolels2[elfaces[i]]-1]++;
              }
          }
      }

    xadj[0] = 0;
    for ( int n = 1; n <= ne; n++ )
      {
        xadj[n] = idxtype(xadj[n-1] + cnt[n-1]); 
      }

    adjacency = new idxtype[xadj[ne]];
    cnt = 0;

    for (int face = 0; face < nfaces; face++)
      {
        int e1, e2;
        e1 = facevolels1[face];
        e2 = facevolels2[face];
        if ( e2 == -1 ) continue;
        adjacency[ xadj[e1-1] + cnt[e1-1] ] = e2-1;
        adjacency[ xadj[e2-1] + cnt[e2-1] ] = e1-1;
        cnt[e1-1]++;
        cnt[e2-1]++;
      }

    for ( int el = 0; el < ne; el++ )
      {
        FlatArray<idxtype> array ( cnt[el], &adjacency[ xadj[el] ] );
        BubbleSort(array);
      }

    Timer timermetis("Metis itself");
    timermetis.Start();

#ifdef METIS4
    METIS_PartGraphKway ( &ne, xadj, adjacency, v_weights, e_weights, &weightflag, 
                          &numflag, &nparts, options, &edgecut, part );
#else
    cout << "currently not supported (metis5), B" << endl;
#endif


    timermetis.Stop();

    Array<int> nodesinpart(ntasks);

    vol_partition.SetSize(ne);
    for ( int el = 1; el <= ne; el++ )
      {
        // Element & volel = VolumeElement(el);
        nodesinpart = 0;

        // VolumeElement(el).SetPartition(part[el-1 ] + 1);
        vol_partition[el-1] = part[el-1 ] + 1;
      }

    /*    
    for ( int i=1; i<=ne; i++)
      {
        neloc[ VolumeElement(i).GetPartition() ] ++;
      }
    */

    delete [] xadj;
    delete [] part;
    delete [] adjacency;
#else
    cout << "partdualmesh not available" << endl;
#endif
#endif
    
  }





  void Mesh :: PartDualHybridMesh2D ( ) 
  {
#ifdef METIS
    idxtype ne = GetNSE();
    int nv = GetNV();

    Array<idxtype> xadj(ne+1);
    Array<idxtype> adjacency(ne*4);

    // first, build the vertex 2 element table:
    Array<int, PointIndex> cnt(nv);
    cnt = 0;
    for (auto el : SurfaceElements())
      for (int j = 0; j < el.GetNP(); j++)
        cnt[ el[j] ] ++;
    
    DynamicTable<SurfaceElementIndex, PointIndex> vert2els(nv);
    for (SurfaceElementIndex sei : SurfaceElements().Range())
      for (int j = 0; j < (*this)[sei].GetNP(); j++)
        vert2els.Add ((*this)[sei][j], sei);
    

    // find all neighbour elements
    int cntnb = 0;
    Array<SurfaceElementIndex, SurfaceElementIndex> marks(ne);   // to visit each neighbour just once
    marks = SurfaceElementIndex::INVALID;
    for (SurfaceElementIndex sei : T_Range<SurfaceElementIndex>(ne))
      {
        xadj[sei.Nr0()] = cntnb;
        for (int j = 0; j < (*this)[sei].GetNP(); j++)
          {
            PointIndex vnr = (*this)[sei][j];

            // all elements with at least one common vertex
            for (int k = 0; k < vert2els[vnr].Size(); k++)   
              {
                SurfaceElementIndex sei2 = vert2els[vnr][k];
                if (sei == sei2) continue;
                if (marks[sei2] == sei) continue;
                
                // neighbour, if two common vertices
                int common = 0;
                for (int m1 = 0; m1 < (*this)[sei].GetNP(); m1++)
                  for (int m2 = 0; m2 < (*this)[sei2].GetNP(); m2++)
                    if ( (*this)[sei][m1] == (*this)[sei2][m2])
                      common++;
                
                if (common >= 2)
                  {
                    marks[sei2] = sei;     // mark as visited
                    adjacency[cntnb++] = sei2.Nr0();
                  }
              }
          }
      }
    xadj[ne] = cntnb;

    idxtype *v_weights = NULL, *e_weights = NULL;

    // int numflag = 0;
    idxtype nparts = ntasks - 1;

    idxtype edgecut;
    Array<idxtype> part(ne);

    for ( int el = 0; el < ne; el++ )
      BubbleSort (adjacency.Range (xadj[el], xadj[el+1]));

#ifdef METIS4   
    idxtype weightflag = 0;
    int options[5];
    options[0] = 0;
    METIS_PartGraphKway ( &ne, &xadj[0], &adjacency[0], v_weights, e_weights, &weightflag, 
                          &numflag, &nparts, options, &edgecut, &part[0] );
#else
    idx_t ncon = 1;
    METIS_PartGraphKway ( &ne, &ncon, &xadj[0], &adjacency[0], 
                          v_weights, NULL, e_weights, 
                          &nparts, 
                          NULL, NULL, NULL,
                          &edgecut, &part[0] );
#endif


    surf_partition.SetSize(ne);
    for (SurfaceElementIndex sei : T_Range<SurfaceElementIndex>(ne))
      // (*this) [sei].SetPartition (part[sei]+1);
      surf_partition[sei] = part[sei.Nr0()]+1;
#else
    cout << "partdualmesh not available" << endl;
#endif

  }



}

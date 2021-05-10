#ifdef PARALLEL

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

/*
namespace ngcore {
  template <> struct MPI_typetrait<netgen::PointIndex> {
    static MPI_Datatype MPIType () { return MPI_INT; } };  
}
*/

namespace ngcore
{

  /** An MPI-Package for a Surface element **/
  class SurfPointPackage
  {
  public:
    int num;     // point numebr
    int trignum; // STL geo info
    double u, v; // OCC geo info
    SurfPointPackage () { ; }
    SurfPointPackage & operator = (const SurfPointPackage & other) {
      num = other.num;
      trignum = other.trignum;
      u = other.u;
      v = other.v;
      return *this;
    }
  }; // class SurfPointPackage

  template<> struct MPI_typetrait<SurfPointPackage> {
    static MPI_Datatype MPIType () {
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  int block_len[2] = { 2, 2 };
	  MPI_Aint displs[3] = { 0, 2*sizeof(int) };
	  MPI_Datatype types[2] = { MPI_INT, MPI_DOUBLE };
	  MPI_Type_create_struct(2, block_len, displs, types, &MPI_T);
	  MPI_Type_commit(&MPI_T);
	}
      return MPI_T;
    }
  }; // struct MPI_typetrait<SurfPointPackage>


  class SelPackage
  {
  public:
    int sei;
    int index;
    int np;
    /** we send too much here, especially in 2d! **/
    SurfPointPackage points[ELEMENT2D_MAXPOINTS];
    SelPackage () { ; }
    SelPackage (const netgen::Mesh & mesh, netgen::SurfaceElementIndex _sei)
    {
      const netgen::Element2d & el = mesh[_sei];
      sei = _sei;
      index = el.GetIndex();
      np = el.GetNP();
      for (int k : Range(1, np+1)) {
	auto & pnt = points[k-1];;
	pnt.num = el.PNum(k);
	pnt.trignum = el.GeomInfoPi(k).trignum;
	pnt.u = el.GeomInfoPi(k).u;
	pnt.v = el.GeomInfoPi(k).v;
      }
      /** otherwise, we use uninitialized values **/
      for (int k : Range(np, ELEMENT2D_MAXPOINTS)) {
	points[k].num = -1;
	points[k].trignum = -1;
	points[k].u = -1;
	points[k].v = -1;
      }
    }
    void Unpack (netgen::Element2d & el) const {
      	el.SetIndex(index);
	for (int k : Range(1, np + 1)) {
	  auto & pnt = points[k-1];
	  el.PNum(k) = pnt.num;
	  el.GeomInfoPi(k).trignum = pnt.trignum;
	  el.GeomInfoPi(k).u = pnt.u;
	  el.GeomInfoPi(k).v = pnt.v;
	}
    }
    SelPackage & operator = (const SelPackage & other) {
      sei = other.sei;
      index = other.index;
      np = other.np;
      for (int k : Range(ELEMENT2D_MAXPOINTS))
	{ points[k] = other.points[k]; }
      return *this;
    }
  }; // class SelPackage

  template<> struct MPI_typetrait<SelPackage> {
    static MPI_Datatype MPIType () {
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  int block_len[2] = { 3, ELEMENT2D_MAXPOINTS };
	  MPI_Aint displs[3] = { 0, 3*sizeof(int) };
	  MPI_Datatype types[2] = { MPI_INT, GetMPIType<SurfPointPackage>() };
	  MPI_Type_create_struct(2, block_len, displs, types, &MPI_T);
	  MPI_Type_commit(&MPI_T);
	}
      return MPI_T;
    }
  }; // MPI_typetrait<SelPackage>


  class PointElPackage
  {
  public:
    netgen::PointIndex pnum;
    int index;
    PointElPackage () { pnum = -1; index = -1; }
    PointElPackage (const netgen::Element0d & el)
    { pnum = el.pnum; index = el.index; }
  }; // class PointElPackage

  template<> struct MPI_typetrait<PointElPackage> {
    static MPI_Datatype MPIType () {
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  int block_len[2] = { 1, 1 };
	  MPI_Aint displs[3] = { 0, sizeof(netgen::PointIndex) };
	  MPI_Datatype types[2] = { GetMPIType<netgen::PointIndex>(), MPI_INT };
	  MPI_Type_create_struct(2, block_len, displs, types, &MPI_T);
	  MPI_Type_commit(&MPI_T);
	}
      return MPI_T;
    }
  }; // MPI_typetrait<Element0d>


} // namespace ngcore

namespace netgen
{
  /*
  template <>
  inline MPI_Datatype MyGetMPIType<PointIndex> ( )
  { return MPI_INT; }
  */

  void Mesh :: SendRecvMesh ()
  {
    int id = GetCommunicator().Rank();
    int np = GetCommunicator().Size();

    if (np == 1) {
      throw NgException("SendRecvMesh called, but only one rank in communicator!!");
    }
    
    if (id == 0)
      PrintMessage (1, "Send/Receive mesh");

    // Why is this here??
    if (id == 0)
      {
	paralleltop -> SetNV (GetNV());
	paralleltop -> SetNV_Loc2Glob (GetNV());
	paralleltop -> SetNE (GetNE());
	paralleltop -> SetNSegm (GetNSeg());
	paralleltop -> SetNSE (GetNSE());
      }

    if (id == 0)
      SendMesh ();
    else
      ReceiveParallelMesh();

    paralleltop -> UpdateCoarseGrid();
  }







  void Mesh :: SendMesh () const   
  {

    NgMPI_Comm comm = GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    int dim = GetDimension();
    comm.Bcast(dim);

    Array<MPI_Request> sendrequests(8*(ntasks-1));
    sendrequests.SetSize0();
    
    // If the topology is not already updated, we do not need to
    // build edges/faces.
    auto & top = const_cast<MeshTopology&>(GetTopology());
    if(top.NeedsUpdate()) {
      top.SetBuildEdges(false);
      top.SetBuildFaces(false);
      top.Update();
    }
    
    PrintMessage ( 3, "Sending nr of elements");
    
    Array<int> num_els_on_proc(ntasks);
    num_els_on_proc = 0;
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      num_els_on_proc[vol_partition[ei]]++;

    comm.ScatterRoot (num_els_on_proc);

    Table<ElementIndex> els_of_proc (num_els_on_proc);
    num_els_on_proc = 0;
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      {
        auto nr = vol_partition[ei];
        els_of_proc[nr][num_els_on_proc[nr]++] = ei;
      }
    
    PrintMessage ( 3, "Building vertex/proc mapping");

    Array<int> num_sels_on_proc(ntasks);
    num_sels_on_proc = 0;
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      num_sels_on_proc[surf_partition[ei]]++;

    Table<SurfaceElementIndex> sels_of_proc (num_sels_on_proc);
    num_sels_on_proc = 0;
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      {
        auto nr = surf_partition[ei];
        sels_of_proc[nr][num_sels_on_proc[nr]++] = ei;
      }
    

    NgArray<int> num_segs_on_proc(ntasks);
    num_segs_on_proc = 0;
    for (SegmentIndex ei = 0; ei < GetNSeg(); ei++)
      // num_segs_on_proc[(*this)[ei].GetPartition()]++;
      num_segs_on_proc[seg_partition[ei]]++;

    TABLE<SegmentIndex> segs_of_proc (num_segs_on_proc);
    for (SegmentIndex ei = 0; ei < GetNSeg(); ei++)
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
    NgArray<INDEX_2> per_pairs;
    NgArray<INDEX_2> pp2;
    auto & idents = GetIdentifications();
    bool has_periodic = false; 
    for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
      {
	if(idents.GetType(idnr)!=Identifications::PERIODIC) continue;
	has_periodic = true;
	idents.GetPairs(idnr, pp2);
	per_pairs.Append(pp2);
      }
    NgArray<int, PointIndex::BASE> npvs(GetNV());
    npvs = 0;
    for (int k = 0; k < per_pairs.Size(); k++) {
      npvs[per_pairs[k].I1()]++;
      npvs[per_pairs[k].I2()]++;
    }

    /** for each vertex, gives us all identified vertices **/
    TABLE<PointIndex, PointIndex::BASE> per_verts(npvs);
    for (int k = 0; k < per_pairs.Size(); k++) {
      per_verts.Add(per_pairs[k].I1(), per_pairs[k].I2());
      per_verts.Add(per_pairs[k].I2(), per_pairs[k].I1());
    }
    for (int k = PointIndex::BASE; k < GetNV()+PointIndex::BASE; k++) {
      BubbleSort(per_verts[k]);
    }

    /** The same table as per_verts, but TRANSITIVE!! **/
    auto iterate_per_verts_trans = [&](auto f){
      NgArray<int> allvs;
      for (int k = PointIndex::BASE; k < GetNV()+PointIndex::BASE; k++)
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
    TABLE<PointIndex, PointIndex::BASE> per_verts_trans(npvs);
    iterate_per_verts_trans([&](auto k, auto & allvs) {
	for (int j = 0; j<allvs.Size(); j++)
	  per_verts_trans.Add(k, allvs[j]);
      });
    for (int k = PointIndex::BASE; k < GetNV()+PointIndex::BASE; k++) {
      BubbleSort(per_verts_trans[k]);
    }

    /** Now we build the vertex-data to send to the workers. **/
    NgArray<int, PointIndex::BASE> vert_flag (GetNV());
    NgArray<int, PointIndex::BASE> num_procs_on_vert (GetNV());
    NgArray<int> num_verts_on_proc (ntasks);
    num_verts_on_proc = 0;
    num_procs_on_vert = 0;
    auto iterate_vertices = [&](auto f) {
      vert_flag = -1;
      for (int dest = 1; dest < ntasks; dest++)
	{
          /*
	  FlatArray<ElementIndex> els = els_of_proc[dest];
	  for (int hi = 0; hi < els.Size(); hi++)
	    {
	      const Element & el = (*this) [ els[hi] ];
	      for (int i = 0; i < el.GetNP(); i++)
		f(el[i], dest);
	    }
          */
          for (auto & ei : els_of_proc[dest])
            for (auto pnum : (*this)[ei].PNums())
              f(pnum, dest);
          /*
	  FlatArray<SurfaceElementIndex> sels = sels_of_proc[dest];
	  for (int hi = 0; hi < sels.Size(); hi++)
	    {
	      const Element2d & el = (*this) [ sels[hi] ];
	      for (int i = 0; i < el.GetNP(); i++)
		f(el[i], dest);
	    }
          */
          for (auto & ei : sels_of_proc[dest])
            for (auto pnum : (*this)[ei].PNums())
              f(pnum, dest);

	  NgFlatArray<SegmentIndex> segs = segs_of_proc[dest];
	  for (int hi = 0; hi < segs.Size(); hi++)
	    {
	      const Segment & el = (*this) [segs[hi]];
	      for (int i = 0; i < 2; i++)
		f(el[i], dest);
	    }
	}
    };
    /** count vertices per proc and procs per vertex **/
    iterate_vertices([&](auto vertex, auto dest){
	auto countit = [&] (auto vertex, auto dest) {
	  if (vert_flag[vertex] < dest)
	    {
	      vert_flag[vertex] = dest;
	      num_verts_on_proc[dest]++;
	      num_procs_on_vert[vertex]++;
	      // GetParallelTopology().SetDistantPNum (dest, vertex);
              GetParallelTopology().AddDistantProc (PointIndex(vertex), dest); 
	    }
	};
	countit(vertex, dest);
	auto pers = per_verts_trans[vertex];
	for(int j = 0; j < pers.Size(); j++)
	  countit(pers[j], dest);
      });
    TABLE<PointIndex> verts_of_proc (num_verts_on_proc);
    TABLE<int, PointIndex::BASE> procs_of_vert (num_procs_on_vert);
    TABLE<int, PointIndex::BASE> loc_num_of_vert (num_procs_on_vert);
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
	auto pers = per_verts_trans[vertex];
	for(int j = 0; j < pers.Size(); j++)
	  addit(pers[j], dest);
      });
    /** 
	local vertex numbers on distant procs 
	(I think this was only used for debugging??) 
    **/
    for (int vert = 1; vert <= GetNP(); vert++ )
      {
	NgFlatArray<int> procs = procs_of_vert[vert];
	for (int j = 0; j < procs.Size(); j++)
	  {
	    int dest = procs[j];
	    // !! we also use this as offsets for MPI-type, if this is changed, also change ReceiveParallelMesh
	    verts_of_proc.Add (dest, vert - IndexBASE<T_POINTS::index_type>());
	    loc_num_of_vert.Add (vert, verts_of_proc[dest].Size());
	  }
      }
    PrintMessage ( 3, "Sending Vertices - vertices");

    Array<MPI_Datatype> point_types(ntasks-1);
    for (int dest = 1; dest < ntasks; dest++)
      {
	NgFlatArray<PointIndex> verts = verts_of_proc[dest];
	// sendrequests.Append (MyMPI_ISend (verts, dest, MPI_TAG_MESH+1, comm));
        sendrequests.Append (comm.ISend (FlatArray<PointIndex>(verts), dest, MPI_TAG_MESH+1));

	MPI_Datatype mptype = MeshPoint::MyGetMPIType();

	int numv = verts.Size();

	NgArray<int> blocklen (numv);  
	blocklen = 1;
	
	MPI_Type_indexed (numv, (numv == 0) ? nullptr : &blocklen[0], 
			  (numv == 0) ? nullptr : reinterpret_cast<int*> (&verts[0]), 
			  mptype, &point_types[dest-1]);
	MPI_Type_commit (&point_types[dest-1]);

	MPI_Request request;
	MPI_Isend( points.Data(), 1, point_types[dest-1], dest, MPI_TAG_MESH+1, comm, &request);
	sendrequests.Append (request);
      }


    /**
       Next, we send the identifications themselfs.
       
       Info about periodic identifications sent to each proc is an array of
       integers.
       - maxidentnr
       - type for each identification
       - nr of pairs for each identification (each pair is local!)
       - pairs for each periodic ident (global numbers)
    **/
    PrintMessage ( 3, "Sending Vertices - identifications");
    int maxidentnr = idents.GetMaxNr();
    NgArray<int> ppd_sizes(ntasks);
    ppd_sizes = 1 + 2*maxidentnr;
    for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
      {
	if(idents.GetType(idnr)!=Identifications::PERIODIC) continue;
	idents.GetPairs(idnr, pp2);
	for(int j = 0; j<pp2.Size(); j++)
	  {
	    INDEX_2 & pair = pp2[j];
	    // both are on same procs!
	    auto ps = procs_of_vert[pair.I1()];
	    for (int l = 0; l < ps.Size(); l++)
	      {
		ppd_sizes[ps[l]] += 2;
	      }
	  }
      }
    TABLE<int> pp_data(ppd_sizes);
    for(int dest = 0; dest < ntasks; dest++)
      pp_data.Add(dest, maxidentnr);
    for (int dest = 0; dest < ntasks; dest++)
      {
	for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
	  pp_data.Add(dest, idents.GetType(idnr));
	for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
	  pp_data.Add(dest, 0);
      }
    for (int idnr = 1; idnr < idents.GetMaxNr()+1; idnr++)
      {
	if(idents.GetType(idnr)!=Identifications::PERIODIC) continue;
	idents.GetPairs(idnr, pp2);
	for(int j = 0; j<pp2.Size(); j++)
	  {
	    INDEX_2 & pair = pp2[j];
	    auto ps = procs_of_vert[pair.I1()];
	    for (int l = 0; l < ps.Size(); l++)
	      {
		auto p = ps[l];
		pp_data[p][maxidentnr + idnr]++;
		pp_data.Add(p, pair.I1());
		pp_data.Add(p, pair.I2());
	      }
	  }
      }
    Array<MPI_Request> req_per;
    for(int dest = 1; dest < ntasks; dest++)
      // req_per.Append(MyMPI_ISend(pp_data[dest], dest, MPI_TAG_MESH+1, comm));
      req_per.Append(comm.ISend(FlatArray<int>(pp_data[dest]), dest, MPI_TAG_MESH+1));
    MyMPI_WaitAll(req_per);

    PrintMessage ( 3, "Sending Vertices - distprocs");

    Array<int> num_distpnums(ntasks);
    num_distpnums = 0;
    
    for (int vert = 1; vert <= GetNP(); vert++)
      {
	FlatArray<int> procs = procs_of_vert[vert];
	for (auto p : procs)
	  num_distpnums[p] += 3 * (procs.Size()-1);
      }

    DynamicTable<int> distpnums (num_distpnums);

    for (int vert = 1; vert <= GetNP(); vert++)
      {
	NgFlatArray<int> procs = procs_of_vert[vert];
	for (int j = 0; j < procs.Size(); j++)
	  for (int k = 0; k < procs.Size(); k++)
	    if (j != k)
	      {
		distpnums.Add (procs[j], loc_num_of_vert[vert][j]);
		distpnums.Add (procs[j], procs_of_vert[vert][k]);
		distpnums.Add (procs[j], loc_num_of_vert[vert][k]);
	      }
      }
    
    for ( int dest = 1; dest < ntasks; dest ++ )
      sendrequests.Append (comm.ISend (distpnums[dest], dest, MPI_TAG_MESH+1));



    PrintMessage ( 3, "Sending elements" );

    Array<int> elarraysize (ntasks);
    elarraysize = 0;
    for ( int ei = 1; ei <= GetNE(); ei++)
      {
	const Element & el = VolumeElement (ei);
	// int dest = el.GetPartition();
        int dest = vol_partition[ei-1];
	elarraysize[dest] += 3 + el.GetNP();
      }

    DynamicTable<int> elementarrays(elarraysize);

    for (int ei = 1; ei <= GetNE(); ei++)
      {
	const Element & el = VolumeElement (ei);
	// int dest = el.GetPartition();
        int dest = vol_partition[ei-1];
        
	elementarrays.Add (dest, ei);
	elementarrays.Add (dest, el.GetIndex());
	elementarrays.Add (dest, el.GetNP());
	for (int i = 0; i < el.GetNP(); i++)
	  elementarrays.Add (dest, el[i]);
      }

    for (int dest = 1; dest < ntasks; dest ++ )
      // sendrequests.Append (MyMPI_ISend (elementarrays[dest], dest, MPI_TAG_MESH+2, comm));
      sendrequests.Append (comm.ISend (elementarrays[dest], dest, MPI_TAG_MESH+2));


    PrintMessage ( 3, "Sending Face Descriptors" );

    Array<double> fddata (6 * GetNFD());
    for (int fdi = 1; fdi <= GetNFD(); fdi++)
      {
	fddata[6*fdi-6] = GetFaceDescriptor(fdi).SurfNr();
	fddata[6*fdi-5] = GetFaceDescriptor(fdi).DomainIn();	
	fddata[6*fdi-4] = GetFaceDescriptor(fdi).DomainOut();
	fddata[6*fdi-3] = GetFaceDescriptor(fdi).BCProperty();
	fddata[6*fdi-2] = GetFaceDescriptor(fdi).domin_singular;
	fddata[6*fdi-1] = GetFaceDescriptor(fdi).domout_singular;
	
      }
    for (int dest = 1; dest < ntasks; dest++)
      sendrequests.Append (comm.ISend (fddata, dest, MPI_TAG_MESH+3));
    
    /** Surface Elements **/

    PrintMessage ( 3, "Sending Surface elements" );
    // build sel-identification
    size_t nse = GetNSE();
    NgArray<SurfaceElementIndex> ided_sel(nse);
    ided_sel = -1;
    bool has_ided_sels = false;
    if(GetNE() && has_periodic) //we can only have identified surf-els if we have vol-els (right?)
      {
	Array<SurfaceElementIndex> os1, os2;
	for(SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
	    if(ided_sel[sei]!=-1) continue;
	    const Element2d & sel = (*this)[sei];
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
	    const Element2d & sel2 = (*this)[sei];
	    auto points2 = sel2.PNums();
	    has_ided_sels = true;
	    ided_sel[sei] = os1[0];
	    ided_sel[os1[0]] = sei;
	  }
      }
    // build sel data to send
    auto iterate_sels = [&](auto f) {
      for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++ )
	{
	  const Element2d & sel = (*this)[sei];
	  // int dest = (*this)[sei].GetPartition();
          int dest = surf_partition[sei];
	  f(sei, sel, dest);
	  if(ided_sel[sei]!=-1)
	    {
	      // int dest2 = (*this)[ided_sel[sei]].GetPartition();
              int dest2 = surf_partition[ided_sel[sei]];
	      f(sei, sel, dest2);
	    }
	}      
    };
    Array <int> nlocsel(ntasks), bufsize(ntasks);
    nlocsel = 0;
    bufsize = 0;
    iterate_sels([&](SurfaceElementIndex sei, const Element2d & sel, int dest){
	nlocsel[dest]++;
	bufsize[dest]++;
      });
    DynamicTable<SelPackage> selbuf(bufsize);
    iterate_sels([&](SurfaceElementIndex sei, const auto & sel, int dest) {
	selbuf.Add (dest, SelPackage(*this, sei));
      });
    // distribute sel data
    for (int dest = 1; dest < ntasks; dest++)
      sendrequests.Append (comm.ISend(selbuf[dest], dest, MPI_TAG_MESH+4));
    

    /** Segments **/
    PrintMessage ( 3, "Sending Edge Segments");
    auto iterate_segs1 = [&](auto f) {
      NgArray<SegmentIndex> osegs1, osegs2, osegs_both;
      NgArray<int> type1, type2;
      for(SegmentIndex segi = 0; segi < GetNSeg(); segi++)
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
    NgArray<int> per_seg_size(GetNSeg());
    per_seg_size = 0;
    iterate_segs1([&](SegmentIndex segi1, SegmentIndex segi2)
		  { per_seg_size[segi1]++; });
    TABLE<SegmentIndex> per_seg(per_seg_size);
    iterate_segs1([&](SegmentIndex segi1, SegmentIndex segi2)
		  { per_seg.Add(segi1, segi2); });
    // make per_seg transitive
    auto iterate_per_seg_trans = [&](auto f){
      NgArray<SegmentIndex> allsegs;
      for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
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
    iterate_per_seg_trans([&](SegmentIndex segi, NgArray<SegmentIndex> & segs){
	for (int j = 0; j < segs.Size(); j++)
	  per_seg_size[segi] = segs.Size();
      });
    TABLE<SegmentIndex> per_seg_trans(per_seg_size);
    iterate_per_seg_trans([&](SegmentIndex segi, NgArray<SegmentIndex> & segs){
	for (int j = 0; j < segs.Size(); j++)
	  per_seg_trans.Add(segi, segs[j]);
      });
    // build segment data
    NgArray<int> dests;
    auto iterate_segs2 = [&](auto f)
      {
	for (SegmentIndex segi = 0; segi<GetNSeg(); segi++)
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
    NgArray<int> nloc_seg(ntasks);
    // bufsize = 1; //was originally this - why??
    bufsize = 0;
    nloc_seg = 0;
    iterate_segs2([&](auto segi, const auto & seg, int dest)
		  {
		    nloc_seg[dest]++;
		    bufsize[dest] += 14;
		  });
    DynamicTable<double> segm_buf(bufsize);
    iterate_segs2([&](auto segi, const auto & seg, int dest)
		  {
		    segm_buf.Add (dest, segi);
		    segm_buf.Add (dest, seg.si);
		    segm_buf.Add (dest, seg.pnums[0]);
		    segm_buf.Add (dest, seg.pnums[1]);
		    segm_buf.Add (dest, seg.geominfo[0].trignum);
		    segm_buf.Add (dest, seg.geominfo[1].trignum);
		    segm_buf.Add (dest, seg.surfnr1);
		    segm_buf.Add (dest, seg.surfnr2);
		    segm_buf.Add (dest, seg.edgenr);
		    segm_buf.Add (dest, seg.epgeominfo[0].dist);
		    segm_buf.Add (dest, seg.epgeominfo[1].edgenr);
		    segm_buf.Add (dest, seg.epgeominfo[1].dist);
		    segm_buf.Add (dest, seg.singedge_right);
		    segm_buf.Add (dest, seg.singedge_left);
		  });
    // distrubute segment data
    for (int dest = 1; dest < ntasks; dest++)
      sendrequests.Append (comm.ISend(segm_buf[dest], dest, MPI_TAG_MESH+5));

    /** Point-Elements **/
    PrintMessage ( 3, "Point-Elements ...");

    auto iterate_zdes = [&](auto f) {
      for (auto k : Range(pointelements)) {
	auto & el = pointelements[k];
	PointElPackage pack(el);
	auto dests = procs_of_vert[el.pnum];
	for (auto dest : dests)
	  { f(pack, dest); }
      }
    };

    bufsize = 0;
    iterate_zdes([&](const auto & pack, auto dest) { bufsize[dest]++; });
    DynamicTable<PointElPackage> zde_buf(bufsize); // zero dim elements
    iterate_zdes([&](const auto & pack, auto dest) { zde_buf.Add(dest, pack); });

    for (int dest = 1; dest < ntasks; dest++)
      { sendrequests.Append (comm.ISend(zde_buf[dest], dest, MPI_TAG_MESH+6)); }

    PrintMessage ( 3, "now wait ...");

    MyMPI_WaitAll (sendrequests);

    // clean up MPI-datatypes we allocated earlier
    for (auto t : point_types)
      { MPI_Type_free(&t); }

    paralleltop -> SetNV_Loc2Glob (0);
    paralleltop -> SetNV (0);
    paralleltop -> EnumeratePointsGlobally();
    PrintMessage ( 3, "Sending names");

    sendrequests.SetSize(3*ntasks);
    /** Send bc/mat/cd*-names **/
    // nr of names
    ArrayMem<int,4> nnames{0,0,0,0};
    nnames[0] = materials.Size();
    nnames[1] = bcnames.Size();
    nnames[2] = GetNCD2Names();
    nnames[3] = GetNCD3Names();
    int tot_nn = nnames[0] + nnames[1] + nnames[2] + nnames[3];
    for( int k = 1; k < ntasks; k++)
      sendrequests[k] = comm.ISend(nnames, k, MPI_TAG_MESH+7);
      // (void) MPI_Isend(nnames, 4, MPI_INT, k, MPI_TAG_MESH+6, comm, &sendrequests[k]);
    auto iterate_names = [&](auto func) {
      for (int k = 0; k < nnames[0]; k++) func(materials[k]);
      for (int k = 0; k < nnames[1]; k++) func(bcnames[k]);
      for (int k = 0; k < nnames[2]; k++) func(cd2names[k]);
      for (int k = 0; k < nnames[3]; k++) func(cd3names[k]);
    };
    // sizes of names
    NgArray<int> name_sizes(tot_nn);
    tot_nn = 0;
    iterate_names([&](auto ptr) { name_sizes[tot_nn++] = (ptr==NULL) ? 0 : ptr->size(); });
    for( int k = 1; k < ntasks; k++)
      (void) MPI_Isend(&name_sizes[0], tot_nn, MPI_INT, k, MPI_TAG_MESH+7, comm, &sendrequests[ntasks+k]);
    // names
    int strs = 0;
    iterate_names([&](auto ptr) { strs += (ptr==NULL) ? 0 : ptr->size(); });
    NgArray<char> compiled_names(strs);
    strs = 0;
    iterate_names([&](auto ptr) {
	if (ptr==NULL) return;
	auto& name = *ptr;
	for (int j=0; j < name.size(); j++) compiled_names[strs++] = name[j];
      });
    for( int k = 1; k < ntasks; k++)
      (void) MPI_Isend(&(compiled_names[0]), strs, MPI_CHAR, k, MPI_TAG_MESH+7, comm, &sendrequests[2*ntasks+k]);

    PrintMessage ( 3, "wait for names");

    MyMPI_WaitAll (sendrequests);
    
    comm.Barrier();

    PrintMessage( 3, "Clean up local memory");

    auto & self = const_cast<Mesh&>(*this);
    self.points = T_POINTS(0);
    self.surfelements = Array<Element2d>(0);
    self.volelements = Array<Element>(0);
    self.segments = Array<Segment>(0);
    self.pointelements = Array<Element0d>(0);
    self.lockedpoints = Array<PointIndex>(0);
    auto cleanup_ptr = [](auto & ptr) {
      if (ptr != nullptr) {
	delete ptr;
	ptr = nullptr;
      }
    };
    /*
    cleanup_ptr(self.boundaryedges);
    cleanup_ptr(self.segmentht);
    cleanup_ptr(self.surfelementht);
    */
    self.boundaryedges = nullptr;
    self.segmentht = nullptr;
    self.surfelementht = nullptr;
    
    self.openelements = NgArray<Element2d>(0);
    self.opensegments = NgArray<Segment>(0);
    self.numvertices = 0;
    self.mlbetweennodes = NgArray<PointIndices<2>,PointIndex::BASE> (0);
    self.mlparentelement = NgArray<int>(0);
    self.mlparentsurfaceelement = NgArray<int>(0);
    self.curvedelems = make_unique<CurvedElements> (self);
    self.clusters = make_unique<AnisotropicClusters> (self);
    self.ident = make_unique<Identifications> (self);
    self.topology = MeshTopology(*this);
    self.topology.Update();
    self.BuildElementSearchTree();
    
    // const_cast<Mesh&>(*this).DeleteMesh();

    // paralleltop -> SetNV (0);
    // paralleltop->EnumeratePointsGlobally();
    
    PrintMessage( 3, "send mesh complete");
  }








  // workers receive the mesh from the master
  void Mesh :: ReceiveParallelMesh ( )
  {
    int timer = NgProfiler::CreateTimer ("ReceiveParallelMesh");
    int timer_pts = NgProfiler::CreateTimer ("Receive points");
    int timer_els = NgProfiler::CreateTimer ("Receive elements");
    int timer_sels = NgProfiler::CreateTimer ("Receive surface elements");
    NgProfiler::RegionTimer reg(timer);

    NgMPI_Comm comm = GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();
    
    int dim;
    comm.Bcast(dim);
    SetDimension(dim);
    
    // Receive number of local elements
    int nelloc;
    comm.Scatter (nelloc);
    paralleltop -> SetNE (nelloc);
    
    // receive vertices
    NgProfiler::StartTimer (timer_pts);

    Array<int> verts;
    comm.Recv (verts, 0, MPI_TAG_MESH+1);

    int numvert = verts.Size();
    paralleltop -> SetNV (numvert);
    paralleltop -> SetNV_Loc2Glob (numvert);
    
    // INDEX_CLOSED_HASHTABLE<int> glob2loc_vert_ht (3*numvert+1);
    INDEX_HASHTABLE<int> glob2loc_vert_ht (3*numvert+1);

    for (int vert = 0; vert < numvert; vert++)
      {
	int globvert = verts[vert] + IndexBASE<T_POINTS::index_type>();
        // paralleltop->SetLoc2Glob_Vert ( vert+1, globvert  );
        paralleltop->L2G (PointIndex(vert+PointIndex::BASE)) = globvert;
	glob2loc_vert_ht.Set (globvert, vert+1);
      }
    
    for (int i = 0; i < numvert; i++)
      AddPoint (netgen::Point<3> (0,0,0));
    
    MPI_Datatype mptype = MeshPoint::MyGetMPIType();
    MPI_Status status;
    MPI_Recv( points.Data(), numvert, mptype, 0, MPI_TAG_MESH+1, comm, &status);

    Array<int> pp_data;
    comm.Recv(pp_data, 0, MPI_TAG_MESH+1);

    int maxidentnr = pp_data[0];
    auto & idents = GetIdentifications();
    for (int idnr = 1; idnr < maxidentnr+1; idnr++)
      idents.SetType(idnr, (Identifications::ID_TYPE)pp_data[idnr]);

    int offset = 2*maxidentnr+1;
    for(int idnr = 1; idnr < maxidentnr+1; idnr++)
      {
    	int npairs = pp_data[maxidentnr+idnr];
    	NgFlatArray<int> pairdata(2*npairs, &pp_data[offset]);
    	offset += 2*npairs;
	for (int k = 0; k<npairs; k++) {
	  PointIndex loc1 = glob2loc_vert_ht.Get(pairdata[2*k]);
	  PointIndex loc2 = glob2loc_vert_ht.Get(pairdata[2*k+1]);
	  idents.Add(loc1, loc2, idnr);
	}
      }
    
    Array<int> dist_pnums; 
    comm.Recv (dist_pnums, 0, MPI_TAG_MESH+1);
    
    for (int hi = 0; hi < dist_pnums.Size(); hi += 3)
      paralleltop ->
	// SetDistantPNum (dist_pnums[hi+1], dist_pnums[hi]); // , dist_pnums[hi+2]);
        AddDistantProc (PointIndex(dist_pnums[hi]), dist_pnums[hi+1]);
    
    NgProfiler::StopTimer (timer_pts);
    *testout << "got " << numvert << " vertices" << endl;

    
    {
      Array<int> elarray;
      comm.Recv (elarray, 0, MPI_TAG_MESH+2);
      
      NgProfiler::RegionTimer reg(timer_els);

      for (int ind = 0, elnum = 1; ind < elarray.Size(); elnum++)
	{
	  paralleltop->SetLoc2Glob_VolEl ( elnum,  elarray[ind++]);

          int index = elarray[ind++];
          Element el(elarray[ind++]);          
	  el.SetIndex(index);
	  
	  for ( int j = 0; j < el.GetNP(); j++)
	    el[j] = glob2loc_vert_ht.Get (elarray[ind++]); 
	  
	  AddVolumeElement (el);
	}
    }

    {
      Array<double> fddata;
      comm.Recv (fddata, 0, MPI_TAG_MESH+3);
      for (int i = 0; i < fddata.Size(); i += 6)
	{
	  int faceind = AddFaceDescriptor 
	    (FaceDescriptor(int(fddata[i]), int(fddata[i+1]), int(fddata[i+2]), 0));
	  GetFaceDescriptor(faceind).SetBCProperty (int(fddata[i+3]));
	  GetFaceDescriptor(faceind).domin_singular = fddata[i+4];
	  GetFaceDescriptor(faceind).domout_singular = fddata[i+5];
	}
    }

    {
      NgProfiler::RegionTimer reg(timer_sels);
      Array<SelPackage> selbuf;

      comm.Recv ( selbuf, 0, MPI_TAG_MESH+4);
      
      int nlocsel = selbuf.Size();
      paralleltop -> SetNSE ( nlocsel );
      
      int sel = 0;
      for (auto k : Range(selbuf)) {
	auto & pack = selbuf[k];
	Element2d el(pack.np);
	pack.Unpack(el);
	/** map global point numbers to local ones **/
	for (int k : Range(1, 1+el.GetNP()))
	  { el.PNum(k) = glob2loc_vert_ht.Get(el.PNum(k)); }
	paralleltop->SetLoc2Glob_SurfEl (sel+1, pack.sei);
	AddSurfaceElement (el);
	sel++;
      }
    }
    


    {
      NgArray<double> segmbuf;
      MyMPI_Recv ( segmbuf, 0, MPI_TAG_MESH+5, comm);

      Segment seg;
      int globsegi;
      int ii = 0;
      int segi = 1;
      int nsegloc = int ( segmbuf.Size() / 14 ) ;
      paralleltop -> SetNSegm ( nsegloc );

      while ( ii < segmbuf.Size() )
	{
	  globsegi = int (segmbuf[ii++]);
	  seg.si = int (segmbuf[ii++]);
	  
	  seg.pnums[0] = glob2loc_vert_ht.Get (int(segmbuf[ii++]));
	  seg.pnums[1] = glob2loc_vert_ht.Get (int(segmbuf[ii++]));
	  seg.geominfo[0].trignum = int( segmbuf[ii++] );
	  seg.geominfo[1].trignum = int ( segmbuf[ii++]);
	  seg.surfnr1 = int ( segmbuf[ii++]);
	  seg.surfnr2 = int ( segmbuf[ii++]);
	  seg.edgenr = int ( segmbuf[ii++]);
	  seg.epgeominfo[0].dist = segmbuf[ii++];
	  seg.epgeominfo[1].edgenr = int (segmbuf[ii++]);
	  seg.epgeominfo[1].dist = segmbuf[ii++];
	  
	  seg.singedge_left = segmbuf[ii++];
	  seg.singedge_right = segmbuf[ii++];
	  
	  seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;
	  
	  seg.domin = seg.surfnr1;
	  seg.domout = seg.surfnr2;
	  if ( seg.pnums[0] >0 && seg.pnums[1] > 0 )
	    {
	      paralleltop-> SetLoc2Glob_Segm ( segi,  globsegi );
	      
	      AddSegment (seg);
	      segi++;
	    }
	}
    }

    { /** 0d-Elements **/
      Array<PointElPackage> zdes;
      comm.Recv ( zdes, 0, MPI_TAG_MESH+6);
      pointelements.SetSize(zdes.Size());
      for (auto k : Range(pointelements)) {
	auto & el = pointelements[k];
	el.pnum = glob2loc_vert_ht.Get(zdes[k].pnum);
	el.index = zdes[k].index;
      }
    }

    // paralleltop -> SetNV_Loc2Glob (0);
    paralleltop -> EnumeratePointsGlobally();
    /** Recv bc-names **/
    ArrayMem<int,4> nnames{0,0,0,0};
    // MPI_Recv(nnames, 4, MPI_INT, 0, MPI_TAG_MESH+6, comm, MPI_STATUS_IGNORE);
    comm.Recv(nnames, 0, MPI_TAG_MESH+7);
    // cout << "nnames = " << FlatArray(nnames) << endl;
    materials.SetSize(nnames[0]);
    bcnames.SetSize(nnames[1]);
    cd2names.SetSize(nnames[2]);
    cd3names.SetSize(nnames[3]);

    int tot_nn = nnames[0] + nnames[1] + nnames[2] + nnames[3];
    NgArray<int> name_sizes(tot_nn);
    MPI_Recv(&name_sizes[0], tot_nn, MPI_INT, 0, MPI_TAG_MESH+7, comm, MPI_STATUS_IGNORE);
    int tot_size = 0;
    for (int k = 0; k < tot_nn; k++) tot_size += name_sizes[k];
    
    NgArray<char> compiled_names(tot_size);
    MPI_Recv(&(compiled_names[0]), tot_size, MPI_CHAR, 0, MPI_TAG_MESH+7, comm, MPI_STATUS_IGNORE);

    tot_nn = tot_size = 0;
    auto write_names = [&] (auto & array) {
      for (int k = 0; k < array.Size(); k++) {
	int s = name_sizes[tot_nn];
	array[k] = new string(&compiled_names[tot_size], s);
	tot_nn++;
	tot_size += s;
      }
    };
    write_names(materials);
    write_names(bcnames);
    write_names(cd2names);
    write_names(cd3names);
    
    comm.Barrier();

    int timerloc = NgProfiler::CreateTimer ("Update local mesh");
    int timerloc2 = NgProfiler::CreateTimer ("CalcSurfacesOfNode");

    NgProfiler::RegionTimer regloc(timerloc);
    stringstream str;
    str << "p" << id << ": got " << GetNE() << " elements and " 
	 << GetNSE() << " surface elements";
    PrintMessage(2, str.str());
    // cout << str.str() << endl;
    // PrintMessage (2, "Got ", GetNE(), " elements and ", GetNSE(), " surface elements");
    // PrintMessage (2, "Got ", GetNSE(), " surface elements");

    NgProfiler::StartTimer (timerloc2);

    CalcSurfacesOfNode ();

    NgProfiler::StopTimer (timerloc2);

    topology.Update();
    clusters -> Update();

    // paralleltop -> UpdateCoarseGrid();
    // paralleltop->EnumeratePointsGlobally();
    SetNextMajorTimeStamp();
  }
  


  
  
  // distribute the mesh to the worker processors
  // call it only for the master !
  void Mesh :: Distribute ()
  {
    NgMPI_Comm comm = GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    if (id != 0 || ntasks == 1 ) return;

#ifdef METIS
    ParallelMetis ();
#else
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      (*this)[ei].SetPartition(ntasks * ei/GetNE() + 1);
#endif

    /*
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
      */
    
    // MyMPI_SendCmd ("mesh");
    SendRecvMesh (); 
  }
  

#ifdef METIS5
  void Mesh :: ParallelMetis ( )  
  {
    PrintMessage (3, "call metis 5 ...");

    int timer = NgProfiler::CreateTimer ("Mesh::Partition");
    NgProfiler::RegionTimer reg(timer);

    idx_t ne = GetNE() + GetNSE() + GetNSeg();
    idx_t nn = GetNP();

    NgArray<idx_t> eptr, eind;
    for (int i = 0; i < GetNE(); i++)
      {
	eptr.Append (eind.Size());
	const Element & el = VolumeElement(i+1);
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSE(); i++)
      {
	eptr.Append (eind.Size());
	const Element2d & el = SurfaceElement(i+1);
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSeg(); i++)
      {
	eptr.Append (eind.Size());
	const Segment & el = LineSegment(i+1);
	eind.Append (el[0]-1);
	eind.Append (el[1]-1);
      }
    eptr.Append (eind.Size());
    NgArray<idx_t> epart(ne), npart(nn);

    idxtype nparts = GetCommunicator().Size()-1;

    vol_partition.SetSize(GetNE());
    surf_partition.SetSize(GetNSE());
    seg_partition.SetSize(GetNSeg());
    if (nparts == 1)
      {
        for (int i = 0; i < GetNE(); i++)
          // VolumeElement(i+1).SetPartition(1);
          vol_partition[i]= 1;
        for (int i = 0; i < GetNSE(); i++)
          // SurfaceElement(i+1).SetPartition(1);
          surf_partition[i] = 1;
        for (int i = 0; i < GetNSeg(); i++)
          // LineSegment(i+1).SetPartition(1);
          seg_partition[i] = 1;
      }

    else
      
      {

        idxtype edgecut;
        
        idxtype ncommon = 3;
        METIS_PartMeshDual (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &ncommon, &nparts,
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
          vol_partition[i]= epart[i] + 1;
        for (int i = 0; i < GetNSE(); i++)
          // SurfaceElement(i+1).SetPartition(epart[i+GetNE()] + 1);
          surf_partition[i] = epart[i+GetNE()] + 1;
        for (int i = 0; i < GetNSeg(); i++)
          // LineSegment(i+1).SetPartition(epart[i+GetNE()+GetNSE()] + 1);
          seg_partition[i] = epart[i+GetNE()+GetNSE()] + 1;
      }
    
        
    // surface elements attached to volume elements
    NgArray<bool, PointIndex::BASE> boundarypoints (GetNP());
    boundarypoints = false;

    if(GetDimension() == 3)
      for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	{
	  const Element2d & el = (*this)[sei];
	  for (int j = 0; j < el.GetNP(); j++)
	    boundarypoints[el[j]] = true;
	}
    else
      for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
	{
	  const Segment & seg = (*this)[segi];
	  for (int j = 0; j < 2; j++)
	    boundarypoints[seg[j]] = true;
	}

    
    // Build Pnt2Element table, boundary points only
    NgArray<int, PointIndex::BASE> cnt(GetNP());
    cnt = 0;

    auto loop_els_2d = [&](auto f) {
      for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	{
	  const Element2d & el = (*this)[sei];
	  for (int j = 0; j < el.GetNP(); j++) {
	    f(el[j], sei);
	  }
	}
    };
    auto loop_els_3d = [&](auto f) {
      for (ElementIndex ei = 0; ei < GetNE(); ei++)
	{
	  const Element & el = (*this)[ei];
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

    
    loop_els([&](auto vertex, int index)
	{
	  if(boundarypoints[vertex])
	    cnt[vertex]++;
	});
    TABLE<int, PointIndex::BASE> pnt2el(cnt);
    loop_els([&](auto vertex, int index)
	{
	  if(boundarypoints[vertex])
	    pnt2el.Add(vertex, index);
	});


    if (GetDimension() == 3)
      {
	for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
	    Element2d & sel = (*this)[sei];
	    PointIndex pi1 = sel[0];
	    // NgFlatArray<ElementIndex> els = pnt2el[pi1];
	    NgFlatArray<int> els = pnt2el[pi1];
	    
	    // sel.SetPartition (-1);
            surf_partition[sei] = -1;
	    
	    for (int j = 0; j < els.Size(); j++)
	      {
		const Element & el = (*this)[ElementIndex(els[j])];
		
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
                    surf_partition[sei] = vol_partition[ElementIndex(els[j])];
		    break;
		  }
	      }
	    // if (sel.GetPartition() == -1)
            if (surf_partition[sei] == -1)
	      cerr << "no volume element found" << endl;
	  }


	for (SegmentIndex si = 0; si < GetNSeg(); si++)
	  {
	    Segment & sel = (*this)[si];
	    PointIndex pi1 = sel[0];
	    NgFlatArray<int> els = pnt2el[pi1];
	    
	    // sel.SetPartition (-1);
            seg_partition[si] = -1;
	    
	    for (int j = 0; j < els.Size(); j++)
	      {
		const Element & el = (*this)[ElementIndex(els[j])];
		
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
                    seg_partition[si] = vol_partition[ElementIndex(els[j])];
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
	for (SegmentIndex segi = 0; segi < GetNSeg(); segi++)
	  {
	    Segment & seg = (*this)[segi];
	    // seg.SetPartition(-1);
            seg_partition[segi] = -1;
	    PointIndex pi1 = seg[0];

	    NgFlatArray<int> sels = pnt2el[pi1];
	    for (int j = 0; j < sels.Size(); j++)
	      {
		SurfaceElementIndex sei = sels[j];
		Element2d & se = (*this)[sei];
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
  void Mesh :: Distribute (NgArray<int> & volume_weights , NgArray<int>  & surface_weights, NgArray<int>  & segment_weights)
  {
    NgMPI_Comm comm = GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    if (id != 0 || ntasks == 1 ) return;

#ifdef METIS
    ParallelMetis (volume_weights, surface_weights, segment_weights);
#else
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      (*this)[ei].SetPartition(ntasks * ei/GetNE() + 1);
#endif

    /*
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
      */
    
    // MyMPI_SendCmd ("mesh");
    SendRecvMesh (); 
  }
  

#ifdef METIS5
  void Mesh :: ParallelMetis (NgArray<int> & volume_weights , NgArray<int> & surface_weights, NgArray<int> & segment_weights)  
  {
    PrintMessage (3, "call metis 5 with weights ...");
    
    // cout << "segment_weights " << segment_weights << endl;
    // cout << "surface_weights " << surface_weights << endl;
    // cout << "volume_weights " << volume_weights << endl;

    int timer = NgProfiler::CreateTimer ("Mesh::Partition");
    NgProfiler::RegionTimer reg(timer);

    idx_t ne = GetNE() + GetNSE() + GetNSeg();
    idx_t nn = GetNP();
    
    NgArray<idx_t> eptr, eind , nwgt;
    for (int i = 0; i < GetNE(); i++)
      {
	eptr.Append (eind.Size());
	
	const Element & el = VolumeElement(i+1);
	
	int ind = el.GetIndex();	
	if (volume_weights.Size()<ind)
	    nwgt.Append(0);
	else
	    nwgt.Append (volume_weights[ind -1]);
	
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSE(); i++)
      {
	eptr.Append (eind.Size());
	const Element2d & el = SurfaceElement(i+1);
	
	
	int ind = el.GetIndex(); 
	ind = GetFaceDescriptor(ind).BCProperty();
	if (surface_weights.Size()<ind)
	    nwgt.Append(0);
	else
	    nwgt.Append (surface_weights[ind -1]);

	
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSeg(); i++)
      {
	eptr.Append (eind.Size());
	
	const Segment & el = LineSegment(i+1);	
	
	int ind = el.si;
	if (segment_weights.Size()<ind)
	    nwgt.Append(0);
	else
	    nwgt.Append (segment_weights[ind -1]);
	
	eind.Append (el[0]);
	eind.Append (el[1]);
      }
      
    eptr.Append (eind.Size());
    NgArray<idx_t> epart(ne), npart(nn);

    idxtype nparts = GetCommunicator().Size()-1;
    vol_partition.SetSize(GetNE());
    surf_partition.SetSize(GetNSE());
    seg_partition.SetSize(GetNSeg());
    
    if (nparts == 1)
      {
        for (int i = 0; i < GetNE(); i++)
          // VolumeElement(i+1).SetPartition(1);
          vol_partition[i] = 1;
        for (int i = 0; i < GetNSE(); i++)
          // SurfaceElement(i+1).SetPartition(1);
          surf_partition[i] = 1;
        for (int i = 0; i < GetNSeg(); i++)
          // LineSegment(i+1).SetPartition(1);
          seg_partition[i] = 1;
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
      vol_partition[i] = epart[i] + 1;
    for (int i = 0; i < GetNSE(); i++)
      // SurfaceElement(i+1).SetPartition(epart[i+GetNE()] + 1);
      surf_partition[i] = epart[i+GetNE()] + 1;
    for (int i = 0; i < GetNSeg(); i++)
      // LineSegment(i+1).SetPartition(epart[i+GetNE()+GetNSE()] + 1);
      seg_partition[i] = epart[i+GetNE()+GetNSE()] + 1;
  }
#endif 



//===========================================================================================









#ifdef METIS4
  void Mesh :: ParallelMetis ( )  
  {
    int timer = NgProfiler::CreateTimer ("Mesh::Partition");
    NgProfiler::RegionTimer reg(timer);

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
	NgArray<idxtype> elmnts(ne*npe);
	
	int etype;
	if (elementtype == TET)
	  etype = 2;
	else if (elementtype == HEX)
	  etype = 3;
	
    
	for (int i=1; i<=ne; i++)
	  for (int j=1; j<=npe; j++)
	    elmnts[(i-1)*npe+(j-1)] = VolumeElement(i).PNum(j)-1;
	
	int numflag = 0;
	int nparts = ntasks-1;
	int ncommon = 3;
	int edgecut;
	NgArray<idxtype> epart(ne), npart(nn);
	
	//     if ( ntasks == 1 ) 
	//       {
	// 	(*this) = *mastermesh;
	// 	nparts = 4;	   
	// 	metis :: METIS_PartMeshDual (&ne, &nn, elmnts, &etype, &numflag, &nparts,
	// 				     &edgecut, epart, npart);
	// 	cout << "done" << endl;
	
	// 	cout << "edge-cut: " << edgecut << ", balance: " << metis :: ComputeElementBalance(ne, nparts, epart) << endl;
	
	// 	for (int i=1; i<=ne; i++)
	// 	  {
	// 	    mastermesh->VolumeElement(i).SetPartition(epart[i-1]);
	// 	  }
	
	// 	return;
	//       }
	
	
	int timermetis = NgProfiler::CreateTimer ("Metis itself");
	NgProfiler::StartTimer (timermetis);
	
#ifdef METIS4
	cout << "call metis(4)_PartMeshDual ... " << flush;
	METIS_PartMeshDual (&ne, &nn, &elmnts[0], &etype, &numflag, &nparts,
			    &edgecut, &epart[0], &npart[0]);
#else
	cout << "call metis(5)_PartMeshDual ... " << endl;
	// idx_t options[METIS_NOPTIONS];
	
	NgArray<idx_t> eptr(ne+1);
	for (int j = 0; j < ne+1; j++)
	  eptr[j] = 4*j;
	
	METIS_PartMeshDual (&ne, &nn, &eptr[0], &elmnts[0], NULL, NULL, &ncommon, &nparts,
			    NULL, NULL,
			    &edgecut, &epart[0], &npart[0]);
#endif
	
	NgProfiler::StopTimer (timermetis);
	
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
    

    for (int sei = 1; sei <= GetNSE(); sei++ )
      {
	int ei1, ei2;
	GetTopology().GetSurface2VolumeElement (sei, ei1, ei2);
	Element2d & sel = SurfaceElement (sei);

        for (int j = 0; j < 2; j++)
          {
            int ei = (j == 0) ? ei1 : ei2;
            if ( ei > 0 && ei <= GetNE() )
              {
		sel.SetPartition (VolumeElement(ei).GetPartition());
		break;
	      }
	  }	
      }
    
  }
#endif


  void Mesh :: PartHybridMesh () 
  {
#ifdef METIS
    int ne = GetNE();
    
    int nn = GetNP();
    int nedges = topology.GetNEdges();

    idxtype  *xadj, * adjacency, *v_weights = NULL, *e_weights = NULL;

    int weightflag = 0;
    int numflag = 0;
    int nparts = ntasks - 1;

    int options[5];
    options[0] = 0;
    int edgecut;
    idxtype * part;

    xadj = new idxtype[nn+1];
    part = new idxtype[nn];

    NgArray<int> cnt(nn+1);
    cnt = 0;

    for ( int edge = 1; edge <= nedges; edge++ )
      {
	int v1, v2;
	topology.GetEdgeVertices ( edge, v1, v2);
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

    for ( int edge = 1; edge <= nedges; edge++ )
      {
	int v1, v2;
	topology.GetEdgeVertices ( edge, v1, v2);
	adjacency[ xadj[v1-1] + cnt[v1-1] ] = v2-1;
	adjacency[ xadj[v2-1] + cnt[v2-1] ] = v1-1;
	cnt[v1-1]++;
	cnt[v2-1]++;
      }

    for ( int vert = 0; vert < nn; vert++ )
      {
	NgFlatArray<idxtype> array ( cnt[vert], &adjacency[ xadj[vert] ] );
	BubbleSort(array);
      }

#ifdef METIS4
    METIS_PartGraphKway ( &nn, xadj, adjacency, v_weights, e_weights, &weightflag, 
			  &numflag, &nparts, options, &edgecut, part );
#else
    cout << "currently not supported (metis5), A" << endl;
#endif

    NgArray<int> nodesinpart(ntasks);
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


  void Mesh :: PartDualHybridMesh ( ) // NgArray<int> & neloc ) 
  {
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

    NgArray<int, 0> facevolels1(nfaces), facevolels2(nfaces);
    facevolels1 = -1;
    facevolels2 = -1;

    NgArray<int, 0> elfaces;
    xadj = new idxtype[ne+1];
    part = new idxtype[ne];

    NgArray<int, 0> cnt(ne+1);
    cnt = 0;

    for ( int el=1; el <= ne; el++ )
      {
	Element volel = VolumeElement(el);
	topology.GetElementFaces(el, elfaces);
	for ( int i = 0; i < elfaces.Size(); i++ )
	  {
	    if ( facevolels1[elfaces[i]-1] == -1 )
	      facevolels1[elfaces[i]-1] = el;
	    else
	      {
		facevolels2[elfaces[i]-1] = el;
		cnt[facevolels1[elfaces[i]-1]-1]++;
		cnt[facevolels2[elfaces[i]-1]-1]++;
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

    for ( int face = 1; face <= nfaces; face++ )
      {
	int e1, e2;
	e1 = facevolels1[face-1];
	e2 = facevolels2[face-1];
	if ( e2 == -1 ) continue;
	adjacency[ xadj[e1-1] + cnt[e1-1] ] = e2-1;
	adjacency[ xadj[e2-1] + cnt[e2-1] ] = e1-1;
	cnt[e1-1]++;
	cnt[e2-1]++;
      }

    for ( int el = 0; el < ne; el++ )
      {
	NgFlatArray<idxtype> array ( cnt[el], &adjacency[ xadj[el] ] );
	BubbleSort(array);
      }

    int timermetis = NgProfiler::CreateTimer ("Metis itself");
    NgProfiler::StartTimer (timermetis);

#ifdef METIS4
    METIS_PartGraphKway ( &ne, xadj, adjacency, v_weights, e_weights, &weightflag, 
			  &numflag, &nparts, options, &edgecut, part );
#else
    cout << "currently not supported (metis5), B" << endl;
#endif


    NgProfiler::StopTimer (timermetis);

    NgArray<int> nodesinpart(ntasks);

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

  }





  void Mesh :: PartDualHybridMesh2D ( ) 
  {
#ifdef METIS
    idxtype ne = GetNSE();
    int nv = GetNV();

    NgArray<idxtype> xadj(ne+1);
    NgArray<idxtype> adjacency(ne*4);

    // first, build the vertex 2 element table:
    NgArray<int, PointIndex::BASE> cnt(nv);
    cnt = 0;
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      for (int j = 0; j < (*this)[sei].GetNP(); j++)
	cnt[ (*this)[sei][j] ] ++;
    
    TABLE<SurfaceElementIndex, PointIndex::BASE> vert2els(cnt);
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      for (int j = 0; j < (*this)[sei].GetNP(); j++)
	vert2els.Add ((*this)[sei][j], sei);
    

    // find all neighbour elements
    int cntnb = 0;
    NgArray<int> marks(ne);   // to visit each neighbour just once
    marks = -1;
    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      {
	xadj[sei] = cntnb;
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
		    adjacency[cntnb++] = sei2;
		  }
	      }
	  }
      }
    xadj[ne] = cntnb;

    idxtype *v_weights = NULL, *e_weights = NULL;

    idxtype weightflag = 0;
    // int numflag = 0;
    idxtype nparts = ntasks - 1;

    idxtype edgecut;
    NgArray<idxtype> part(ne);

    for ( int el = 0; el < ne; el++ )
      BubbleSort (adjacency.Range (xadj[el], xadj[el+1]));

#ifdef METIS4	
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
    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      // (*this) [sei].SetPartition (part[sei]+1);
      surf_partition[sei] = part[sei]+1;
#else
    cout << "partdualmesh not available" << endl;
#endif

  }



}



#endif

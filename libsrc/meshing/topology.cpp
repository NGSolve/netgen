#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{
  using ngcore::ParallelForRange;
  using ngcore::ParallelFor;
  using ngcore::INT;
  using ngcore::TasksPerThread;

  /*
  template <class T>
  void QuickSortRec (NgFlatArray<T> data,
		     int left, int right)
  {
    int i = left;
    int j = right;
    T midval = data[(left+right)/2];
  
    do
      {
	while (data[i] < midval) i++;
	while (midval < data[j]) j--;
      
	if (i <= j)
	  {
	    Swap (data[i], data[j]);
	    i++; j--;
	  }
      }
    while (i <= j);
    if (left < j) QuickSortRec (data, left, j);
    if (i < right) QuickSortRec (data, i, right);
  }

  template <class T>
  void QuickSort (NgFlatArray<T> data)
  {
    if (data.Size() > 1)
      QuickSortRec (data, 0, data.Size()-1);
  }
  */



  
  MeshTopology ::  MeshTopology (const Mesh & amesh)
    : mesh(&amesh)
  {
    buildedges = static_buildedges;
    buildfaces = static_buildfaces;
    buildvertex2element = static_buildvertex2element;
    timestamp = -1;
  }

  MeshTopology :: ~MeshTopology () { ;  }

  bool MeshTopology :: NeedsUpdate() const
  { return (timestamp <= mesh->GetTimeStamp()); }

  
  void MeshTopology :: EnableTable (string name, bool set)
  {
    if (name == "edges")
      SetBuildEdges(set);
    else if (name == "faces")
      SetBuildFaces(set);
    else if (name == "parentedges")
      SetBuildParentEdges(set);
    else if (name == "parentfaces")
      SetBuildParentFaces(set);
    else
      throw Exception ("nothing known about table "+name +"\n"
                       "known are 'edges', 'faces', 'parentedges', 'parentfaces'");
  }

  bool MeshTopology :: static_buildedges = true; 
  bool MeshTopology :: static_buildfaces = true; 
  bool MeshTopology :: static_buildvertex2element = true; 
  
  void MeshTopology :: EnableTableStatic (string name, bool set)
  {
    if (name == "edges")
      static_buildedges = set;
    else if (name == "faces")
      static_buildfaces = set;      
    else if (name == "vertex2element")
      static_buildvertex2element = set;      
    else
      throw Exception ("nothing known about table "+name +"\n"
                       "known are 'edges', 'faces', 'vertex2element'");
  }

  
  template <typename FUNC>
  void LoopOverEdges (const Mesh & mesh, MeshTopology & top, PointIndex v,
                      FUNC func)
  {
    for (ElementIndex elnr : top.GetVertexElements(v))
      {
        const Element & el = mesh[elnr];

        auto eledges = MeshTopology::GetEdges (el.GetType());
        for (int k = 0; k < eledges.Size(); k++)
          {
            INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);

            int edgedir = (edge.I1() > edge.I2());
            if (edgedir) swap (edge.I1(), edge.I2());
            if (edge.I1() != v) continue;

            func (edge, elnr, k, 3);
          }      
      }
    
    for (SurfaceElementIndex elnr : top.GetVertexSurfaceElements(v))
      {
        const Element2d & el = mesh[elnr];

        auto eledges = MeshTopology::GetEdges (el.GetType());
        for (int k = 0; k < eledges.Size(); k++)
          {
            INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);

            int edgedir = (edge.I1() > edge.I2());
            if (edgedir) swap (edge.I1(), edge.I2());
            
            if (edge.I1() != v) continue;

            func (edge, elnr, k, 2);
          }        
      }
    
    for (SegmentIndex elnr : top.GetVertexSegments(v))
      {
        const Segment & el = mesh[elnr];
        INDEX_2 edge(el[0], el[1]);
        int edgedir = (edge.I1() > edge.I2());
        if (edgedir) swap (edge.I1(), edge.I2());
        
        edge.Sort();
        if (edge.I1() != v) continue;
        
        func (edge, elnr, 0, 1);
      }
  }
  
  template <typename FUNC>
  void LoopOverFaces (const Mesh & mesh, MeshTopology & top, PointIndex v,
                      FUNC func)
  {
    for (ElementIndex elnr : top.GetVertexElements(v))
      {
        const Element & el = mesh[elnr];
	
        int nelfaces = MeshTopology::GetNFaces (el.GetType());
        const ELEMENT_FACE * elfaces = MeshTopology::GetFaces0 (el.GetType());
	
        for (int j = 0; j < nelfaces; j++)
          if (elfaces[j][3] < 0)
            
            { // triangle
              INDEX_4 face(el[elfaces[j][0]], el[elfaces[j][1]], 
                           el[elfaces[j][2]], 0);
              
              int facedir = 0;
              if (face.I1() > face.I2())
                { swap (face.I1(), face.I2()); facedir += 1; }
              if (face.I2() > face.I3())
                { swap (face.I2(), face.I3()); facedir += 2; }
              if (face.I1() > face.I2())
                { swap (face.I1(), face.I2()); facedir += 4; }

              if (face.I1() != v) continue;

              func (face, elnr, j, true);
            }
        /*
              if (pass == 1)
                {
                  if (!vert2face.Used (face))
                    {
                      nfa++;
                      vert2face.Set (face, nfa);
                      INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
                      face2vert.Append (hface);
                    }
                }
              else
                {
                  int facenum = vert2face.Get(face);
                  faces[elnr][j].fnr = facenum-1;
                  faces[elnr][j].forient = facedir;
                }
        */
          else
            {
              // quad
              // int facenum;
              INDEX_4 face4(el[elfaces[j][0]], el[elfaces[j][1]],
                            el[elfaces[j][2]], el[elfaces[j][3]]);
              
              int facedir = 0;
              if (min2 (face4.I1(), face4.I2()) > 
                  min2 (face4.I4(), face4.I3())) 
                {  // z - flip
                  facedir += 1; 
                  swap (face4.I1(), face4.I4());
                  swap (face4.I2(), face4.I3());
                }
              if (min2 (face4.I1(), face4.I4()) >
                  min2 (face4.I2(), face4.I3())) 
                {  // x - flip
                  facedir += 2; 
                  swap (face4.I1(), face4.I2());
                  swap (face4.I3(), face4.I4());
                }
              if (face4.I2() > face4.I4())
                {  // diagonal flip
                  facedir += 4; 
                  swap (face4.I2(), face4.I4());
                }
              
              if (face4.I1() != v) continue;
              
              func(face4, elnr, j, true);
                /*
              INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
              
              
              if (vert2face.Used (face))
                {
                  facenum = vert2face.Get(face);
                }
              else
                {
                  if (pass == 2) cout << "hier in pass 2" << endl;
                  nfa++;
                  vert2face.Set (face, nfa);
                  facenum = nfa;
                  
                  INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I4());
                  face2vert.Append (hface);
                }
              
              faces[elnr][j].fnr = facenum-1;
                          faces[elnr][j].forient = facedir;
			}
                */
                }
      }

    
    for (SurfaceElementIndex elnr : top.GetVertexSurfaceElements(v))          
      {
        const Element2d & el = mesh[elnr];
        
        const ELEMENT_FACE * elfaces = MeshTopology::GetFaces1 (el.GetType());
	
        if (elfaces[0][3] == 0)
          
          { // triangle
            
            // int facenum;
            int facedir;
            
            INDEX_4 face(el.PNum(elfaces[0][0]),
                         el.PNum(elfaces[0][1]),
                         el.PNum(elfaces[0][2]),0);
            
            facedir = 0;
            if (face.I1() > face.I2())
              {
                swap (face.I1(), face.I2());
                facedir += 1;
              }
            if (face.I2() > face.I3())
              {
                swap (face.I2(), face.I3());
                facedir += 2;
              }
            if (face.I1() > face.I2())
              {
                swap (face.I1(), face.I2());
                facedir += 4;
              }
            
            if (face.I1() != v) continue;
            
            func(face, elnr, 0, false);
            /*
              if (vert2face.Used (face))
              facenum = vert2face.Get(face);
              else
              {
              nfa++;
              vert2face.Set (face, nfa);
              facenum = nfa;
              
              INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
              face2vert.Append (hface);
              }
                
              surffaces[elnr].fnr = facenum-1;
              surffaces[elnr].forient = facedir;
            */
          }
        
        else
          
          {
            // quad
            // int facenum;
            int facedir;
            
            INDEX_4 face4(el.PNum(elfaces[0][0]),
                          el.PNum(elfaces[0][1]),
                          el.PNum(elfaces[0][2]),
                          el.PNum(elfaces[0][3]));
            
            facedir = 0;
            if (min2 (face4.I1(), face4.I2()) > 
                min2 (face4.I4(), face4.I3())) 
              {  // z - orientation
                facedir += 1; 
                swap (face4.I1(), face4.I4());
                swap (face4.I2(), face4.I3());
              }
            if (min2 (face4.I1(), face4.I4()) >
                min2 (face4.I2(), face4.I3())) 
              {  // x - orientation
                facedir += 2; 
                swap (face4.I1(), face4.I2());
                swap (face4.I3(), face4.I4());
              }
            if (face4.I2() > face4.I4())
              { 
                facedir += 4; 
                swap (face4.I2(), face4.I4());
              }
            
            if (face4.I1() != v) continue;
            func(face4, elnr, 0, false);
            /*
              INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
		
              if (vert2face.Used (face))
              facenum = vert2face.Get(face);
              else
              {
              nfa++;
              vert2face.Set (face, nfa);
              facenum = nfa;
                    
              INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I4());
              face2vert.Append (hface);
              }
                
              surffaces[elnr].fnr = facenum-1;
              surffaces[elnr].forient = facedir;
              }
            */
          }
            
      }
  }
  
  
  void MeshTopology :: Update (NgTaskManager tm_unused, NgTracer tracer)
  {
    static Timer timer("Topology::Update");
    static Timer timer_tables("Build vertex to element table");
    RegionTimer reg (timer);

#ifdef PARALLEL
    // ParallelMeshTopology & paralleltop = mesh.GetParallelTopology();
#endif

    auto id = this->mesh->GetCommunicator().Rank();
    auto ntasks = this->mesh->GetCommunicator().Size();
  
    if (timestamp > mesh->GetTimeStamp()) return;
  
    int ne = mesh->GetNE();
    int nse = mesh->GetNSE();
    int nseg = mesh->GetNSeg();
    int np = mesh->GetNP();
    int nv = mesh->GetNV(); 

    if (id == 0)
      PrintMessage (3, "Update mesh topology");

    (*testout) << " UPDATE MESH TOPOLOGY " << endl; 
    (*testout) << "ne   = " << ne << endl;
    (*testout) << "nse  = " << nse << endl;
    (*testout) << "nseg = " << nseg << endl;
    (*testout) << "np   = " << np << endl;
    (*testout) << "nv   = " << nv << endl;

    (*tracer) ("Topology::Update setup tables", false);
    NgArray<int,PointIndex::BASE> cnt(nv);

    /*
      generate:
      vertex to element 
      vertex to surface element
      vertex to segment 
    */

    if (buildvertex2element)
      {
        timer_tables.Start();
        vert2element = mesh->CreatePoint2ElementTable();
        vert2surfelement = mesh->CreatePoint2SurfaceElementTable(0);

        vert2segment = ngcore::CreateSortedTable<SegmentIndex, PointIndex>( mesh->LineSegments().Range(),
                                                                            [&](auto & table, SegmentIndex segi)
                                                                            {
                                                                              const Segment & seg = (*mesh)[segi];
                                                                              table.Add (seg[0], segi);
                                                                              table.Add (seg[1], segi);
                                                                            }, np);
        
        vert2pointelement = ngcore::CreateSortedTable<int, PointIndex>( mesh->pointelements.Range(),
                                                                        [&](auto & table, int pei)
                                                                        {
                                                                          const Element0d & pointel = mesh->pointelements[pei];
                                                                          table.Add(pointel.pnum, pei);
                                                                        }, np);
        timer_tables.Stop();
      }


    (*tracer) ("Topology::Update setup tables", true);

    
    if (buildedges)
      {
        static Timer timer1("topology::buildedges");
        RegionTimer reg1(timer1);
	
	if (id == 0)
	  PrintMessage (5, "Update edges ");
      
	edges.SetSize(ne);
	surfedges.SetSize(nse); 
	segedges.SetSize(nseg);

        /*
	for (int i = 0; i < ne; i++)
	  for (int j = 0; j < 12; j++)
	    edges[i][j].nr = -1;
	for (int i = 0; i < nse; i++)
	  for (int j = 0; j < 4; j++)
	    surfedges[i][j].nr = -1;
        */
        ParallelFor (ne, [this](auto i)
                     {
                       for (auto & e : edges[i])
                         e = -1;
                     });
	ParallelFor (nse, [this](auto i)
                     {
                       for (auto & e : surfedges[i])
                         e = -1;
                     });


        
	// keep existing edges
	cnt = 0;
	for (int i = 0; i < edge2vert.Size(); i++)
	  cnt[edge2vert[i][0]]++;
	TABLE<int,PointIndex::BASE> vert2edge (cnt);
	for (int i = 0; i < edge2vert.Size(); i++)
	  vert2edge.AddSave (edge2vert[i][0], i);

	// ensure all coarse grid and intermediate level edges
	cnt = 0;
	// for (int i = mesh->mlbetweennodes.Begin(); i < mesh->mlbetweennodes.End(); i++)
        for (int i : mesh->mlbetweennodes.Range())
	  {
	    INDEX_2 parents = Sort (mesh->mlbetweennodes[i]);
	    if (parents[0] >= PointIndex::BASE) cnt[parents[0]]++;
	  }
	TABLE<int,PointIndex::BASE> vert2vertcoarse (cnt);
	// for (int i = mesh->mlbetweennodes.Begin(); i < mesh->mlbetweennodes.End(); i++)
        for (int i : mesh->mlbetweennodes.Range())
	  {
	    INDEX_2 parents = Sort (mesh->mlbetweennodes[i]);
	    if (parents[0] >= PointIndex::BASE) vert2vertcoarse.AddSave (parents[0], parents[1]);
	  }



	int max_edge_on_vertex = 0;
	for (int i = PointIndex::BASE; i < nv+PointIndex::BASE; i++)
	  {
	    int onv = vert2edge[i].Size() + vert2vertcoarse[i].Size() +
              4*(vert2element)[i].Size() + 2*(vert2surfelement)[i].Size() + (vert2segment)[i].Size();
	    max_edge_on_vertex = max (onv, max_edge_on_vertex);
	  }

        
        // count edges associated with vertices
        cnt = 0;

        ParallelForRange
          (mesh->GetNV(), // Points().Size(),
           [&] (IntRange r)
           {
             auto begin = r.First();
             auto end = r.Next();
             // INDEX_CLOSED_HASHTABLE<int> v2eht(2*max_edge_on_vertex+10);
             ngcore::ClosedHashTable<int, int> v2eht(2*max_edge_on_vertex+10);
             for (PointIndex v = begin+PointIndex::BASE;
                  v < end+PointIndex::BASE; v++)
               {
                 v2eht.DeleteData();
                 for (int ednr : vert2edge[v])
                   {
                     int v2 = edge2vert[ednr][1];
                     v2eht.Set (v2, ednr);
                   }

                 size_t usedold = v2eht.UsedElements();
                 
                 for (int v2 : vert2vertcoarse[v])
                   v2eht.Set (v2, 33);   // some value                   
                 
                 LoopOverEdges (*mesh, *this, v,
                                [&] (INDEX_2 edge, int elnr, int loc_edge, int element_dim)
                                {
                                  v2eht.Set (edge[1], 33); // something                                  
                                });
                 
                 cnt[v] = v2eht.UsedElements()-usedold;
               }
           }, TasksPerThread(4) );

        // accumulate number of edges
        int ned = edge2vert.Size();

        for (size_t v : cnt.Range())
          {
            auto hv = cnt[v];
            cnt[v] = ned;
            ned += hv;
          }
        edge2vert.SetSize(ned);
        edge2segment.SetSize(ned);
        edge2segment = -1;

        // INDEX_CLOSED_HASHTABLE<int> v2eht(2*max_edge_on_vertex+10);
	// NgArray<int> vertex2;
	// for (PointIndex v = PointIndex::BASE; v < nv+PointIndex::BASE; v++)

        ParallelForRange
          (mesh->GetNV(), // Points().Size(),
           [&] (IntRange r)
           {
             auto begin = r.First();
             auto end = r.Next();
             // INDEX_CLOSED_HASHTABLE<int> v2eht(2*max_edge_on_vertex+10);
             ngcore::ClosedHashTable<int, int> v2eht(2*max_edge_on_vertex+10);

             Array<int> vertex2;
             for (PointIndex v = begin+PointIndex::BASE;
                  v < end+PointIndex::BASE; v++)
               {
                 int ned = cnt[v];
                 v2eht.DeleteData();            
                 vertex2.SetSize0 ();
                 
                 for (int ednr : vert2edge[v])
                   {
                     int v2 = edge2vert[ednr][1];
                     v2eht.Set (v2, ednr);
                   }
                 
                 for (int v2 : vert2vertcoarse[v])
                   if (!v2eht.Used(v2))
                     {
                       v2eht.Set (v2, 33);   // some value
                       vertex2.Append (v2);
                     }
                 
                 LoopOverEdges (*mesh, *this, v,
                                [&](INDEX_2 edge, int elnr, int loc_edge, int element_dim)
                                {
                                  size_t pos;
                                  if (v2eht.PositionCreate(edge[1], pos))
                                    {
                                      vertex2.Append (edge[1]);
                                      v2eht.SetData (pos, 33);
                                    }
                                  /*
                                  if (!v2eht.Used(edge.I2()))
                                    {
                                      vertex2.Append (edge.I2());
                                      v2eht.Set (edge.I2(), 33); 
                                    }
                                  */
                                });
                 
                 QuickSort (vertex2);

                 /*
                 for (int j = 0; j < vertex2.Size(); j++)
                   {
                     v2eht.Set (vertex2[j], ned);
                     edge2vert[ned] = { v, vertex2[j] };
                     ned++;
                   }
                 */
                 for (auto v2 : vertex2)
                   {
                     v2eht.Set (v2, ned);
                     edge2vert[ned] = { v, v2 };
                     ned++;
                   }
                 
                 LoopOverEdges (*mesh, *this, v,
                                [&](INDEX_2 edge, int elnr, int loc_edge, int element_dim)
                                {
                                  int edgenum = v2eht.Get(edge[1]);
                                  switch (element_dim)
                                    {
                                    case 3:
                                      edges[elnr][loc_edge] = edgenum;
                                      break;
                                    case 2:
                                      surfedges[elnr][loc_edge] = edgenum;
                                      break;
                                    case 1:
                                      segedges[elnr] = edgenum;
                                      edge2segment[edgenum] = elnr;
                                      break;
                                    }
                                });
               }
           }, TasksPerThread(4) );


        if (build_parent_edges)
        {
          static Timer t("build_hierarchy"); RegionTimer reg(t);
          cnt = 0;
          for (auto verts : edge2vert) cnt[verts[0]]++;
          TABLE<int,PointIndex::BASE> vert2edge (cnt);
          for (auto i : edge2vert.Range())
            vert2edge.AddSave (edge2vert[i][0], i);

          // build edge hierarchy:
          parent_edges.SetSize (ned);
          parent_edges = { -1, { -1, -1, -1 } };

          for (size_t i = 0; i < ned; i++)
          {
            auto verts = edge2vert[i];  // 2 vertices of edge

            if (verts[0] >= mesh->mlbetweennodes.Size()+PointIndex::BASE ||
                verts[1] >= mesh->mlbetweennodes.Size()+PointIndex::BASE)
              continue;

            auto pa0 = mesh->mlbetweennodes[verts[0]]; // two parent vertices of v0
            auto pa1 = mesh->mlbetweennodes[verts[1]]; // two parent vertices of v1

            // both vertices are on coarsest mesh
            if (!pa0[0].IsValid() && !pa1[0].IsValid())
              continue;

            int issplitedge = 0;
            if (pa0[0] == verts[1] || pa0[1] == verts[1])
              issplitedge = 1;
            if (pa1[0] == verts[0] || pa1[1] == verts[0])
              issplitedge = 2;

            if (issplitedge)
            {
              // cout << "split edge " << endl;
              // edge is obtained by splitting one edge into two parts:
              auto paedge = issplitedge == 1 ? pa0 : pa1;

              if (paedge[0] > paedge[1]) 
                Swap (paedge[0], paedge[1]);

              for (int ednr : vert2edge[paedge[0]])
                if (auto cverts = edge2vert[ednr]; cverts[1] == paedge[1])
                {
                  int orient = (paedge[0] == verts[0] || paedge[1] == verts[1]) ? 1 : 0;
                  parent_edges[i] = { orient, { ednr, -1, -1 } };			  
                }
            }
            else
            {
              bool bisect_edge = false;
              // edge is splitting edge in middle of triangle:
              for (int j = 1; j <= 2; j++)
              {
                INT<2> paedge1, paedge2, paedge3;
                int orient_inner = 0;
                if (j == 1)
                {
                  paedge1 = INT<2> (pa0[0], verts[1]);
                  paedge2 = INT<2> (pa0[1], verts[1]);
                  paedge3 = INT<2> (pa0[0], pa0[1]);
                  orient_inner = 0;
                }
                else
                {
                  paedge1 = INT<2> (pa1[0], verts[0]);
                  paedge2 = INT<2> (pa1[1], verts[0]);
                  paedge3 = INT<2> (pa1[0], pa1[1]);
                  orient_inner = 1;
                }
                if (paedge1[0] > paedge1[1]) 
                  Swap (paedge1[0], paedge1[1]);
                if (paedge2[0] > paedge2[1]) 
                  Swap (paedge2[0], paedge2[1]);
                if (paedge3[0] > paedge3[1]) 
                  Swap (paedge3[0], paedge3[1]);

                // if first vertex number is -1, then don't try to find entry in node2edge hash table
                if ( paedge1[0] == PointIndex::BASE-1 || paedge2[0] == PointIndex::BASE-1 )
                  continue;

                int paedgenr1=-1, paedgenr2=-1, paedgenr3=-1, orient1 = 0, orient2 = 0;
                for (int ednr : vert2edge[paedge1[0]])
                  if (auto cverts = edge2vert[ednr]; cverts[1] == paedge1[1])
                  {
                    paedgenr1 = ednr;
                    orient1 = (paedge1[0] == verts[0] || paedge1[1] == verts[1]) ? 1 : 0;
                  }
                for (int ednr : vert2edge[paedge2[0]])
                  if (auto cverts = edge2vert[ednr]; cverts[1] == paedge2[1])
                  {
                    paedgenr2 = ednr;
                    orient2 = (paedge2[0] == verts[0] || paedge2[1] == verts[1]) ? 1 : 0;
                  }

                for (int ednr : vert2edge[paedge3[0]])
                  if (auto cverts = edge2vert[ednr]; cverts[1] == paedge3[1])
                    paedgenr3 = ednr;

                if (paedgenr1 != -1 && paedgenr2 != -1){
                  bisect_edge = true;
                  parent_edges[i] = { orient1+2*orient2+4*orient_inner, { paedgenr1, paedgenr2, paedgenr3 } };
                }
              }

              if (!bisect_edge) // not a bisect edge (then a red edge)
              {
                INT<2> paedge1, paedge2, paedge3;
                int orient1 = 0, orient2 = 0, orient3=0;
                int orient_inner = 0;
                paedge1 = INT<2> (pa0[0], pa0[1]);
                paedge2 = INT<2> (pa1[0], pa1[1]);
                // find common vertex and the third pa edge
                if (pa0[0]==pa1[0]){// 00
                  //orient1 = 0; 
                  orient2 = 1; 
                  if (pa0[1]<pa1[1]){
                    orient3 = 1;
                    paedge3 = INT<2> (pa0[1], pa1[1]);
                  }else{
                    //orient3 = 0;
                    paedge3 = INT<2> (pa1[1], pa0[1]);
                  }
                }
                else if (pa0[0]==pa1[1]){//01
                  //orient1 = 0; 
                  //orient2 = 0; 
                  if (pa0[1]<pa1[0]){
                    orient3 = 1;
                    paedge3 = INT<2> (pa0[1], pa1[0]);
                  }else{
                    //orient3 = 0;
                    paedge3 = INT<2> (pa1[0], pa0[1]);
                  }
                }
                else if (pa0[1]==pa1[0]){//10
                  orient1 = 1; 
                  orient2 = 1; 
                  if (pa0[0]<pa1[1]){
                    orient3 = 1;
                    paedge3 = INT<2> (pa0[0], pa1[1]);
                  }else{
                    //orient3 = 0;
                    paedge3 = INT<2> (pa1[1], pa0[0]);
                  }
                }
                else if (pa0[1]==pa1[1]){//11
                  orient1 = 1; 
                  //orient2 = 0; 
                  if (pa0[0]<pa1[0]){
                    orient3 = 1;
                    paedge3 = INT<2> (pa0[0], pa1[0]);
                  }else{
                    //orient3 = 0;
                    paedge3 = INT<2> (pa1[0], pa0[0]);
                  }
                }

                int paedgenr1=-1, paedgenr2=-1, paedgenr3=-1;
                for (int ednr : vert2edge[paedge1[0]])
                  if (auto cverts = edge2vert[ednr]; cverts[1] == paedge1[1])
                    paedgenr1 = ednr;
                for (int ednr : vert2edge[paedge2[0]])
                  if (auto cverts = edge2vert[ednr]; cverts[1] == paedge2[1])
                    paedgenr2 = ednr;

                for (int ednr : vert2edge[paedge3[0]])
                  if (auto cverts = edge2vert[ednr]; cverts[1] == paedge3[1])
                    paedgenr3 = ednr;

                parent_edges[i] = { 8+orient1+2*orient2+4*orient3, { paedgenr1, paedgenr2, paedgenr3 } };

                //cout <<8+orient1+2*orient2+4*orient3  <<":"<<paedgenr1 <<", "<< paedgenr2 << ", "<< paedgenr3 << endl;
              }


              // TODO: quad edges
              /*
                 if (parentedges[i][0] == -1)
                 {
              // quad split
              if (pa1[0] != pa2[0] && 
              pa1[0] != pa2[1] && 
              pa1[1] != pa2[0] && 
              pa1[1] != pa2[1])
              for (int j = 1; j <= 2; j++)
              {
              INT<2> paedge1, paedge2;
              if (j == 1)
              {
              paedge1 = INT<2> (pa1[0], pa2[0]);
              paedge2 = INT<2> (pa1[1], pa2[1]);
              }
              else
              {
              paedge1 = INT<2> (pa1[0], pa2[1]);
              paedge2 = INT<2> (pa1[1], pa2[0]);
              }

              int paedgenr1 = 0, paedgenr2 = 0;
              int orient1 = 1, orient2 = 1;

              if (paedge1[0] > paedge1[1]) 
              {
              Swap (paedge1[0], paedge1[1]);
              orient1 = 0;
              }
              if (paedge2[0] > paedge2[1]) 
              {
              Swap (paedge2[0], paedge2[1]);
              orient2 = 0;
              }

              if ( paedge1[0] == -1 || paedge2[0] == -1 )
              continue;

              if (node2edge.Used (paedge1) && node2edge.Used (paedge2))
              {
              paedgenr1 = node2edge.Get (paedge1);
              paedgenr2 = node2edge.Get (paedge2);
              parentedges[i][0] = 2 * paedgenr1 + orient1;	      
              parentedges[i][1] = 2 * paedgenr2 + orient2;	      
              }
              }
              }

              if (parentedges[i][0] == -1)
              {
              // triangle split into quad+trig (from anisotropic pyramids)
              for (int j = 0; j < 2; j++)
              for (int k = 0; k < 2; k++)
              {
              INT<2> paedge (pa1[1-j], pa2[1-k]);
              int orientpa = 1;
              if (paedge[0] > paedge[1]) 
              {
              Swap (paedge[0], paedge[1]);
              orientpa = 0;
              }	
              if (pa1[j] == pa2[k] && node2edge.Used(paedge))
              {
              int paedgenr = node2edge.Get (paedge);
              parentedges[i][0] = 2 * paedgenr + orientpa;
              }
              }

              }
              */
            }

          }

          /*
             for (int i : Range(parent_edges))
             {
             auto [info, nrs] = parent_edges[i];
             cout << "edge " << i << " has " << info << ", nrs = " << nrs[0] << " " << nrs[1] << endl;
             }
             */
        }
      }
    

    // edge hashtable:: needed for getting parent faces	
    ngcore::ClosedHashTable<INT<2>, int> v2e(nv);
    if (build_parent_faces)
      for (auto i : Range(edge2vert))
        {
          auto edge = edge2vert[i];
          INT<2> e2(edge[0], edge[1]);
          e2.Sort();
          v2e[e2] = i;
        }

    
    // generate faces
    if (buildfaces) 
      {
	static Timer timer2("topology::buildfaces");
	// static int timer2a = NgProfiler::CreateTimer ("topology::buildfacesa");
	// static int timer2b = NgProfiler::CreateTimer ("topology::buildfacesb");
        // static int timer2b1 = NgProfiler::CreateTimer ("topology::buildfacesb1");
	// static int timer2c = NgProfiler::CreateTimer ("topology::buildfacesc");
	RegionTimer reg2 (timer2);

	if (id == 0)
	  PrintMessage (5, "Update faces ");

        // NgProfiler::StartTimer (timer2a);

	faces.SetSize(ne);
	surffaces.SetSize(nse);
  

	cnt = 0;
	for (int i = 0; i < face2vert.Size(); i++)
	  cnt[face2vert[i][0]]++;
	TABLE<int,PointIndex::BASE> vert2oldface(cnt);
	for (int i = 0; i < face2vert.Size(); i++)
	  vert2oldface.AddSave (face2vert[i][0], i);

        // find all potential intermediate faces
        Array<INT<3>> intermediate_faces;
        if (build_parent_faces)
          {
            for (ElementIndex ei = 0; ei < ne; ei++)
              for (int i = 0; i < 4; i++)
                {
                  Element2d face;
                  // cout << "element: " << (*mesh)[ei].PNums() << endl;
                  (*mesh)[ei].GetFace(i+1, face);
                  // cout << "face " << face.PNums() << endl;
                  INT<3,PointIndex> f3 = { face[0], face[1], face[2] };
                  for (int j = 0; j < 3; j++)
                    {
                      PointIndex v = f3[j];
                      if (v >= mesh->mlbetweennodes.Size()+PointIndex::BASE)
                        continue;

                      auto pa = mesh->mlbetweennodes[v];
                      for (int k = 0; k < 2; k++)
                        if (f3.Contains(pa[k]))
                          {
                            PointIndex v0 = pa[k]; // also in face
                            PointIndex v1 = pa[1-k];
                            PointIndex v2 = f3[0]+f3[1]+f3[2] - v - v0;
                            // if there is an edge connecting v1 and v2, accept
                            // the new face
                            INT<2> parentedge(v1, v2);
                            parentedge.Sort();
                            if (v2e.Used(parentedge)){ 
                              INT<3> cf3 = { v0, v1, v2 };
                              cf3.Sort();
                              // cout << "intermediate: " << cf3 << " of " << f3 << endl;
                              intermediate_faces.Append (cf3);
                            }
                          }
                    }
                }
          }
	cnt = 0;
	for (int i = 0; i < intermediate_faces.Size(); i++)
	  cnt[intermediate_faces[i][0]]++;
	TABLE<int,PointIndex::BASE> vert2intermediate(cnt);
	for (int i = 0; i < intermediate_faces.Size(); i++)
	  vert2intermediate.AddSave (intermediate_faces[i][0], i);
        // cout << "vert2intermediate = " << endl << vert2intermediate << endl;

        
	for (int elnr = 0; elnr < ne; elnr++)
	  for (int j = 0; j < 6; j++)
	    faces[elnr][j] = -1;
	

	int max_face_on_vertex = 0;
	for (int i = PointIndex::BASE; i < nv+PointIndex::BASE; i++)
	  {
	    int onv = vert2oldface[i].Size() + vert2element[i].Size() + vert2surfelement[i].Size();
	    max_face_on_vertex = max (onv, max_face_on_vertex);
	  }
	



        // NgProfiler::StopTimer (timer2a);
        // NgProfiler::StartTimer (timer2b);

        // INDEX_3_CLOSED_HASHTABLE<int> vert2face(2*max_face_on_vertex+10);         

	int oldnfa = face2vert.Size();

        // count faces associated with vertices
        cnt = 0;
        // for (auto v : mesh.Points().Range())
        // NgProfiler::StartTimer (timer2b1);
        ParallelForRange
          (mesh->GetNV(), // Points().Size(),
           [&] (IntRange r)
            {
              // auto begin = r.First();
              // auto end = r.Next();
              // INDEX_3_CLOSED_HASHTABLE<int> vert2face(2*max_face_on_vertex+10);
              ClosedHashTable<INDEX_3, int> vert2face(2*max_face_on_vertex+10); 
              // for (PointIndex v = begin+PointIndex::BASE;
              // v < end+PointIndex::BASE; v++)
              for (PointIndex v : r+PointIndex::BASE)                
                {
                  vert2face.DeleteData();
                  
                  for (int j = 0; j < vert2oldface[v].Size(); j++)
                    {
                      int fnr = vert2oldface[v][j];
                      INDEX_3 face (face2vert[fnr][0],
                                    face2vert[fnr][1],
                                    face2vert[fnr][2]);
                      vert2face.Set (face, 33);  // something
                    }
                  int cnti = 0;

                  for (int j = 0; j < vert2intermediate[v].Size(); j++)
                    {
                      int fnr = vert2intermediate[v][j];
                      INDEX_3 face (intermediate_faces[fnr][0],
                                    intermediate_faces[fnr][1],
                                    intermediate_faces[fnr][2]);
                      face.Sort();
                      if (!vert2face.Used(face))
                        {
                          cnti++;
                          vert2face.Set (face, 33); // something
                        }
                    }
                  LoopOverFaces (*mesh, *this, v,
                                 [&] (INDEX_4 i4, int elnr, int j, bool volume)
                                 {
                                   INDEX_3 face(i4[0], i4[1], i4[2]);
                                   if (!vert2face.Used (face))
                                     {
                                       cnti++;
                                       vert2face.Set (face, 33); // something
                                     }
                                 });
                  cnt[v] = cnti;
                }
            }, TasksPerThread(4) );
        // NgProfiler::StopTimer (timer2b1);
        
        // accumulate number of faces
        int nfa = oldnfa;
        // for (auto v : Range(mesh->GetNV())) // Points().Range())
        // for (size_t v = 0; v < mesh->GetNV(); v++)
        for (auto v : cnt.Range())
          {
            auto hv = cnt[v];
            cnt[v] = nfa;
            nfa += hv;
          }
        face2vert.SetSize(nfa);
        

        ParallelForRange
          (mesh->GetNV(),
           [&] (IntRange r)
            {
              // auto begin = r.First();
              // auto end = r.Next();
              // INDEX_3_CLOSED_HASHTABLE<int> vert2face(2*max_face_on_vertex+10);
              ClosedHashTable<INDEX_3, int> vert2face(2*max_face_on_vertex+10);
              /*
              for (PointIndex v = begin+PointIndex::BASE;
                   v < end+PointIndex::BASE; v++)
              */
              for (PointIndex v : r+PointIndex::BASE)
                {
                  int first_fa = cnt[v];
                  int nfa = first_fa;
                  vert2face.DeleteData();
                  
                  for (int j = 0; j < vert2oldface[v].Size(); j++)
                    {
                      int fnr = vert2oldface[v][j];
                      INDEX_3 face (face2vert[fnr][0], 
                                    face2vert[fnr][1],
                                    face2vert[fnr][2]);
                      vert2face.Set (face, fnr);
                    }
                  
                  for (int j = 0; j < vert2intermediate[v].Size(); j++)
                    {
                      int fnr = vert2intermediate[v][j];
                      INDEX_3 face (intermediate_faces[fnr][0],
                                    intermediate_faces[fnr][1],
                                    intermediate_faces[fnr][2]);
                      face.Sort();
                      /*
                      if (!vert2face.Used(face))
                        {
                          face2vert[nfa] = { face[0], face[1], face[2], 0 }; // i4;
                          vert2face.Set (face, nfa);
                          nfa++;
                        }
                      */
                      size_t pos;
                      if (vert2face.PositionCreate(face, pos))
                        {
                          face2vert[nfa] = { face[0], face[1], face[2], 0 }; // i4;
                          vert2face.SetData (pos, face, nfa);
                          nfa++;
                        }
                      
                    }
                  
                  LoopOverFaces (*mesh, *this, v,
                                 [&] (INDEX_4 i4, int elnr, int j, bool volume)
                                 {
                                   INDEX_3 face(i4.I1(), i4.I2(), i4.I3());
                                   /*
                                   if (!vert2face.Used (face))
                                     {
                                       face2vert[nfa] = { i4[0], i4[1], i4[2], i4[3] }; // i4;
                                       vert2face.Set (face, nfa);
                                       nfa++;
                                     }
                                   */
                                   size_t pos;
                                   if (vert2face.PositionCreate(face, pos))
                                     {
                                       face2vert[nfa] = { i4[0], i4[1], i4[2], i4[3] }; // i4;
                                       vert2face.SetData (pos, face, nfa);
                                       nfa++;
                                     }
                                 });
                  
                  
                  QuickSort (face2vert.Range(first_fa, nfa));
                  
                  for (int j = first_fa; j < nfa; j++)
                    {
                      if (face2vert[j][0] == v)
                        {
                          INDEX_3 face (face2vert[j][0], 
                                        face2vert[j][1], 
                                        face2vert[j][2]);
                          vert2face.Set (face, j);
                        }
                      else
                        break;
                    }
                  
                  
                  LoopOverFaces (*mesh, *this, v,
                                 [&] (INDEX_4 i4, int elnr, int j, bool volume)
                                 {
                                   INDEX_3 face(i4.I1(), i4.I2(), i4.I3());
                                   int facenum = vert2face.Get(face);
                                   if (volume)
                                     faces[elnr][j] = facenum;
                                   else
                                     surffaces[elnr] = facenum;
                                 });
                }
            }, TasksPerThread(4) );
        

	// *testout << "face2vert = " << endl << face2vert << endl;

        // NgProfiler::StopTimer (timer2b);
        // NgProfiler::StartTimer (timer2c);


	face2surfel.SetSize (nfa);
	face2surfel = 0;
	for (int i = 1; i <= nse; i++)
	  face2surfel.Elem(GetSurfaceElementFace(i)) = i;

	/*
	  cout << "build table complete" << endl;

	  cout << "faces = " << endl;

	  cout << "face2vert = " << endl << face2vert << endl;
	  cout << "surffaces = " << endl << surffaces << endl;
	  cout << "face2surfel = " << endl << face2surfel << endl;
	*/

	
	surf2volelement.SetSize (nse);
	for (int i = 1; i <= nse; i++)
	  {
	    surf2volelement.Elem(i)[0] = 0;
	    surf2volelement.Elem(i)[1] = 0;
	  }
        (*tracer) ("Topology::Update build surf2vol", false);        
	// for (int i = 0; i < ne; i++)
        ParallelFor (ne, [this](auto i)
                     {
                       for (int j = 0; j < 6; j++)
                         {
                           // int fnum = (faces.Get(i)[j]+7) / 8;
                           int fnum = faces[i][j]+1;
                           if (fnum > 0 && face2surfel.Elem(fnum))
                             {
                               int sel = face2surfel.Elem(fnum);
                               surf2volelement.Elem(sel)[1] = 
                                 surf2volelement.Elem(sel)[0];
                               surf2volelement.Elem(sel)[0] = i+1;
                             }
                         }});
        (*tracer) ("Topology::Update build surf2vol", true);        

	face2vert.SetAllocSize (face2vert.Size());

	// face table complete


#ifdef PARALLEL
	// (*testout) << " RESET Paralleltop" << endl;
	// paralleltop.Reset ();
#endif

        (*tracer) ("Topology::Update count face_els", false);
	NgArray<short int> face_els(nfa), face_surfels(nfa);
	face_els = 0;
	face_surfels = 0;

        ParallelForRange
          (ne,
           [&] (IntRange r)
            {
              /*
              NgArray<int> hfaces;              
              for (ElementIndex ei : r)
              {
                  GetElementFaces (ei+1, hfaces);
                  for (auto f : hfaces)
                    AsAtomic(face_els[f-1])++;
                }
              */
              for (ElementIndex ei : r)
                for (auto f : GetFaces(ei))
                  AsAtomic(face_els[f])++;
              
            }, TasksPerThread(4));
	for (int i = 1; i <= nse; i++)
	  face_surfels[GetSurfaceElementFace (i)-1]++;
        (*tracer) ("Topology::Update count face_els", true);


	if (ne)
	  {
	    int cnt_err = 0;
	    for (int i = 0; i < nfa; i++)
	      {
		/*
		  (*testout) << "face " << i << " has " << int(face_els[i]) << " els, " 
		  << int(face_surfels[i]) << " surfels, tot = "
		  << face_els[i] + face_surfels[i] << endl; 
		*/
		if (face_els[i] + face_surfels[i] == 1)
		  {
		    cnt_err++;
#ifdef PARALLEL
		    if ( ntasks > 1 )
		      {
			continue;
			// if ( !paralleltop.DoCoarseUpdate() ) continue;
		      }
		    else
#endif
		      {
			(*testout) << "illegal face : " << i << endl;
			(*testout) << "points = "
                                   << face2vert[i][0] << ","
                                   << face2vert[i][1] << ","
                                   << face2vert[i][2] << ","
                                   << face2vert[i][3] 
                                   << endl;
			(*testout) << "pos = ";
			for (int j = 0; j < 4; j++)
                          if (face2vert[i][j] >= 1)
                            (*testout) << (*mesh)[(PointIndex)face2vert[i][j]] << " ";
			(*testout) << endl;

			FlatArray<ElementIndex> vertels = GetVertexElements (face2vert[i][0]);
			for (int k = 0; k < vertels.Size(); k++)
			  {
			    int elfaces[10], orient[10];
			    int nf = GetElementFaces (vertels[k]+1, elfaces, orient);
			    for (int l = 0; l < nf; l++)
			      if (elfaces[l] == i)
				{
				  // (*testout) << "is face of element " << vertels[k] << endl;
			    
				  if (mesh->coarsemesh && mesh->hpelements->Size() == mesh->GetNE() )
				    {
				      const HPRefElement & hpref_el =
					(*mesh->hpelements) [ (*mesh)[vertels[k]].hp_elnr];
				      (*testout) << "coarse eleme = " << hpref_el.coarse_elnr << endl;
				    }

				}
			  }
		      }
		  }
	      }

	    if (cnt_err && ntasks == 1)
	      cout << IM(5) << cnt_err << " elements are not matching !!!" << endl;
	  }
        // NgProfiler::StopTimer (timer2c);


        if (build_parent_faces)
          {
            // tets only
            if (id == 0)
              PrintMessage (5, "build face hierarchy");

            // cout << "f2v = " << face2vert << endl;
            
            ngcore::ClosedHashTable<INT<3>, int> v2f(nv);
            for (auto i : Range(face2vert))
              {
                auto face = face2vert[i];
                INT<3> f3(face[0], face[1], face[2]);
                f3.Sort();
                v2f[f3] = i;
              }

            // cout << "v2f:" << endl << v2f << endl;
           
            parent_faces.SetSize (nfa);
            parent_faces = { -1, { -1, -1, -1, -1 } };

            for (auto i : Range(nfa))
              {
                INT<3,PointIndex> f3(face2vert[i][0], face2vert[i][1], face2vert[i][2]);


                // face on coarses level ?
                bool all_vert_coarse = true;
                for (int k = 0; k < 3; k++)
                  {
                    PointIndex vb = f3[k]; 
                    if (vb >= mesh->mlbetweennodes.Size()+PointIndex::BASE)
                      continue;
                    auto parents = mesh->mlbetweennodes[vb];
                    if (parents[0] >= PointIndex::BASE)
                      all_vert_coarse = false;
                  }
                if (all_vert_coarse) continue;

                
                
                // find a vertex, such that one of its parent is a trig vertex
                bool issplit = false;
                for (int k = 0; k < 3; k++)
                  {
                    PointIndex vb = f3[k]; // assume vb as the new bisect vert
                    if (vb >= mesh->mlbetweennodes.Size()+PointIndex::BASE)
                      continue;
                    auto parents = mesh->mlbetweennodes[vb];
                    
                    // is face part of one parent face (boundary-face) ?
                    for (int j = 0; j < 2; j++)
                      {
                        if (f3.Contains(parents[j]))
                          {
                            PointIndex v0 = parents[j];
                            PointIndex v1 = parents[1-j];
                            
                            // the third one, on the tip
                            PointIndex v2 = f3[0]+f3[1]+f3[2] - v0 - vb;
                            
                            // if there is an edge connecting v1 and v2, accept
                            // the new face
                            INT<2> parentedge(v1, v2);
                            parentedge.Sort();
                            if (v2e.Used(parentedge)){ 
                              INT<3> parentverts(v0, v1, v2);
                              parentverts.Sort();

                              int classnr = 0;
                              if (v2 > vb) { Swap (v2, vb); classnr += 1; }                            
                              if (v0 > v1) { Swap (v0, v1); classnr += 2; }
                              if (v1 > v2) { Swap (v1, v2); classnr += 4; }
                              if (v0 > v1) { Swap (v0, v1); classnr += 8; }

                              if (v2f.Used(parentverts))
                              {
                                int pafacenr = v2f[parentverts];
                                // cout << "parent-face = " << pafacenr << endl;
                                parent_faces[i] = { classnr, { pafacenr, -1, -1, -1 } };
                              }
                              else
                              {
                                cout << "missing parent face: " << parentverts << endl;
                              }
                              issplit=true;
                              break;
                            }
                          }
                        }
                  }

                /*
                // is face a new face (bisect-face) ?
                if (!issplit)
                  for (int k = 0; k < 3; k++)
                    {
                      PointIndex vb = f3[k]; // assume vb as the new bisect vert
                      if (vb >= mesh->mlbetweennodes.Size()+PointIndex::BASE)
                        continue;
                      auto parents = mesh->mlbetweennodes[vb];

                      PointIndex v0 = parents[0];
                      PointIndex v1 = parents[1];
                      PointIndex v2 = f3[(k+1)%3];
                      PointIndex v3 = f3[(k+2)%3];
                      INT<3> parentedge1(v0, v2);
                      parentedge1.Sort();
                      INT<3> parentedge2(v0, v3);
                      parentedge2.Sort();
                      INT<3> parentedge3(v1, v2);
                      parentedge3.Sort();
                      INT<3> parentedge4(v1, v3);
                      parentedge4.Sort();
                      // if edges [v0,v2], [v0, v3], [v1,v2], [v1,v3] exists
                      // then vb is the bisecting edge
                      if (v2e.Used(parentedge1) && v2e.Used(parentedge2) 
                          && v2e.Used(parentedge3) && v2e.Used(parentedge4) 
                         ){
                        int classnr;
                        if (k==2){// vb is the largest vert: 6 cases
                          // by default v0 < v1, v2 < v3
                          if (v1 < v2) classnr = 0;
                          else if (v1 < v3 && v0 < v2) classnr = 1;
                          else if (v0 < v2) classnr = 2;
                          else if (v1 < v3) classnr = 3;
                          else if (v0 < v3) classnr = 4;
                          else classnr = 5;
                        }else if (k==1){// vb is the second largest vert: 3 cases
                          // by default v0 < v1, v3 < v2
                          if (v1 < v3) classnr = 6;
                          else if (v0 < v3) classnr = 7;
                          else classnr = 8;
                        }else {// vb is the third largest vert: 1 case
                          // by default v0 < v1 < vb <  v2 < v3
                          classnr=9;
                        }
                        INT<3> parentverts1(v0, v2, v3);
                        parentverts1.Sort();
                        INT<3> parentverts2(v1, v2, v3);
                        parentverts2.Sort();
                        INT<3> parentverts3(v0, v1, v2);
                        parentverts3.Sort();
                        INT<3> parentverts4(v0, v1, v3);
                        parentverts4.Sort();
                        int pafacenr1=-1, pafacenr2=-1, pafacenr3=-1, pafacenr4=-1;
                        if (v2f.Used(parentverts1))
                        {
                          pafacenr1 = v2f[parentverts1];
                          // cout << "parent-face1 = " << pafacenr1<< endl ;
                        }
                        if (v2f.Used(parentverts2))
                        {
                          pafacenr2 = v2f[parentverts2];
                          // cout << "parent-face2 = " << pafacenr2<< endl  ;
                        }
                        if (v2f.Used(parentverts3))
                        {
                          pafacenr3 = v2f[parentverts3];
                          // cout << "parent-face3 = " << pafacenr3<< endl  ;
                        }
                        if (v2f.Used(parentverts4))
                        {
                          pafacenr4 = v2f[parentverts4];
                          // cout << "parent-face4 = " << pafacenr4<< endl  ;
                        }

                        if (k == 0 || k == 2)
                          parent_faces[i] = { classnr, { pafacenr2, pafacenr1,
                                                         pafacenr4, pafacenr3} };
                        else
                          parent_faces[i] = { classnr, { pafacenr2, pafacenr1,
                                                         pafacenr3, pafacenr4} };
                        break;
                      }
                    }
                */

                // is face a new face (bisect-face) ?
                if (!issplit)
                  for (int k = 0; k < 3; k++)
                    {
                      PointIndex vb = f3[k]; // assume vb as the new bisect vert
                      if (vb >= mesh->mlbetweennodes.Size()+PointIndex::BASE)
                        continue;
                      auto parents = mesh->mlbetweennodes[vb];

                      PointIndex v0 = parents[0];
                      PointIndex v1 = parents[1];
                      PointIndex v2 = f3[(k+1)%3];
                      PointIndex v3 = f3[(k+2)%3];
                      INT<2> parentedge1(v0, v2);
                      parentedge1.Sort();
                      INT<2> parentedge2(v0, v3);
                      parentedge2.Sort();
                      INT<2> parentedge3(v1, v2);
                      parentedge3.Sort();
                      INT<2> parentedge4(v1, v3);
                      parentedge4.Sort();

                      // if edges [v0,v2], [v0, v3], [v1,v2], [v1,v3] exists
                      // then vb is the bisecting edge
                      if (v2e.Used(parentedge1) && v2e.Used(parentedge2) 
                          && v2e.Used(parentedge3) && v2e.Used(parentedge4))
                        {
                          int verts[5] = { v0, v1, v2, v3, vb };
                          /*
                          cout << "verts5: ";
                          for (int j = 0; j < 5; j++)
                            cout << verts[j] << " ";
                          */
                          // classify permutation of verts
                          int classnr = 0;
                          for (int j = 0; j < 4; j++)
                            {
                              int maxk = 0;
                              for (int k = 0; k < 5-j; k++)
                                if (verts[k] > verts[maxk]) maxk = k;
                              // compress
                              for (int k = maxk; k < 4-j; k++)
                                verts[k] = verts[k+1];
                              classnr = maxk + (5-j) * classnr;
                            }
                          // cout << "classnr = " << classnr << endl;

                          INT<3> parentverts1(v1, v2, v3);
                          parentverts1.Sort();
                          INT<3> parentverts2(v0, v2, v3);
                          parentverts2.Sort();
                          INT<3> parentverts3(v0, v1, v3);
                          parentverts3.Sort();
                          INT<3> parentverts4(v0, v1, v2);
                          parentverts4.Sort();
                          
                          if (!v2f.Used(parentverts1) || !v2f.Used(parentverts2) ||
                              !v2f.Used(parentverts3) || !v2f.Used(parentverts4))
                            {
                              cout << "all edges are used, but not faces ????" << endl;
                              continue;
                            }
                                        
                          int pafacenr1 = v2f[parentverts1];
                          int pafacenr2 = v2f[parentverts2];
                          int pafacenr3 = v2f[parentverts3];
                          int pafacenr4 = v2f[parentverts4];

                          
                          parent_faces[i] = { classnr, { pafacenr1, pafacenr2,
                                                         pafacenr3, pafacenr4} };
                          
                          break;
                        }
                    }

                auto [info, nrs] = parent_faces[i];
                if (nrs[0] == -1){
                  // hacking for tet red refinements
                  PointIndex v0 = f3[0];
                  auto pa0 = mesh->mlbetweennodes[v0];
                  auto pa1 = mesh->mlbetweennodes[f3[1]];
                  auto pa2 = mesh->mlbetweennodes[f3[2]];
                  // v0 is a coarse vertex ==> f3 is a boundary face
                  if (v0==pa1[0] || v0==pa1[1]){
                    if (pa1[0]==v0){// type 0: bottom left corner
                      INT<3> parentverts(v0, pa1[1], pa2[1]);
                      int pafacenr = v2f[parentverts];
                      parent_faces[i] = { 16, { pafacenr, -1, -1, -1} };
                      //cout << "f "<<i<<":pf "<< pafacenr<< "A" <<endl;
                    }else if (pa2[0]==v0) {// type 1: bottom right corner
                      INT<3> parentverts(pa1[0], v0, pa2[1]);
                      int pafacenr = v2f[parentverts];
                      parent_faces[i] = { 17, { pafacenr, -1, -1, -1} };
                      //cout << "f "<<i<<":pf "<< pafacenr<< "B" <<endl;
                    }else if (pa1[1]==v0){// type 2: top left corner
                      INT<3> parentverts(pa1[0], pa2[0], v0);
                      int pafacenr = v2f[parentverts];
                      parent_faces[i] = { 18, { pafacenr, -1, -1, -1} };
                      //cout << "f "<<i<<":pf "<< pafacenr<< "C" <<endl;
                    }else{
                      cout << "************************** unhandled parent-face case **********************" << endl;
                    }
                  }
                  else{// all vertices are on fine level [fff]
                    // Here we only work with boundary fff face
                    if (pa0[0]==pa1[0] && pa0[1]==pa2[0] && pa1[1]==pa2[1]){//type 3 bdry face
                      INT<3> parentverts(pa0[0], pa0[1], pa1[1]);
                      int pafacenr = v2f[parentverts];
                      parent_faces[i] = { 19, { pafacenr, -1, -1, -1} };
                      //cout << "f "<<i<<":pf "<< pafacenr<< "D" <<endl;
                    }else{// this is an interior face FIXME 
                      parent_faces[i] = { 20, { -1, -1, -1, -1} };
                      //cout << "face "<< i << ":"<< f3 <<" is an int face"<< endl;
                    }
                  }
                }
              }
          }
      }
    

#ifdef PARALLEL
    if (id != 0)  
      {
	// if ( paralleltop.DoCoarseUpdate() )
	// paralleltop.UpdateCoarseGrid();
      }
#endif
 
 
  
    /* 
       for (i = 1; i <= ne; i++)
       {
       (*testout) << "Element " << i << endl;
       (*testout) << "PNums " << endl; 
       for( int l=1;l<=8;l++) *testout << mesh.VolumeElement(i).PNum(l) << "\t"; 
       *testout << endl; 
       (*testout) << "edges: " << endl;
       for (j = 0; j < 9; j++)
       (*testout) << edges.Elem(i)[j] << " ";
       (*testout) << "faces: " << endl;
       for (j = 0; j < 6; j++)m
       (*testout) << faces.Elem(i)[j] << " ";
       }

       for (i = 1; i <= nse; i++)
       {
       (*testout) << "SElement " << i << endl;
       (*testout) << "PNums " << endl; 
       for( int l=1;l<=4;l++) *testout << mesh.SurfaceElement(i).PNum(l) << "\t"; 
       *testout << endl; 
       }
    */
    timestamp = NextTimeStamp();
  }

  



  const Point3d * MeshTopology :: GetVertices (ELEMENT_TYPE et)
  {
    static Point3d segm_points [] = 
      { Point3d (1, 0, 0),
	Point3d (0, 0, 0) };
  
    static Point3d trig_points [] = 
      { Point3d ( 1, 0, 0 ),
	Point3d ( 0, 1, 0 ),
	Point3d ( 0, 0, 0 ) };

    static Point3d quad_points [] = 
      { Point3d ( 0, 0, 0 ),
	Point3d ( 1, 0, 0 ),
	Point3d ( 1, 1, 0 ),
	Point3d ( 0, 1, 0 ) };

    static Point3d tet_points [] = 
      { Point3d ( 1, 0, 0 ),
	Point3d ( 0, 1, 0 ),
	Point3d ( 0, 0, 1 ),
	Point3d ( 0, 0, 0 ) };

    static Point3d pyramid_points [] =
      {
	Point3d ( 0, 0, 0 ),
	Point3d ( 1, 0, 0 ),
	Point3d ( 1, 1, 0 ),
	Point3d ( 0, 1, 0 ),
	Point3d ( 0, 0, 1-1e-7 ),
      };    
  
    static Point3d prism_points[] = 
      {
	Point3d ( 1, 0, 0 ),
	Point3d ( 0, 1, 0 ),
	Point3d ( 0, 0, 0 ),
	Point3d ( 1, 0, 1 ),
	Point3d ( 0, 1, 1 ),
	Point3d ( 0, 0, 1 )
      };


    static Point3d hex_points [] = 
      { Point3d ( 0, 0, 0 ),
	Point3d ( 1, 0, 0 ),
	Point3d ( 1, 1, 0 ),
	Point3d ( 0, 1, 0 ),
	Point3d ( 0, 0, 1 ),
	Point3d ( 1, 0, 1 ),
	Point3d ( 1, 1, 1 ),
	Point3d ( 0, 1, 1 ) };


    switch (et)
      {
      case SEGMENT:
      case SEGMENT3:
	return segm_points;

      case TRIG:
      case TRIG6:
	return trig_points;

      case QUAD:
      case QUAD6:
      case QUAD8:
	return quad_points;

      case TET:
      case TET10:
	return tet_points;

      case PYRAMID:
	return pyramid_points;

      case PRISM:
      case PRISM12:
	return prism_points;

      case HEX:
	return hex_points;
      default:
	cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
      }
    return 0;
  }








  void MeshTopology :: GetElementEdges (int elnr, NgArray<int> & eledges) const
  {
    int ned = GetNEdges (mesh->VolumeElement(elnr).GetType());
    eledges.SetSize (ned);
    for (int i = 0; i < ned; i++)
      eledges[i] = edges.Get(elnr)[i]+1;
      // eledges[i] = abs (edges.Get(elnr)[i]);
  }
  void MeshTopology :: GetElementFaces (int elnr, NgArray<int> & elfaces, bool withorientation) const
  {
    int nfa = GetNFaces (mesh->VolumeElement(elnr).GetType());
    elfaces.SetSize (nfa);

    for (auto i : Range(nfa))
      elfaces[i] = faces.Get(elnr)[i]+1;
    
    if(withorientation)
    {
        for(auto & face : elfaces)
        {
            auto v = face2vert[face-1];
            if(v[3]!=0)
                cerr << "GetElementFaces with orientation currently not supported for quads" << endl;

            int classnr = 0;
            if (v[0] > v[1]) { classnr++; }
            if (v[1] > v[2]) { classnr++; }
            if (v[2] > v[0]) { classnr++; }

            if(classnr==1)
                face = -face;
        }
    }
  }

  void MeshTopology :: GetElementEdgeOrientations (int elnr, NgArray<int> & eorient) const
  {
    int ned = GetNEdges (mesh->VolumeElement(elnr).GetType());
    eorient.SetSize (ned);
    for (int i = 1; i <= ned; i++)
      // eorient.Elem(i) = (edges.Get(elnr)[i-1] > 0) ? 1 : -1;
      // eorient.Elem(i) = (edges.Get(elnr)[i-1].orient) ? -1 : 1;
      eorient.Elem(i) = GetElementEdgeOrientation (elnr, i-1) ? -1 : 1;
  }

  void MeshTopology :: GetElementFaceOrientations (int elnr, NgArray<int> & forient) const
  {
    int nfa = GetNFaces (mesh->VolumeElement(elnr).GetType());
    forient.SetSize (nfa);
    for (int i = 1; i <= nfa; i++)
      // forient.Elem(i) = faces.Get(elnr)[i-1].forient;
      // forient.Elem(i) = (faces.Get(elnr)[i-1]-1) % 8;
      forient.Elem(i) = GetElementFaceOrientation(elnr, i-1);
  }


  
  int MeshTopology :: GetElementEdges (int elnr, int * eledges, int * orient) const
  {
    //  int ned = GetNEdges (mesh.VolumeElement(elnr).GetType());
    
    if (mesh->GetDimension()==3 || 1)
      {
        if (orient)
	  {
	    for (int i = 0; i < 12; i++)
	      {
                /*
		if (!edges.Get(elnr)[i]) return i;
		eledges[i] = abs (edges.Get(elnr)[i]);
		orient[i] = (edges.Get(elnr)[i] > 0 ) ? 1 : -1;
                */
                if (edges.Get(elnr)[i] == -1) return i;
                eledges[i] = edges.Get(elnr)[i]+1;
		// orient[i] = edges.Get(elnr)[i].orient ? -1 : 1;
                orient[i] = GetElementEdgeOrientation(elnr, i) ? -1 : 1;
	      }
	  }
	else
	  {
	    for (int i = 0; i < 12; i++)
	      {
		// if (!edges.Get(elnr)[i]) return i;
		// eledges[i] = abs (edges.Get(elnr)[i]);
                if (edges.Get(elnr)[i] == -1) return i;
                eledges[i] = edges.Get(elnr)[i]+1;

	      }
	  }
	return 12;
      }
    else
      {
	throw NgException("rethink implementation");
	/*
	  if (orient)
	  {
	  for (i = 0; i < 4; i++)
	  {
	  if (!surfedges.Get(elnr)[i]) return i;
	  eledges[i] = abs (surfedges.Get(elnr)[i]);
	  orient[i] = (surfedges.Get(elnr)[i] > 0 ) ? 1 : -1;
	  }
	  }
	  else
	  {
	  if (!surfedges.Get(elnr)[i]) return i;
	  for (i = 0; i < 4; i++)
	  eledges[i] = abs (surfedges.Get(elnr)[i]);
	  }
	*/
	return 4;
	//      return GetSurfaceElementEdges (elnr, eledges, orient);
      }
  }

  int MeshTopology :: GetElementFaces (int elnr, int * elfaces, int * orient) const
  {
    //  int nfa = GetNFaces (mesh.VolumeElement(elnr).GetType());
    if (orient)
      {
	for (int i = 0; i < 6; i++)
	  {
            /*
	    if (!faces.Get(elnr)[i]) return i;
	    elfaces[i] = (faces.Get(elnr)[i]-1) / 8 + 1;
	    orient[i] = (faces.Get(elnr)[i]-1) % 8;
            */
	    if (faces.Get(elnr)[i] == -1) return i;
	    elfaces[i] = faces.Get(elnr)[i]+1;
	    // orient[i] = faces.Get(elnr)[i].forient;
            orient[i] = GetElementFaceOrientation (elnr, i);
	  }
      }
    else
      {
	for (int i = 0; i < 6; i++)
	  {
	    // if (!faces.Get(elnr)[i]) return i;
	    // elfaces[i] = (faces.Get(elnr)[i]-1) / 8 + 1;
	    if (faces.Get(elnr)[i] == -1) return i;
	    elfaces[i] = faces.Get(elnr)[i]+1;
	  }
      }
    return 6;
  }

  
  void MeshTopology :: GetSurfaceElementEdges (int elnr, NgArray<int> & eledges) const
  {
    int ned = GetNEdges (mesh->SurfaceElement(elnr).GetType());
    eledges.SetSize (ned);
    for (int i = 0; i < ned; i++)
      eledges[i] = surfedges.Get(elnr)[i]+1;
  }

  void MeshTopology :: GetEdges (SurfaceElementIndex elnr, NgArray<int> & eledges) const
  {
    int ned = GetNEdges ( (*mesh)[elnr].GetType());
    eledges.SetSize (ned);
    for (int i = 0; i < ned; i++)
      eledges[i] = surfedges[elnr][i];
  }

  /*
  FlatArray<T_EDGE> MeshTopology :: GetEdges (SurfaceElementIndex elnr) const
  {
    return FlatArray<T_EDGE>(GetNEdges ( (*mesh)[elnr].GetType()), &surfedges[elnr][0]);
  }

  FlatArray<T_EDGE> MeshTopology :: GetEdges (ElementIndex elnr) const
  {
    return FlatArray<T_EDGE>(GetNEdges ( (*mesh)[elnr].GetType()), &edges[elnr][0]);
  }

  FlatArray<T_FACE> MeshTopology :: GetFaces (ElementIndex elnr) const
  {
    return FlatArray<T_FACE>(GetNFaces ( (*mesh)[elnr].GetType()), &faces[elnr][0]);
  }
  */
  
  
  int MeshTopology :: GetSurfaceElementFace (int elnr) const
  {
    return surffaces.Get(elnr)+1;
  }
  
  /*
  int MeshTopology :: GetFace (SurfaceElementIndex elnr) const
  {
    return surffaces[elnr].fnr;
  }
  */


  void MeshTopology :: 
  GetSurfaceElementEdgeOrientations (int elnr, NgArray<int> & eorient) const
  {
    int ned = GetNEdges (mesh->SurfaceElement(elnr).GetType());
    eorient.SetSize (ned);
    for (int i = 0; i < ned; i++)
      // eorient[i] = (surfedges.Get(elnr)[i] > 0) ? 1 : -1;
      // eorient[i] = (surfedges.Get(elnr)[i].orient) ? -1 : 1;
      eorient[i] = GetSurfaceElementEdgeOrientation(elnr, i) ? -1 : 1;
  }

  int MeshTopology :: GetSurfaceElementFaceOrientation (int elnr) const
  {
    // return (surffaces.Get(elnr)-1) % 8;
    // return surffaces.Get(elnr).forient;
    return GetSurfaceElementFaceOrientation2(elnr);
  }

  int MeshTopology :: GetSurfaceElementEdges (int elnr, int * eledges, int * orient) const
  {
    int i;
    if (mesh->GetDimension() == 3 || 1)
      {
	if (orient)
	  {
	    for (i = 0; i < 4; i++)
	      {
                /*
		if (!surfedges.Get(elnr)[i]) return i;
		eledges[i] = abs (surfedges.Get(elnr)[i]);
		orient[i] = (surfedges.Get(elnr)[i] > 0 ) ? 1 : -1;
                */
		if (surfedges.Get(elnr)[i] == -1) return i;
		eledges[i] = surfedges.Get(elnr)[i]+1;
		// orient[i] = (surfedges.Get(elnr)[i].orient) ? -1 : 1;
                // orient[i] = GetSurfaceElementEdgeOrientation(elnr, i) ? -1 : 1;
                orient[i] = 1;

	      }
	  }
	else
	  {
	    for (i = 0; i < 4; i++)
	      {
                /*
		if (!surfedges.Get(elnr)[i]) return i;
		eledges[i] = abs (surfedges.Get(elnr)[i]);
                */
		if (surfedges.Get(elnr)[i] == -1) return i;
		eledges[i] = surfedges.Get(elnr)[i]+1;
	      }
	  }
	return 4;
      }
    else
      {
        /*
	eledges[0] = abs (segedges.Get(elnr));
	if (orient)
	  orient[0] = segedges.Get(elnr) > 0 ? 1 : -1;
        */
	eledges[0] = segedges.Get(elnr)+1;
	if (orient)
	  // orient[0] = segedges.Get(elnr).orient ? -1 : 1;
          // orient[0] = GetSegmentEdgeOrientation(elnr) ? -1 : 1;
          orient[0] = 1;
      }
    return 1;
  }


  int MeshTopology :: GetElementEdgeOrientation (int elnr, int locedgenr) const
  {
    const Element & el = mesh->VolumeElement (elnr);
    const ELEMENT_EDGE * eledges = MeshTopology::GetEdges0 (el.GetType());    

    int k = locedgenr;
    INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);
    int edgedir = (edge.I1() > edge.I2());
    return edgedir;
  }

  
  int MeshTopology :: GetElementFaceOrientation (int elnr, int locfacenr) const
  {
    const Element & el = mesh->VolumeElement (elnr);
    
    const ELEMENT_FACE * elfaces = MeshTopology::GetFaces0 (el.GetType());

    int j = locfacenr;
    if (elfaces[j][3] < 0)
      { // triangle
        INDEX_4 face(el[elfaces[j][0]], el[elfaces[j][1]], 
                     el[elfaces[j][2]], 0);
        
        int facedir = 0;
        if (face.I1() > face.I2())
          { swap (face.I1(), face.I2()); facedir += 1; }
        if (face.I2() > face.I3())
          { swap (face.I2(), face.I3()); facedir += 2; }
        if (face.I1() > face.I2())
          { swap (face.I1(), face.I2()); facedir += 4; }

        return facedir;
      }
    else
      {
        // quad
        // int facenum;
        INDEX_4 face4(el[elfaces[j][0]], el[elfaces[j][1]],
                      el[elfaces[j][2]], el[elfaces[j][3]]);
        
        int facedir = 0;
        if (min2 (face4.I1(), face4.I2()) > 
            min2 (face4.I4(), face4.I3())) 
          {  // z - flip
            facedir += 1; 
            swap (face4.I1(), face4.I4());
            swap (face4.I2(), face4.I3());
          }
        if (min2 (face4.I1(), face4.I4()) >
            min2 (face4.I2(), face4.I3())) 
          {  // x - flip
            facedir += 2; 
            swap (face4.I1(), face4.I2());
            swap (face4.I3(), face4.I4());
          }
        if (face4.I2() > face4.I4())
          {  // diagonal flip
            facedir += 4; 
            swap (face4.I2(), face4.I4());
          }
        
        return facedir;
      }
  }        
    

  
  int MeshTopology :: GetSurfaceElementEdgeOrientation (int elnr, int locedgenr) const
  {
    const Element2d & el = mesh->SurfaceElement (elnr);
    const ELEMENT_EDGE * eledges = MeshTopology::GetEdges0 (el.GetType());    

    int k = locedgenr;
    INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);
    int edgedir = (edge.I1() > edge.I2());
    return edgedir;
  }
  
  int MeshTopology :: GetSurfaceElementFaceOrientation2 (int elnr) const
  {
    const Element2d & el = mesh->SurfaceElement (elnr);
    
    const ELEMENT_FACE * elfaces = MeshTopology::GetFaces0 (el.GetType());

    int j = 0;
    if (elfaces[j][3] < 0)
      { // triangle
        INDEX_4 face(el[elfaces[j][0]], el[elfaces[j][1]], 
                     el[elfaces[j][2]], 0);
        
        int facedir = 0;
        if (face.I1() > face.I2())
          { swap (face.I1(), face.I2()); facedir += 1; }
        if (face.I2() > face.I3())
          { swap (face.I2(), face.I3()); facedir += 2; }
        if (face.I1() > face.I2())
          { swap (face.I1(), face.I2()); facedir += 4; }

        return facedir;
      }
    else
      {
        // quad
        // int facenum;
        INDEX_4 face4(el[elfaces[j][0]], el[elfaces[j][1]],
                      el[elfaces[j][2]], el[elfaces[j][3]]);
        
        int facedir = 0;
        if (min2 (face4.I1(), face4.I2()) > 
            min2 (face4.I4(), face4.I3())) 
          {  // z - flip
            facedir += 1; 
            swap (face4.I1(), face4.I4());
            swap (face4.I2(), face4.I3());
          }
        if (min2 (face4.I1(), face4.I4()) >
            min2 (face4.I2(), face4.I3())) 
          {  // x - flip
            facedir += 2; 
            swap (face4.I1(), face4.I2());
            swap (face4.I3(), face4.I4());
          }
        if (face4.I2() > face4.I4())
          {  // diagonal flip
            facedir += 4; 
            swap (face4.I2(), face4.I4());
          }
        
        return facedir;
      }
  }

  void MeshTopology :: GetSegmentEdge (int segnr, int & enr, int & orient) const
  {
    enr = segedges.Get(segnr)+1;
    orient = GetSegmentEdgeOrientation(segnr);
  }

  
  int MeshTopology :: GetSegmentEdgeOrientation (int elnr) const
  {
    const Segment & el = mesh->LineSegment (elnr);
    const ELEMENT_EDGE * eledges = MeshTopology::GetEdges0 (el.GetType());    

    int k = 0;
    INDEX_2 edge(el[eledges[k][0]], el[eledges[k][1]]);
    int edgedir = (edge.I1() > edge.I2());
    return edgedir;
  }


  
  void MeshTopology :: GetFaceVertices (int fnr, NgArray<int> & vertices) const
  {
    vertices.SetSize(4);
    for (int i = 0; i < 4; i++)
      vertices[i] = face2vert[fnr-1][i];
    if (vertices[3] == 0)
      vertices.SetSize(3);
  }

  void MeshTopology :: GetFaceVertices (int fnr, int * vertices) const
  {
    for (int i = 0; i <= 3; i++)
      vertices[i] = face2vert[fnr-1][i];
  }


  void MeshTopology :: GetEdgeVertices (int ednr, int & v1, int & v2) const
  {
    // cout << "id = " << id << "getedgevertices, ednr = " << ednr << ", ned = " << edge2vert.Size() << "&v1 = " << &v1 << endl;
    if (ednr < 1 || ednr > edge2vert.Size())
      cerr << "illegal edge nr: " << ednr << ", numedges = " << edge2vert.Size() 
	   << " id = " << id 
	   << endl;
    v1 = edge2vert[ednr-1][0];
    v2 = edge2vert[ednr-1][1];
  }

  void MeshTopology :: GetEdgeVertices (int ednr, PointIndex & v1, PointIndex & v2) const
  {
    v1 = edge2vert[ednr-1][0];
    v2 = edge2vert[ednr-1][1];
  }


  void MeshTopology :: GetFaceEdges (int fnr, NgArray<int> & fedges, bool withorientation) const
  {
    NgArrayMem<int,4> pi(4);
    NgArrayMem<int,12> eledges;
  
    fedges.SetSize (0);
    GetFaceVertices(fnr, pi);

    // Sort Edges according to global vertex numbers 
    // e1 = fmax, f2 
    // e2 = fmax, f1 
    // e3 = op e1(f2,f3) 
    // e4 = op e2(f1,f3) 

    /*  NgArrayMem<int,4> fp; 
	fp[0] = pi[0]; 
	for(int k=1;k<pi.Size();k++) 
	if(fp[k]>fp[0]) swap(fp[k],fp[0]); 
  
	fp[1] = fp[0]+ */ 
  

    //  GetVertexElements (pi[0], els);
    FlatArray<ElementIndex> els = GetVertexElements (pi[0]);

    // find one element having all vertices of the face
    for (int i = 0; i < els.Size(); i++)
      {
	const Element & el = (*mesh)[els[i]];
	int nref_faces = GetNFaces (el.GetType());
	const ELEMENT_FACE * ref_faces = GetFaces1 (el.GetType());
	int nfa_ref_edges = GetNEdges (GetFaceType(fnr));
      
	int cntv = 0,fa=-1; 
	for(int m=0;m<nref_faces;m++)
	  { 
	    cntv=0;
	    for(int j=0;j<nfa_ref_edges && ref_faces[m][j]>0;j++)
	      for(int k=0;k<pi.Size();k++)
		{
		  if(el[ref_faces[m][j]-1] == pi[k])
		    cntv++;
		}
	    if (cntv == pi.Size())
	      {
		fa=m;
		break;
	      }
	  }
     
	if(fa>=0)
	  {
	    const ELEMENT_EDGE * fa_ref_edges = GetEdges1 (GetFaceType(fnr)); 
	    fedges.SetSize(nfa_ref_edges);
	    GetElementEdges (els[i]+1, eledges);
	  
	    for (int j = 0; j < eledges.Size(); j++)
	      {
		int vi1, vi2;
		GetEdgeVertices (eledges[j], vi1, vi2);
	    
		bool has1 = 0;
		bool has2 = 0;
		for (int k = 0; k < pi.Size(); k++)
		  {
		    if (vi1 == pi[k]) has1 = 1;
		    if (vi2 == pi[k]) has2 = 1;
		  
		  }
	      
		if (has1 && has2) // eledges[j] is on face 
		  {
		    // fedges.Append (eledges[j]);
		    for(int k=0;k<nfa_ref_edges;k++)
		      {
			int w1 = el[ref_faces[fa][fa_ref_edges[k][0]-1]-1]; 
			int w2 = el[ref_faces[fa][fa_ref_edges[k][1]-1]-1]; 

			if(withorientation)
			  {
			    if(w1==vi1 && w2==vi2)
			      fedges[k] = eledges[j];
			    if(w1==vi2 && w2==vi1)
			      fedges[k] = -eledges[j];
			  }
			else
			  if((w1==vi1 && w2==vi2) || (w1==vi2 && w2==vi1))
			    fedges[k] = eledges[j];
		      }
		  }
	      }
	  
	    // *testout << " Face " << fnr << endl; 
	    // *testout << " GetFaceEdges " << fedges << endl;
	  
	    return;
	  }
      }   

    int surfel = GetFace2SurfaceElement(fnr);
    if (surfel != 0)
      {
	GetSurfaceElementEdges (surfel, fedges);
	return;
      }
  }


  void MeshTopology :: GetVertexElements (int vnr, Array<ElementIndex> & elements) const
  {
    if (vert2element.Size())
      elements = vert2element[vnr];
  }

  void MeshTopology :: GetVertexSurfaceElements( int vnr, 
						 Array<SurfaceElementIndex> & elements ) const
  {
    if (vert2surfelement.Size())
      elements = vert2surfelement[vnr];
  }


  int MeshTopology :: GetVerticesEdge ( int v1, int v2 ) const
  {
    // NgArray<int> elementedges;
    // Array<ElementIndex> elements_v1;
    // GetVertexElements ( v1, elements_v1);
    auto elements_v1 = GetVertexElements ( v1 );
    int edv1, edv2;

    for ( int i = 0; i < elements_v1.Size(); i++ )
      {
	// GetElementEdges( elements_v1[i]+1, elementedges );
        auto elementedges = GetEdges(ElementIndex(elements_v1[i]));
	for ( int ed = 0; ed < elementedges.Size(); ed ++)
	  {
	    GetEdgeVertices( elementedges[ed]+1, edv1, edv2 );
	    if ( ( edv1 == v1 && edv2 == v2 ) || ( edv1 == v2 && edv2 == v1 ) )
	      return elementedges[ed];
	  }
      }

    return -1;
  }



  void MeshTopology :: 
  GetSegmentVolumeElements ( int segnr, NgArray<ElementIndex> & volels ) const
  {
    /*
    int v1, v2;
    // GetEdgeVertices ( GetSegmentEdge (segnr), v1, v2 );
    GetEdgeVertices ( GetEdge (segnr-1)+1, v1, v2 );
    */
    auto [v1,v2] = GetEdgeVertices ( GetEdge (segnr-1) );
    auto volels1 = GetVertexElements ( v1 );
    auto volels2 = GetVertexElements ( v2 );
    volels.SetSize(0);

    for ( auto volel1 : volels1 )
      if ( volels2.Contains( volel1 ) )
	volels.Append ( volel1 );
  }

  void MeshTopology :: 
  GetSegmentSurfaceElements (int segnr, NgArray<SurfaceElementIndex> & els) const
  {
    int v1, v2;
    // GetEdgeVertices ( GetSegmentEdge (segnr), v1, v2 );
    GetEdgeVertices ( GetEdge (segnr-1)+1, v1, v2 );
    auto els1 = GetVertexSurfaceElements ( v1 );
    auto els2 = GetVertexSurfaceElements ( v2 );
    els.SetSize(0);

    for ( auto el1 : els1 )
      if ( els2.Contains( el1 ) )
	els.Append ( el1 );
  }




}

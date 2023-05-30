#include <mystdlib.h>

#include "meshing.hpp"

namespace netgen
{

  AnisotropicClusters ::  AnisotropicClusters (const Mesh & amesh)
    : mesh(amesh)
  {
    ;
  }

  AnisotropicClusters ::  ~AnisotropicClusters ()
  {
    ;
  }
  
  void AnisotropicClusters ::  Update()
  {
    static Timer timer("clusters");
    // static int timer1 = NgProfiler::CreateTimer ("clusters1");
    // static int timer2 = NgProfiler::CreateTimer ("clusters2");
    // static int timer3 = NgProfiler::CreateTimer ("clusters3");
    RegionTimer reg (timer);

    const MeshTopology & top = mesh.GetTopology();

    auto id = this->mesh.GetCommunicator().Rank();
    // auto ntasks = this->mesh.GetCommunicator().Size();

    bool hasedges = top.HasEdges();
    bool hasfaces = top.HasFaces();

    if (!hasedges || !hasfaces) return;

    if (id == 0)
      PrintMessage (3, "Update clusters");

    nv = mesh.GetNV();
    ned = top.GetNEdges();
    nfa = top.GetNFaces();
    ne = mesh.GetNE();
    int nse = mesh.GetNSE();

    cluster_reps.SetSize (nv+ned+nfa+ne);
    cluster_reps = -1;

    NgArray<int> llist (nv+ned+nfa+ne);
    llist = 0;
  
    NgArray<int> nnums, ednums, fanums;
    int changed;

    // NgProfiler::StartTimer(timer1);    

  
    /*
    for (int i = 1; i <= ne; i++)
      {
	const Element & el = mesh.VolumeElement(i);
	ELEMENT_TYPE typ = el.GetType();
      
	top.GetElementEdges (i, ednums);
	top.GetElementFaces (i, fanums);
      
	int elnv = top.GetNVertices (typ);
	int elned = ednums.Size();
	int elnfa = fanums.Size();
	  
	nnums.SetSize(elnv+elned+elnfa+1);
	for (int j = 1; j <= elnv; j++)
	  nnums.Elem(j) = el.PNum(j);
	for (int j = 1; j <= elned; j++)
	  nnums.Elem(elnv+j) = nv+ednums.Elem(j);
	for (int j = 1; j <= elnfa; j++)
	  nnums.Elem(elnv+elned+j) = nv+ned+fanums.Elem(j);
	nnums.Elem(elnv+elned+elnfa+1) = nv+ned+nfa+i;

	for (int j = 0; j < nnums.Size(); j++)
	  cluster_reps.Elem(nnums[j]) = nnums[j];
      }
    */
    ngcore::ParallelForRange
      (mesh.VolumeElements().Range(),
       [&] (auto myrange)
       {
         NgArray<int> nnums; // , ednums, fanums;
         for (int i_ : myrange)
           {
             int i = i_+1;
             const Element & el = mesh.VolumeElement(i);
             ELEMENT_TYPE typ = el.GetType();
             
             // top.GetElementEdges (i, ednums);
             auto ednums = top.GetEdges (ElementIndex(i_));
             // top.GetElementFaces (i, fanums);
             auto fanums = top.GetFaces (ElementIndex(i_));
             
             int elnv = top.GetNVertices (typ);
             int elned = ednums.Size();
             int elnfa = fanums.Size();
             
             nnums.SetSize(elnv+elned+elnfa+1);
             for (int j = 0; j < elnv; j++)
               nnums[j] = el[j]+1-PointIndex::BASE;
             for (int j = 0; j < elned; j++)
               nnums[elnv+j] = nv+ednums[j]+1;
             for (int j = 0; j < elnfa; j++)
               nnums[elnv+elned+j] = nv+ned+fanums[j]+1;
             nnums[elnv+elned+elnfa] = nv+ned+nfa+i;

             for (int j = 0; j < nnums.Size(); j++)
               cluster_reps.Elem(nnums[j]) = nnums[j];
           }
       }, ngcore::TasksPerThread(4));
    
    // NgProfiler::StopTimer(timer1);
    // NgProfiler::StartTimer(timer2);      
    /*
    for (int i = 1; i <= nse; i++)
      {
	const Element2d & el = mesh.SurfaceElement(i);
	ELEMENT_TYPE typ = el.GetType();
      
	top.GetSurfaceElementEdges (i, ednums);
	int fanum = top.GetSurfaceElementFace (i);
      
	int elnv = top.GetNVertices (typ);
	int elned = ednums.Size();
	  
	nnums.SetSize(elnv+elned+1);
	for (int j = 1; j <= elnv; j++)
	  nnums.Elem(j) = el.PNum(j)+1-PointIndex::BASE;
	for (int j = 1; j <= elned; j++)
	  nnums.Elem(elnv+j) = nv+ednums.Elem(j);
	nnums.Elem(elnv+elned+1) = fanum;

	for (int j = 0; j < nnums.Size(); j++)
	  cluster_reps.Elem(nnums[j]) = nnums[j];
      }
    */
    ngcore::ParallelForRange
      (mesh.SurfaceElements().Range(),
       [&] (auto myrange)
       {
         NgArrayMem<int,9> nnums; // , ednums;
         for (int i_ : myrange)
           {
             int i = i_+1;
             const Element2d & el = mesh.SurfaceElement(i);
             ELEMENT_TYPE typ = el.GetType();
             
             // top.GetSurfaceElementEdges (i, ednums);
             auto ednums = top.GetEdges (SurfaceElementIndex(i_));
             // cout << "ednums = " << ednums << endl;
             
             int fanum = top.GetSurfaceElementFace (i);             
             
             int elnv = top.GetNVertices (typ);
             int elned = ednums.Size();
             
             nnums.SetSize(elnv+elned+1);
             for (int j = 0; j < elnv; j++)
               nnums[j] = el[j]+1-PointIndex::BASE;
             for (int j = 0; j < elned; j++)
               nnums[elnv+j] = nv+ednums[j]+1;
             nnums[elnv+elned] = fanum;             
             
             for (int j = 0; j < nnums.Size(); j++)
               cluster_reps.Elem(nnums[j]) = nnums[j];
           }
       }, ngcore::TasksPerThread(4));

    
    // NgProfiler::StopTimer(timer2);
    // NgProfiler::StartTimer(timer3);      

    
    static const int hex_cluster[] =
      { 
	1, 2, 3, 4, 1, 2, 3, 4, 
	5, 6, 7, 8, 5, 6, 7, 8, 1, 2, 3, 4, 
	9, 9, 5, 8, 6, 7, 
	9
      };

    static const int prism_cluster[] =
      { 
	1, 2, 3, 1, 2, 3,
	4, 5, 6, 4, 5, 6, 3, 1, 2,
	7, 7, 4, 5, 6, 
	7
      };

    static const int pyramid_cluster[] =
      { 
	1, 2, 2, 1, 3,
	4, 2, 1, 4, 6, 5, 5, 6,
	7, 5, 7, 6, 4, 
	7
      };
    static const int tet_cluster14[] =
      { 1, 2, 3, 1,   1, 4, 5, 4, 5, 6,   7, 5, 4, 7, 7 };
  
    static const int tet_cluster12[] =
      { 1, 1, 2, 3,   4, 4, 5, 1, 6, 6,   7, 7, 4, 6, 7 };

    static const int tet_cluster13[] =
      { 1, 2, 1, 3,   4, 6, 4, 5, 1, 5,   7, 4, 7, 5, 7 };

    static const int tet_cluster23[] =
      { 2, 1, 1, 3, 6, 5, 5, 4, 4, 1, 5, 7, 7, 4, 7 };

    static const int tet_cluster24[] =
      { 2, 1, 3, 1, 4, 1, 5, 4, 6, 5, 5, 7, 4, 7, 7 };

    static const int tet_cluster34[] =
      { 2, 3, 1, 1, 4, 5, 1, 6, 4, 5, 5, 4, 7, 7, 7 };

    int cnt = 0;

    do
      {
        static Timer t("update cluster, identify");
        RegionTimer rtr(t);
	cnt++;
	changed = 0;
      
	for (int i = 1; i <= ne; i++)
	  {
	    const Element & el = mesh.VolumeElement(i);
	    ELEMENT_TYPE typ = el.GetType();
	  	  
	    const int * clustertab = NULL;
	    switch (typ)
	      {
	      case PRISM:
	      case PRISM12:
		clustertab = prism_cluster;
		break;
	      case HEX: 
		clustertab = hex_cluster; 
		break; 
	      case PYRAMID:
		clustertab = pyramid_cluster;
		break;
	      case TET:
	      case TET10:
		if (cluster_reps.Get(el.PNum(1)+1-PointIndex::BASE) == 
		    cluster_reps.Get(el.PNum(2)+1-PointIndex::BASE))
		  clustertab = tet_cluster12;
		else if (cluster_reps.Get(el.PNum(1)+1-PointIndex::BASE) == 
			 cluster_reps.Get(el.PNum(3)+1-PointIndex::BASE))
		  clustertab = tet_cluster13;
		else if (cluster_reps.Get(el.PNum(1)+1-PointIndex::BASE) == 
			 cluster_reps.Get(el.PNum(4)+1-PointIndex::BASE))
		  clustertab = tet_cluster14;
		else if (cluster_reps.Get(el.PNum(2)+1-PointIndex::BASE) == 
			 cluster_reps.Get(el.PNum(3)+1-PointIndex::BASE))
		  clustertab = tet_cluster23;
		else if (cluster_reps.Get(el.PNum(2)+1-PointIndex::BASE) == 
			 cluster_reps.Get(el.PNum(4)+1-PointIndex::BASE))
		  clustertab = tet_cluster24;
		else if (cluster_reps.Get(el.PNum(3)+1-PointIndex::BASE) == 
			 cluster_reps.Get(el.PNum(4)+1-PointIndex::BASE))
		  clustertab = tet_cluster34;

		else
		  clustertab = NULL;
		break;
	      default:
		clustertab = NULL;
	      }
	  
	    if (clustertab)
              {
                // top.GetElementEdges (i, ednums);
                // top.GetElementFaces (i, fanums);
                auto ednums = top.GetEdges (ElementIndex(i-1));
                auto fanums = top.GetFaces (ElementIndex(i-1));                
                
                int elnv = top.GetNVertices (typ);
                int elned = ednums.Size();
                int elnfa = fanums.Size();
                
                nnums.SetSize(elnv+elned+elnfa+1);
                for (int j = 0; j < elnv; j++)
                  nnums[j] = el[j]+1-PointIndex::BASE;
                for (int j = 0; j < elned; j++)
                  nnums[elnv+j] = nv+ednums[j]+1;
                for (int j = 0; j < elnfa; j++)
                  nnums[elnv+elned+j] = nv+ned+fanums[j]+1;
                nnums[elnv+elned+elnfa] = nv+ned+nfa+i;
                
                
	      for (int j = 0; j < nnums.Size(); j++)
		for (int k = 0; k < j; k++)
		  if (clustertab[j] == clustertab[k])
		    {
		      int jj = nnums[j];
		      int kk = nnums[k];

		      if (cluster_reps.Get(kk) < cluster_reps.Get(jj))
			swap (jj,kk);

		      if (cluster_reps.Get(jj) < cluster_reps.Get(kk))
			{
			  /*
			  cluster_reps.Elem(kk) = cluster_reps.Get(jj);
			  changed = 1;
			  */
			  
			  int rep  = cluster_reps.Get(jj);
			  int next = cluster_reps.Get(kk);
			  do
			    {
			      int cur = next;
			      next = llist.Elem(next);
				
			      cluster_reps.Elem(cur) = rep;
			      llist.Elem(cur) = llist.Elem(rep);
			      llist.Elem(rep) = cur;
			    }
			  while (next);
			  changed = 1;
			}
		    }
              }

	    /*
	      if (clustertab)
	      {
	      if (typ == PYRAMID)
	      (*testout) << "pyramid";
	      else if (typ == PRISM || typ == PRISM12)
	      (*testout) << "prism";
	      else if (typ == TET || typ == TET10)
	      (*testout) << "tet";
	      else
	      (*testout) << "unknown type" << endl;
		
	      (*testout) << ", nnums  = ";
	      for (j = 0; j < nnums.Size(); j++)
	      (*testout) << "node " << j << " = " << nnums[j] << ", rep = "
	      << cluster_reps.Get(nnums[j]) << endl;
	      }
	    */
	  }
      }
    while (changed);
    // NgProfiler::StopTimer(timer3);
    /*
      (*testout) << "cluster reps:" << endl;
      for (i = 1; i <= cluster_reps.Size(); i++)
      {
      (*testout) << i << ": ";
      if (i <= nv)
      (*testout) << "v" << i << " ";
      else if (i <= nv+ned)
      (*testout) << "e" << i-nv << " ";
      else if (i <= nv+ned+nfa)
      (*testout) << "f" << i-nv-ned << " ";
      else
      (*testout) << "c" << i-nv-ned-nfa << " ";
      (*testout) << cluster_reps.Get(i) << endl;
      }
    */
  }
}

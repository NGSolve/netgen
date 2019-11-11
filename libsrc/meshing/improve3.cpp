#include <mystdlib.h>
#include <algorithm>

#include <core/taskmanager.hpp>
#include <core/logging.hpp>

#include "meshing.hpp"

#ifdef SOLIDGEOM
#include <csg.hpp>
#endif
#include <opti.hpp>

namespace netgen
{


// Calc badness of new element where pi1 and pi2 are replaced by pnew
double CalcBadReplacePoints (const Mesh::T_POINTS & points, const MeshingParameters & mp, const Element & elem, double h, PointIndex &pi1, PointIndex &pi2, MeshPoint &pnew)
  {
    if (elem.GetType() != TET) return 0;

    MeshPoint* p[] = {&points[elem[0]], &points[elem[1]], &points[elem[2]], &points[elem[3]]};

    for (auto i : Range(4))
        if(elem[i]==pi1 || elem[i]==pi2) p[i] = &pnew;

    return CalcTetBadness (*p[0], *p[1], *p[2], *p[3], h, mp);
  }

/*
  Combine two points to one.
  Set new point into the center, if both are
  inner points.
  Connect inner point to boundary point, if one
  point is inner point.
*/
double MeshOptimize3d :: CombineImproveEdge (Mesh & mesh,
                            const MeshingParameters & mp,
                            Table<ElementIndex, PointIndex> & elements_of_point,
                            Array<double> & elerrs,
                            PointIndex pi0, PointIndex pi1,
                            FlatArray<bool, PointIndex> is_point_removed,
                            bool check_only)
{
  if (pi1 < pi0) Swap (pi0, pi1);
  if(is_point_removed[pi0] || is_point_removed[pi1]) return false;

  MeshPoint p0 = mesh[pi0];
  MeshPoint p1 = mesh[pi1];

  if (p1.Type() != INNERPOINT)
      return false;

  ArrayMem<ElementIndex, 50> has_one_point;
  ArrayMem<ElementIndex, 50> has_both_points;

  for (auto ei : elements_of_point[pi0] )
  {
      Element & elem = mesh[ei];
      if (elem.IsDeleted()) return false;

      if (elem[0] == pi1 || elem[1] == pi1 || elem[2] == pi1 || elem[3] == pi1)
      {
          if(!has_both_points.Contains(ei))
            has_both_points.Append (ei);
      }
      else
      {
          if(!has_one_point.Contains(ei))
              has_one_point.Append (ei);
      }
  }

  for (auto ei : elements_of_point[pi1] )
  {
      Element & elem = mesh[ei];
      if (elem.IsDeleted()) return false;

      if (elem[0] == pi0 || elem[1] == pi0 || elem[2] == pi0 || elem[3] == pi0)
      {
          ;
      }
      else
      {
          if(!has_one_point.Contains(ei))
            has_one_point.Append (ei);
      }
  }

  double badness_old = 0.0;
  for (auto ei : has_one_point)
      badness_old += elerrs[ei];
  for (auto ei : has_both_points)
      badness_old += elerrs[ei];

  MeshPoint pnew = p0;
  if (p0.Type() == INNERPOINT)
      pnew = Center (p0, p1);

  ArrayMem<double, 50> one_point_badness(has_one_point.Size());

  double badness_new = 0;
  for (auto i : Range(has_one_point))
  {
      const Element & elem = mesh[has_one_point[i]];
      double badness = CalcBadReplacePoints (mesh.Points(), mp, elem, 0, pi0, pi1, pnew);
      badness_new += badness;
      one_point_badness[i] = badness;
  }

  // Check if changed tets are topologically legal
  if (p0.Type() != INNERPOINT)
  {
      for (auto ei : has_one_point)
      {
          Element elem = mesh[ei];
          int l;
          for (int l = 0; l < 4; l++)
              if (elem[l] == pi1)
              {
                  elem[l] = pi0;
                  break;
              }

          elem.flags.illegal_valid = 0;
          if (!mesh.LegalTet(elem))
              badness_new += 1e4;
      }
  }

  double d_badness = badness_new / has_one_point.Size() - badness_old / (has_one_point.Size()+has_both_points.Size());

  // Do the actual combine operation
  if (d_badness < 0.0 && !check_only)
  {
      is_point_removed[pi1] = true;
      mesh[pi0] = pnew;

      for (auto ei : elements_of_point[pi1])
      {
          Element & elem = mesh[ei];
          if (elem.IsDeleted()) continue;

          for (int l = 0; l < elem.GetNP(); l++)
              if (elem[l] == pi1)
                  elem[l] = pi0;

          elem.flags.illegal_valid = 0;
          if (!mesh.LegalTet (elem))
              (*testout) << "illegal tet " << ei << endl;
      }

      for (auto i : Range(has_one_point))
          elerrs[has_one_point[i]] = one_point_badness[i];

      for (auto ei : has_both_points)
      {
          mesh[ei].flags.illegal_valid = 0;
          mesh[ei].Delete();
      }
  }

  return d_badness;
}

void MeshOptimize3d :: CombineImproveSequential (Mesh & mesh,
				       OPTIMIZEGOAL goal)
{
  static Timer t("MeshOptimize3d::CombineImprove"); RegionTimer reg(t);
  
  int np = mesh.GetNP();
  int ne = mesh.GetNE();

  TABLE<ElementIndex, PointIndex::BASE> elementsonnode(np); 
  NgArray<ElementIndex> hasonepi, hasbothpi;

  NgArray<double> oneperr;
  NgArray<double> elerrs (ne);

  PrintMessage (3, "CombineImprove");
  (*testout)  << "Start CombineImprove" << "\n";

  //  mesh.CalcSurfacesOfNode ();
  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Combine Improve";


  double totalbad = 0;
  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      if(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh[ei].GetIndex())
	continue;
      double elerr = CalcBad (mesh.Points(), mesh[ei], 0);
      totalbad += elerr;
      elerrs[ei] = elerr;
    }

  if (goal == OPT_QUALITY)
    {
      totalbad = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << totalbad << endl;
      PrintMessage (5, "Total badness = ", totalbad);
    }

  for (ElementIndex ei = 0; ei < ne; ei++)
    if (!mesh[ei].IsDeleted())
      for (int j = 0; j < mesh[ei].GetNP(); j++)
	elementsonnode.Add (mesh[ei][j], ei);
  
  INDEX_2_HASHTABLE<int> edgetested (np+1);

  int cnt = 0;

  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      if(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
      	continue;
      if (multithread.terminate)
	break;
      
      multithread.percent = 100.0 * (ei+1) / ne;

      if (mesh.ElementType(ei) == FIXEDELEMENT)
	continue;

      for (int j = 0; j < 6; j++)
	{
	  Element & elemi = mesh[ei];
	  if (elemi.IsDeleted()) continue;
	  
	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  PointIndex pi1 = elemi[tetedges[j][0]];
	  PointIndex pi2 = elemi[tetedges[j][1]];

	  if (pi2 < pi1) Swap (pi1, pi2);
	  
	  INDEX_2 si2 (pi1, pi2);
	  si2.Sort();
	  
	  if (edgetested.Used (si2)) continue;
	  edgetested.Set (si2, 1);


	  // hasonepoint.SetSize(0);
	  //	  hasbothpoints.SetSize(0);
	  hasonepi.SetSize(0);
	  hasbothpi.SetSize(0);

	  NgFlatArray<ElementIndex> row1 = elementsonnode[pi1];
	  for (int k = 0; k < row1.Size(); k++)
	    {
	      Element & elem = mesh[row1[k]];
	      if (elem.IsDeleted()) continue;

	      if (elem[0] == pi2 || elem[1] == pi2 ||
		  elem[2] == pi2 || elem[3] == pi2)
		{
		  hasbothpi.Append (row1[k]);
		}
	      else
		{
		  hasonepi.Append (row1[k]);
		}
	    } 
	  
	  NgFlatArray<ElementIndex> row2 = elementsonnode[pi2];	  
	  for (int k = 0; k < row2.Size(); k++)
	    {
	      Element & elem = mesh[row2[k]];
	      if (elem.IsDeleted()) continue;

	      if (elem[0] == pi1 || elem[1] == pi1 ||
		  elem[2] == pi1 || elem[3] == pi1)
		;
	      else
		{
		  hasonepi.Append (row2[k]);
		}
	    } 
	  
	  double bad1 = 0;
	  for (int k = 0; k < hasonepi.Size(); k++)
	    bad1 += elerrs[hasonepi[k]];
	  for (int k = 0; k < hasbothpi.Size(); k++)
	    bad1 += elerrs[hasbothpi[k]];
	  
	  MeshPoint p1 = mesh[pi1];
	  MeshPoint p2 = mesh[pi2];
	  

	  // if (mesh.PointType(pi2) != INNERPOINT)
	  if (p2.Type() != INNERPOINT)
	    continue;
	  
	  MeshPoint pnew;
	  // if (mesh.PointType(pi1) != INNERPOINT)
	  if (p1.Type() != INNERPOINT)
	    pnew = p1;
	  else
	    pnew = Center (p1, p2);
	    
	  mesh[pi1] = pnew;
	  mesh[pi2] = pnew;

	  oneperr.SetSize (hasonepi.Size());

	  double bad2 = 0;
	  for (int k = 0; k < hasonepi.Size(); k++)
	    {
	      const Element & elem = mesh[hasonepi[k]];
	      double err = CalcBad (mesh.Points(), elem, 0);
	      // CalcTetBadness (mesh[elem[0]], mesh[elem[1]],  
	      // mesh[elem[2]], mesh[elem[3]], 0, mparam);
	      bad2 += err;
	      oneperr[k] = err;
	    }
	  
	  mesh[pi1] = p1;
	  mesh[pi2] = p2;

	  // if (mesh.PointType(pi1) != INNERPOINT)
	  if (p1.Type() != INNERPOINT)
	    {
	      for (int k = 0; k < hasonepi.Size(); k++)
		{
		  Element & elem = mesh[hasonepi[k]];
		  int l;
		  for (l = 0; l < 4; l++)
		    if (elem[l] == pi2)
		      {
			elem[l] = pi1;
			break;
		      }
		      
		  elem.flags.illegal_valid = 0;
		  if (!mesh.LegalTet(elem))
		    bad2 += 1e4;
		  
		  if (l < 4)
		    {
		      elem.flags.illegal_valid = 0;
		      elem[l] = pi2;
		    }
		}
	    }

	  if (bad2 / hasonepi.Size()  < 
	      bad1 / (hasonepi.Size()+hasbothpi.Size()))
	    {
	      mesh[pi1] = pnew;
	      cnt++;

	      NgFlatArray<ElementIndex> row = elementsonnode[pi2];
	      for (int k = 0; k < row.Size(); k++)
		{
		  Element & elem = mesh[row[k]];
		  if (elem.IsDeleted()) continue;

		  elementsonnode.Add (pi1, row[k]);
		  for (int l = 0; l < elem.GetNP(); l++)
		    if (elem[l] == pi2)
		      elem[l] = pi1;
		  
		  elem.flags.illegal_valid = 0;
		  if (!mesh.LegalTet (elem))
		    (*testout) << "illegal tet " << elementsonnode[pi2][k] << endl;
		}

	      for (int k = 0; k < hasonepi.Size(); k++)
		elerrs[hasonepi[k]] = oneperr[k];
	      
	      for (int k = 0; k < hasbothpi.Size(); k++)
		{
		  mesh[hasbothpi[k]].flags.illegal_valid = 0;
		  mesh[hasbothpi[k]].Delete();
		}
	    }
	}
    }

  mesh.Compress();
  mesh.MarkIllegalElements();

  PrintMessage (5, cnt, " elements combined");
  (*testout) << "CombineImprove done" << "\n";

  totalbad = 0;
  for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
    if(!(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex()))
      totalbad += CalcBad (mesh.Points(), mesh[ei], 0);

  if (goal == OPT_QUALITY)
    {
      totalbad = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << totalbad << endl;

      int cntill = 0;
      for (ElementIndex ei = 0; ei < ne; ei++)
	if(!(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex()))
	  if (!mesh.LegalTet (mesh[ei]))
	    cntill++;

      PrintMessage (5, cntill, " illegal tets");
    }
  multithread.task = savetask;
} 

void MeshOptimize3d :: CombineImprove (Mesh & mesh,
				       OPTIMIZEGOAL goal)
{
  static Timer t("MeshOptimize3d::CombineImprove"); RegionTimer reg(t);
  static Timer topt("Optimize");
  static Timer tsearch("Search");
  static Timer tbuild_elements_table("Build elements table");
  static Timer tbad("CalcBad");

  // return CombineImproveSequential(mesh, goal);

  mesh.BoundaryEdge (1,2); // ensure the boundary-elements table is built
  
  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int ntasks = 4*ngcore::TaskManager::GetNumThreads();

  Array<double> elerrs (ne);
  Array<bool, PointIndex> is_point_removed (np);
  is_point_removed = false;

  PrintMessage (3, "CombineImprove");
  (*testout)  << "Start CombineImprove" << "\n";

  //  mesh.CalcSurfacesOfNode ();
  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Combine Improve";


  tbad.Start();
  double totalbad = 0.0;
  ParallelForRange(Range(ne), [&] (auto myrange)
  {
    double totalbad_local = 0.0;
    for (ElementIndex ei : myrange)
    {
      if(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh[ei].GetIndex())
	continue;
      double elerr = CalcBad (mesh.Points(), mesh[ei], 0);
      totalbad_local += elerr;
      elerrs[ei] = elerr;
    }
    AtomicAdd(totalbad, totalbad_local);
  }, ntasks);
  tbad.Stop();

  if (goal == OPT_QUALITY)
    {
      totalbad = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << totalbad << endl;
    }

  auto elementsonnode = mesh.CreatePoint2ElementTable();

  Array<std::tuple<PointIndex,PointIndex>> edges;
  BuildEdgeList(mesh, elementsonnode, edges);

  // Find edges with improvement
  Array<std::tuple<double, int>> combine_candidate_edges(edges.Size());
  std::atomic<int> improvement_counter(0);

  tsearch.Start();
  ParallelForRange(Range(edges), [&] (auto myrange)
  {
    for(auto i : myrange)
    {
      auto [p0,p1] = edges[i];
      double d_badness = CombineImproveEdge (mesh, mp, elementsonnode, elerrs, p0, p1, is_point_removed, true);
      if(d_badness<0.0)
      {
        int index = improvement_counter++;
        combine_candidate_edges[index] = make_tuple(d_badness, i);
      }
    }
  }, ntasks);
  tsearch.Stop();

  auto edges_with_improvement = combine_candidate_edges.Part(0, improvement_counter.load());

  QuickSort(edges_with_improvement);
  PrintMessage(5, edges.Size(), " edges");
  PrintMessage(5, edges_with_improvement.Size(), " edges with improvement");

  // Apply actual optimizations
  topt.Start();
  int cnt = 0;
  for(auto [d_badness, ei] : edges_with_improvement)
  {
      auto [p0,p1] = edges[ei];
      if (CombineImproveEdge (mesh, mp, elementsonnode, elerrs, p0, p1, is_point_removed, false) < 0.0)
        cnt++;
  }
  topt.Stop();

  mesh.Compress();
  mesh.MarkIllegalElements();

  PrintMessage (5, cnt, " elements combined");
  (*testout) << "CombineImprove done" << "\n";

  if (goal == OPT_QUALITY)
    {
      totalbad = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << totalbad << endl;

      int cntill = 0;
      for (ElementIndex ei = 0; ei < ne; ei++)
	if(!(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex()))
	  if (!mesh.LegalTet (mesh[ei]))
	    cntill++;

      PrintMessage (5, cntill, " illegal tets");
    }
  multithread.task = savetask;
} 



double MeshOptimize3d :: SplitImproveEdge (Mesh & mesh, OPTIMIZEGOAL goal, Table<ElementIndex,PointIndex> & elementsonnode, Array<double> &elerrs, NgArray<INDEX_3> &locfaces, double badmax, PointIndex pi1, PointIndex pi2, PointIndex ptmp, bool check_only)
{
  double d_badness = 0.0;
  int cnt = 0;

  ArrayMem<ElementIndex, 20> hasbothpoints;

  if (mesh.BoundaryEdge (pi1, pi2)) return 0.0;

  for (ElementIndex ei : elementsonnode[pi1])
    {
      Element & el = mesh[ei];

      if(el.IsDeleted()) return 0.0;
      if (mesh[ei].GetType() != TET) return 0.0;

      bool has1 = el.PNums().Contains(pi1);
      bool has2 = el.PNums().Contains(pi2);

      if (has1 && has2)
          if (!hasbothpoints.Contains (ei))
              hasbothpoints.Append (ei);
    }

  if(mp.only3D_domain_nr)
      for(auto ei : hasbothpoints)
          if(mp.only3D_domain_nr != mesh[ei].GetIndex())
              return 0.0;

  if (goal == OPT_LEGAL)
    {
      bool all_tets_legal = true;
      for(auto ei : hasbothpoints)
          if( !mesh.LegalTet (mesh[ei]) || elerrs[ei] > 1e3)
              all_tets_legal = false;
      if(all_tets_legal)
          return 0.0;
    }

  double bad1 = 0.0;
  double bad1_max = 0.0;
  for (ElementIndex ei : hasbothpoints)
    {
      double bad = elerrs[ei];
      bad1 += bad;
      bad1_max = max(bad1_max, bad);
    }

  if(bad1_max < 100.0)
      return 0.0;

  bool puretet = 1;
  for (ElementIndex ei : hasbothpoints)
      if (mesh[ei].GetType() != TET)
          puretet = 0;
  if (!puretet) return 0.0;

  Point3d p1 = mesh[pi1];
  Point3d p2 = mesh[pi2];

  locfaces.SetSize(0);
  for (ElementIndex ei : hasbothpoints)
    {
      const Element & el = mesh[ei];

      for (int l = 0; l < 4; l++)
          if (el[l] == pi1 || el[l] == pi2)
            {
              INDEX_3 i3;
              Element2d face(TRIG);
              el.GetFace (l+1, face);
              for (int kk = 1; kk <= 3; kk++)
                  i3.I(kk) = face.PNum(kk);
              locfaces.Append (i3);
            }
    }

  PointFunction1 pf (mesh.Points(), locfaces, mp, -1);
  OptiParameters par;
  par.maxit_linsearch = 50;
  par.maxit_bfgs = 20;

  Point3d pnew = Center (p1, p2);
  Vector px(3);
  px(0) = pnew.X();
  px(1) = pnew.Y();
  px(2) = pnew.Z();

  if (bad1_max > 0.1 * badmax)
    {
      int pok = pf.Func (px) < 1e10;
      if (!pok)
          pok = FindInnerPoint (mesh.Points(), locfaces, pnew);

      if(pok)
        {
          px(0) = pnew.X();
          px(1) = pnew.Y();
          px(2) = pnew.Z();
          BFGS (px, pf, par);
          pnew.X() = px(0);
          pnew.Y() = px(1);
          pnew.Z() = px(2);
        }
    }

  double bad2 = pf.Func (px);

  mesh[ptmp] = Point<3>(pnew);

  for (int k = 0; k < hasbothpoints.Size(); k++)
    {
      Element & oldel = mesh[hasbothpoints[k]];
      Element newel1 = oldel;
      Element newel2 = oldel;

      oldel.flags.illegal_valid = 0;
      newel1.flags.illegal_valid = 0;
      newel2.flags.illegal_valid = 0;

      for (int l = 0; l < 4; l++)
        {
          if (newel1[l] == pi2) newel1[l] = ptmp;
          if (newel2[l] == pi1) newel2[l] = ptmp;
        }

      if (!mesh.LegalTet (oldel)) bad1 += 1e6;
      if (!mesh.LegalTet (newel1)) bad2 += 1e6;
      if (!mesh.LegalTet (newel2)) bad2 += 1e6;
    }

  d_badness = bad2-bad1;
  if(check_only)
      return d_badness;

  if (d_badness<0.0)
    {
      cnt++;

      PointIndex pinew = mesh.AddPoint (pnew);

      for (ElementIndex ei : hasbothpoints)
        {
          Element & oldel = mesh[ei];
          Element newel1 = oldel;
          Element newel2 = oldel;

          oldel.flags.illegal_valid = 0;
          oldel.Delete();

          newel1.flags.illegal_valid = 0;
          newel2.flags.illegal_valid = 0;

          for (int l = 0; l < 4; l++)
            {
              if (newel1[l] == pi2) newel1[l] = pinew;
              if (newel2[l] == pi1) newel2[l] = pinew;
            }

          mesh.AddVolumeElement (newel1);
          mesh.AddVolumeElement (newel2);
        }
    }
  return d_badness;
}

void MeshOptimize3d :: SplitImprove (Mesh & mesh,
				     OPTIMIZEGOAL goal)
{
  static Timer t("MeshOptimize3d::SplitImprove"); RegionTimer reg(t);
  static Timer topt("Optimize");
  static Timer tsearch("Search");

  // return SplitImproveSequential(mesh, goal);

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  double bad = 0.0;
  double badmax = 0.0;

  auto elementsonnode = mesh.CreatePoint2ElementTable();

  Array<double> elerrs(ne);

  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Split Improve";

  PrintMessage (3, "SplitImprove");
  (*testout)  << "start SplitImprove" << "\n";
  mesh.BoundaryEdge (1,2); // ensure the boundary-elements table is built

  ParallelFor( mesh.VolumeElements().Range(), [&] (ElementIndex ei) NETGEN_LAMBDA_INLINE
    {
      if(mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
          return;

      elerrs[ei] = CalcBad (mesh.Points(), mesh[ei], 0);
      bad += elerrs[ei];
      AtomicMax(badmax, elerrs[ei]);
    });

  if (goal == OPT_QUALITY)
    {
      bad = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad << endl;
    }

  Array<std::tuple<PointIndex,PointIndex>> edges;
  BuildEdgeList(mesh, elementsonnode, edges);

  // Find edges with improvement
  Array<std::tuple<double, int>> candidate_edges(edges.Size());
  std::atomic<int> improvement_counter(0);
  auto ptmp = mesh.AddPoint( {0,0,0} );

  tsearch.Start();
  ParallelForRange(Range(edges), [&] (auto myrange)
  {
    NgArray<INDEX_3> locfaces;

    for(auto i : myrange)
    {
      auto [p0,p1] = edges[i];
      double d_badness = SplitImproveEdge (mesh, goal, elementsonnode, elerrs, locfaces, badmax, p0, p1, ptmp, true);
      if(d_badness<0.0)
      {
        int index = improvement_counter++;
        candidate_edges[index] = make_tuple(d_badness, i);
      }
    }
  }, ngcore::TasksPerThread(4));
  tsearch.Stop();

  auto edges_with_improvement = candidate_edges.Part(0, improvement_counter.load());

  QuickSort(edges_with_improvement);
  PrintMessage(5, edges.Size(), " edges");
  PrintMessage(5, edges_with_improvement.Size(), " edges with improvement");

  // Apply actual optimizations
  topt.Start();
  int cnt = 0;
  NgArray<INDEX_3> locfaces;
  for(auto [d_badness, ei] : edges_with_improvement)
  {
      auto [p0,p1] = edges[ei];
      if (SplitImproveEdge (mesh, goal, elementsonnode, elerrs, locfaces, badmax, p0, p1, ptmp, false) < 0.0)
        cnt++;
  }
  topt.Stop();
  mesh.Compress();
  PrintMessage (5, cnt, " splits performed");
  (*testout) << "Splitt - Improve done" << "\n";

  if (goal == OPT_QUALITY)
    {
      bad = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad << endl;

      int cntill = 0;
      ne = mesh.GetNE();
      for (ElementIndex ei = 0; ei < ne; ei++)
        if (!mesh.LegalTet (mesh[ei]))
          cntill++;
      //      cout << cntill << " illegal tets" << endl;
    }

  multithread.task = savetask;
}


/*
  Mesh improvement by edge splitting.
  If mesh quality is improved by inserting a node into an inner edge,
  the edge is split into two parts.
*/
void MeshOptimize3d :: SplitImproveSequential (Mesh & mesh,
				     OPTIMIZEGOAL goal)
{
  static Timer t("MeshOptimize3d::SplitImprove"); RegionTimer reg(t);
  static Timer tloop("MeshOptimize3d::SplitImprove loop");
  
  double bad1, bad2, badmax, badlimit;
  int cnt = 0;
  int np = mesh.GetNP();
  int ne = mesh.GetNE();

  auto elementsonnode = mesh.CreatePoint2ElementTable();

  NgArray<ElementIndex> hasbothpoints;

  NgBitArray origpoint(np+1), boundp(np+1);  // big enough for 0 and 1-based
  origpoint.Set();

  NgArray<double> elerrs(ne);
  NgBitArray illegaltet(ne);
  illegaltet.Clear();

  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Split Improve";

  PrintMessage (3, "SplitImprove");
  (*testout)  << "start SplitImprove" << "\n";

  NgArray<INDEX_3> locfaces;

  INDEX_2_HASHTABLE<int> edgetested (np);

  bad1 = 0;
  badmax = 0;
  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      if(mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
      	continue;
      
      elerrs[ei] = CalcBad (mesh.Points(), mesh[ei], 0);
      bad1 += elerrs[ei];
      if (elerrs[ei] > badmax) badmax = elerrs[ei];
    }

  PrintMessage (5, "badmax = ", badmax);
  badlimit = 0.5 * badmax;

  boundp.Clear();
  for (auto & el : mesh.SurfaceElements())
    for (PointIndex pi : el.PNums())
      boundp.Set (pi);
  
  if (goal == OPT_QUALITY)
    {
      bad1 = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad1 << endl;
    }

  mesh.MarkIllegalElements();
  if (goal == OPT_QUALITY || goal == OPT_LEGAL)
    {
      int cntill = 0;
      for (ElementIndex ei = 0; ei < ne; ei++)
	{
          //	  if (!LegalTet (volelements.Get(i)))
          if (mesh[ei].flags.illegal)
            {
              cntill++;
              illegaltet.Set (ei);
            }
	}
    }

  tloop.Start();
  for (ElementIndex ei : mesh.VolumeElements().Range())    
    {
      Element & elem = mesh[ei];
      
      if(mp.only3D_domain_nr && mp.only3D_domain_nr != elem.GetIndex())
      	continue;
      if (multithread.terminate)
	break;

      multithread.percent = 100.0 * (ei+1) / ne;

      bool ltestmode = 0;

      if (elerrs[ei] < badlimit && !illegaltet.Test(ei)) continue;

      if ((goal == OPT_LEGAL) &&
	  !illegaltet.Test(ei) &&
	  CalcBad (mesh.Points(), elem, 0) < 1e3) 
	continue;

      if (ltestmode)
	{
	  (*testout) << "test el " << ei << endl;
	  for (int j = 0; j < 4; j++)
	    (*testout) << elem[j] << " ";
	  (*testout) << endl;
	}

      for (int j = 0; j < 6; j++)
	{
	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  PointIndex pi1 = elem[tetedges[j][0]];
	  PointIndex pi2 = elem[tetedges[j][1]];

	  if (pi2 < pi1) Swap (pi1, pi2);
	  if (pi2 >= elementsonnode.Size()+PointIndex::BASE) continue; // old number of points

	  if (!origpoint.Test(pi1) || !origpoint.Test(pi2))
	    continue;
  
	  INDEX_2 i2(pi1, pi2);
	  i2.Sort();

	  if (mesh.BoundaryEdge (pi1, pi2)) continue;

	  if (edgetested.Used (i2) && !illegaltet.Test(ei)) continue;
	  edgetested.Set (i2, 1);

	  hasbothpoints.SetSize (0);
          /*
	  for (int k = 1; k <= elementsonnode.EntrySize(pi1); k++)
	    {
	      ElementIndex elnr = elementsonnode.Get(pi1, k);
          */
          for (ElementIndex ei : elementsonnode[pi1])
            {
	      Element & el = mesh[ei];
              bool has1 = el.PNums().Contains(pi1);
              bool has2 = el.PNums().Contains(pi2);

	      if (has1 && has2) 
		if (!hasbothpoints.Contains (ei))
		  hasbothpoints.Append (ei);
	    }
	  
	  bad1 = 0;

          for (ElementIndex ei : hasbothpoints)
	    bad1 += CalcBad (mesh.Points(), mesh[ei], 0);            
          
	  bool puretet = 1;
	  for (ElementIndex ei : hasbothpoints)
	    if (mesh[ei].GetType() != TET)
	      puretet = 0;
	  if (!puretet) continue;

	  Point3d p1 = mesh[pi1];
	  Point3d p2 = mesh[pi2];

	  /*
	    pnew = Center (p1, p2);
	  
	    points.Elem(pi1) = pnew;
	    bad2 = 0;
	    for (k = 1; k <= hasbothpoints.Size(); k++)
	    bad2 += CalcBad (points, 
	    volelements.Get(hasbothpoints.Get(k)), 0);

	    points.Elem(pi1) = p1;
	    points.Elem(pi2) = pnew;
	      
	    for (k = 1; k <= hasbothpoints.Size(); k++)
	    bad2 += CalcBad (points, 
	    volelements.Get(hasbothpoints.Get(k)), 0);
	    points.Elem(pi2) = p2;
	  */

	  locfaces.SetSize (0);
	  for (ElementIndex ei : hasbothpoints)
	    {
	      const Element & el = mesh[ei];
              
	      for (int l = 0; l < 4; l++)
		if (el[l] == pi1 || el[l] == pi2)
		  {
		    INDEX_3 i3;
		    Element2d face(TRIG);
		    el.GetFace (l+1, face);
		    for (int kk = 1; kk <= 3; kk++)
		      i3.I(kk) = face.PNum(kk);
		    locfaces.Append (i3);
		  }
	    }

	  PointFunction1 pf (mesh.Points(), locfaces, mp, -1);
	  OptiParameters par;
	  par.maxit_linsearch = 50;
	  par.maxit_bfgs = 20;

	  Point3d pnew = Center (p1, p2);
	  Vector px(3);
	  px(0) = pnew.X();
	  px(1) = pnew.Y();
	  px(2) = pnew.Z();

	  if (elerrs[ei] > 0.1 * badmax)
	    BFGS (px, pf, par);

	  bad2 = pf.Func (px);

	  pnew.X() = px(0);
	  pnew.Y() = px(1);
	  pnew.Z() = px(2);

	  PointIndex hpinew = mesh.AddPoint (pnew);
	  //	  ptyps.Append (INNERPOINT);

	  for (int k = 0; k < hasbothpoints.Size(); k++)
	    {
	      Element & oldel = mesh[hasbothpoints[k]];
	      Element newel1 = oldel;
	      Element newel2 = oldel;

	      oldel.flags.illegal_valid = 0;
	      newel1.flags.illegal_valid = 0;
	      newel2.flags.illegal_valid = 0;
	      
	      for (int l = 0; l < 4; l++)
		{
		  if (newel1[l] == pi2) newel1[l] = hpinew;
		  if (newel2[l] == pi1) newel2[l] = hpinew;
		}
	      
	      if (!mesh.LegalTet (oldel)) bad1 += 1e6;
	      if (!mesh.LegalTet (newel1)) bad2 += 1e6;
	      if (!mesh.LegalTet (newel2)) bad2 += 1e6;
	    }	  
	  
	  // mesh.PointTypes().DeleteLast();
	  mesh.Points().DeleteLast();

	  if (bad2 < bad1) 
	    /*	      (bad1 > 1e4 && boundp.Test(pi1) && boundp.Test(pi2)) */ 
	    {
	      cnt++;

	      PointIndex pinew = mesh.AddPoint (pnew);
	      
	      for (ElementIndex ei : hasbothpoints)
		{
		  Element & oldel = mesh[ei];
		  Element newel = oldel;
                  
		  newel.flags.illegal_valid = 0;
		  oldel.flags.illegal_valid = 0;

		  for (int l = 0; l < 4; l++)
		    {
		      origpoint.Clear (oldel[l]);

		      if (oldel[l] == pi2) oldel[l] = pinew;
		      if (newel[l] == pi1) newel[l] = pinew;
		    }
		  mesh.AddVolumeElement (newel);
		}
	      
	      j = 10;  // end j-loop
	    }
	}
    }
  tloop.Stop();

  mesh.Compress();
  PrintMessage (5, cnt, " splits performed");

  (*testout) << "Splitt - Improve done" << "\n";

  if (goal == OPT_QUALITY)
    {
      bad1 = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad1 << endl;

      int cntill = 0;
      ne = mesh.GetNE();
      for (ElementIndex ei = 0; ei < ne; ei++)
        if (!mesh.LegalTet (mesh[ei]))
          cntill++;
      //      cout << cntill << " illegal tets" << endl;
    }

  multithread.task = savetask;
}

void MeshOptimize3d :: SwapImproveSequential (Mesh & mesh, OPTIMIZEGOAL goal,
				    const NgBitArray * working_elements)
{
  static Timer t("MeshOptimize3d::SwapImprove"); RegionTimer reg(t);
  static Timer tloop("MeshOptimize3d::SwapImprove loop");
  
  PointIndex pi3(PointIndex::INVALID), pi4(PointIndex::INVALID),
    pi5(PointIndex::INVALID), pi6(PointIndex::INVALID);
  int cnt = 0;

  Element el21(TET), el22(TET), el31(TET), el32(TET), el33(TET);
  Element el1(TET), el2(TET), el3(TET), el4(TET);
  Element el1b(TET), el2b(TET), el3b(TET), el4b(TET);
  
  double bad1, bad2, bad3;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  
  // contains at least all elements at node
  TABLE<ElementIndex,PointIndex::BASE> elementsonnode(np);

  NgArray<ElementIndex> hasbothpoints;

  PrintMessage (3, "SwapImprove ");
  (*testout) << "\n" << "Start SwapImprove" << endl;

  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Swap Improve";
  
  //  mesh.CalcSurfacesOfNode ();
    
  INDEX_3_HASHTABLE<int> faces(mesh.GetNOpenElements()/3 + 2);
  if (goal == OPT_CONFORM)
    {
      for (int i = 1; i <= mesh.GetNOpenElements(); i++)
	{
	  const Element2d & hel = mesh.OpenElement(i);
	  INDEX_3 face(hel[0], hel[1], hel[2]);
	  face.Sort();
	  faces.Set (face, 1);
	}
    }
  
  // Calculate total badness
  if (goal == OPT_QUALITY)
    {
      bad1 = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad1 << endl;
    }
  
  // find elements on node
  for (ElementIndex ei = 0; ei < ne; ei++)
    for (PointIndex pi : mesh[ei].PNums())
      elementsonnode.Add (pi, ei);
    /*
    for (int j = 0; j < mesh[ei].GetNP(); j++)
      elementsonnode.Add (mesh[ei][j], ei);
    */


  // INDEX_2_HASHTABLE<int> edgeused(2 * ne + 5);
  INDEX_2_CLOSED_HASHTABLE<int> edgeused(12 * ne + 5);

  tloop.Start();
  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      if (multithread.terminate)
	break;
      
      if (mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
      	continue;
      
      multithread.percent = 100.0 * (ei+1) / ne;

      if ((mesh.ElementType(ei)) == FIXEDELEMENT)
	continue;

      if(working_elements && 
	 ei < working_elements->Size() &&
	 !working_elements->Test(ei))
	continue;

      if (mesh[ei].IsDeleted())
	continue;

      if ((goal == OPT_LEGAL) && 
	  mesh.LegalTet (mesh[ei]) &&
	  CalcBad (mesh.Points(), mesh[ei], 0) < 1e3)
	continue;

      //      int onlybedges = 1;

      for (int j = 0; j < 6; j++)
	{
	  // loop over edges

	  const Element & elemi = mesh[ei];
	  if (elemi.IsDeleted()) continue;

	  int mattyp = elemi.GetIndex();
	  
	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  PointIndex pi1 = elemi[tetedges[j][0]];
	  PointIndex pi2 = elemi[tetedges[j][1]];

	  if (pi2 < pi1) Swap (pi1, pi2);
	 	  
	  if (mesh.BoundaryEdge (pi1, pi2)) continue;


	  INDEX_2 i2 (pi1, pi2);
	  i2.Sort();
	  if (edgeused.Used(i2)) continue;
	  edgeused.Set (i2, 1);
	  
	  hasbothpoints.SetSize (0);
	  for (int k = 0; k < elementsonnode[pi1].Size(); k++)
	    {
	      bool has1 = 0, has2 = 0;
	      ElementIndex elnr = elementsonnode[pi1][k];
	      const Element & elem = mesh[elnr];
	      
	      if (elem.IsDeleted()) continue;
	      
	      for (int l = 0; l < elem.GetNP(); l++)
		{
		  if (elem[l] == pi1) has1 = 1;
		  if (elem[l] == pi2) has2 = 1;
		}

	      if (has1 && has2) 
		{ // only once
                  if (hasbothpoints.Contains (elnr))
                    has1 = false;
                  
		  if (has1)
                    {
                      hasbothpoints.Append (elnr);
                    }
		}
	    }
	  
	  bool puretet = true;
	  for (ElementIndex ei : hasbothpoints)
	    if (mesh[ei].GetType () != TET)
	      puretet = false;
          
	  if (!puretet) continue;

	  int nsuround = hasbothpoints.Size();

	  if ( nsuround == 3 )
	    {
	      Element & elem = mesh[hasbothpoints[0]];
	      for (int l = 0; l < 4; l++)
		if (elem[l] != pi1 && elem[l] != pi2)
		  {
		    pi4 = pi3;
		    pi3 = elem[l];
		  }
	      
	      el31[0] = pi1;
	      el31[1] = pi2;
	      el31[2] = pi3;
	      el31[3] = pi4;
	      el31.SetIndex (mattyp);
	      
	      if (WrongOrientation (mesh.Points(), el31))
		{
		  Swap (pi3, pi4);
		  el31[2] = pi3;
		  el31[3] = pi4;
		}
	      
	      pi5.Invalidate();
	      for (int k = 0; k < 3; k++)   // JS, 201212
		{
		  const Element & elemk = mesh[hasbothpoints[k]];
		  bool has1 = false;
		  for (int l = 0; l < 4; l++)
		    if (elemk[l] == pi4)
		      has1 = true;
		  if (has1)
		    {
		      for (int l = 0; l < 4; l++)
			if (elemk[l] != pi1 && elemk[l] != pi2 && elemk[l] != pi4)
			  pi5 = elemk[l];
		    }
		}

	      if (!pi5.IsValid())
                throw NgException("Illegal state observed in SwapImprove");

	      
	      el32[0] = pi1;  
	      el32[1] = pi2;  
	      el32[2] = pi4;  
	      el32[3] = pi5;  
	      el32.SetIndex (mattyp);
	      
	      el33[0] = pi1;  
	      el33[1] = pi2;  
	      el33[2] = pi5;  
	      el33[3] = pi3;  
	      el33.SetIndex (mattyp);
	      
	      elementsonnode.Add (pi4, hasbothpoints[1]);
	      elementsonnode.Add (pi3, hasbothpoints[2]);
	      
	      bad1 = CalcBad (mesh.Points(), el31, 0) + 
		CalcBad (mesh.Points(), el32, 0) +
		CalcBad (mesh.Points(), el33, 0);

	      el31.flags.illegal_valid = 0;
	      el32.flags.illegal_valid = 0;
	      el33.flags.illegal_valid = 0;

	      if (!mesh.LegalTet(el31) ||
		  !mesh.LegalTet(el32) ||
		  !mesh.LegalTet(el33))
		bad1 += 1e4;
	      
	      el21[0] = pi3;
	      el21[1] = pi4;
	      el21[2] = pi5;
	      el21[3] = pi2;
	      el21.SetIndex (mattyp);
	      
	      el22[0] = pi5;
	      el22[1] = pi4;
	      el22[2] = pi3;
	      el22[3] = pi1;
	      el22.SetIndex (mattyp);	      

	      bad2 = CalcBad (mesh.Points(), el21, 0) + 
		CalcBad (mesh.Points(), el22, 0);

	      el21.flags.illegal_valid = 0;
	      el22.flags.illegal_valid = 0;

	      if (!mesh.LegalTet(el21) ||
		  !mesh.LegalTet(el22))
		bad2 += 1e4;


	      if (goal == OPT_CONFORM && bad2 < 1e4)
		{
		  INDEX_3 face(pi3, pi4, pi5);
		  face.Sort();
		  if (faces.Used(face))
		    {
		      // (*testout) << "3->2 swap, could improve conformity, bad1 = " << bad1
		      //				 << ", bad2 = " << bad2 << endl;
		      if (bad2 < 1e4)
			bad1 = 2 * bad2;
		    }
		  /*
		    else
		    {
		    INDEX_2 hi1(pi3, pi4);
		    hi1.Sort();
		    INDEX_2 hi2(pi3, pi5);
		    hi2.Sort();
		    INDEX_2 hi3(pi4, pi5);
		    hi3.Sort();

		    if (boundaryedges->Used (hi1) ||
		    boundaryedges->Used (hi2) ||
		    boundaryedges->Used (hi3) )
		    bad1 = 2 * bad2;
		    }
		  */
		}

	      if (bad2 < bad1)
		{
		  //		  (*mycout) << "3->2 " << flush;
		  //		  (*testout) << "3->2 conversion" << endl;
		  cnt++;
		  

		  /*
		  (*testout) << "3->2 swap, old els = " << endl
			     << mesh[hasbothpoints[0]] << endl
			     << mesh[hasbothpoints[1]] << endl
			     << mesh[hasbothpoints[2]] << endl
			     << "new els = " << endl
			     << el21 << endl
			     << el22 << endl;
		  */
                  
		  el21.flags.illegal_valid = 0;
		  el22.flags.illegal_valid = 0;
		  mesh[hasbothpoints[0]] = el21;
		  mesh[hasbothpoints[1]] = el22;
		  for (int l = 0; l < 4; l++)
		    mesh[hasbothpoints[2]][l].Invalidate();
		  mesh[hasbothpoints[2]].Delete();

		  for (int k = 0; k < 2; k++)
		    for (int l = 0; l < 4; l++)
		      elementsonnode.Add (mesh[hasbothpoints[k]][l], hasbothpoints[k]);
		}
	    }	  
          
	  if (nsuround == 4)
	    {
	      const Element & elem1 = mesh[hasbothpoints[0]];
	      for (int l = 0; l < 4; l++)
		if (elem1[l] != pi1 && elem1[l] != pi2)
		  {
		    pi4 = pi3;
		    pi3 = elem1[l];
		  }
	      
	      el1[0] = pi1; el1[1] = pi2;
	      el1[2] = pi3; el1[3] = pi4;
	      el1.SetIndex (mattyp);

	      if (WrongOrientation (mesh.Points(), el1))
		{
		  Swap (pi3, pi4);
		  el1[2] = pi3;
		  el1[3] = pi4;
		}
	      
	      pi5.Invalidate();
	      for (int k = 0; k < 4; k++)
		{
		  const Element & elem = mesh[hasbothpoints[k]];
                  bool has1 = elem.PNums().Contains(pi4);
		  if (has1)
		    {
		      for (int l = 0; l < 4; l++)
			if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi4)
			  pi5 = elem[l];
		    }
		}
	      
	      pi6.Invalidate();
	      for (int k = 0; k < 4; k++)
		{
		  const Element & elem = mesh[hasbothpoints[k]];
                  bool has1 = elem.PNums().Contains(pi3);
		  if (has1)
		    {
		      for (int l = 0; l < 4; l++)
			if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi3)
			  pi6 = elem[l];
		    }
		}


	      /*
	      INDEX_2 i22(pi3, pi5);
	      i22.Sort();
	      INDEX_2 i23(pi4, pi6);
	      i23.Sort();
	      */
	      
	      el1[0] = pi1; el1[1] = pi2;  
	      el1[2] = pi3; el1[3] = pi4; 
	      el1.SetIndex (mattyp);
	      
	      el2[0] = pi1; el2[1] = pi2;  
	      el2[2] = pi4; el2[3] = pi5;  
	      el2.SetIndex (mattyp);
	      
	      el3[0] = pi1; el3[1] = pi2;  
	      el3[2] = pi5; el3[3] = pi6;  
	      el3.SetIndex (mattyp);
	      
	      el4[0] = pi1; el4[1] = pi2;  
	      el4[2] = pi6; el4[3] = pi3;  
	      el4.SetIndex (mattyp);
	      
	      //        elementsonnode.Add (pi4, hasbothpoints.Elem(2));
	      //        elementsonnode.Add (pi3, hasbothpoints.Elem(3));
	      
	      bad1 = CalcBad (mesh.Points(), el1, 0) + 
		CalcBad (mesh.Points(), el2, 0) +
		CalcBad (mesh.Points(), el3, 0) + 
		CalcBad (mesh.Points(), el4, 0);
	      

	      el1.flags.illegal_valid = 0;
	      el2.flags.illegal_valid = 0;
	      el3.flags.illegal_valid = 0;
	      el4.flags.illegal_valid = 0;


	      if (goal != OPT_CONFORM)
		{
		  if (!mesh.LegalTet(el1) ||
		      !mesh.LegalTet(el2) ||
		      !mesh.LegalTet(el3) ||
		      !mesh.LegalTet(el4))
		    bad1 += 1e4;
		}
	      
	      el1[0] = pi3; el1[1] = pi5;  
	      el1[2] = pi2; el1[3] = pi4; 
	      el1.SetIndex (mattyp);
	 
	      el2[0] = pi3; el2[1] = pi5;  
	      el2[2] = pi4; el2[3] = pi1;  
	      el2.SetIndex (mattyp);
	      
	      el3[0] = pi3; el3[1] = pi5;  
	      el3[2] = pi1; el3[3] = pi6;  
	      el3.SetIndex (mattyp);
	      
	      el4[0] = pi3; el4[1] = pi5;  
	      el4[2] = pi6; el4[3] = pi2;  	
	      el4.SetIndex (mattyp);
      
	      bad2 = CalcBad (mesh.Points(), el1, 0) + 
		CalcBad (mesh.Points(), el2, 0) +
		CalcBad (mesh.Points(), el3, 0) + 
		CalcBad (mesh.Points(), el4, 0);

	      el1.flags.illegal_valid = 0;
	      el2.flags.illegal_valid = 0;
	      el3.flags.illegal_valid = 0;
	      el4.flags.illegal_valid = 0;

	      if (goal != OPT_CONFORM)
		{
		  if (!mesh.LegalTet(el1) ||
		      !mesh.LegalTet(el2) ||
		      !mesh.LegalTet(el3) ||
		      !mesh.LegalTet(el4))
		    bad2 += 1e4;
		}

	      
	      el1b[0] = pi4; el1b[1] = pi6;  
	      el1b[2] = pi3; el1b[3] = pi2; 
	      el1b.SetIndex (mattyp);
	      
	      el2b[0] = pi4; el2b[1] = pi6;  
	      el2b[2] = pi2; el2b[3] = pi5;  
	      el2b.SetIndex (mattyp);
	      
	      el3b[0] = pi4; el3b[1] = pi6;  
	      el3b[2] = pi5; el3b[3] = pi1;  
	      el3b.SetIndex (mattyp);
	      
	      el4b[0] = pi4; el4b[1] = pi6;  
	      el4b[2] = pi1; el4b[3] = pi3;  
	      el4b.SetIndex (mattyp);
	      
	      bad3 = CalcBad (mesh.Points(), el1b, 0) + 
		CalcBad (mesh.Points(), el2b, 0) +
		CalcBad (mesh.Points(), el3b, 0) +
		CalcBad (mesh.Points(), el4b, 0);
	      
	      el1b.flags.illegal_valid = 0;
	      el2b.flags.illegal_valid = 0;
	      el3b.flags.illegal_valid = 0;
	      el4b.flags.illegal_valid = 0;

	      if (goal != OPT_CONFORM)
		{
		  if (!mesh.LegalTet(el1b) ||
		      !mesh.LegalTet(el2b) ||
		      !mesh.LegalTet(el3b) ||
		      !mesh.LegalTet(el4b))
		    bad3 += 1e4;
		}


	      /*
	      int swap2 = (bad2 < bad1) && (bad2 < bad3);
	      int swap3 = !swap2 && (bad3 < bad1);
	      
	      if ( ((bad2 < 10 * bad1) || 
		    (bad2 < 1e6)) && mesh.BoundaryEdge (pi3, pi5))
		swap2 = 1;
	      else  if ( ((bad3 < 10 * bad1) || 
			  (bad3 < 1e6)) && mesh.BoundaryEdge (pi4, pi6))
		{
		  swap3 = 1;
		  swap2 = 0;
		}
	      */
	      bool swap2, swap3;

	      if (goal != OPT_CONFORM)
		{
		  swap2 = (bad2 < bad1) && (bad2 < bad3);
		  swap3 = !swap2 && (bad3 < bad1);
		}
	      else
		{
		  if (mesh.BoundaryEdge (pi3, pi5)) bad2 /= 1e6;
		  if (mesh.BoundaryEdge (pi4, pi6)) bad3 /= 1e6;

		  swap2 = (bad2 < bad1) && (bad2 < bad3);
		  swap3 = !swap2 && (bad3 < bad1);
		}
		

	      if (swap2 || swap3)
		{
		  // (*mycout) << "4->4 " << flush;
		  cnt++;
		  //		  (*testout) << "4->4 conversion" << "\n";
		  /*
		    (*testout) << "bad1 = " << bad1 
		    << " bad2 = " << bad2
		    << " bad3 = " << bad3 << "\n";
		  
		    (*testout) << "Points: " << pi1 << " " << pi2 << " " << pi3 
		    << " " << pi4 << " " << pi5 << " " << pi6 << "\n";
		    (*testout) << "Elements: " 
		    << hasbothpoints.Get(1) << "  "
		    << hasbothpoints.Get(2) << "  "
		    << hasbothpoints.Get(3) << "  "
		    << hasbothpoints.Get(4) << "  " << "\n";
		  */

		  /*
		    {
		    int i1, j1;
		    for (i1 = 1; i1 <= 4; i1++)
		    {
		    for (j1 = 1; j1 <= 4; j1++)
		    (*testout) << volelements.Get(hasbothpoints.Get(i1)).PNum(j1)
		    << "  ";
		    (*testout) << "\n";
		    }
		    }
		  */
		}
	      

	      if (swap2)
		{
		  //		  (*mycout) << "bad1 = " << bad1 << " bad2 = " << bad2 << "\n";


		  /*
		  (*testout) << "4->4 swap A, old els = " << endl
			     << mesh[hasbothpoints[0]] << endl
			     << mesh[hasbothpoints[1]] << endl
			     << mesh[hasbothpoints[2]] << endl
			     << mesh[hasbothpoints[3]] << endl
			     << "new els = " << endl
			     << el1 << endl
			     << el2 << endl
			     << el3 << endl
			     << el4 << endl;
		  */



		  el1.flags.illegal_valid = 0;
		  el2.flags.illegal_valid = 0;
		  el3.flags.illegal_valid = 0;
		  el4.flags.illegal_valid = 0;
		  
		  mesh[hasbothpoints[0]] = el1;
		  mesh[hasbothpoints[1]] = el2;
		  mesh[hasbothpoints[2]] = el3;
		  mesh[hasbothpoints[3]] = el4;

		  for (int k = 0; k < 4; k++)
		    for (int l = 0; l < 4; l++)
		      elementsonnode.Add (mesh[hasbothpoints[k]][l], hasbothpoints[k]);
		}
	      else if (swap3)
		{
		  // (*mycout) << "bad1 = " << bad1 << " bad3 = " << bad3 << "\n";
		  el1b.flags.illegal_valid = 0;
		  el2b.flags.illegal_valid = 0;
		  el3b.flags.illegal_valid = 0;
		  el4b.flags.illegal_valid = 0;


		  /*
		  (*testout) << "4->4 swap A, old els = " << endl
			     << mesh[hasbothpoints[0]] << endl
			     << mesh[hasbothpoints[1]] << endl
			     << mesh[hasbothpoints[2]] << endl
			     << mesh[hasbothpoints[3]] << endl
			     << "new els = " << endl
			     << el1b << endl
			     << el2b << endl
			     << el3b << endl
			     << el4b << endl;
		  */


		  mesh[hasbothpoints[0]] = el1b;
		  mesh[hasbothpoints[1]] = el2b;
		  mesh[hasbothpoints[2]] = el3b;
		  mesh[hasbothpoints[3]] = el4b;

		  for (int k = 0; k < 4; k++)
		    for (int l = 0; l < 4; l++)
		      elementsonnode.Add (mesh[hasbothpoints[k]][l], hasbothpoints[k]);
		}
	    }

          // if (goal == OPT_QUALITY)
	  if (nsuround >= 5) 
	    {
	      Element hel(TET);

	      NgArrayMem<PointIndex, 50> suroundpts(nsuround);
	      NgArrayMem<bool, 50> tetused(nsuround);

	      Element & elem = mesh[hasbothpoints[0]];

	      for (int l = 0; l < 4; l++)
		if (elem[l] != pi1 && elem[l] != pi2)
		  {
		    pi4 = pi3;
		    pi3 = elem[l];
		  }

	      hel[0] = pi1;
	      hel[1] = pi2;
	      hel[2] = pi3;
	      hel[3] = pi4;
	      hel.SetIndex (mattyp);
	      
	      if (WrongOrientation (mesh.Points(), hel))
		{
		  Swap (pi3, pi4);
		  hel[2] = pi3;
		  hel[3] = pi4;
		}

	      
	      // suroundpts.SetSize (nsuround);
              suroundpts = PointIndex::INVALID;
	      suroundpts[0] = pi3;
	      suroundpts[1] = pi4;

	      tetused = false;
	      tetused[0] = true;

	      for (int l = 2; l < nsuround; l++)
		{
		  PointIndex oldpi = suroundpts[l-1];
		  PointIndex newpi;
                  newpi.Invalidate();

		  for (int k = 0; k < nsuround && !newpi.IsValid(); k++)
		    if (!tetused[k])
		      {
			const Element & nel = mesh[hasbothpoints[k]];
			for (int k2 = 0; k2 < 4 && !newpi.IsValid(); k2++)
			  if (nel[k2] == oldpi)
			    {
			      newpi = 
				nel[0] + nel[1] + nel[2] + nel[3] 
				- pi1 - pi2 - oldpi;

			      tetused[k] = true; 
			      suroundpts[l] = newpi;
			    }
			}
		}
	      
	      
	      bad1 = 0;
	      for (int k = 0; k < nsuround; k++)
		{
		  hel[0] = pi1;
		  hel[1] = pi2;
		  hel[2] = suroundpts[k];
		  hel[3] = suroundpts[(k+1) % nsuround];
		  hel.SetIndex (mattyp);

		  bad1 += CalcBad (mesh.Points(), hel, 0);
		}

	      //  (*testout) << "nsuround = " << nsuround << " bad1 = " << bad1 << endl;


	      int bestl = -1;
	      int confface = -1;
	      int confedge = -1;
	      double badopt = bad1;

	      for (int l = 0; l < nsuround; l++)
		{
		  bad2 = 0;

		  for (int k = l+1; k <= nsuround + l - 2; k++)
		    {
		      hel[0] = suroundpts[l];
		      hel[1] = suroundpts[k % nsuround];
		      hel[2] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi2;

		      bad2 += CalcBad (mesh.Points(), hel, 0);
		      hel.flags.illegal_valid = 0;
		      if (!mesh.LegalTet(hel)) bad2 += 1e4;

		      hel[2] = suroundpts[k % nsuround];
		      hel[1] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi1;

		      bad2 += CalcBad (mesh.Points(), hel, 0);

		      hel.flags.illegal_valid = 0;
		      if (!mesh.LegalTet(hel)) bad2 += 1e4;
		    }
		  // (*testout) << "bad2," << l << " = " << bad2 << endl;
		  
		  if ( bad2 < badopt )
		    {
		      bestl = l;
		      badopt = bad2;
		    }
		  
		  
		  if (goal == OPT_CONFORM)
		       // (bad2 <= 100 * bad1 || bad2 <= 1e6))
		    {
		      bool nottoobad =
			(bad2 <= bad1) ||
			(bad2 <= 100 * bad1 && bad2 <= 1e18) ||
			(bad2 <= 1e8);
		      
		      for (int k = l+1; k <= nsuround + l - 2; k++)
			{
			  INDEX_3 hi3(suroundpts[l],
				      suroundpts[k % nsuround],
				      suroundpts[(k+1) % nsuround]);
			  hi3.Sort();
			  if (faces.Used(hi3))
			    {
			      // (*testout) << "could improve face conformity, bad1 = " << bad1
			      // << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
			      if (nottoobad)
				confface = l;
			    }
			}

		      for (int k = l+2; k <= nsuround+l-2; k++)
			{
			  if (mesh.BoundaryEdge (suroundpts[l],
						 suroundpts[k % nsuround]))
			    {
			      /*
			      *testout << "could improve edge conformity, bad1 = " << bad1
				   << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
			      */
			      if (nottoobad)
				confedge = l;
			    }
			}
		    }
		}
	      
	      if (confedge != -1)
		bestl = confedge;
	      if (confface != -1)
		bestl = confface;

	      if (bestl != -1)
		{
		  // (*mycout) << nsuround << "->" << 2 * (nsuround-2) << " " << flush;
		  cnt++;
		  
		  for (int k = bestl+1; k <= nsuround + bestl - 2; k++)
		    {
		      int k1;

		      hel[0] = suroundpts[bestl];
		      hel[1] = suroundpts[k % nsuround];
		      hel[2] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi2;
		      hel.flags.illegal_valid = 0;

		      /*
		      (*testout) << nsuround << "-swap, new el,top = "
				 << hel << endl;
		      */
		      mesh.AddVolumeElement (hel);

		      for (k1 = 0; k1 < 4; k1++)
			elementsonnode.Add (hel[k1], mesh.GetNE()-1);
		      

		      hel[2] = suroundpts[k % nsuround];
		      hel[1] = suroundpts[(k+1) % nsuround];
		      hel[3] = pi1;

		      /*
		      (*testout) << nsuround << "-swap, new el,bot = "
				 << hel << endl;
		      */

		      mesh.AddVolumeElement (hel);

		      for (k1 = 0; k1 < 4; k1++)
			elementsonnode.Add (hel[k1], mesh.GetNE()-1);
		    }		  
		  
		  for (int k = 0; k < nsuround; k++)
		    {
		      Element & rel = mesh[hasbothpoints[k]];
		      /*
		      (*testout) << nsuround << "-swap, old el = "
				 << rel << endl;
		      */
		      rel.Delete();
		      for (int k1 = 0; k1 < 4; k1++)
			rel[k1].Invalidate();
		    }
		}
	    }
	}

      /*
	if (onlybedges)
	{
	(*testout) << "bad tet: " 
	<< volelements.Get(i)[0] 
	<< volelements.Get(i)[1] 
	<< volelements.Get(i)[2] 
	<< volelements.Get(i)[3] << "\n";

	if (!mesh.LegalTet (volelements.Get(i)))
	cerr << "Illegal tet" << "\n";
	}
      */
    }
  //  (*mycout) << endl;
  tloop.Stop();
  /*  
      cout << "edgeused: ";
      edgeused.PrintMemInfo(cout);
  */
  PrintMessage (5, cnt, " swaps performed");





  mesh.Compress ();

  /*
  if (goal == OPT_QUALITY)
    {
      bad1 = CalcTotalBad (mesh.Points(), mesh.VolumeElements());
      //      (*testout) << "Total badness = " << bad1 << endl;
    }
  */

  /*
    for (i = 1; i <= GetNE(); i++)
    if (volelements.Get(i)[0])
    if (!mesh.LegalTet (volelements.Get(i)))
    {
    cout << "detected illegal tet, 2" << endl;
    (*testout) << "detected illegal tet1: " << i << endl;
    }
  */

  multithread.task = savetask;
}


double MeshOptimize3d :: SwapImproveEdge (Mesh & mesh, OPTIMIZEGOAL goal,
        const NgBitArray * working_elements,
        Table<ElementIndex, PointIndex> & elementsonnode,
        INDEX_3_HASHTABLE<int> & faces,
        PointIndex pi1, PointIndex pi2, bool check_only)
{
  PointIndex pi3(PointIndex::INVALID), pi4(PointIndex::INVALID),
             pi5(PointIndex::INVALID), pi6(PointIndex::INVALID);

  double bad1, bad2, bad3;

  Element el21(TET), el22(TET), el31(TET), el32(TET), el33(TET);
  Element el1(TET), el2(TET), el3(TET), el4(TET);
  Element el1b(TET), el2b(TET), el3b(TET), el4b(TET);
  ArrayMem<ElementIndex, 20> hasbothpoints;

  double d_badness = 0.0;
  if (pi2 < pi1) Swap (pi1, pi2);

  if (mesh.BoundaryEdge (pi1, pi2)) return 0.0;


  hasbothpoints.SetSize (0);
  for (ElementIndex elnr : elementsonnode[pi1])
    {
      bool has1 = 0, has2 = 0;
      const Element & elem = mesh[elnr];

      if (elem.IsDeleted()) return 0.0;

      for (int l = 0; l < elem.GetNP(); l++)
        {
          if (elem[l] == pi1) has1 = 1;
          if (elem[l] == pi2) has2 = 1;
        }

      if (has1 && has2)
        { // only once
          if (hasbothpoints.Contains (elnr))
              has1 = false;

          if (has1)
            {
              hasbothpoints.Append (elnr);
            }
        }
    }

  for (ElementIndex ei : hasbothpoints)
    {
      if (mesh[ei].GetType () != TET)
          return 0.0;

      if (mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(ei).GetIndex())
          return 0.0;


      if ((mesh.ElementType(ei)) == FIXEDELEMENT)
          return 0.0;

      if(working_elements &&
              ei < working_elements->Size() &&
              !working_elements->Test(ei))
          return 0.0;

      if (mesh[ei].IsDeleted())
          return 0.0;

      if ((goal == OPT_LEGAL) &&
              mesh.LegalTet (mesh[ei]) &&
              CalcBad (mesh.Points(), mesh[ei], 0) < 1e3)
          return 0.0;
    }

  int nsuround = hasbothpoints.Size();
  int mattyp = mesh[hasbothpoints[0]].GetIndex();

  if ( nsuround == 3 )
    {
      Element & elem = mesh[hasbothpoints[0]];
      for (int l = 0; l < 4; l++)
          if (elem[l] != pi1 && elem[l] != pi2)
            {
              pi4 = pi3;
              pi3 = elem[l];
            }

      el31[0] = pi1;
      el31[1] = pi2;
      el31[2] = pi3;
      el31[3] = pi4;
      el31.SetIndex (mattyp);

      if (WrongOrientation (mesh.Points(), el31))
        {
          Swap (pi3, pi4);
          el31[2] = pi3;
          el31[3] = pi4;
        }

      pi5.Invalidate();
      for (int k = 0; k < 3; k++)   // JS, 201212
        {
          const Element & elemk = mesh[hasbothpoints[k]];
          bool has1 = false;
          for (int l = 0; l < 4; l++)
              if (elemk[l] == pi4)
                  has1 = true;
          if (has1)
            {
              for (int l = 0; l < 4; l++)
                  if (elemk[l] != pi1 && elemk[l] != pi2 && elemk[l] != pi4)
                      pi5 = elemk[l];
            }
        }

      if (!pi5.IsValid())
          throw NgException("Illegal state observed in SwapImprove");


      el32[0] = pi1;
      el32[1] = pi2;
      el32[2] = pi4;
      el32[3] = pi5;
      el32.SetIndex (mattyp);

      el33[0] = pi1;
      el33[1] = pi2;
      el33[2] = pi5;
      el33[3] = pi3;
      el33.SetIndex (mattyp);

      bad1 = CalcBad (mesh.Points(), el31, 0) +
          CalcBad (mesh.Points(), el32, 0) +
          CalcBad (mesh.Points(), el33, 0);

      el31.flags.illegal_valid = 0;
      el32.flags.illegal_valid = 0;
      el33.flags.illegal_valid = 0;

      if (!mesh.LegalTet(el31) ||
              !mesh.LegalTet(el32) ||
              !mesh.LegalTet(el33))
          bad1 += 1e4;

      el21[0] = pi3;
      el21[1] = pi4;
      el21[2] = pi5;
      el21[3] = pi2;
      el21.SetIndex (mattyp);

      el22[0] = pi5;
      el22[1] = pi4;
      el22[2] = pi3;
      el22[3] = pi1;
      el22.SetIndex (mattyp);

      bad2 = CalcBad (mesh.Points(), el21, 0) +
          CalcBad (mesh.Points(), el22, 0);

      el21.flags.illegal_valid = 0;
      el22.flags.illegal_valid = 0;

      if (!mesh.LegalTet(el21) ||
              !mesh.LegalTet(el22))
          bad2 += 1e4;


      if (goal == OPT_CONFORM && bad2 < 1e4)
        {
          INDEX_3 face(pi3, pi4, pi5);
          face.Sort();
          if (faces.Used(face))
            {
              // (*testout) << "3->2 swap, could improve conformity, bad1 = " << bad1
              //				 << ", bad2 = " << bad2 << endl;
              if (bad2 < 1e4)
                  bad1 = 2 * bad2;
            }
        }

      if (bad2 < bad1)
        {
          //		  (*mycout) << "3->2 " << flush;
          //		  (*testout) << "3->2 conversion" << endl;
          d_badness = bad2-bad1;
          if(check_only)
              return d_badness;


          /*
             (*testout) << "3->2 swap, old els = " << endl
             << mesh[hasbothpoints[0]] << endl
             << mesh[hasbothpoints[1]] << endl
             << mesh[hasbothpoints[2]] << endl
             << "new els = " << endl
             << el21 << endl
             << el22 << endl;
             */

          mesh[hasbothpoints[0]].Delete();
          mesh[hasbothpoints[1]].Delete();
          mesh[hasbothpoints[2]].Delete();

          el21.flags.illegal_valid = 0;
          el22.flags.illegal_valid = 0;
          mesh.AddVolumeElement(el21);
          mesh.AddVolumeElement(el22);
        }
    }

  if (nsuround == 4)
    {
      const Element & elem1 = mesh[hasbothpoints[0]];
      for (int l = 0; l < 4; l++)
          if (elem1[l] != pi1 && elem1[l] != pi2)
            {
              pi4 = pi3;
              pi3 = elem1[l];
            }

      el1[0] = pi1; el1[1] = pi2;
      el1[2] = pi3; el1[3] = pi4;
      el1.SetIndex (mattyp);

      if (WrongOrientation (mesh.Points(), el1))
        {
          Swap (pi3, pi4);
          el1[2] = pi3;
          el1[3] = pi4;
        }

      pi5.Invalidate();
      for (int k = 0; k < 4; k++)
        {
          const Element & elem = mesh[hasbothpoints[k]];
          bool has1 = elem.PNums().Contains(pi4);
          if (has1)
            {
              for (int l = 0; l < 4; l++)
                  if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi4)
                      pi5 = elem[l];
            }
        }

      pi6.Invalidate();
      for (int k = 0; k < 4; k++)
        {
          const Element & elem = mesh[hasbothpoints[k]];
          bool has1 = elem.PNums().Contains(pi3);
          if (has1)
            {
              for (int l = 0; l < 4; l++)
                  if (elem[l] != pi1 && elem[l] != pi2 && elem[l] != pi3)
                      pi6 = elem[l];
            }
        }

      el1[0] = pi1; el1[1] = pi2;
      el1[2] = pi3; el1[3] = pi4;
      el1.SetIndex (mattyp);

      el2[0] = pi1; el2[1] = pi2;
      el2[2] = pi4; el2[3] = pi5;
      el2.SetIndex (mattyp);

      el3[0] = pi1; el3[1] = pi2;
      el3[2] = pi5; el3[3] = pi6;
      el3.SetIndex (mattyp);

      el4[0] = pi1; el4[1] = pi2;
      el4[2] = pi6; el4[3] = pi3;
      el4.SetIndex (mattyp);

      bad1 = CalcBad (mesh.Points(), el1, 0) +
          CalcBad (mesh.Points(), el2, 0) +
          CalcBad (mesh.Points(), el3, 0) +
          CalcBad (mesh.Points(), el4, 0);


      el1.flags.illegal_valid = 0;
      el2.flags.illegal_valid = 0;
      el3.flags.illegal_valid = 0;
      el4.flags.illegal_valid = 0;


      if (goal != OPT_CONFORM)
        {
          if (!mesh.LegalTet(el1) ||
                  !mesh.LegalTet(el2) ||
                  !mesh.LegalTet(el3) ||
                  !mesh.LegalTet(el4))
              bad1 += 1e4;
        }

      el1[0] = pi3; el1[1] = pi5;
      el1[2] = pi2; el1[3] = pi4;
      el1.SetIndex (mattyp);

      el2[0] = pi3; el2[1] = pi5;
      el2[2] = pi4; el2[3] = pi1;
      el2.SetIndex (mattyp);

      el3[0] = pi3; el3[1] = pi5;
      el3[2] = pi1; el3[3] = pi6;
      el3.SetIndex (mattyp);

      el4[0] = pi3; el4[1] = pi5;
      el4[2] = pi6; el4[3] = pi2;  	
      el4.SetIndex (mattyp);

      bad2 = CalcBad (mesh.Points(), el1, 0) +
          CalcBad (mesh.Points(), el2, 0) +
          CalcBad (mesh.Points(), el3, 0) +
          CalcBad (mesh.Points(), el4, 0);

      el1.flags.illegal_valid = 0;
      el2.flags.illegal_valid = 0;
      el3.flags.illegal_valid = 0;
      el4.flags.illegal_valid = 0;

      if (goal != OPT_CONFORM)
        {
          if (!mesh.LegalTet(el1) ||
                  !mesh.LegalTet(el2) ||
                  !mesh.LegalTet(el3) ||
                  !mesh.LegalTet(el4))
              bad2 += 1e4;
        }


      el1b[0] = pi4; el1b[1] = pi6;
      el1b[2] = pi3; el1b[3] = pi2;
      el1b.SetIndex (mattyp);

      el2b[0] = pi4; el2b[1] = pi6;
      el2b[2] = pi2; el2b[3] = pi5;
      el2b.SetIndex (mattyp);

      el3b[0] = pi4; el3b[1] = pi6;
      el3b[2] = pi5; el3b[3] = pi1;
      el3b.SetIndex (mattyp);

      el4b[0] = pi4; el4b[1] = pi6;
      el4b[2] = pi1; el4b[3] = pi3;
      el4b.SetIndex (mattyp);

      bad3 = CalcBad (mesh.Points(), el1b, 0) +
          CalcBad (mesh.Points(), el2b, 0) +
          CalcBad (mesh.Points(), el3b, 0) +
          CalcBad (mesh.Points(), el4b, 0);

      el1b.flags.illegal_valid = 0;
      el2b.flags.illegal_valid = 0;
      el3b.flags.illegal_valid = 0;
      el4b.flags.illegal_valid = 0;

      if (goal != OPT_CONFORM)
        {
          if (!mesh.LegalTet(el1b) ||
                  !mesh.LegalTet(el2b) ||
                  !mesh.LegalTet(el3b) ||
                  !mesh.LegalTet(el4b))
              bad3 += 1e4;
        }

      bool swap2, swap3;

      if (goal != OPT_CONFORM)
        {
          swap2 = (bad2 < bad1) && (bad2 < bad3);
          swap3 = !swap2 && (bad3 < bad1);
        }
      else
        {
          if (mesh.BoundaryEdge (pi3, pi5)) bad2 /= 1e6;
          if (mesh.BoundaryEdge (pi4, pi6)) bad3 /= 1e6;

          swap2 = (bad2 < bad1) && (bad2 < bad3);
          swap3 = !swap2 && (bad3 < bad1);
        }

      d_badness = swap2 ? bad2-bad1 : bad3-bad1;
      if(check_only)
          return d_badness;

      if (swap2)
        {
          for (auto i : IntRange(4))
              mesh[hasbothpoints[i]].Delete();

          el1.flags.illegal_valid = 0;
          el2.flags.illegal_valid = 0;
          el3.flags.illegal_valid = 0;
          el4.flags.illegal_valid = 0;
          mesh.AddVolumeElement (el1);
          mesh.AddVolumeElement (el2);
          mesh.AddVolumeElement (el3);
          mesh.AddVolumeElement (el4);
        }
      else if (swap3)
        {
          for (auto i : IntRange(4))
              mesh[hasbothpoints[i]].Delete();

          el1b.flags.illegal_valid = 0;
          el2b.flags.illegal_valid = 0;
          el3b.flags.illegal_valid = 0;
          el4b.flags.illegal_valid = 0;
          mesh.AddVolumeElement (el1b);
          mesh.AddVolumeElement (el2b);
          mesh.AddVolumeElement (el3b);
          mesh.AddVolumeElement (el4b);
        }
    }

  // if (goal == OPT_QUALITY)
  if (nsuround >= 5)
    {
      Element hel(TET);

      NgArrayMem<PointIndex, 50> suroundpts(nsuround);
      NgArrayMem<bool, 50> tetused(nsuround);

      Element & elem = mesh[hasbothpoints[0]];

      for (int l = 0; l < 4; l++)
          if (elem[l] != pi1 && elem[l] != pi2)
            {
              pi4 = pi3;
              pi3 = elem[l];
            }

      hel[0] = pi1;
      hel[1] = pi2;
      hel[2] = pi3;
      hel[3] = pi4;
      hel.SetIndex (mattyp);

      if (WrongOrientation (mesh.Points(), hel))
        {
          Swap (pi3, pi4);
          hel[2] = pi3;
          hel[3] = pi4;
        }


      // suroundpts.SetSize (nsuround);
      suroundpts = PointIndex::INVALID;
      suroundpts[0] = pi3;
      suroundpts[1] = pi4;

      tetused = false;
      tetused[0] = true;

      for (int l = 2; l < nsuround; l++)
        {
          PointIndex oldpi = suroundpts[l-1];
          PointIndex newpi;
          newpi.Invalidate();

          for (int k = 0; k < nsuround && !newpi.IsValid(); k++)
              if (!tetused[k])
                {
                  const Element & nel = mesh[hasbothpoints[k]];
                  for (int k2 = 0; k2 < 4 && !newpi.IsValid(); k2++)
                      if (nel[k2] == oldpi)
                        {
                          newpi =
                              nel[0] + nel[1] + nel[2] + nel[3]
                              - pi1 - pi2 - oldpi;

                          tetused[k] = true;
                          suroundpts[l] = newpi;
                        }
                }
        }


      bad1 = 0;
      for (int k = 0; k < nsuround; k++)
        {
          hel[0] = pi1;
          hel[1] = pi2;
          hel[2] = suroundpts[k];
          hel[3] = suroundpts[(k+1) % nsuround];
          hel.SetIndex (mattyp);

          bad1 += CalcBad (mesh.Points(), hel, 0);
        }

      //  (*testout) << "nsuround = " << nsuround << " bad1 = " << bad1 << endl;


      int bestl = -1;
      int confface = -1;
      int confedge = -1;
      double badopt = bad1;

      for (int l = 0; l < nsuround; l++)
        {
          bad2 = 0;

          for (int k = l+1; k <= nsuround + l - 2; k++)
            {
              hel[0] = suroundpts[l];
              hel[1] = suroundpts[k % nsuround];
              hel[2] = suroundpts[(k+1) % nsuround];
              hel[3] = pi2;

              bad2 += CalcBad (mesh.Points(), hel, 0);
              hel.flags.illegal_valid = 0;
              if (!mesh.LegalTet(hel)) bad2 += 1e4;

              hel[2] = suroundpts[k % nsuround];
              hel[1] = suroundpts[(k+1) % nsuround];
              hel[3] = pi1;

              bad2 += CalcBad (mesh.Points(), hel, 0);

              hel.flags.illegal_valid = 0;
              if (!mesh.LegalTet(hel)) bad2 += 1e4;
            }
          // (*testout) << "bad2," << l << " = " << bad2 << endl;

          if ( bad2 < badopt )
            {
              bestl = l;
              badopt = bad2;
            }


          if (goal == OPT_CONFORM)
              // (bad2 <= 100 * bad1 || bad2 <= 1e6))
            {
              bool nottoobad =
                  (bad2 <= bad1) ||
                  (bad2 <= 100 * bad1 && bad2 <= 1e18) ||
                  (bad2 <= 1e8);

              for (int k = l+1; k <= nsuround + l - 2; k++)
                {
                  INDEX_3 hi3(suroundpts[l],
                          suroundpts[k % nsuround],
                          suroundpts[(k+1) % nsuround]);
                  hi3.Sort();
                  if (faces.Used(hi3))
                    {
                      // (*testout) << "could improve face conformity, bad1 = " << bad1
                      // << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
                      if (nottoobad)
                          confface = l;
                    }
                }

              for (int k = l+2; k <= nsuround+l-2; k++)
                {
                  if (mesh.BoundaryEdge (suroundpts[l],
                              suroundpts[k % nsuround]))
                    {
                      /*
                       *testout << "could improve edge conformity, bad1 = " << bad1
                       << ", bad 2 = " << bad2 << ", nottoobad = " << nottoobad << endl;
                       */
                      if (nottoobad)
                          confedge = l;
                    }
                }
            }
        }

      if (confedge != -1)
          bestl = confedge;
      if (confface != -1)
          bestl = confface;

      if (bestl != -1)
        {
          // (*mycout) << nsuround << "->" << 2 * (nsuround-2) << " " << flush;
          d_badness = badopt-bad1;
          if(check_only)
              return d_badness;

          for (int k = bestl+1; k <= nsuround + bestl - 2; k++)
            {
              int k1;

              hel[0] = suroundpts[bestl];
              hel[1] = suroundpts[k % nsuround];
              hel[2] = suroundpts[(k+1) % nsuround];
              hel[3] = pi2;
              hel.flags.illegal_valid = 0;

              /*
                 (*testout) << nsuround << "-swap, new el,top = "
                 << hel << endl;
                 */
              mesh.AddVolumeElement (hel);

              hel[2] = suroundpts[k % nsuround];
              hel[1] = suroundpts[(k+1) % nsuround];
              hel[3] = pi1;

              /*
                 (*testout) << nsuround << "-swap, new el,bot = "
                 << hel << endl;
                 */

              mesh.AddVolumeElement (hel);
            }

          for (int k = 0; k < nsuround; k++)
            {
              Element & rel = mesh[hasbothpoints[k]];
              /*
                 (*testout) << nsuround << "-swap, old el = "
                 << rel << endl;
                 */
              rel.Delete();
              for (int k1 = 0; k1 < 4; k1++)
                  rel[k1].Invalidate();
            }
        }
    }
  return d_badness;
}

void MeshOptimize3d :: SwapImprove (Mesh & mesh, OPTIMIZEGOAL goal,
				    const NgBitArray * working_elements)
{
  static Timer t("MeshOptimize3d::SwapImprove"); RegionTimer reg(t);
  static Timer tloop("MeshOptimize3d::SwapImprove loop");

  // return SwapImproveSequential(mesh, goal, working_elements);

  int cnt = 0;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();

  mesh.BoundaryEdge (1,2); // ensure the boundary-elements table is built

  auto elementsonnode = mesh.CreatePoint2ElementTable();

  NgArray<ElementIndex> hasbothpoints;

  PrintMessage (3, "SwapImprove ");
  (*testout) << "\n" << "Start SwapImprove" << endl;

  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Swap Improve";

  INDEX_3_HASHTABLE<int> faces(mesh.GetNOpenElements()/3 + 2);
  if (goal == OPT_CONFORM)
    {
      for (int i = 1; i <= mesh.GetNOpenElements(); i++)
	{
	  const Element2d & hel = mesh.OpenElement(i);
	  INDEX_3 face(hel[0], hel[1], hel[2]);
	  face.Sort();
	  faces.Set (face, 1);
	}
    }

  // Calculate total badness
  if (goal == OPT_QUALITY)
    {
      double bad1 = mesh.CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad1 << endl;
    }

  Array<std::tuple<PointIndex,PointIndex>> edges;
  BuildEdgeList(mesh, elementsonnode, edges);

  Array<std::tuple<double, int>> candidate_edges(edges.Size());
  std::atomic<int> improvement_counter(0);

  tloop.Start();

  ParallelForRange(Range(edges), [&] (auto myrange)
  {
    for(auto i : myrange)
    {
      if (multithread.terminate)
        break;

      auto [pi0, pi1] = edges[i];
      double d_badness = SwapImproveEdge (mesh, goal, working_elements, elementsonnode, faces, pi0, pi1, true);
      if(d_badness<0.0)
      {
        int index = improvement_counter++;
        candidate_edges[index] = make_tuple(d_badness, i);
      }
    }
  });

  auto edges_with_improvement = candidate_edges.Part(0, improvement_counter.load());
  QuickSort(edges_with_improvement);

  for(auto [d_badness, ei] : edges_with_improvement)
  {
      auto [pi0,pi1] = edges[ei];
      if(SwapImproveEdge (mesh, goal, working_elements, elementsonnode, faces, pi0, pi1, false) < 0.0)
          cnt++;
  }

  tloop.Stop();

  PrintMessage (5, cnt, " swaps performed");

  mesh.Compress ();

  multithread.task = savetask;
}
  





void MeshOptimize3d :: SwapImproveSurface (Mesh & mesh, OPTIMIZEGOAL goal,
					   const NgBitArray * working_elements,
					   const NgArray< NgArray<int,PointIndex::BASE>* > * idmaps)
{
  NgArray< NgArray<int,PointIndex::BASE>* > locidmaps;
  const NgArray< NgArray<int,PointIndex::BASE>* > * used_idmaps;

  if(idmaps)
    used_idmaps = idmaps;
  else
    {
      used_idmaps = &locidmaps;
      
      for(int i=1; i<=mesh.GetIdentifications().GetMaxNr(); i++)
	{
	  if(mesh.GetIdentifications().GetType(i) == Identifications::PERIODIC)
	    {
	      locidmaps.Append(new NgArray<int,PointIndex::BASE>);
	      mesh.GetIdentifications().GetMap(i,*locidmaps.Last(),true);
	    }
	}
    }


  PointIndex pi1, pi2, pi3, pi4, pi5, pi6;
  PointIndex pi1other, pi2other;
  int cnt = 0;

  //double bad1, bad2, bad3, sbad;
  double bad1, sbad;
  double h;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();

  int mattype, othermattype;

  
  // contains at least all elements at node
  TABLE<ElementIndex,PointIndex::BASE> elementsonnode(np);
  TABLE<SurfaceElementIndex,PointIndex::BASE> surfaceelementsonnode(np);
  TABLE<int,PointIndex::BASE> surfaceindicesonnode(np);

  NgArray<ElementIndex> hasbothpoints;
  NgArray<ElementIndex> hasbothpointsother;

  PrintMessage (3, "SwapImproveSurface ");
  (*testout) << "\n" << "Start SwapImproveSurface" << endl;

  const char * savetask = multithread.task;
  multithread.task = "Swap Improve Surface";
    
      
  
  // find elements on node
  for (ElementIndex ei = 0; ei < ne; ei++)
    for (int j = 0; j < mesh[ei].GetNP(); j++)
      elementsonnode.Add (mesh[ei][j], ei);

  for (SurfaceElementIndex sei = 0; sei < nse; sei++)
    for(int j=0; j<mesh[sei].GetNP(); j++)
      {
	surfaceelementsonnode.Add(mesh[sei][j], sei);
	if(!surfaceindicesonnode[mesh[sei][j]].Contains(mesh[sei].GetIndex()))
	  surfaceindicesonnode.Add(mesh[sei][j],mesh[sei].GetIndex());
      }

  bool periodic;
  int idnum(-1);

  // INDEX_2_HASHTABLE<int> edgeused(2 * ne + 5);
  INDEX_2_CLOSED_HASHTABLE<int> edgeused(12 * ne + 5);

  for (ElementIndex ei = 0; ei < ne; ei++)
    {
      if (multithread.terminate)
	break;
      
      multithread.percent = 100.0 * (ei+1) / ne;

      if (mesh.ElementType(ei) == FIXEDELEMENT)
	continue;
      
      if(working_elements && 
	 ei < working_elements->Size() &&
	 !working_elements->Test(ei))
	continue;

      if (mesh[ei].IsDeleted())
	continue;

      if ((goal == OPT_LEGAL) && 
	  mesh.LegalTet (mesh[ei]) &&
	  CalcBad (mesh.Points(), mesh[ei], 0) < 1e3)
	continue;

      const Element & elemi = mesh[ei];
      //Element elemi = mesh[ei];
      if (elemi.IsDeleted()) continue;


      mattype = elemi.GetIndex();

      bool swapped = false;

      for (int j = 0; !swapped && j < 6; j++)
	{
	  // loop over edges

	  
	  static const int tetedges[6][2] =
	    { { 0, 1 }, { 0, 2 }, { 0, 3 },
	      { 1, 2 }, { 1, 3 }, { 2, 3 } };

	  pi1 = elemi[tetedges[j][0]];
	  pi2 = elemi[tetedges[j][1]];

	  
	  if (pi2 < pi1)
	    Swap (pi1, pi2);
	    	  
	  
	  bool found = false;
	  for(int k=0; !found && k<used_idmaps->Size(); k++)
	    {
	      if(pi2 < (*used_idmaps)[k]->Size() + PointIndex::BASE)
		{
		  pi1other = (*(*used_idmaps)[k])[pi1];
		  pi2other = (*(*used_idmaps)[k])[pi2];
		  found = (pi1other != 0 && pi2other != 0 && pi1other != pi1 && pi2other != pi2);
		  if(found)
		    idnum = k;
		}
	    }
	  if(found)
	    periodic = true;
	  else
	    {
	      periodic = false;
	      pi1other = pi1; pi2other = pi2;
	    }


	 	  
	  if (!mesh.BoundaryEdge (pi1, pi2) ||
	      mesh.IsSegment(pi1, pi2)) continue;

	  othermattype = -1;

	  
	  INDEX_2 i2 (pi1, pi2);
	  i2.Sort();
	  if (edgeused.Used(i2)) continue;
	  edgeused.Set (i2, 1);
	  if(periodic)
	    {
	      i2.I1() = pi1other;
	      i2.I2() = pi2other;
	      i2.Sort();
	      edgeused.Set(i2,1);
	    }
	  
	  
	  hasbothpoints.SetSize (0);
	  hasbothpointsother.SetSize (0);
	  for (int k = 0; k < elementsonnode[pi1].Size(); k++)
	    {
	      bool has1 = false, has2 = false;
	      ElementIndex elnr = elementsonnode[pi1][k];
	      const Element & elem = mesh[elnr];
	      
	      if (elem.IsDeleted()) continue;
	      
	      for (int l = 0; l < elem.GetNP(); l++)
		{
		  if (elem[l] == pi1) has1 = true;
		  if (elem[l] == pi2) has2 = true;
		}

	      if (has1 && has2) 
		{ 
		  if(othermattype == -1 && elem.GetIndex() != mattype)
		    othermattype = elem.GetIndex();

		  if(elem.GetIndex() == mattype)
		    {
		      // only once
		      for (int l = 0; l < hasbothpoints.Size(); l++)
			if (hasbothpoints[l] == elnr)
			  has1 = 0;
		      
		      if (has1)
			hasbothpoints.Append (elnr);
		    }
		  else if(elem.GetIndex() == othermattype)
		    {
		      // only once
		      for (int l = 0; l < hasbothpointsother.Size(); l++)
			if (hasbothpointsother[l] == elnr)
			  has1 = 0;
		      
		      if (has1)
			hasbothpointsother.Append (elnr);
		    }
		  else
		    {
		      cout << "problem with domain indices" << endl;
		      (*testout) << "problem: mattype = " << mattype << ", othermattype = " << othermattype 
				 << " elem " << elem << " mt " << elem.GetIndex() << endl
				 << " pi1 " << pi1 << " pi2 " << pi2 << endl;
		      (*testout) << "hasbothpoints:" << endl;
		      for(int ii=0; ii < hasbothpoints.Size(); ii++)
			(*testout) << mesh[hasbothpoints[ii]] << endl;
		      (*testout) << "hasbothpointsother:" << endl;
		      for(int ii=0; ii < hasbothpointsother.Size(); ii++)
			(*testout) << mesh[hasbothpointsother[ii]] << endl;
		    }
		}
	    }

	  if(hasbothpointsother.Size() > 0 && periodic)
	    throw NgException("SwapImproveSurface: Assumption about interface/periodicity wrong!");

	  if(periodic)
	    {
	      for (int k = 0; k < elementsonnode[pi1other].Size(); k++)
		{
		  bool has1 = false, has2 = false;
		  ElementIndex elnr = elementsonnode[pi1other][k];
		  const Element & elem = mesh[elnr];
	      
		  if (elem.IsDeleted()) continue;
	      
		  for (int l = 0; l < elem.GetNP(); l++)
		    {
		      if (elem[l] == pi1other) has1 = true;
		      if (elem[l] == pi2other) has2 = true;
		    }
		  
		  if (has1 && has2) 
		    { 
		      if(othermattype == -1)
			othermattype = elem.GetIndex();

		      // only once
		      for (int l = 0; l < hasbothpointsother.Size(); l++)
			if (hasbothpointsother[l] == elnr)
			  has1 = 0;
		      
		      if (has1)
			hasbothpointsother.Append (elnr);
		    }
		}
	    }


	  //for(k=0; k<hasbothpoints.Size(); k++)
	  //  (*testout) << "hasbothpoints["<<k<<"]: " << mesh[hasbothpoints[k]] << endl;

	  
	  SurfaceElementIndex sel1=-1,sel2=-1;
	  SurfaceElementIndex sel1other=-1,sel2other=-1;
	  for(int k = 0; k < surfaceelementsonnode[pi1].Size(); k++)
	    {
	      bool has1 = false, has2 = false;
	      SurfaceElementIndex elnr = surfaceelementsonnode[pi1][k];
	      const Element2d & elem = mesh[elnr];

	      if (elem.IsDeleted()) continue;

	      for (int l = 0; l < elem.GetNP(); l++)
		{
		  if (elem[l] == pi1) has1 = true;
		  if (elem[l] == pi2) has2 = true;
		}

	      if(has1 && has2 && elnr != sel2)
		{
		  sel1 = sel2;
		  sel2 = elnr;
		}
	    }

	  if(periodic)
	    {
	      for(int k = 0; k < surfaceelementsonnode[pi1other].Size(); k++)
		{
		  bool has1 = false, has2 = false;
		  SurfaceElementIndex elnr = surfaceelementsonnode[pi1other][k];
		  const Element2d & elem = mesh[elnr];

		  if (elem.IsDeleted()) continue;

		  for (int l = 0; l < elem.GetNP(); l++)
		    {
		      if (elem[l] == pi1other) has1 = true;
		      if (elem[l] == pi2other) has2 = true;
		    }

		  if(has1 && has2 && elnr != sel2other)
		    {
		      sel1other = sel2other;
		      sel2other = elnr;
		    }
		}
	    }
	  else
	    {
	      sel1other = sel1; sel2other = sel2;
	    }

	  //(*testout) << "sel1 " << sel1 << " sel2 " << sel2 << " el " << mesh[sel1] << " resp. " << mesh[sel2] << endl;

	  PointIndex sp1(0), sp2(0);
	  PointIndex sp1other, sp2other;
	  for(int l=0; l<mesh[sel1].GetNP(); l++)
	    if(mesh[sel1][l] != pi1 && mesh[sel1][l] != pi2)
	      sp1 = mesh[sel1][l];
	  for(int l=0; l<mesh[sel2].GetNP(); l++)
	    if(mesh[sel2][l] != pi1 && mesh[sel2][l] != pi2)
	      sp2 = mesh[sel2][l];

	  if(periodic)
	    {
	      sp1other = (*(*used_idmaps)[idnum])[sp1];
	      sp2other = (*(*used_idmaps)[idnum])[sp2];

	      bool change = false;
	      for(int l=0; !change && l<mesh[sel1other].GetNP(); l++)
		change = (sp2other == mesh[sel1other][l]);
	      
	      if(change)
		{
		  SurfaceElementIndex aux = sel1other;
		  sel1other = sel2other;
		  sel2other = aux;
		}

	    }
	  else
	    {
	      sp1other = sp1; sp2other = sp2;
	    }
	  
	  Vec<3> v1 = mesh[sp1]-mesh[pi1],
	    v2 = mesh[sp2]-mesh[pi1],
	    v3 = mesh[sp1]-mesh[pi2],
	    v4 = mesh[sp2]-mesh[pi2];
	  double vol = 0.5*(Cross(v1,v2).Length() + Cross(v3,v4).Length());
	  h = sqrt(vol);
	  h = 0;

	  sbad = CalcTriangleBadness (mesh[pi1],mesh[pi2],mesh[sp1],0,0) + 
	    CalcTriangleBadness (mesh[pi2],mesh[pi1],mesh[sp2],0,0);
	  


	  bool puretet = true;
	  for (int k = 0; puretet && k < hasbothpoints.Size(); k++)
	    if (mesh[hasbothpoints[k]].GetType () != TET)
	      puretet = false;
	  for (int k = 0; puretet && k < hasbothpointsother.Size(); k++)
	    if (mesh[hasbothpointsother[k]].GetType () != TET)
	      puretet = false;
	  if (!puretet)
	    continue;

	  int nsuround = hasbothpoints.Size();
	  int nsuroundother = hasbothpointsother.Size();

	  NgArray < int > outerpoints(nsuround+1);
	  outerpoints[0] = sp1;

	  for(int i=0; i<nsuround; i++)
	    {
	      bool done = false;
	      for(int jj=i; !done && jj<hasbothpoints.Size(); jj++)
		{
		  for(int k=0; !done && k<4; k++)
		    if(mesh[hasbothpoints[jj]][k] == outerpoints[i])
		      {
			done = true;
			for(int l=0; l<4; l++)
			  if(mesh[hasbothpoints[jj]][l] != pi1 &&
			     mesh[hasbothpoints[jj]][l] != pi2 &&
			     mesh[hasbothpoints[jj]][l] != outerpoints[i])
			    outerpoints[i+1] = mesh[hasbothpoints[jj]][l];
		      }
		  if(done)
		    {
		      ElementIndex aux = hasbothpoints[i];
		      hasbothpoints[i] = hasbothpoints[jj];
		      hasbothpoints[jj] = aux;
		    }
		}
	    }
	  if(outerpoints[nsuround] != sp2)
	    {
	      cerr << "OJE OJE OJE" << endl;
	      (*testout) << "OJE OJE OJE" << endl;
	      (*testout) << "hasbothpoints: " << endl;
	      for(int ii=0; ii < hasbothpoints.Size(); ii++)
		{
		  (*testout) << mesh[hasbothpoints[ii]] << endl;
		  for(int jj=0; jj<mesh[hasbothpoints[ii]].GetNP(); jj++)
		    if(mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] > 0)
		      (*testout) << mesh[hasbothpoints[ii]][jj] << " between "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] << " and "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][1] << endl;
		}
	      (*testout) << "outerpoints: " << outerpoints << endl;
	      (*testout) << "sel1 " << mesh[sel1] << endl
			 << "sel2 " << mesh[sel2] << endl;
	      for(int ii=0; ii<3; ii++)
		{
		  if(mesh.mlbetweennodes[mesh[sel1][ii]][0] > 0)
		    (*testout) << mesh[sel1][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][1] << endl;
		  if(mesh.mlbetweennodes[mesh[sel2][ii]][0] > 0)
		    (*testout) << mesh[sel2][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][1] << endl;
		}
	    }

	  
	  NgArray < int > outerpointsother;

	  if(nsuroundother > 0)
	    {
	      outerpointsother.SetSize(nsuroundother+1);
	      outerpointsother[0] = sp2other;
	    }

	  for(int i=0; i<nsuroundother; i++)
	    {
	      bool done = false;
	      for(int jj=i; !done && jj<hasbothpointsother.Size(); jj++)
		{
		  for(int k=0; !done && k<4; k++)
		    if(mesh[hasbothpointsother[jj]][k] == outerpointsother[i])
		      {
			done = true;
			for(int l=0; l<4; l++)
			  if(mesh[hasbothpointsother[jj]][l] != pi1other &&
			     mesh[hasbothpointsother[jj]][l] != pi2other &&
			     mesh[hasbothpointsother[jj]][l] != outerpointsother[i])
			    outerpointsother[i+1] = mesh[hasbothpointsother[jj]][l];
		      }
		  if(done)
		    {
		      ElementIndex aux = hasbothpointsother[i];
		      hasbothpointsother[i] = hasbothpointsother[jj];
		      hasbothpointsother[jj] = aux;
		    }
		}
	    }
	  if(nsuroundother > 0 && outerpointsother[nsuroundother] != sp1other)
	    {
	      cerr << "OJE OJE OJE (other)" << endl;
	      (*testout) << "OJE OJE OJE (other)" << endl;
	      (*testout) << "pi1 " << pi1 << " pi2 " << pi2 << " sp1 " << sp1 << " sp2 " << sp2 << endl;
	      (*testout) << "hasbothpoints: " << endl;
	      for(int ii=0; ii < hasbothpoints.Size(); ii++)
		{
		  (*testout) << mesh[hasbothpoints[ii]] << endl;
		  for(int jj=0; jj<mesh[hasbothpoints[ii]].GetNP(); jj++)
		    if(mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] > 0)
		      (*testout) << mesh[hasbothpoints[ii]][jj] << " between "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][0] << " and "
				 << mesh.mlbetweennodes[mesh[hasbothpoints[ii]][jj]][1] << endl;
		}
	      (*testout) << "outerpoints: " << outerpoints << endl;
	      (*testout) << "sel1 " << mesh[sel1] << endl
			 << "sel2 " << mesh[sel2] << endl;
	      for(int ii=0; ii<3; ii++)
		{
		  if(mesh.mlbetweennodes[mesh[sel1][ii]][0] > 0)
		    (*testout) << mesh[sel1][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel1][ii]][1] << endl;
		  if(mesh.mlbetweennodes[mesh[sel2][ii]][0] > 0)
		    (*testout) << mesh[sel2][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel2][ii]][1] << endl;
		}
		  
	      (*testout) << "pi1other " << pi1other << " pi2other " << pi2other << " sp1other " << sp1other << " sp2other " << sp2other << endl;
	      (*testout) << "hasbothpointsother: " << endl;
	      for(int ii=0; ii < hasbothpointsother.Size(); ii++)
		{
		  (*testout) << mesh[hasbothpointsother[ii]] << endl;
		  for(int jj=0; jj<mesh[hasbothpointsother[ii]].GetNP(); jj++)
		    if(mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][0] > 0)
		      (*testout) << mesh[hasbothpointsother[ii]][jj] << " between "
				 << mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][0] << " and "
				 << mesh.mlbetweennodes[mesh[hasbothpointsother[ii]][jj]][1] << endl;
		}
	      (*testout) << "outerpoints: " << outerpointsother << endl;
	      (*testout) << "sel1other " << mesh[sel1other] << endl
			 << "sel2other " << mesh[sel2other] << endl;
	      for(int ii=0; ii<3; ii++)
		{
		  if(mesh.mlbetweennodes[mesh[sel1other][ii]][0] > 0)
		    (*testout) << mesh[sel1other][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel1other][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel1other][ii]][1] << endl;
		  if(mesh.mlbetweennodes[mesh[sel2other][ii]][0] > 0)
		    (*testout) << mesh[sel2other][ii] << " between "
			       << mesh.mlbetweennodes[mesh[sel2other][ii]][0] << " and "
			       << mesh.mlbetweennodes[mesh[sel2other][ii]][1] << endl;
		}
	    }

	  bad1=0;
	  for(int i=0; i<hasbothpoints.Size(); i++)
	    bad1 += CalcBad(mesh.Points(), mesh[hasbothpoints[i]],h);
	  for(int i=0; i<hasbothpointsother.Size(); i++)
	    bad1 += CalcBad(mesh.Points(), mesh[hasbothpointsother[i]],h);
	  bad1 /= double(hasbothpoints.Size() + hasbothpointsother.Size());

	  
	  int startpoints,startpointsother;


	  if(outerpoints.Size() == 3)
	    startpoints = 1;
	  else if(outerpoints.Size() == 4)
	    startpoints = 2;
	  else
	    startpoints = outerpoints.Size();
	  
	  if(outerpointsother.Size() == 3)
	    startpointsother = 1;
	  else if(outerpointsother.Size() == 4)
	    startpointsother = 2;
	  else
	    startpointsother = outerpointsother.Size();
	  

	  NgArray < NgArray < Element* > * > newelts(startpoints);
	  NgArray < NgArray < Element* > * > neweltsother(startpointsother);

	  double minbad = 1e50, minbadother = 1e50, currbad;
	  int minpos = -1, minposother = -1;

	  //(*testout) << "pi1 " << pi1 << " pi2 " << pi2 << " outerpoints " << outerpoints << endl;

	  for(int i=0; i<startpoints; i++)
	    {
	      newelts[i] = new NgArray <Element*>(2*(nsuround-1));
	      
	      for(int jj=0; jj<nsuround-1; jj++)
		{
		  (*newelts[i])[2*jj] = new Element(TET);
		  (*newelts[i])[2*jj+1] = new Element(TET);
		  Element & newel1 = *((*newelts[i])[2*jj]);
		  Element & newel2 = *((*newelts[i])[2*jj+1]);

		  newel1[0] = pi1;
		  newel1[1] = outerpoints[i];
		  newel1[2] = outerpoints[(i+jj+1)%outerpoints.Size()];
		  newel1[3] = outerpoints[(i+jj+2)%outerpoints.Size()];

		  newel2[0] = pi2;
		  newel2[1] = outerpoints[i];
		  newel2[2] = outerpoints[(i+jj+2)%outerpoints.Size()];
		  newel2[3] = outerpoints[(i+jj+1)%outerpoints.Size()];
		  

		  //(*testout) << "j " << j << " newel1 " << newel1[0] << " "<< newel1[1] << " "<< newel1[2] << " "<< newel1[3] << endl
		  //     << " newel2 " << newel2[0] << " "<< newel2[1] << " "<< newel2[2] << " "<< newel2[3] << endl;
		  
		  newel1.SetIndex(mattype);
		  newel2.SetIndex(mattype);

		}

	      bool wrongorientation = true;
	      for(int jj = 0; wrongorientation && jj<newelts[i]->Size(); jj++)
		wrongorientation = wrongorientation && WrongOrientation(mesh.Points(), *(*newelts[i])[jj]);
	      
	      currbad = 0;

	      for(int jj=0; jj<newelts[i]->Size(); jj++)
		{
		  if(wrongorientation)
		    Swap((*(*newelts[i])[jj])[2],(*(*newelts[i])[jj])[3]);


		  // not two new faces on same surface
		  NgArray<int> face_index;
		  for(int k = 0; k<surfaceindicesonnode[(*(*newelts[i])[jj])[0]].Size(); k++)
		    face_index.Append(surfaceindicesonnode[(*(*newelts[i])[jj])[0]][k]);

		  for(int k=1; k<4; k++)
		    {
		      for(int l=0; l<face_index.Size(); l++)
			{
			  if(face_index[l] != -1 && 
			     !(surfaceindicesonnode[(*(*newelts[i])[jj])[k]].Contains(face_index[l])))
			    face_index[l] = -1;
			}

		    }
		      
		  for(int k=0; k<face_index.Size(); k++)
		    if(face_index[k] != -1)
		      currbad += 1e12;


		  currbad += CalcBad(mesh.Points(),*(*newelts[i])[jj],h);


		}  

	      //currbad /= double(newelts[i]->Size());
		    


	      if(currbad < minbad)
		{
		  minbad = currbad;
		  minpos = i;
		}

	    }

	  if(startpointsother == 0)
	    minbadother = 0;

	  for(int i=0; i<startpointsother; i++)
	    {
	      neweltsother[i] = new NgArray <Element*>(2*(nsuroundother));
	      
	      for(int jj=0; jj<nsuroundother; jj++)
		{
		  (*neweltsother[i])[2*jj] = new Element(TET);
		  (*neweltsother[i])[2*jj+1] = new Element(TET);
		  Element & newel1 = *((*neweltsother[i])[2*jj]);
		  Element & newel2 = *((*neweltsother[i])[2*jj+1]);

		  newel1[0] = pi1other;
		  newel1[1] = outerpointsother[i];
		  newel1[2] = outerpointsother[(i+jj+1)%outerpointsother.Size()];
		  newel1[3] = outerpointsother[(i+jj+2)%outerpointsother.Size()];

		  newel2[0] = pi2other;
		  newel2[1] = outerpointsother[i];
		  newel2[2] = outerpointsother[(i+jj+2)%outerpointsother.Size()];
		  newel2[3] = outerpointsother[(i+jj+1)%outerpointsother.Size()];
		  

		  //(*testout) << "j " << j << " newel1 " << newel1[0] << " "<< newel1[1] << " "<< newel1[2] << " "<< newel1[3] << endl
		  //	     << " newel2 " << newel2[0] << " "<< newel2[1] << " "<< newel2[2] << " "<< newel2[3] << endl;
		  
		  newel1.SetIndex(othermattype);
		  newel2.SetIndex(othermattype);

		}

	      bool wrongorientation = true;
	      for(int jj = 0; wrongorientation && jj<neweltsother[i]->Size(); jj++)
		wrongorientation = wrongorientation && WrongOrientation(mesh.Points(), *(*neweltsother[i])[jj]);
	      
	      currbad = 0;

	      for(int jj=0; jj<neweltsother[i]->Size(); jj++)
		{
		  if(wrongorientation)
		    Swap((*(*neweltsother[i])[jj])[2],(*(*neweltsother[i])[jj])[3]);

		  currbad += CalcBad(mesh.Points(),*(*neweltsother[i])[jj],h);
		}  

	      //currbad /= double(neweltsother[i]->Size());
		    


	      if(currbad < minbadother)
		{
		  minbadother = currbad;
		  minposother = i;
		}

	    }

	  //(*testout) << "minbad " << minbad << " bad1 " << bad1 << endl;

	  
	  double sbadnew = CalcTriangleBadness (mesh[pi1],mesh[sp2],mesh[sp1],0,0) + 
	    CalcTriangleBadness (mesh[pi2],mesh[sp1],mesh[sp2],0,0);
	  

	  int denom = newelts[minpos]->Size();
	  if(minposother >= 0)
	    denom += neweltsother[minposother]->Size();
	  

	  if((minbad+minbadother)/double(denom) < bad1 && 
	     sbadnew < sbad)
	    {
	      cnt++;

	      swapped = true;


	      int start1 = -1;
	      for(int l=0; l<3; l++)
		if(mesh[sel1][l] == pi1)
		  start1 = l;
	      if(mesh[sel1][(start1+1)%3] == pi2)
		{
		  mesh[sel1][0] = pi1;
		  mesh[sel1][1] = sp2;
		  mesh[sel1][2] = sp1;
		  mesh[sel2][0] = pi2;
		  mesh[sel2][1] = sp1;
		  mesh[sel2][2] = sp2;
		}
	      else
		{
		  mesh[sel1][0] = pi2;
		  mesh[sel1][1] = sp2;
		  mesh[sel1][2] = sp1;
		  mesh[sel2][0] = pi1;
		  mesh[sel2][1] = sp1;
		  mesh[sel2][2] = sp2;
		}
	      //(*testout) << "changed surface element " << sel1 << " to " << mesh[sel1] << ", " << sel2 << " to " << mesh[sel2] << endl;

	      for(int l=0; l<3; l++)
		{
		  surfaceelementsonnode.Add(mesh[sel1][l],sel1);
		  surfaceelementsonnode.Add(mesh[sel2][l],sel2);
		}
	      


	      if(periodic)
		{
		  start1 = -1;
		  for(int l=0; l<3; l++)
		    if(mesh[sel1other][l] == pi1other)
		      start1 = l;
		  


		  //(*testout) << "changed surface elements " << mesh[sel1other] << " and " << mesh[sel2other] << endl;
		  if(mesh[sel1other][(start1+1)%3] == pi2other)
		    {
		      mesh[sel1other][0] = pi1other;
		      mesh[sel1other][1] = sp2other;
		      mesh[sel1other][2] = sp1other;
		      mesh[sel2other][0] = pi2other;
		      mesh[sel2other][1] = sp1other;
		      mesh[sel2other][2] = sp2other;
		      //(*testout) << "       with rule 1" << endl;
		    }
		  else
		    {
		      mesh[sel1other][0] = pi2other;
		      mesh[sel1other][1] = sp2other;
		      mesh[sel1other][2] = sp1other;
		      mesh[sel2other][0] = pi1other;
		      mesh[sel2other][1] = sp1other;
		      mesh[sel2other][2] = sp2other;
		      //(*testout) << "       with rule 2" << endl;
		    }
		  //(*testout) << "         to " << mesh[sel1other] << " and " << mesh[sel2other] << endl;
		  
		  //(*testout) << "  and surface element " << sel1other << " to " << mesh[sel1other] << ", " << sel2other << " to " << mesh[sel2other] << endl;

		  for(int l=0; l<3; l++)
		    {
		      surfaceelementsonnode.Add(mesh[sel1other][l],sel1other);
		      surfaceelementsonnode.Add(mesh[sel2other][l],sel2other);
		    }
		}




	      for(int i=0; i<hasbothpoints.Size(); i++)
		{
		  mesh[hasbothpoints[i]] = *(*newelts[minpos])[i];

		  for(int l=0; l<4; l++)
		    elementsonnode.Add((*(*newelts[minpos])[i])[l],hasbothpoints[i]);
		}

	      for(int i=hasbothpoints.Size(); i<(*newelts[minpos]).Size(); i++)
		{
		  ElementIndex ni = mesh.AddVolumeElement(*(*newelts[minpos])[i]);
		  
		  for(int l=0; l<4; l++)
		    elementsonnode.Add((*(*newelts[minpos])[i])[l],ni);
		}

	      if(hasbothpointsother.Size() > 0)
		{
		  for(int i=0; i<hasbothpointsother.Size(); i++)
		    {
		      mesh[hasbothpointsother[i]] = *(*neweltsother[minposother])[i];
		      for(int l=0; l<4; l++)
			elementsonnode.Add((*(*neweltsother[minposother])[i])[l],hasbothpointsother[i]);
		    }
		  
		  for(int i=hasbothpointsother.Size(); i<(*neweltsother[minposother]).Size(); i++)
		    {
		      ElementIndex ni = mesh.AddVolumeElement(*(*neweltsother[minposother])[i]);
		      for(int l=0; l<4; l++)
			elementsonnode.Add((*(*neweltsother[minposother])[i])[l],ni);
		    }
		}

	      

	    }

	  for(int i=0; i<newelts.Size(); i++)
	    {
	      for(int jj=0; jj<newelts[i]->Size(); jj++)
		delete (*newelts[i])[jj];
	      delete newelts[i];
	    }

	  for(int i=0; i<neweltsother.Size(); i++)
	    {
	      for(int jj=0; jj<neweltsother[i]->Size(); jj++)
		delete (*neweltsother[i])[jj];
	      delete neweltsother[i];
	    }
	
	}
    }

  PrintMessage (5, cnt, " swaps performed");


  for(int i=0; i<locidmaps.Size(); i++)
    delete locidmaps[i];


  mesh.Compress ();

  multithread.task = savetask;
}
  







/*
  2 -> 3 conversion
*/

double MeshOptimize3d :: SwapImprove2 ( Mesh & mesh, OPTIMIZEGOAL goal, ElementIndex eli1, int face,
  Table<ElementIndex, PointIndex> & elementsonnode,
  TABLE<SurfaceElementIndex, PointIndex::BASE> & belementsonnode, bool check_only )
{
  PointIndex pi1, pi2, pi3, pi4, pi5;
  Element el21(TET), el22(TET), el31(TET), el32(TET), el33(TET);
  int j = face;
  double bad1, bad2;
  double d_badness = 0.0;

  Element & elem = mesh[eli1];
  if (elem.IsDeleted()) return 0.0;

  int mattyp = elem.GetIndex();

  switch (j)
  {
    case 0:
      pi1 = elem.PNum(1); pi2 = elem.PNum(2);
      pi3 = elem.PNum(3); pi4 = elem.PNum(4);
      break;
    case 1:
      pi1 = elem.PNum(1); pi2 = elem.PNum(4);
      pi3 = elem.PNum(2); pi4 = elem.PNum(3);
      break;
    case 2:
      pi1 = elem.PNum(1); pi2 = elem.PNum(3);
      pi3 = elem.PNum(4); pi4 = elem.PNum(2);
      break;
    case 3:
      pi1 = elem.PNum(2); pi2 = elem.PNum(4);
      pi3 = elem.PNum(3); pi4 = elem.PNum(1);
      break;
  }


  bool bface = 0;
  for (int k = 0; k < belementsonnode[pi1].Size(); k++)
  {
      const Element2d & bel =
        mesh[belementsonnode[pi1][k]];

      bool bface1 = 1;
      for (int l = 0; l < 3; l++)
          if (bel[l] != pi1 && bel[l] != pi2 && bel[l] != pi3)
          {
              bface1 = 0;
              break;
          }

      if (bface1)
      {
          bface = 1;
          break;
      }
  }

  if (bface) return 0.0;


  FlatArray<ElementIndex> row = elementsonnode[pi1];
  for (int k = 0; k < row.Size(); k++)
  {
      ElementIndex eli2 = row[k];

      if ( eli1 != eli2 )
      {
          Element & elem2 = mesh[eli2];
          if (elem2.IsDeleted()) continue;
          if (elem2.GetType() != TET)
              continue;

          int comnodes=0;
          for (int l = 1; l <= 4; l++)
              if (elem2.PNum(l) == pi1 || elem2.PNum(l) == pi2 ||
                  elem2.PNum(l) == pi3)
              {
                  comnodes++;
              }
              else
              {
                  pi5 = elem2.PNum(l);
              }

          if (comnodes == 3)
          {
              bad1 = CalcBad (mesh.Points(), elem, 0) +
                CalcBad (mesh.Points(), elem2, 0);

              if (!mesh.LegalTet(elem) ||
                  !mesh.LegalTet(elem2))
                  bad1 += 1e4;


              el31.PNum(1) = pi1;
              el31.PNum(2) = pi2;
              el31.PNum(3) = pi5;
              el31.PNum(4) = pi4;
              el31.SetIndex (mattyp);

              el32.PNum(1) = pi2;
              el32.PNum(2) = pi3;
              el32.PNum(3) = pi5;
              el32.PNum(4) = pi4;
              el32.SetIndex (mattyp);

              el33.PNum(1) = pi3;
              el33.PNum(2) = pi1;
              el33.PNum(3) = pi5;
              el33.PNum(4) = pi4;
              el33.SetIndex (mattyp);

              bad2 = CalcBad (mesh.Points(), el31, 0) +
                CalcBad (mesh.Points(), el32, 0) +
                CalcBad (mesh.Points(), el33, 0);


              el31.flags.illegal_valid = 0;
              el32.flags.illegal_valid = 0;
              el33.flags.illegal_valid = 0;

              if (!mesh.LegalTet(el31) ||
                  !mesh.LegalTet(el32) ||
                  !mesh.LegalTet(el33))
                  bad2 += 1e4;


              d_badness = bad2 - bad1;

              if ( ((bad2 < 1e6) || (bad2 < 10 * bad1)) &&
                  mesh.BoundaryEdge (pi4, pi5))
                  d_badness = -1e4;

              if(check_only)
                  return d_badness;

              if (d_badness<0.0)
              {
                  el31.flags.illegal_valid = 0;
                  el32.flags.illegal_valid = 0;
                  el33.flags.illegal_valid = 0;

                  mesh[eli1].Delete();
                  mesh[eli2].Delete();
                  mesh.AddVolumeElement (el31);
                  mesh.AddVolumeElement (el32);
                  mesh.AddVolumeElement (el33);
              }
              return d_badness;
          }
      }
  }
  return d_badness;
}

/*
  2 -> 3 conversion
*/

void MeshOptimize3d :: SwapImprove2Sequential (Mesh & mesh, OPTIMIZEGOAL goal)
{
  static Timer t("MeshOptimize3d::SwapImprove2"); RegionTimer reg(t);
  
  PointIndex pi1, pi2, pi3, pi4, pi5;
  Element el21(TET), el22(TET), el31(TET), el32(TET), el33(TET);

  int cnt = 0;
  double bad1, bad2;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();

  if (goal == OPT_CONFORM) return;

  // contains at least all elements at node
  // TABLE<ElementIndex, PointIndex::BASE> elementsonnode(np); 
  TABLE<SurfaceElementIndex, PointIndex::BASE> belementsonnode(np);

  PrintMessage (3, "SwapImprove2 ");
  (*testout) << "\n" << "Start SwapImprove2" << "\n";
  //  TestOk();


  /*
    CalcSurfacesOfNode ();
    for (i = 1; i <= GetNE(); i++)
    if (volelements.Get(i)[0])
    if (!mesh.LegalTet (volelements.Get(i)))
    {
    cout << "detected illegal tet, 1" << endl;
    (*testout) << "detected illegal tet1: " << i << endl;
    }
  */

  
  // Calculate total badness

  bad1 = mesh.CalcTotalBad (mp);
  (*testout) << "Total badness = " << bad1 << endl;
  cout << "tot bad = " << bad1 << endl;

  auto elementsonnode = mesh.CreatePoint2ElementTable();
  
  for (SurfaceElementIndex sei = 0; sei < nse; sei++)
    for (int j = 0; j < 3; j++)
      belementsonnode.Add (mesh[sei][j], sei);

  for (ElementIndex eli1 = 0; eli1 < ne; eli1++)
    {
      if (multithread.terminate)
	break;

      if (mesh.ElementType (eli1) == FIXEDELEMENT)
	continue;

      if (mesh[eli1].GetType() != TET)
	continue;

      if ((goal == OPT_LEGAL) && 
	  mesh.LegalTet (mesh[eli1]) &&
	  CalcBad (mesh.Points(), mesh[eli1], 0) < 1e3)
	continue;
      
      if(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(eli1).GetIndex())
      	continue;

      // cout << "eli = " << eli1 << endl;
      //      (*testout) << "swapimp2, eli = " << eli1 << "; el = " << mesh[eli1] << endl;

      for (int j = 0; j < 4; j++)
	{
          cnt += SwapImprove2( mesh, goal, eli1, j, elementsonnode, belementsonnode, false);
	}
    }


  PrintMessage (5, cnt, " swaps performed");



  /*
    CalcSurfacesOfNode ();
    for (i = 1; i <= GetNE(); i++)
    if (volelements.Get(i).PNum(1))
    if (!LegalTet (volelements.Get(i)))
    {
    cout << "detected illegal tet, 2" << endl;
    (*testout) << "detected illegal tet2: " << i << endl;
    }
  */


  mesh.Compress();
  bad1 = mesh.CalcTotalBad (mp);
  (*testout) << "Total badness = " << bad1 << endl;
  (*testout) << "swapimprove2 done" << "\n";
  //  (*mycout) << "Vol = " << CalcVolume (points, volelements) << "\n";
}

void MeshOptimize3d :: SwapImprove2 (Mesh & mesh, OPTIMIZEGOAL goal)
{
  static Timer t("MeshOptimize3d::SwapImprove2"); RegionTimer reg(t);

  // return SwapImprove2Sequential(mesh, goal);

  mesh.BoundaryEdge (1,2); // ensure the boundary-elements table is built

  int cnt = 0;
  double bad1, bad2;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();

  if (goal == OPT_CONFORM) return;

  // contains at least all elements at node
  TABLE<SurfaceElementIndex, PointIndex::BASE> belementsonnode(np);

  PrintMessage (3, "SwapImprove2 ");
  (*testout) << "\n" << "Start SwapImprove2" << "\n";

  bad1 = mesh.CalcTotalBad (mp);
  (*testout) << "Total badness = " << bad1 << endl;

  // find elements on node

  auto elementsonnode = mesh.CreatePoint2ElementTable();
  
  for (SurfaceElementIndex sei = 0; sei < nse; sei++)
    for (int j = 0; j < 3; j++)
      belementsonnode.Add (mesh[sei][j], sei);

  int num_threads = ngcore::TaskManager::GetNumThreads();

  Array<std::tuple<double, ElementIndex, int>> faces_with_improvement;
  Array<Array<std::tuple<double, ElementIndex, int>>> faces_with_improvement_threadlocal(num_threads);

  ParallelForRange( Range(ne), [&]( auto myrange )
      {
        int tid = ngcore::TaskManager::GetThreadId();
        auto & my_faces_with_improvement = faces_with_improvement_threadlocal[tid];
        for (ElementIndex eli1 : myrange)
          {
            if (multithread.terminate)
              break;

            if (mesh.ElementType (eli1) == FIXEDELEMENT)
              continue;

            if (mesh[eli1].GetType() != TET)
              continue;

            if ((goal == OPT_LEGAL) &&
                mesh.LegalTet (mesh[eli1]) &&
                CalcBad (mesh.Points(), mesh[eli1], 0) < 1e3)
              continue;

            if(mesh.GetDimension()==3 && mp.only3D_domain_nr && mp.only3D_domain_nr != mesh.VolumeElement(eli1).GetIndex())
              continue;

            for (int j = 0; j < 4; j++)
              {
                double d_badness = SwapImprove2( mesh, goal, eli1, j, elementsonnode, belementsonnode, true);
                if(d_badness<0.0)
                    my_faces_with_improvement.Append( std::make_tuple(d_badness, eli1, j) );
              }
          }
      });

  for (auto & a : faces_with_improvement_threadlocal)
    faces_with_improvement.Append(a);

  QuickSort(faces_with_improvement);

  for (auto [dummy, eli,j] : faces_with_improvement)
      if(SwapImprove2( mesh, goal, eli, j, elementsonnode, belementsonnode, false) < 0.0)
          cnt++;

  PrintMessage (5, cnt, " swaps performed");

  mesh.Compress();
  bad1 = mesh.CalcTotalBad (mp);
  (*testout) << "Total badness = " << bad1 << endl;
  (*testout) << "swapimprove2 done" << "\n";
}


/*
  void Mesh :: SwapImprove2 (OPTIMIZEGOAL goal)
  {
  int i, j;
  int eli1, eli2;
  int mattyp;

  Element el31(4), el32(4), el33(4);
  double bad1, bad2;


  INDEX_3_HASHTABLE<INDEX_2> elsonface (GetNE());

  (*mycout) << "SwapImprove2 " << endl;
  (*testout) << "\n" << "Start SwapImprove2" << "\n";

  // Calculate total badness

  if (goal == OPT_QUALITY)
  {
  double bad1 = CalcTotalBad (points, volelements);
  (*testout) << "Total badness = " << bad1 << endl;
  }

  // find elements on node


  Element2d face;
  for (i = 1; i <= GetNE(); i++)
  if ( (i > eltyps.Size()) || (eltyps.Get(i) != FIXEDELEMENT) )
  {
  const Element & el = VolumeElement(i);
  if (!el.PNum(1)) continue;

  for (j = 1; j <= 4; j++)
  {
  el.GetFace (j, face);
  INDEX_3 i3 (face.PNum(1), face.PNum(2), face.PNum(3));
  i3.Sort();


  int bnr, posnr;
  if (!elsonface.PositionCreate (i3, bnr, posnr))
  {
  INDEX_2 i2;
  elsonface.GetData (bnr, posnr, i3, i2);
  i2.I2() = i;
  elsonface.SetData (bnr, posnr, i3, i2);
  }
  else
  {
  INDEX_2 i2 (i, 0);
  elsonface.SetData (bnr, posnr, i3, i2);
  }

  //  	    if (elsonface.Used (i3))
  //  	      {
  //  		INDEX_2 i2 = elsonface.Get(i3);
  //  		i2.I2() = i;
  //  		elsonface.Set (i3, i2);
  //  	      }
  //  	    else
  //  	      {
  //  		INDEX_2 i2 (i, 0);
  //  		elsonface.Set (i3, i2);
  //  	      }

  }
  }

  NgBitArray original(GetNE());
  original.Set();

  for (i = 1; i <= GetNSE(); i++)
  {
  const Element2d & sface = SurfaceElement(i);
  INDEX_3 i3 (sface.PNum(1), sface.PNum(2), sface.PNum(3));
  i3.Sort();
  INDEX_2 i2(0,0);
  elsonface.Set (i3, i2);
  }


  for (i = 1; i <= elsonface.GetNBags(); i++)
  for (j = 1; j <= elsonface.GetBagSize(i); j++)
  {
  INDEX_3 i3;
  INDEX_2 i2;
  elsonface.GetData (i, j, i3, i2);


  int eli1 = i2.I1();
  int eli2 = i2.I2();

  if (eli1 && eli2 && original.Test(eli1) && original.Test(eli2) )
  {
  Element & elem = volelements.Elem(eli1);
  Element & elem2 = volelements.Elem(eli2);

  int pi1 = i3.I1();
  int pi2 = i3.I2();
  int pi3 = i3.I3();

  int pi4 = elem.PNum(1) + elem.PNum(2) + elem.PNum(3) + elem.PNum(4) - pi1 - pi2 - pi3;
  int pi5 = elem2.PNum(1) + elem2.PNum(2) + elem2.PNum(3) + elem2.PNum(4) - pi1 - pi2 - pi3;






  el31.PNum(1) = pi1;
  el31.PNum(2) = pi2;
  el31.PNum(3) = pi3;
  el31.PNum(4) = pi4;
  el31.SetIndex (mattyp);
	    
  if (WrongOrientation (points, el31))
  swap (pi1, pi2);


  bad1 = CalcBad (points, elem, 0) + 
  CalcBad (points, elem2, 0); 
	    
  //	    if (!LegalTet(elem) || !LegalTet(elem2))
  //	      bad1 += 1e4;

	    
  el31.PNum(1) = pi1;
  el31.PNum(2) = pi2;
  el31.PNum(3) = pi5;
  el31.PNum(4) = pi4;
  el31.SetIndex (mattyp);
	    
  el32.PNum(1) = pi2;
  el32.PNum(2) = pi3;
  el32.PNum(3) = pi5;
  el32.PNum(4) = pi4;
  el32.SetIndex (mattyp);
		      
  el33.PNum(1) = pi3;
  el33.PNum(2) = pi1;
  el33.PNum(3) = pi5;
  el33.PNum(4) = pi4;
  el33.SetIndex (mattyp);
	    
  bad2 = CalcBad (points, el31, 0) + 
  CalcBad (points, el32, 0) +
  CalcBad (points, el33, 0); 
	    
  //	    if (!LegalTet(el31) || !LegalTet(el32) ||
  //		!LegalTet(el33))
  //	      bad2 += 1e4;
	    
	    
  int swap = (bad2 < bad1);

  INDEX_2 hi2b(pi4, pi5);
  hi2b.Sort();
	    
  if ( ((bad2 < 1e6) || (bad2 < 10 * bad1)) &&
  boundaryedges->Used (hi2b) )
  swap = 1;
	    
  if (swap)
  {
  (*mycout) << "2->3 " << flush;
		
  volelements.Elem(eli1) = el31;
  volelements.Elem(eli2) = el32;
  volelements.Append (el33);
		
  original.Clear (eli1);
  original.Clear (eli2);
  }
  }
  }
  
  (*mycout) << endl;

  if (goal == OPT_QUALITY)
  {
  bad1 = CalcTotalBad (points, volelements);
  (*testout) << "Total badness = " << bad1 << endl;
  }

  //  FindOpenElements ();

  (*testout) << "swapimprove2 done" << "\n";
  }

*/
}

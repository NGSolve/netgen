#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <gprim.hpp>

#include <meshing.hpp>


#include "stlgeom.hpp"


namespace netgen
{

static void STLFindEdges (STLGeometry & geom, Mesh & mesh,
                          const MeshingParameters& mparam,
                          const STLParameters& stlparam)
{
  double h = mparam.maxh;

  // mark edge points:
  //int ngp = geom.GetNP();

  geom.RestrictLocalH(mesh, h, stlparam, mparam);
  
  PushStatusF("Mesh Lines");

  NgArray<STLLine*> meshlines;
  NgArray<Point3d> meshpoints;

  PrintMessage(3,"Mesh Lines");

  /*
  cout << geom.GetNLines() << " lines" << endl;
  double totnp = 0;
  for (int i = 1; i <= geom.GetNLines(); i++)
    totnp += geom.GetLine(i)->NP();
  cout << "avg np per line " << totnp/geom.GetNLines() << endl;
  */

  for (int i = 1; i <= geom.GetNLines(); i++)
    {
      meshlines.Append(geom.GetLine(i)->Mesh(geom.GetPoints(), meshpoints, h, mesh)); 
      SetThreadPercent(100.0 * (double)i/(double)geom.GetNLines());
    }

  geom.meshpoints.SetSize(0); //testing
  geom.meshlines.SetSize(0);  //testing
  for (int i = 1; i <= meshpoints.Size(); i++)
    {
      geom.meshpoints.Append(meshpoints.Get(i)); //testing
      mesh.AddPoint(meshpoints.Get(i));
    }
  //(++++++++++++++testing
  for (int i = 1; i <= geom.GetNLines(); i++)
    {
      geom.meshlines.Append(meshlines.Get(i));
    }
  //++++++++++++++testing)

  PrintMessage(7,"feed with edges");

  for (int i = 1; i <= meshlines.Size(); i++)
    {
      STLLine* line = meshlines.Get(i);
      (*testout) << "store line " << i << endl;
      for (int j = 1; j <= line->GetNS(); j++)
	{
	  int p1, p2;
	  
	  line->GetSeg(j, p1, p2);
	  int trig1, trig2, trig1b, trig2b;

	  if (p1 == p2) 
	    cout << "Add Segment, p1 == p2 == " << p1 << endl;

	  // Test auf geschlossener Rand mit 2 Segmenten 
	      
	  if ((j == 2) && (line->GetNS() == 2))
	    {
	      int oldp1, oldp2;
	      line->GetSeg (1, oldp1, oldp2);
	      if (oldp1 == p2 && oldp2 == p1)
		{
		  PrintMessage(7,"MESSAGE: don't use second segment");
		  continue;
		}
	    }


	  //mesh point number
	  //p1 = geom2meshnum.Get(p1); // for unmeshed lines!!!
	  //p2 = geom2meshnum.Get(p2); // for unmeshed lines!!!
	  
	  //left and right trigs
	  trig1 = line->GetLeftTrig(j);
	  trig2 = line->GetRightTrig(j);
	  trig1b = line->GetLeftTrig(j+1);
	  trig2b = line->GetRightTrig(j+1);
	  
	  (*testout) << "j = " << j << ", p1 = " << p1 << ", p2 = " << p2 << endl;
	  (*testout) << "segm-trigs: "
		   << "trig1 = " << trig1
		   << ", trig1b = " << trig1b
		   << ", trig2 = " << trig2
		   << ", trig2b = " << trig2b << endl;

	  if (trig1 <= 0 || trig2 < 0 || trig1b <= 0 || trig2b < 0)
	    {
	      cout << "negative trigs, "
		   << ", trig1 = " << trig1
		   << ", trig1b = " << trig1b
		   << ", trig2 = " << trig2
		   << ", trig2b = " << trig2b << endl;
	    }
	  /*
	  (*testout) << "   trigs p1: " << trig1 << " - " << trig2 << endl;
	  (*testout) << "   trigs p2: " << trig1b << " - " << trig2b << endl;
	  (*testout) << "   charts p1: " << geom.GetChartNr(trig1) << " - " << geom.GetChartNr(trig2) << endl;
	  (*testout) << "   charts p2: " << geom.GetChartNr(trig1b) << " - " << geom.GetChartNr(trig2b) << endl;
	  */
	  Point3d hp, hp2;
	  Segment seg;
	  seg[0] = p1 + PointIndex::BASE-1;
	  seg[1] = p2 + PointIndex::BASE-1;
	  seg.si = geom.GetTriangle(trig1).GetFaceNum();
	  seg.edgenr = i;

	  seg.epgeominfo[0].edgenr = i;
	  seg.epgeominfo[0].dist = line->GetDist(j);
	  seg.epgeominfo[1].edgenr = i;
	  seg.epgeominfo[1].dist = line->GetDist(j+1);
	  /*
	  (*testout) << "seg = " 
		     << "edgenr " << seg.epgeominfo[0].edgenr
		     << " dist " << seg.epgeominfo[0].dist
		     << " edgenr " << seg.epgeominfo[1].edgenr
		     << " dist " << seg.epgeominfo[1].dist << endl;
	  */
	  
	  seg.geominfo[0].trignum = trig1;
	  seg.geominfo[1].trignum = trig1b;

	  /*
	  geom.SelectChartOfTriangle (trig1);
	  hp = hp2 = mesh.Point (seg[0]);
	  seg.geominfo[0].trignum = geom.Project (hp);

	  (*testout) << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[0].trignum << endl;
	  if (Dist (hp, hp2) > 1e-5 || seg.geominfo[0].trignum == 0) 
	    {
	      (*testout) << "PROBLEM" << endl;
	    }

	  geom.SelectChartOfTriangle (trig1b);
	  hp = hp2 = mesh.Point (seg[1]);
	  seg.geominfo[1].trignum = geom.Project (hp);

	  (*testout) << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[1].trignum << endl;
	  if (Dist (hp, hp2) > 1e-5 || seg.geominfo[1].trignum == 0) 
	    {
	      (*testout) << "PROBLEM" << endl;
	    }
	  */


	  if (Dist (mesh.Point(seg[0]), mesh.Point(seg[1])) < 1e-10)
	    {
	      (*testout) << "ERROR: Line segment of length 0" << endl;
	      (*testout) << "pi1, 2 = " << seg[0] << ", " << seg[1] << endl;
	      (*testout) << "p1, 2 = " << mesh.Point(seg[0])
			 << ", " << mesh.Point(seg[1]) << endl;
	      throw NgException ("Line segment of length 0");
	    }
	  
	  mesh.AddSegment (seg);


          if(trig2 != 0)
            {
	  Segment seg2;
	  seg2[0] = p2 + PointIndex::BASE-1;;
	  seg2[1] = p1 + PointIndex::BASE-1;;
	  seg2.si = geom.GetTriangle(trig2).GetFaceNum();

	  seg2.edgenr = i;

	  seg2.epgeominfo[0].edgenr = i;
	  seg2.epgeominfo[0].dist = line->GetDist(j+1);
	  seg2.epgeominfo[1].edgenr = i;
	  seg2.epgeominfo[1].dist = line->GetDist(j);
	  /*
	  (*testout) << "seg = " 
		     << "edgenr " << seg2.epgeominfo[0].edgenr
		     << " dist " << seg2.epgeominfo[0].dist
		     << " edgenr " << seg2.epgeominfo[1].edgenr
		     << " dist " << seg2.epgeominfo[1].dist << endl;
	  */
	  
	  seg2.geominfo[0].trignum = trig2b;
	  seg2.geominfo[1].trignum = trig2;
	  
	  /*
	  geom.SelectChartOfTriangle (trig2);
	  hp = hp2 = mesh.Point (seg[0]);
	  seg2.geominfo[0].trignum = geom.Project (hp);

	  (*testout) << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[0].trignum << endl;
	  if (Dist (hp, hp2) > 1e-5 || seg2.geominfo[0].trignum == 0) 
	    {
	      (*testout) << "Get GeomInfo PROBLEM" << endl;
	    }


	  geom.SelectChartOfTriangle (trig2b);
	  hp = hp2 = mesh.Point (seg[1]);
	  seg2.geominfo[1].trignum = geom.Project (hp);
	  (*testout) << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[1].trignum << endl;
	  if (Dist (hp, hp2) > 1e-5 || seg2.geominfo[1].trignum == 0) 
	    {
	      (*testout) << "Get GeomInfo PROBLEM" << endl;
	    }
	  */	  
	  mesh.AddSegment (seg2);
            }
	}
    }

  PopStatus();
}




void STLSurfaceMeshing1 (STLGeometry & geom, class Mesh & mesh, const MeshingParameters& mparam,
			 int retrynr, const STLParameters& stlparam);


int STLSurfaceMeshing (STLGeometry & geom, class Mesh & mesh, const MeshingParameters& mparam,
                       const STLParameters& stlparam)
{
  PrintFnStart("Do Surface Meshing");

  geom.PrepareSurfaceMeshing();

  if (mesh.GetNSeg() == 0)
    STLFindEdges (geom, mesh, mparam, stlparam);

  int nopen;
  int outercnt = 20;

  for (int i = 1; i <= mesh.GetNSeg(); i++)
    {
      const Segment & seg = mesh.LineSegment (i);
      if (seg.geominfo[0].trignum <= 0 || seg.geominfo[1].trignum <= 0)
	(*testout) << "Problem with segment " << i << ": " << seg << endl;
    }


  do
    {
      outercnt--;
      if (outercnt <= 0)
	return MESHING3_OUTERSTEPSEXCEEDED;
      
      if (multithread.terminate) return MESHING3_TERMINATE;

      mesh.FindOpenSegments();
      nopen = mesh.GetNOpenSegments();

      if (nopen)
	{
	  int trialcnt = 0;
	  while (nopen && trialcnt <= 5)
	    {
	      if (multithread.terminate) { return MESHING3_TERMINATE; }

	      trialcnt++;
	      STLSurfaceMeshing1 (geom, mesh, mparam, trialcnt, stlparam);

	      mesh.FindOpenSegments();
	      nopen = mesh.GetNOpenSegments();

              auto n_illegal_trigs = mesh.FindIllegalTrigs();
              PrintMessage (3, n_illegal_trigs, " illegal triangles");

	      if (nopen)
		{
		  geom.ClearMarkedSegs();
		  for (int i = 1; i <= nopen; i++)
		    {
		      const Segment & seg = mesh.GetOpenSegment (i);
		      geom.AddMarkedSeg(mesh.Point(seg[0]),mesh.Point(seg[1]));
		    }

		  geom.InitMarkedTrigs();
		  for (int i = 1; i <= nopen; i++)
		    {
		      const Segment & seg = mesh.GetOpenSegment (i);
		      geom.SetMarkedTrig(seg.geominfo[0].trignum,1);
		      geom.SetMarkedTrig(seg.geominfo[1].trignum,1);
		    }

		  MeshOptimize2d optmesh(mesh);
		  optmesh.SetFaceIndex (0);
		  optmesh.SetImproveEdges (0);
		  optmesh.SetMetricWeight (0);
		  
		  mesh.CalcSurfacesOfNode();
		  optmesh.EdgeSwapping (0);
		  optmesh.ImproveMesh (mparam);
		}

	      mesh.Compress();
	      mesh.FindOpenSegments();
	      nopen = mesh.GetNOpenSegments();

	      if (trialcnt <= 5 && nopen)
		{
		  mesh.RemoveOneLayerSurfaceElements();

		  if (trialcnt >= 4)
		    {
		      mesh.FindOpenSegments();
		      mesh.RemoveOneLayerSurfaceElements();

		      mesh.FindOpenSegments ();		  
		      nopen = mesh.GetNOpenSegments();
		    }
		}
	    }


	  if (multithread.terminate)
	    return MESHING3_TERMINATE;

	  if (nopen)
	    {
	      
	      PrintMessage(3,"Meshing failed, trying to refine");

	      mesh.FindOpenSegments ();
	      nopen = mesh.GetNOpenSegments();
			  
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();

	      // Open edge-segments will be refined !
	      INDEX_2_HASHTABLE<int> openseght (nopen+1);
	      for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
		{
		  const Segment & seg = mesh.GetOpenSegment (i);
		  INDEX_2 i2(seg[0], seg[1]);
		  i2.Sort();
		  openseght.Set (i2, 1);
		}

	      
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();
	      

	      INDEX_2_HASHTABLE<int> newpht(100);

	      int nsegold = mesh.GetNSeg();
	      for (int i = 1; i <= nsegold; i++)
		{
		  Segment seg = mesh.LineSegment(i);
		  INDEX_2 i2(seg[0], seg[1]);
		  i2.Sort();
		  if (openseght.Used (i2))
		    {
		      // segment will be split
		      PrintMessage(7,"Split segment ", int(seg[0]), "-", int(seg[1]));
	      
		      Segment nseg1, nseg2;
		      EdgePointGeomInfo newgi;
		      
		      const EdgePointGeomInfo & gi1 = seg.epgeominfo[0];
		      const EdgePointGeomInfo & gi2 = seg.epgeominfo[1];
		      
		      newgi.dist = 0.5 * (gi1.dist + gi2.dist);
		      newgi.edgenr = gi1.edgenr;

		      int hi;
		      
		      Point3d newp;
		      int newpi;
		      
		      if (!newpht.Used (i2))
			{
			  newp = geom.GetLine (gi1.edgenr)->
			    GetPointInDist (geom.GetPoints(), newgi.dist, hi);
			  newpi = mesh.AddPoint (newp);
			  newpht.Set (i2, newpi);
			}
		      else
			{
			  newpi = newpht.Get (i2);
			  newp = mesh.Point (newpi);
			}

		      nseg1 = seg;
		      nseg2 = seg;
		      nseg1[1] = newpi;
		      nseg1.epgeominfo[1] = newgi;
		      
		      nseg2[0] = newpi;
		      nseg2.epgeominfo[0] = newgi;
		      
		      mesh.LineSegment(i) = nseg1;
		      mesh.AddSegment (nseg2);
		      
		      mesh.RestrictLocalH (Center (mesh.Point(nseg1[0]),
						   mesh.Point(nseg1[1])),
					   Dist (mesh.Point(nseg1[0]),
						 mesh.Point(nseg1[1])));
		      mesh.RestrictLocalH (Center (mesh.Point(nseg2[0]),
						   mesh.Point(nseg2[1])),
					   Dist (mesh.Point(nseg2[0]),
						 mesh.Point(nseg2[1])));
		    }
		}

	    }

	  nopen = -1;
	}
    
      else

	{
	  PrintMessage(5,"mesh is closed, verifying ...");

	  // no open elements, check wrong elements (intersecting..)



	  PrintMessage(5,"check overlapping");
	  // 	  mesh.FindOpenElements(); // would leed to locked points
          mesh.CheckOverlappingBoundary();
	  // if(mesh.CheckOverlappingBoundary()) ;
          // return MESHING3_BADSURFACEMESH;

	  geom.InitMarkedTrigs();

	  for (int i = 1; i <= mesh.GetNSE(); i++)
	    if (mesh.SurfaceElement(i).BadElement())
	      {
		int trig = mesh.SurfaceElement(i).PNum(1);
		geom.SetMarkedTrig(trig,1);
		PrintMessage(7, "overlapping element, will be removed");
	      }
	  
	  

	  NgArray<Point3d> refpts;
	  NgArray<double> refh;

	  // was commented:

	  for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	    if (mesh[sei].BadElement())
	      {
		for (int j = 1; j <= 3; j++)
		  {
		    refpts.Append (mesh.Point (mesh[sei].PNum(j)));
		    refh.Append (mesh.GetH (refpts.Last()) / 2);
		  }
		mesh.Delete(sei);
	      }
	  	  
	  // delete wrong oriented element
	  for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	    {
	      const Element2d & el = mesh[sei];
              if (el.IsDeleted()) continue;
	      if (!el.PNum(1).IsValid()) continue;

	      Vec3d n = Cross (Vec3d (mesh.Point(el.PNum(1)), 
				      mesh.Point(el.PNum(2))),
			       Vec3d (mesh.Point(el.PNum(1)), 
				      mesh.Point(el.PNum(3))));
	      Vec3d ng = geom.GetTriangle(el.GeomInfoPi(1).trignum).Normal();
	      if (n * ng < 0)
		{
		  refpts.Append (mesh.Point (mesh[sei].PNum(1)));
		  refh.Append (mesh.GetH (refpts.Last()) / 2);
		  mesh.Delete(sei);
		}
	    }
	  // end comments

	  for (int i = 1; i <= refpts.Size(); i++)
	    mesh.RestrictLocalH (refpts.Get(i), refh.Get(i));

	  mesh.RemoveOneLayerSurfaceElements();
          // Open edge-segments will be refined !
	      INDEX_2_HASHTABLE<int> openseght (nopen+1);
	      for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
		{
		  const Segment & seg = mesh.GetOpenSegment (i);
		  INDEX_2 i2(seg[0], seg[1]);
		  i2.Sort();
		  openseght.Set (i2, 1);
		}
          mesh.FindOpenSegments ();
          mesh.RemoveOneLayerSurfaceElements();
          mesh.FindOpenSegments ();
          int nsegold = mesh.GetNSeg();
          INDEX_2_HASHTABLE<int> newpht(100);
          for (int i = 1; i <= nsegold; i++)
            {
              Segment seg = mesh.LineSegment(i);
              INDEX_2 i2(seg[0], seg[1]);
              i2.Sort();
              if (openseght.Used (i2))
                {
                  // segment will be split
                  PrintMessage(7,"Split segment ", int(seg[0]), "-", int(seg[1]));
	      
                  Segment nseg1, nseg2;
                  EdgePointGeomInfo newgi;
		      
                  const EdgePointGeomInfo & gi1 = seg.epgeominfo[0];
                  const EdgePointGeomInfo & gi2 = seg.epgeominfo[1];
		      
                  newgi.dist = 0.5 * (gi1.dist + gi2.dist);
                  newgi.edgenr = gi1.edgenr;

                  int hi;
		      
                  Point3d newp;
                  int newpi;
		      
                  if (!newpht.Used (i2))
                    {
                      newp = geom.GetLine (gi1.edgenr)->
                        GetPointInDist (geom.GetPoints(), newgi.dist, hi);
                      newpi = mesh.AddPoint (newp);
                      newpht.Set (i2, newpi);
                    }
                  else
                    {
                      newpi = newpht.Get (i2);
                      newp = mesh.Point (newpi);
                    }

                  nseg1 = seg;
                  nseg2 = seg;
                  nseg1[1] = newpi;
                  nseg1.epgeominfo[1] = newgi;
		      
                  nseg2[0] = newpi;
                  nseg2.epgeominfo[0] = newgi;
		      
                  mesh.LineSegment(i) = nseg1;
                  mesh.AddSegment (nseg2);
		      
                  mesh.RestrictLocalH (Center (mesh.Point(nseg1[0]),
                                               mesh.Point(nseg1[1])),
                                       Dist (mesh.Point(nseg1[0]),
                                             mesh.Point(nseg1[1])));
                  mesh.RestrictLocalH (Center (mesh.Point(nseg2[0]),
                                               mesh.Point(nseg2[1])),
                                       Dist (mesh.Point(nseg2[0]),
                                             mesh.Point(nseg2[1])));
                }
            }
	  mesh.Compress();
	  
	  mesh.FindOpenSegments ();
	  nopen = mesh.GetNOpenSegments();

	  /*
	  if (!nopen)
	    {
	      // mesh is still ok

	      void STLSurfaceOptimization (STLGeometry & geom,
					   class Mesh & mesh,
					   MeshingParameters & mparam)
	      
	    }
	  */
	}
      
    }
  while (nopen);

  if(mesh.CheckOverlappingBoundary())
    return MESHING3_BADSURFACEMESH;

  mesh.Compress();
  mesh.CalcSurfacesOfNode();

  return MESHING3_OK;
}






void STLSurfaceMeshing1 (STLGeometry & geom,
			 Mesh & mesh,
                         const MeshingParameters& mparam,
			 int retrynr,
                         const STLParameters& stlparam)
{
  static int timer1 = NgProfiler::CreateTimer ("STL surface meshing1");
  static int timer1a = NgProfiler::CreateTimer ("STL surface meshing1a");
  static int timer1b = NgProfiler::CreateTimer ("STL surface meshing1b");
  static int timer1c = NgProfiler::CreateTimer ("STL surface meshing1c");
  static int timer1d = NgProfiler::CreateTimer ("STL surface meshing1d");

  double h = mparam.maxh;

  mesh.FindOpenSegments();
  
  NgArray<int> spiralps(0);
  spiralps.SetSize(0);
  for (int i = 1; i <= geom.GetNP(); i++)
    if (geom.GetSpiralPoint(i)) 
      spiralps.Append(i);
  
  PrintMessage(7,"NO spiralpoints = ", spiralps.Size());
  //int spfound;

  /*
  NgArray<int> meshsp(mesh.GetNP());
  meshsp = 0;
  for (int i = 1; i <= mesh.GetNP(); i++)
    for (int j = 1; j <= spiralps.Size(); j++)
      if (Dist2(geom.GetPoint(spiralps.Get(j)), mesh.Point(i)) < 1e-20) 
	meshsp.Elem(i) = spiralps.Get(j);
  NgArray<PointIndex> imeshsp;
  for (int i = 1; i <= meshsp.Size(); i++)
    if (meshsp.Elem(i)) imeshsp.Append(i);
  */
  NgArray<PointIndex> imeshsp;
  NgArray<int> ispiral_point;
  for (int i = 1; i <= mesh.GetNP(); i++)
    {
      for (int j = 1; j <= spiralps.Size(); j++)
	if (Dist2(geom.GetPoint(spiralps.Get(j)), mesh.Point(i)) < 1e-20) 
	  {
	    imeshsp.Append(i);
	    ispiral_point.Append(spiralps.Get(j));
	    break;
	  }
    }

  double starttime = GetTime ();
  mesh.SurfaceArea().ReCalc();

  // int oldnp = mesh.GetNP();

  NgArray<int,PointIndex::BASE> compress(mesh.GetNP());
  compress = 0;
  NgArray<PointIndex> icompress; 

  NgArray<int, 1> opensegsperface(mesh.GetNFD());
  opensegsperface = 0;
  for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
    opensegsperface[mesh.GetOpenSegment(i).si]++;
  
  TABLE<int, 1> opensegments(mesh.GetNFD());
  for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
    {
      const Segment & seg = mesh.GetOpenSegment (i);
      if (seg.si < 1 || seg.si > mesh.GetNFD())
	cerr << "segment index " << seg.si << " out of range [1, " << mesh.GetNFD() << "]" << endl;
      opensegments.Add (seg.si, i);
    }
  

  for (int fnr = 1; fnr <= mesh.GetNFD(); fnr++)
    {
      if (!opensegsperface[fnr]) continue;
      if (multithread.terminate) return;

      NgProfiler::StartTimer (timer1);
      NgProfiler::StartTimer (timer1a);


      PrintMessage(5,"Meshing surface ", fnr, "/", mesh.GetNFD());
      MeshingSTLSurface meshing (geom, mparam);
      meshing.SetStartTime (starttime);
      
      // compress = 0;
      icompress.SetSize(0); 
      int cntused = 0;

      for (int i = 0; i < imeshsp.Size(); i++)
	{
	  compress[imeshsp[i]] = ++cntused;
	  icompress.Append(imeshsp[i]);
	}

      NgProfiler::StopTimer (timer1a);
      NgProfiler::StartTimer (timer1b);



      /*
      for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
	{
	  const Segment & seg = mesh.GetOpenSegment (i);
	  if (seg.si == fnr)
	    for (int j = 0; j < 2; j++)
	      if (compress[seg[j]] == 0)
		{
		  compress[seg[j]] = ++cntused;
		  icompress.Append(seg[j]);
		}
	}
      */
      NgFlatArray<int> segs = opensegments[fnr];
      for (int hi = 0; hi < segs.Size(); hi++)
	{
	  int i = segs[hi];
	  const Segment & seg = mesh.GetOpenSegment (i);
	  for (int j = 0; j < 2; j++)
	    if (compress[seg[j]] == 0)
	      {
		compress[seg[j]] = ++cntused;
		icompress.Append(seg[j]);
	      }
	}

      NgProfiler::StopTimer (timer1b);
      NgProfiler::StartTimer (timer1c);


      for (int hi = 0; hi < icompress.Size(); hi++)
	{
	  PointIndex pi = icompress[hi];
	  
	  /*
	  // int sppointnum = meshsp.Get(i);
	  int sppointnum = 0;
	  if (hi < ispiral_point.Size())
	    sppointnum = ispiral_point[hi];

	  if (sppointnum)
	    {
	  */
	  if (hi < ispiral_point.Size())
	    {
	      int sppointnum = ispiral_point[hi];

	      MultiPointGeomInfo mgi;
	      
	      int ntrigs = geom.NOTrigsPerPoint(sppointnum);
	      for (int j = 0; j < ntrigs; j++)
		{
		  PointGeomInfo gi;
		  gi.trignum = geom.TrigPerPoint(sppointnum, j+1);
		  mgi.AddPointGeomInfo (gi);
		}
	      
	      // Einfuegen von ConePoint: Point bekommt alle
	      // Dreiecke (werden dann intern kopiert)
	      // Ein Segment zum ConePoint muss vorhanden sein !!!
	      
	      // meshing.AddPoint (mesh.Point(i), i, &mgi);
	      meshing.AddPoint (mesh[pi], pi, &mgi);
	    }
	  else
	    meshing.AddPoint (mesh[pi], pi);
	}

      NgProfiler::StopTimer (timer1c);
      NgProfiler::StartTimer (timer1d);

      /*
        for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
	  {
	    const Segment & seg = mesh.GetOpenSegment (i);
	    if (seg.si == fnr)
	      meshing.AddBoundaryElement (compress[seg[0]], compress[seg[1]], 
					  seg.geominfo[0], seg.geominfo[1]);
	  }
      */


      // NgFlatArray<int> segs = opensegments[fnr];
      for (int hi = 0; hi < segs.Size(); hi++)
	{
	  int i = segs[hi];
	  const Segment & seg = mesh.GetOpenSegment (i);
	  meshing.AddBoundaryElement (compress[seg[0]], compress[seg[1]], 
				      seg.geominfo[0], seg.geominfo[1]);
	}



      NgProfiler::StopTimer (timer1d);

      NgProfiler::StopTimer (timer1);
      
      PrintMessage(3,"start meshing, trialcnt = ", retrynr);
      
      meshing.GenerateMesh (mesh, mparam, h, fnr);  
      
      for (int i = 0; i < icompress.Size(); i++)
	compress[icompress[i]] = 0;
      
      
      mparam.Render();
    }     
  
  // NgProfiler::Print(stdout);
  
  mesh.CalcSurfacesOfNode();
}



void STLSurfaceOptimization (STLGeometry & geom,
			     Mesh & mesh,
			     const MeshingParameters & mparam)
{
  PrintFnStart("optimize STL Surface");

  MeshOptimize2d optmesh(mesh);

  optmesh.SetFaceIndex (0);
  optmesh.SetImproveEdges (0);
  optmesh.SetMetricWeight (mparam.elsizeweight);

  PrintMessage(5,"optimize string = ", mparam.optimize2d, " elsizew = ", mparam.elsizeweight);

  for (int i = 1; i <= mparam.optsteps2d; i++)
    for (size_t j = 1; j <= mparam.optimize2d.length(); j++)
      {
	if (multithread.terminate)
	  break;

	//(*testout) << "optimize, before, step = " << meshparam.optimize2d[j-1] << mesh.Point (3679) << endl;

	mesh.CalcSurfacesOfNode();
	switch (mparam.optimize2d[j-1])
	  {
	  case 's': 
	    {
	      optmesh.EdgeSwapping(0);
	      break;
	    }
	  case 'S': 
	    {
	      optmesh.EdgeSwapping(1);
	      break;
	    }
	  case 'm': 
	    {
	      optmesh.ImproveMesh(mparam);
	      break;
	    }
	  case 'c': 
	    {
	      optmesh.CombineImprove();
	      break;
	    }
	  }
        // while(mesh.CheckOverlappingBoundary())
        //   {
        //     for(const auto & el : mesh.SurfaceElements())
        //       {
        //         if(el.BadElement())
        //           {
        //             cout << "Restrict localh at el nr " << el << endl;
        //             for(const auto& p : el.PNums())
        //               {
        //                 const auto& pnt = mesh[p];
        //                 mesh.RestrictLocalH(pnt, 0.5*mesh.GetH(pnt));
        //               }
        //           }
        //       }
        //     optmesh.SplitImprove();
        //   }
	//(*testout) << "optimize, after, step = " << meshparam.optimize2d[j-1] << mesh.Point (3679) << endl;
      }

  geom.surfaceoptimized = 1;

  mesh.Compress();
  mesh.CalcSurfacesOfNode();


}



MeshingSTLSurface :: MeshingSTLSurface (STLGeometry & ageom,
					const MeshingParameters & mp)
  : Meshing2(ageom, mp, ageom.GetBoundingBox()), geom(ageom)
{
  ;
}

void MeshingSTLSurface :: DefineTransformation (const Point<3> & p1, const Point<3> & p2,
						const PointGeomInfo * geominfo,
						const PointGeomInfo * geominfo2)
{
  transformationtrig = geominfo[0].trignum;
  
  geom.DefineTangentialPlane(p1, p2, transformationtrig);
}

void MeshingSTLSurface :: TransformToPlain (const Point<3> & locpoint, const MultiPointGeomInfo & gi,
					    Point<2> & plainpoint, double h, int & zone)
{
  int trigs[10000];

  if (gi.GetNPGI() >= 9999) 
    {
      PrintError("In Transform to plane: increase size of trigs!!!");
    }

  for (int i = 1; i <= gi.GetNPGI(); i++)
    trigs[i-1] = gi.GetPGI(i).trignum;
  trigs[gi.GetNPGI()] = 0;

  //  int trig = gi.trignum;
  //   (*testout) << "locpoint = " << locpoint;

  geom.ToPlane (locpoint, trigs, plainpoint, h, zone, 1);

  //  geom.ToPlane (locpoint, NULL, plainpoint, h, zone, 1);
  /*
  (*testout) << " plainpoint = " << plainpoint
	     << " h = " << h 
	     << endl;
  */
}

/*
int MeshingSTLSurface :: ComputeLineGeoInfo (const Point3d & p1, const Point3d & p2,
					      int & geoinfosize, void *& geoinfo)
{
  static int geomtrig[2] = { 0, 0 };

  Point3d hp;
  hp = p1;
  geomtrig[0] = geom.Project (hp);

  hp = p2;
  geomtrig[1] = geom.Project (hp);
  
  geoinfosize = sizeof (geomtrig);
  geoinfo = &geomtrig;

  if (geomtrig[0] == 0)
    {
      return 1;
    }
  return 0;
}
*/


int MeshingSTLSurface :: ComputePointGeomInfo (const Point3d & p, PointGeomInfo & gi)
{
  // compute triangle of point,
  // if non-unique: 0

  Point<3> hp = p;
  gi.trignum = geom.Project (hp);

  if (!gi.trignum)
    {
      return 1;
    }

  return 0;
}


int MeshingSTLSurface :: 
ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
			  PointGeomInfo & pgi)
{
  for (int i = 1; i <= mpgi.GetNPGI(); i++)
    if (geom.TrigIsInOC (mpgi.GetPGI(i).trignum, geom.meshchart))
      {
	pgi = mpgi.GetPGI(i);
	return 0;
      }
  /*
  for (i = 0; i < mpgi.cnt; i++)
    {
      //      (*testout) << "d" << endl;
      if (geom.TrigIsInOC (mpgi.mgi[i].trignum, geom.meshchart))
	{
	  pgi = mpgi.mgi[i];
	  return 0;
	}
    }
  */
  PrintMessage(7,"INFORM: no gi on chart");
  pgi.trignum = 1;
  return 1;
}



int MeshingSTLSurface :: 
IsLineVertexOnChart (const Point3d & p1, const Point3d & p2,
		     int endpoint, const PointGeomInfo & gi)
{
  int lineendtrig = gi.trignum;
  return geom.TrigIsInOC (lineendtrig, geom.meshchart);

  // Vec3d baselinenormal = geom.meshtrignv;
  //  Vec3d linenormal = geom.GetTriangleNormal (lineendtrig);
  //  return ( (baselinenormal * linenormal) > cos (30 * (M_PI/180)) );
}

void MeshingSTLSurface :: 
GetChartBoundary (NgArray<Point<2>> & points, 
		  NgArray<Point<3>> & points3d,
		  NgArray<INDEX_2> & lines, double h) const
{
  points.SetSize (0);
  points3d.SetSize (0);
  lines.SetSize (0);
  geom.GetMeshChartBoundary (points, points3d, lines, h);
}




int MeshingSTLSurface :: TransformFromPlain (const Point<2> & plainpoint,
					     Point<3> & locpoint, 
					     PointGeomInfo & gi, 
					     double h)
{
  //return 0, wenn alles OK
  Point<3> hp3d;
  int res = geom.FromPlane (plainpoint, hp3d, h);
  locpoint = hp3d;
  ComputePointGeomInfo (locpoint, gi);
  return res;
}


int MeshingSTLSurface :: 
BelongsToActiveChart (const Point3d & p, 
		      const PointGeomInfo & gi)
{
  return (geom.TrigIsInOC(gi.trignum, geom.meshchart) != 0);
}



double MeshingSTLSurface :: CalcLocalH (const Point<3> & p, double gh) const
{
  return gh;
}

double MeshingSTLSurface :: Area () const
{
  return geom.Area();
}

}

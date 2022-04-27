#include <mystdlib.h>
#include "meshing.hpp"
#include "visual_interface.hpp"

namespace netgen
{
  static void glrender (int wait)
  {
    //  cout << "plot adfront" << endl;

    if (multithread.drawing)
      {
        // vssurfacemeshing.DrawScene();
	Render ();

	if (wait || multithread.testmode)
	  {
	    multithread.pause = 1;
	  }
	while (multithread.pause);
      }
  }

  ostream& operator << (ostream& ost, const MultiPointGeomInfo& mpgi)
  {
    for(auto i : Range(mpgi.GetNPGI()))
      {
        ost << "gi[" << i << "] = " << mpgi.GetPGI(i+1) << endl;
      }
    return ost;
  }

  static Array<unique_ptr<netrule>> global_trig_rules;
  static Array<unique_ptr<netrule>> global_quad_rules;

  
  Meshing2 :: Meshing2 (const NetgenGeometry& ageo,
                        const MeshingParameters & mp,
                        const Box<3> & aboundingbox)
    : adfront(aboundingbox), boundingbox(aboundingbox), geo(ageo)
  {
    static Timer t("Mesing2::Meshing2"); RegionTimer r(t);

    auto & globalrules = mp.quad ? global_quad_rules : global_trig_rules;

    {
      static mutex mut;
      lock_guard<mutex> guard(mut);
      if (!globalrules.Size())
        {
          LoadRules (NULL, mp.quad);
          for (auto & rule : rules)
            globalrules.Append (make_unique<netrule>(*rule));
          rules.SetSize(0);
        }
      /*
      else
        {
          for (auto i : globalrules.Range())
            rules.Append (globalrules[i].get());
        }
      */
    }
    for (auto i : globalrules.Range())
      rules.Append (make_unique<netrule>(*globalrules[i]));

    // LoadRules ("rules/quad.rls");
    // LoadRules ("rules/triangle.rls");


    
    // adfront = new AdFront2(boundingbox);
    starttime = GetTime();

    maxarea = -1;
  }


  Meshing2 :: ~Meshing2 ()
  { ; } 

  void Meshing2 :: AddPoint (const Point3d & p, PointIndex globind, 
			     MultiPointGeomInfo * mgi,
			     bool pointonsurface)
  {
    //(*testout) << "add point " << globind << endl;
    adfront.AddPoint (p, globind, mgi, pointonsurface);
  }

  void Meshing2 :: AddBoundaryElement (int i1, int i2,
				       const PointGeomInfo & gi1, const PointGeomInfo & gi2)
  {
    //    (*testout) << "add line " << i1 << " - " << i2 << endl;
    if (!gi1.trignum || !gi2.trignum)
      {
	PrintSysError ("addboundaryelement: illegal geominfo");
      }
    adfront.AddLine (i1-1, i2-1, gi1, gi2);
  }



  void Meshing2 :: StartMesh ()
  {
    foundmap.SetSize (rules.Size());
    canuse.SetSize (rules.Size());
    ruleused.SetSize (rules.Size());

    foundmap = 0;
    canuse = 0;
    ruleused = 0;

    // cntelem = 0;
    // trials = 0;
  }

  void Meshing2 :: EndMesh ()
  {
    for (int i = 0; i < ruleused.Size(); i++)
      (*testout) << setw(4) << ruleused[i]
		 << " times used rule " << rules[i] -> Name() << endl;
  }

  void Meshing2 :: SetStartTime (double astarttime)
  {
    starttime = astarttime;
  }

  
  void Meshing2 :: SetMaxArea (double amaxarea)
  {
    maxarea = amaxarea;
  }


  double Meshing2 :: CalcLocalH (const Point<3> & /* p */, double gh) const
  {
    return gh;
  }

  // should be class variables !!(?)
  // static Vec3d ex, ey;
  // static Point3d globp1;

  void Meshing2 :: DefineTransformation (const Point<3> & ap1,
                                         const Point<3> & ap2,
					 const PointGeomInfo * gi1,
					 const PointGeomInfo * gi2)
  {
    p1 = ap1;
    p2 = ap2;
    auto n1 = geo.GetNormal(gi1->trignum, p1, gi1);
    auto n2 = geo.GetNormal(gi2->trignum, p2, gi2);

    ez = 0.5 * (n1+n2);
    ez.Normalize();
    ex = (p2-p1).Normalize();
    ez -= (ez*ex)*ex;
    ez.Normalize();
    ey = Cross(ez, ex);
  }

  void Meshing2 :: TransformToPlain (const Point<3> & locpoint, 
				     const MultiPointGeomInfo & geominfo,
				     Point<2> & plainpoint, double h, int & zone)
  {
    auto& gi = geominfo.GetPGI(1);
    auto n = geo.GetNormal(gi.trignum, locpoint, &gi);
    auto p1p = locpoint - p1;
    plainpoint(0) = (p1p * ex) / h;
    plainpoint(1) = (p1p * ey) / h;

    if(n*ez < 0)
      zone = -1;
    else
      zone = 0;
  }

  int Meshing2 :: TransformFromPlain (const Point<2> & plainpoint,
				      Point<3> & locpoint, 
				      PointGeomInfo & gi, 
				      double h)
  {
    locpoint = p1 + (h*plainpoint(0)) * ex + (h* plainpoint(1)) * ey;
    if (!geo.ProjectPointGI(gi.trignum, locpoint, gi))
      gi = geo.ProjectPoint(gi.trignum, locpoint);
    return 0;
  }


  int Meshing2 :: BelongsToActiveChart (const Point3d & p, 
					const PointGeomInfo & gi)
  {
    return 1;
  }


  int Meshing2 :: ComputePointGeomInfo (const Point3d & p, PointGeomInfo & gi)
  {
    gi.trignum = 1;
    return 0;
  }


  int Meshing2 :: ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
					    PointGeomInfo & pgi)
  {
    pgi = mpgi.GetPGI(1);
    return 0;
  }



  int Meshing2 :: 
  IsLineVertexOnChart (const Point3d & p1, const Point3d & p2,
		       int endpoint, const PointGeomInfo & geominfo)
  {
    return 1;
  }

  void Meshing2 ::
  GetChartBoundary (NgArray<Point<2>> & points, 
		    NgArray<Point<3>> & points3d, 
		    NgArray<INDEX_2> & lines, double h) const
  {
    points.SetSize (0);
    points3d.SetSize (0);
    lines.SetSize (0);
  }

  double Meshing2 :: Area () const
  {
    return -1;
  }





  MESHING2_RESULT Meshing2 :: GenerateMesh (Mesh & mesh, const MeshingParameters & mp, double gh, int facenr, int layer)
  {
    static Timer timer("surface meshing"); RegionTimer reg(timer);

    static int timer1 = NgProfiler::CreateTimer ("surface meshing1");
    static int timer2 = NgProfiler::CreateTimer ("surface meshing2");
    static int timer3 = NgProfiler::CreateTimer ("surface meshing3");

    static int ts1 = NgProfiler::CreateTimer ("surface meshing start 1");
    static int ts2 = NgProfiler::CreateTimer ("surface meshing start 2");
    static int ts3 = NgProfiler::CreateTimer ("surface meshing start 3");


    NgProfiler::StartTimer (ts1);

    NgArray<int> pindex, lindex;
    NgArray<int> delpoints, dellines;

    NgArray<PointGeomInfo> upgeominfo;  // unique info
    NgArray<MultiPointGeomInfo> mpgeominfo;  // multiple info

    NgArray<Element2d> locelements;

    int z1, z2, oldnp(-1);
    bool found;
    int rulenr(-1);

    const PointGeomInfo * blgeominfo1;
    const PointGeomInfo * blgeominfo2;

    bool morerisc;
    bool debugflag;

    // double h;

    auto locpointsptr = make_shared<NgArray<Point<3>>>();
    auto& locpoints = *locpointsptr;
    NgArray<int> legalpoints;
    auto plainpointsptr = make_shared<NgArray<Point<2>>>();
    auto& plainpoints = *plainpointsptr;
    NgArray<int> plainzones;
    auto loclinesptr = make_shared<NgArray<INDEX_2>>();
    auto &loclines = *loclinesptr;
    int cntelem = 0, trials = 0, nfaces = 0;
    int oldnl = 0;

    UpdateVisSurfaceMeshData(oldnl, locpointsptr, loclinesptr, plainpointsptr);

    int qualclass;



    // test for 3d overlaps
    BoxTree<3> surfeltree (boundingbox.PMin(),
                           boundingbox.PMax());

    NgArray<int> intersecttrias;
    NgArray<Point3d> critpoints;

    // test for doubled edges
    //INDEX_2_HASHTABLE<int> doubleedge(300000);


    testmode = 0;

    StartMesh();

    NgArray<Point<2>> chartboundpoints;
    NgArray<Point<3>> chartboundpoints3d;
    NgArray<INDEX_2> chartboundlines;

    // illegal points: points with more then 50 elements per node
    int maxlegalpoint(-1), maxlegalline(-1);
    NgArray<int,PointIndex::BASE> trigsonnode;
    NgArray<int,PointIndex::BASE> illegalpoint;

    trigsonnode.SetSize (mesh.GetNP());
    illegalpoint.SetSize (mesh.GetNP());

    trigsonnode = 0;
    illegalpoint = 0;
  
    double totalarea = Area ();
    double meshedarea = 0;


    // search tree for surface elements:
    /*
    for (sei = 0; sei < mesh.GetNSE(); sei++)
      {
	const Element2d & sel = mesh[sei];

	if (sel.IsDeleted()) continue;

	if (sel.GetIndex() == facenr)
	  {
	    Box<3> box;
	    box.Set ( mesh[sel[0]] );
	    box.Add ( mesh[sel[1]] );
	    box.Add ( mesh[sel[2]] );
	    surfeltree.Insert (box, sei);
	  }
      }
    */
    Array<SurfaceElementIndex> seia;
    mesh.GetSurfaceElementsOfFace (facenr, seia);
    for (int i = 0; i < seia.Size(); i++)
      {
	const Element2d & sel = mesh[seia[i]];

	if (sel.IsDeleted()) continue;

	Box<3> box;
	box.Set ( mesh[sel[0]] );
	box.Add ( mesh[sel[1]] );
	box.Add ( mesh[sel[2]] );
	surfeltree.Insert (box, seia[i]);
      }

    NgProfiler::StopTimer (ts1);
    NgProfiler::StartTimer (ts2);

    if (totalarea > 0 || maxarea > 0)
      meshedarea = mesh.SurfaceArea();
      /*
      for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	{
	  const Element2d & sel = mesh[sei];
	  if (sel.IsDeleted()) continue;
	
	  double trigarea = Cross ( mesh[sel[1]]-mesh[sel[0]],
				    mesh[sel[2]]-mesh[sel[0]] ).Length() / 2;
	  
	  
	  if (sel.GetNP() == 4)
	    trigarea += Cross (Vec3d (mesh.Point (sel.PNum(1)),
				      mesh.Point (sel.PNum(3))),
			       Vec3d (mesh.Point (sel.PNum(1)),
				      mesh.Point (sel.PNum(4)))).Length() / 2;;
	  meshedarea += trigarea;
	}
      */
      // cout << "meshedarea = " << meshedarea << " =?= " 
      // << mesh.SurfaceArea() << endl;

    NgProfiler::StopTimer (ts2);
    NgProfiler::StartTimer (ts3);

    const char * savetask = multithread.task;
    multithread.task = "Surface meshing";

    adfront.SetStartFront ();


    int plotnexttrial = 999;

    double meshedarea_before = meshedarea;

    NgProfiler::StopTimer (ts3);

    static Timer tloop("surfacemeshing mainloop");
    // static Timer tgetlocals("surfacemeshing getlocals");
    {
      RegionTimer rloop(tloop);
    while (!adfront.Empty() && !multithread.terminate)
      {
	NgProfiler::RegionTimer reg1 (timer1);

	if (multithread.terminate)
	  throw NgException ("Meshing stopped");

	// known for STL meshing
	if (totalarea > 0)
	  multithread.percent = 100 * meshedarea / totalarea;
	/*
	  else
	  multithread.percent = 0;
	*/

	locpoints.SetSize0();
	loclines.SetSize0();
	pindex.SetSize0();
	lindex.SetSize0();
	delpoints.SetSize0();
	dellines.SetSize0();
	locelements.SetSize0();



	// plot statistics
	if (trials > plotnexttrial)
	  {
	    PrintMessage (5, 
			  "faces = ", nfaces,
			  " trials = ", trials,
			  " elements = ", mesh.GetNSE(),
			  " els/sec = ",
			  (mesh.GetNSE() / (GetTime() - starttime + 0.0001)));
	    plotnexttrial += 1000;
	  }


	// unique-pgi, multi-pgi
	upgeominfo.SetSize0();
	mpgeominfo.SetSize0();


	nfaces = adfront.GetNFL();
	trials ++;
    

	if (trials % 1000 == 0)
	  {
	    (*testout) << "\n";
	    for (int i = 1; i <= canuse.Size(); i++)
	      {
		(*testout) << foundmap.Get(i) << "/" 
			   << canuse.Get(i) << "/"
			   << ruleused.Get(i) << " map/can/use rule " << rules[i-1]->Name() << "\n";
	      }
	    (*testout) << "\n";
	  }

        Point<3> p1, p2;
	int baselineindex = adfront.SelectBaseLine (p1, p2, blgeominfo1, blgeominfo2, qualclass);

	found = 1;

	double his = Dist (p1, p2);

	Point<3> pmid = Center (p1, p2);
	double hshould = CalcLocalH (pmid, mesh.GetH (pmid, layer));
	if (gh < hshould) hshould = gh;

	mesh.RestrictLocalH (pmid, hshould);

	double h = hshould;

	double hinner = (3 + qualclass) * max2 (his, hshould);

        // tgetlocals.Start();
	adfront.GetLocals (baselineindex, locpoints, mpgeominfo, loclines, 
			     pindex, lindex, 2*hinner);
        // tgetlocals.Stop();

	NgProfiler::RegionTimer reg2 (timer2);

	//(*testout) << "h for locals: " << 2*hinner << endl;
	

	//(*testout) << "locpoints " << locpoints << endl;

	if (qualclass > mp.giveuptol2d)
	  {
	    PrintMessage (3, "give up with qualclass ", qualclass);
	    PrintMessage (3, "number of frontlines = ", adfront.GetNFL());
	    // throw NgException ("Give up 2d meshing");
	    break;
	  }

	/*
	if (found && qualclass > 60)
	  {
	    found = 0;
	  }
	*/
	//      morerisc = ((qualclass > 20) && (qualclass % 2 == 1));
	//      morerisc = 1;
	morerisc = 0;


	PointIndex gpi1 = adfront.GetGlobalIndex (pindex.Get(loclines[0].I1()));
	PointIndex gpi2 = adfront.GetGlobalIndex (pindex.Get(loclines[0].I2()));


	debugflag = 
	  ( 
	   debugparam.haltsegment &&
	   ( ((debugparam.haltsegmentp1 == gpi1) && (debugparam.haltsegmentp2 == gpi2)) || 
	     ((debugparam.haltsegmentp1 == gpi2) && (debugparam.haltsegmentp2 == gpi1))) 
	   ) 
	  ||
	  (
	   debugparam.haltnode &&
	   ( (debugparam.haltsegmentp1 == gpi1) || (debugparam.haltsegmentp2 == gpi1))
	   );
	
      
	if (debugparam.haltface && debugparam.haltfacenr == facenr)
	  {
	    debugflag = 1;
	    cout << "set debugflag" << endl;
	  }
	
	if (debugparam.haltlargequalclass && qualclass == 50)
	  debugflag = 1;

	// problem recognition !
	if (found && 
	    (gpi1 < illegalpoint.Size()+PointIndex::BASE) && 
	    (gpi2 < illegalpoint.Size()+PointIndex::BASE) )
	  {
	    if (illegalpoint[gpi1] || illegalpoint[gpi2])
	      found = 0;
	  }


	// Point2d p12d, p22d;

	if (found)
	  {
	    oldnp = locpoints.Size();
	    oldnl = loclines.Size();

            UpdateVisSurfaceMeshData(oldnl);
	  
	    if (debugflag)
	      (*testout) << "define new transformation" << endl;

	    DefineTransformation (p1, p2, blgeominfo1, blgeominfo2);
	  
	    plainpoints.SetSize (locpoints.Size());
	    plainzones.SetSize (locpoints.Size());

	    // (*testout) << endl;

	    if (debugflag)
	      {
		*testout << "3d->2d transformation" << endl;
		*testout << "3d points: " << endl << locpoints << endl;
	      }


	    for (size_t i = 0; i < locpoints.Size(); i++)
              {
                Point<2> pp;
                TransformToPlain (locpoints[i], mpgeominfo[i],
                                  pp, h, plainzones[i]);
                plainpoints[i] = pp;
              }
            
            /*
	    for (int i = 1; i <= locpoints.Size(); i++)
	      {
		// (*testout) << "pindex(i) = " << pindex[i-1] << endl;
		TransformToPlain (locpoints.Get(i), 
				  mpgeominfo.Get(i),
				  plainpoints.Elem(i), h, plainzones.Elem(i));
		//		(*testout) << mpgeominfo.Get(i).GetPGI(1).u << " " << mpgeominfo.Get(i).GetPGI(1).v << " ";
		//		(*testout) << plainpoints.Get(i).X() << " " << plainpoints.Get(i).Y() << endl;
		//(*testout) << "transform " << locpoints.Get(i) << " to " << plainpoints.Get(i).X() << " " << plainpoints.Get(i).Y() << endl;
	      }
            */
	    //	    (*testout) << endl << endl << endl;


	    if (debugflag)
	      *testout << "2d points: " << endl << plainpoints << endl;


	    // p12d = plainpoints.Get(1);
	    // p22d = plainpoints.Get(2);

	    /*
	    // last idea on friday
	    plainzones.Elem(1) = 0;
	    plainzones.Elem(2) = 0;
	    */


	    /*
	    // old netgen:
	    for (i = 2; i <= loclines.Size(); i++)  // don't remove first line
	    {
	    z1 = plainzones.Get(loclines.Get(i).I1());
	    z2 = plainzones.Get(loclines.Get(i).I2());
	      
	    if (z1 && z2 && (z1 != z2) || (z1 == -1) || (z2 == -1) )
	    {
	    loclines.DeleteElement(i);
	    lindex.DeleteElement(i);
	    oldnl--;
	    i--;
	    }
	    }

	    // 	  for (i = 1; i <= plainpoints.Size(); i++)
	    // 	    if (plainzones.Elem(i) == -1)
	    // 	      plainpoints.Elem(i) = Point2d (1e4, 1e4);
	    */
	  

	  
	    for (int i = 2; i <= loclines.Size(); i++)  // don't remove first line
	      {
		// (*testout) << "loclines(i) = " << loclines.Get(i).I1() << " - " << loclines.Get(i).I2() << endl;
		z1 = plainzones.Get(loclines.Get(i).I1());
		z2 = plainzones.Get(loclines.Get(i).I2());
	      
	      
		// one inner point, one outer
		if ( (z1 >= 0) != (z2 >= 0))
		  {
		    int innerp = (z1 >= 0) ? 1 : 2;
		    if (IsLineVertexOnChart (locpoints.Get(loclines.Get(i).I1()),
					     locpoints.Get(loclines.Get(i).I2()),
					     innerp,
					     adfront.GetLineGeomInfo (lindex.Get(i), innerp)))
		      // pgeominfo.Get(loclines.Get(i).I(innerp))))
		      {		

			if (!morerisc)
			  {
			    // use one end of line
			    int pini = loclines.Get(i).I(innerp);
			    int pouti = loclines.Get(i).I(3-innerp);
			  
			    const auto& pin = plainpoints.Get(pini);
			    const auto& pout = plainpoints.Get(pouti);
			    auto v = pout - pin;
			    double len = v.Length();
			    if (len <= 1e-6)
			      (*testout) << "WARNING(js): inner-outer: short vector" << endl;
			    else
			      v /= len;
			  
			    /*
			    // don't elongate line towards base-line !!
			    if (Vec2d (pin, p12d) * v > 0 && 
			    Vec2d (pin, p22d) * v > 0)
			    v *= -1;  
			    */

			    Point<2> newpout = pin + 1000. * v;
			    newpout = pout;

			  
			    plainpoints.Append (newpout);
			    auto pout3d = locpoints.Get(pouti);
			    locpoints.Append (pout3d);

			    plainzones.Append (0);
			    pindex.Append (-1);
			    oldnp++;
			    loclines.Elem(i).I(3-innerp) = oldnp;
			  }
			else
			  plainzones.Elem(loclines.Get(i).I(3-innerp)) = 0;
			

			//		  (*testout) << "inner - outer correction" << endl;
		      }
		    else
		      {
			// remove line
			loclines.DeleteElement(i);
			lindex.DeleteElement(i);
			oldnl--;
			i--;
		      }			
		  }
	      
		else if ( (z1 > 0 && z2 > 0 && (z1 != z2)) || ((z1 < 0) && (z2 < 0)) )
		  {
		    loclines.DeleteElement(i);
		    lindex.DeleteElement(i);
		    oldnl--;
		    i--;
		  }
	      }
	  




	    legalpoints.SetSize(plainpoints.Size());
            legalpoints = 1;
            /*
	    for (int i = 1; i <= legalpoints.Size(); i++)
	      legalpoints.Elem(i) = 1;
            */
            
	    double avy = 0;
	    for (size_t i = 0; i < plainpoints.Size(); i++)
	      avy += plainpoints[i][1];
	    avy *= 1./plainpoints.Size();
		

	    for (auto i : Range(plainpoints))
	      {
		if (plainzones[i] < 0)
		  {
		    plainpoints[i] = {1e4, 1e4};
		    legalpoints[i] = 0;
		  }
		if (pindex[i] == -1)
		  {
		    legalpoints[i] = 0;
		  }
		    

		if (plainpoints[i][1] < -1e-10*avy) // changed
		  {
		    legalpoints[i] = 0;
		  }
	      }
	    /*
	      for (i = 3; i <= plainpoints.Size(); i++)
	      if (sqr (plainpoints.Get(i).X()) + sqr (plainpoints.Get(i).Y())
	      > sqr (2 + 0.2 * qualclass))
	      legalpoints.Elem(i) = 0;
	    */  

	    /*	  
		 int clp = 0;
		 for (i = 1; i <= plainpoints.Size(); i++)
		 if (legalpoints.Get(i))
		 clp++;
		 (*testout) << "legalpts: " << clp << "/" << plainpoints.Size() << endl; 

		 // sort legal/illegal lines
		 int lastleg = 2;
		 int firstilleg = oldnl;

		 while (lastleg < firstilleg)
		 {
		 while (legalpoints.Get(loclines.Get(lastleg).I1()) &&
		 legalpoints.Get(loclines.Get(lastleg).I2()) &&
		 lastleg < firstilleg)
		 lastleg++;
		 while ( ( !legalpoints.Get(loclines.Get(firstilleg).I1()) ||
		 !legalpoints.Get(loclines.Get(firstilleg).I2())) &&
		 lastleg < firstilleg)
		 firstilleg--;
	      
		 if (lastleg < firstilleg)
		 {
		 swap (loclines.Elem(lastleg), loclines.Elem(firstilleg));
		 swap (lindex.Elem(lastleg), lindex.Elem(firstilleg));
		 }
		 }

		 (*testout) << "leglines " << lastleg << "/" << oldnl << endl;
	    */
	

	    GetChartBoundary (chartboundpoints, 
			      chartboundpoints3d,
			      chartboundlines, h);

	    oldnp = plainpoints.Size();

	    maxlegalpoint = locpoints.Size();
	    maxlegalline = loclines.Size();



	    if (mp.checkchartboundary)
	      {
		for (int i = 1; i <= chartboundpoints.Size(); i++)
		  {
                    pindex.Append(-1);
		    plainpoints.Append (chartboundpoints.Get(i));
		    locpoints.Append (chartboundpoints3d.Get(i));
		    legalpoints.Append (0);
		  }
	      

		for (int i = 1; i <= chartboundlines.Size(); i++)
		  {
		    INDEX_2 line (chartboundlines.Get(i).I1()+oldnp,
				  chartboundlines.Get(i).I2()+oldnp);
		    loclines.Append (line);
		    //	      (*testout) << "line: " << line.I1() << "-" << line.I2() << endl;
		  }
	      }

	    oldnl = loclines.Size();
	    oldnp = plainpoints.Size();
	  }


	/*
	  if (qualclass > 100)
	  {
	  multithread.drawing = 1;
	  glrender(1);
	  cout << "qualclass 100, nfl = " << adfront.GetNFL() << endl;
	  }
	*/

	if (found)
	  {
            // static Timer t("ApplyRules");
            // RegionTimer r(t);
	    rulenr = ApplyRules (plainpoints, legalpoints, maxlegalpoint,
				 loclines, maxlegalline, locelements,
				 dellines, qualclass, mp);

	    //	    (*testout) << "Rule Nr = " << rulenr << endl;
	    if (!rulenr)
	      {
		found = 0;
		if ( debugflag || debugparam.haltnosuccess )
		  PrintWarning ("no rule found");
	      }
	  }
      
	NgProfiler::RegionTimer reg3 (timer3);


	for (int i = 1; i <= locelements.Size() && found; i++)
	  {
	    const Element2d & el = locelements.Get(i);

	    for (int j = 1; j <= el.GetNP(); j++)
	      if (el.PNum(j) <= oldnp && pindex.Get(el.PNum(j)) == -1)
		{
		  found = 0;
		  PrintSysError ("meshing2, index missing");
		}
	  }


	if (found)
	  {
	    locpoints.SetSize (plainpoints.Size());
	    upgeominfo.SetSize(locpoints.Size());

	    for (int i = oldnp+1; i <= plainpoints.Size(); i++)
	      {
                Point<3> locp;
                upgeominfo.Elem(i) = *blgeominfo1;
		int err =
		  TransformFromPlain (plainpoints.Elem(i), locp, 
				      upgeominfo.Elem(i), h);
                locpoints.Elem(i) = locp;

		if (err)
		  {
		    found = 0;

		    if ( debugflag || debugparam.haltnosuccess )
		      PrintSysError ("meshing2, Backtransformation failed");

		    break;
		  }
	      }
	  }
	  

	//      for (i = 1; i <= oldnl; i++)
	//        adfront.ResetClass (lindex[i]);


	/*
	  double violateminh;
	  if (qualclass <= 10)
	  violateminh = 3;
	  else
	  violateminh = 3 * qualclass;

	  if (uselocalh && found) //  && qualclass <= 10)
	  {
	  for (i = 1; i <= locelements.Size(); i++)
	  {
	  Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
	  Point3d pmax = pmin;
	  for (j = 2; j <= 3; j++)
	  {
	  const Point3d & hp = 
	  locpoints.Get(locelements.Get(i).PNum(j));
	  pmin.SetToMin (hp);
	  pmax.SetToMax (hp);
	  }
	  double minh = mesh.GetMinH (pmin, pmax);
	  if (h > violateminh * minh)
	  {
	  found = 0;
	  loclines.SetSize (oldnl);
	  locpoints.SetSize (oldnp);
	  }
	  }
	  }
	*/


	if (found) 
	  {
	    double violateminh = 3 + 0.1 * sqr (qualclass);
	    double minh = 1e8;
	    double newedgemaxh = 0;
	    for (int i = oldnl+1; i <= loclines.Size(); i++)
	      {
		double eh = Dist (locpoints.Get(loclines.Get(i).I1()),
				  locpoints.Get(loclines.Get(i).I2()));

		// Markus (brute force method to avoid bad elements on geometries like \_/ )
		//if(eh > 4.*mesh.GetH(locpoints.Get(loclines.Get(i).I1()))) found = 0;
		//if(eh > 4.*mesh.GetH(locpoints.Get(loclines.Get(i).I2()))) found = 0;
		// Markus end

		if (eh > newedgemaxh)
		  newedgemaxh = eh;
	      }

	    for (int i = 1; i <= locelements.Size(); i++)
	      {
		Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
		Point3d pmax = pmin;
		for (int j = 2; j <= locelements.Get(i).GetNP(); j++)
		  {
		    const Point3d & hp = 
		      locpoints.Get(locelements.Get(i).PNum(j));
		    pmin.SetToMin (hp);
		    pmax.SetToMax (hp);
		  }
		double eh = mesh.GetMinH (pmin, pmax);
		if (eh < minh)
		  minh = eh;
	      }

	    for (int i = 1; i <= locelements.Size(); i++)
	      for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		if (Dist2 (locpoints.Get(locelements.Get(i).PNum(j)), pmid) > hinner*hinner)
		  found = 0;

	    //	  cout << "violate = " << newedgemaxh / minh << endl;
	    static double maxviolate = 0;
	    if (newedgemaxh / minh > maxviolate)
	      {
		maxviolate = newedgemaxh / minh;
		//	      cout << "max minhviolate = " << maxviolate << endl;
	      }


	    if (newedgemaxh > violateminh * minh)
	      {
		found = 0;
		loclines.SetSize (oldnl);
		locpoints.SetSize (oldnp);

		if ( debugflag || debugparam.haltnosuccess )
		  PrintSysError ("meshing2, maxh too large");


	      }
	  }



	/*
	// test good ComputeLineGeoInfo
	if (found)
	{
	// is line on chart ?
	for (i = oldnl+1; i <= loclines.Size(); i++)
	{
	int gisize;
	void *geominfo;

	if (ComputeLineGeoInfo (locpoints.Get(loclines.Get(i).I1()),
	locpoints.Get(loclines.Get(i).I2()),
	gisize, geominfo))
	found = 0;
	}
	}
	*/


	// changed for OCC meshing
	if (found)
	  {
	    // take geominfo from dellines
	    // upgeominfo.SetSize(locpoints.Size());

	    /*
	      for (i = 1; i <= dellines.Size(); i++)
	      for (j = 1; j <= 2; j++)
	      {
	      upgeominfo.Elem(loclines.Get(dellines.Get(i)).I(j)) =
	      adfront.GetLineGeomInfo (lindex.Get(dellines.Get(i)), j);
	      }
	    */


	    for (int i = 1; i <= locelements.Size(); i++)
	      for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		{
		  int pi = locelements.Get(i).PNum(j);
		  if (pi <= oldnp)
		    {
		    
		      if (ChooseChartPointGeomInfo (mpgeominfo.Get(pi), upgeominfo.Elem(pi)))
			{
			  // cannot select, compute new one
			  PrintWarning ("calc point geominfo instead of using");
			  if (ComputePointGeomInfo (locpoints.Get(pi), upgeominfo.Elem(pi)))
			    {
			      found = 0;
			      PrintSysError ("meshing2d, geominfo failed");
			    }
			}
		    }
		}

	    /*
	    // use upgeominfo from ProjectFromPlane
	    for (i = oldnp+1; i <= locpoints.Size(); i++)
	    {
	    if (ComputePointGeomInfo (locpoints.Get(i), upgeominfo.Elem(i)))
	    {
	    found = 0;
	    if ( debugflag || debugparam.haltnosuccess )
	    PrintSysError ("meshing2d, compute geominfo failed");
	    }
	    }
	    */
	  }


	if (found && mp.checkoverlap)
	  {
	    // cout << "checkoverlap" << endl;
	    // test for overlaps
	  
	    Point3d hullmin(1e10, 1e10, 1e10);
	    Point3d hullmax(-1e10, -1e10, -1e10);
	  
	    for (int i = 1; i <= locelements.Size(); i++)
	      for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		{
		  const Point3d & p = locpoints.Get(locelements.Get(i).PNum(j));
		  hullmin.SetToMin (p);
		  hullmax.SetToMax (p);
		}
	    hullmin += Vec3d (-his, -his, -his);
	    hullmax += Vec3d ( his,  his,  his);

	    surfeltree.GetIntersecting (hullmin, hullmax, intersecttrias);

	    critpoints.SetSize (0);
	    for (int i = oldnp+1; i <= locpoints.Size(); i++)
	      critpoints.Append (locpoints.Get(i));

	    for (int i = 1; i <= locelements.Size(); i++)
	      {
		const Element2d & tri = locelements.Get(i);
		if (tri.GetNP() == 3)
		  {
		    const Point3d & tp1 = locpoints.Get(tri.PNum(1));
		    const Point3d & tp2 = locpoints.Get(tri.PNum(2));
		    const Point3d & tp3 = locpoints.Get(tri.PNum(3));
		  
		    Vec3d tv1 (tp1, tp2);
		    Vec3d tv2 (tp1, tp3);
		  
		    double lam1, lam2;
		    for (lam1 = 0.2; lam1 <= 0.8; lam1 += 0.2)
		      for (lam2 = 0.2; lam2 + lam1 <= 0.8; lam2 += 0.2)
			{
			  Point3d hp = tp1 + lam1 * tv1 + lam2 * tv2;
			  critpoints.Append (hp);
			}
		  }
		else if (tri.GetNP() == 4)
		  {
		    const Point3d & tp1 = locpoints.Get(tri.PNum(1));
		    const Point3d & tp2 = locpoints.Get(tri.PNum(2));
		    const Point3d & tp3 = locpoints.Get(tri.PNum(3));
		    const Point3d & tp4 = locpoints.Get(tri.PNum(4));
		  
		    double l1, l2;
		    for (l1 = 0.1; l1 <= 0.9; l1 += 0.1)
		      for (l2 = 0.1; l2 <= 0.9; l2 += 0.1)
			{
			  Point3d hp;
			  hp.X() = 
			    (1-l1)*(1-l2) * tp1.X() +
			    l1*(1-l2) * tp2.X() +
			    l1*l2 * tp3.X() +
			    (1-l1)*l2 * tp4.X();
			  hp.Y() = 
			    (1-l1)*(1-l2) * tp1.Y() +
			    l1*(1-l2) * tp2.Y() +
			    l1*l2 * tp3.Y() +
			    (1-l1)*l2 * tp4.Y();
			  hp.Z() = 
			    (1-l1)*(1-l2) * tp1.Z() +
			    l1*(1-l2) * tp2.Z() +
			    l1*l2 * tp3.Z() +
			    (1-l1)*l2 * tp4.Z();


			  critpoints.Append (hp);
			}
		  }
	      }
	    /*
	      for (i = oldnl+1; i <= loclines.Size(); i++)
	      {
	      Point3d hp = locpoints.Get(loclines.Get(i).I1());
	      Vec3d hv(hp, locpoints.Get(loclines.Get(i).I2()));
	      int ncp = 2;
	      for (j = 1; j <= ncp; j++)
	      critpoints.Append ( hp + (double(j)/(ncp+1)) * hv);
	      }
	    */


	    /*
	      for (i = oldnp+1; i <= locpoints.Size(); i++)
	      {
	      const Point3d & p = locpoints.Get(i);
	    */


	    for (int i = 1; i <= critpoints.Size(); i++)
	      {
		const Point3d & p = critpoints.Get(i);
		 
		for (int jj = 0; jj < intersecttrias.Size(); jj++)
		  {
		    // int j = intersecttrias.Get(jj);
		    // const Element2d & el = mesh.SurfaceElement(j);

		    SurfaceElementIndex j = intersecttrias[jj];
		    const Element2d & el = mesh[j];

		    int ntrig = (el.GetNP() == 3) ? 1 : 2;

		    int jl;
		    for (jl = 1; jl <= ntrig; jl++)
		      {
			Point3d tp1, tp2, tp3;

			if (jl == 1)
			  {
			    tp1 = mesh.Point(el.PNum(1));
			    tp2 = mesh.Point(el.PNum(2));
			    tp3 = mesh.Point(el.PNum(3));
			  }
			else
			  {
			    tp1 = mesh.Point(el.PNum(1));
			    tp2 = mesh.Point(el.PNum(3));
			    tp3 = mesh.Point(el.PNum(4));
			  }

			int onchart = 0;
			for (int k = 1; k <= el.GetNP(); k++)
			  if (BelongsToActiveChart (mesh.Point(el.PNum(k)),
						    el.GeomInfoPi(k)))
			    onchart = 1;
			if (!onchart)
			  continue;
		      
			Vec3d e1(tp1, tp2);
			Vec3d e2(tp1, tp3);
			Vec3d n = Cross (e1, e2);
			n /= n.Length();
			double lam1, lam2, lam3;
			lam3 = n * Vec3d (tp1, p);
			LocalCoordinates (e1, e2, Vec3d (tp1, p), lam1, lam2);
		      
			if (fabs (lam3) < 0.1 * hshould && 
			    lam1 > 0 && lam2 > 0 && (lam1 + lam2) < 1)
			  {
#ifdef DEVELOP
			    cout << "overlap" << endl;
			    (*testout) << "overlap:" << endl
				       << "tri = " << tp1 << "-" << tp2 << "-" << tp3 << endl
				       << "point = " << p << endl
				       << "lam1, 2 = " << lam1 << ", " << lam2 << endl
				       << "lam3 = " << lam3 << endl;
			  
			    //		      cout << "overlap !!!" << endl;
#endif
			    for (int k = 1; k <= 5; k++)
			      adfront.IncrementClass (lindex.Get(1));

			    found = 0;
			  
			    if ( debugflag || debugparam.haltnosuccess )
			      PrintWarning ("overlapping");
			  
			  
			    if (debugparam.haltoverlap)
			      {
				debugflag = 1;
			      }
			  
			    /*
			      multithread.drawing = 1;
			      glrender(1);
			    */
			  }
		      }
		  }
	      }
	  }


	if (found)
	  {
	    // check, whether new front line already exists

	    for (int i = oldnl+1; i <= loclines.Size(); i++)
	      {
		int nlgpi1 = loclines.Get(i).I1();
		int nlgpi2 = loclines.Get(i).I2();
		if (nlgpi1 <= pindex.Size() && nlgpi2 <= pindex.Size())
		  {
		    nlgpi1 = adfront.GetGlobalIndex (pindex.Get(nlgpi1));
		    nlgpi2 = adfront.GetGlobalIndex (pindex.Get(nlgpi2));

		    int exval = adfront.ExistsLine (nlgpi1, nlgpi2);
		    if (exval)
		      {
			cout << "ERROR: new line exits, val = " << exval << endl;
			(*testout) << "ERROR: new line exits, val = " << exval << endl;
			found = 0;


			if (debugparam.haltexistingline)
			  debugflag = 1;

		      }
		  }
	      }
	  
	  }


	/*
	  if (found)
	  {
	  // check, whether new triangles insert edges twice
	  for (i = 1; i <= locelements.Size(); i++)
	  for (j = 1; j <= 3; j++)
	  {
	  int tpi1 = locelements.Get(i).PNumMod (j);
	  int tpi2 = locelements.Get(i).PNumMod (j+1);
	  if (tpi1 <= pindex.Size() && tpi2 <= pindex.Size())
	  {
	  tpi1 = adfront.GetGlobalIndex (pindex.Get(tpi1));
	  tpi2 = adfront.GetGlobalIndex (pindex.Get(tpi2));

	  if (doubleedge.Used (INDEX_2(tpi1, tpi2)))
	  {
	  if (debugparam.haltexistingline)
	  debugflag = 1;
	  cerr << "ERROR Insert edge "
	  << tpi1 << " - " << tpi2 << " twice !!!" << endl;
	  found = 0;
	  }
	  doubleedge.Set (INDEX_2(tpi1, tpi2), 1);
	  }
	  }
	  }
	*/


	if (found)
	  {
	    // everything is ok, perform mesh update

	    ruleused.Elem(rulenr)++;


	    pindex.SetSize(locpoints.Size());
	      
	    for (int i = oldnp+1; i <= locpoints.Size(); i++)
	      {
		PointIndex globind = mesh.AddPoint (locpoints.Get(i));
		pindex.Elem(i) = adfront.AddPoint (locpoints.Get(i), globind);
	      }
	      
	    for (int i = oldnl+1; i <= loclines.Size(); i++)
	      {
		/*
		  for (j = 1; j <= locpoints.Size(); j++)
		  {
		  (*testout) << j << ": " << locpoints.Get(j) << endl;
		  }
		*/
	      
		/*
		  ComputeLineGeoInfo (locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()),
		  gisize, geominfo);
		*/		  

		if (pindex.Get(loclines.Get(i).I1()) == -1 || 
		    pindex.Get(loclines.Get(i).I2()) == -1)
		  {
		    (*testout) << "pindex is 0" << endl;
		  }

		if (!upgeominfo.Get(loclines.Get(i).I1()).trignum || 
		    !upgeominfo.Get(loclines.Get(i).I2()).trignum)
		  {
		    cout << "new el: illegal geominfo" << endl;
		  }

		adfront.AddLine (pindex.Get(loclines.Get(i).I1()),
				    pindex.Get(loclines.Get(i).I2()),
				    upgeominfo.Get(loclines.Get(i).I1()),
				    upgeominfo.Get(loclines.Get(i).I2()));
	      }
	    for (int i = 1; i <= locelements.Size(); i++)
	      {
		Element2d mtri(locelements.Get(i).GetNP());
		mtri = locelements.Get(i);
		mtri.SetIndex (facenr);


		// compute triangle geominfo:
		//	      (*testout) << "triggeominfo: ";
		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		    mtri.GeomInfoPi(j) = upgeominfo.Get(locelements.Get(i).PNum(j));
		    //		  (*testout) << mtri.GeomInfoPi(j).trignum << " ";
		  }
		//	      (*testout) << endl;

		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		    mtri.PNum(j) = 
		      locelements.Elem(i).PNum(j) =
		      adfront.GetGlobalIndex (pindex.Get(locelements.Get(i).PNum(j)));
		  }
	      
		
	      
	      
		mesh.AddSurfaceElement (mtri);
		cntelem++;
		//	      cout << "elements: " << cntelem << endl;


	      
		Box<3> box;
		box.Set (mesh[mtri[0]]);
		box.Add (mesh[mtri[1]]);
		box.Add (mesh[mtri[2]]);
		surfeltree.Insert (box, mesh.GetNSE()-1);

		const Point3d & sep1 = mesh.Point (mtri.PNum(1));
		const Point3d & sep2 = mesh.Point (mtri.PNum(2));
		const Point3d & sep3 = mesh.Point (mtri.PNum(3));

		double trigarea = Cross (Vec3d (sep1, sep2), 
					 Vec3d (sep1, sep3)).Length() / 2;

		if (mtri.GetNP() == 4)
		  {
		    const Point3d & sep4 = mesh.Point (mtri.PNum(4));
		    trigarea += Cross (Vec3d (sep1, sep3), 
				       Vec3d (sep1, sep4)).Length() / 2;
		  }

		meshedarea += trigarea;

		if(maxarea > 0 && meshedarea-meshedarea_before > maxarea)
		  {
		    cerr << "meshed area = " << meshedarea-meshedarea_before << endl
			 << "maximal area = " << maxarea << endl
			 << "GIVING UP" << endl;
		    return MESHING2_GIVEUP;
		  }
	      


		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		    int gpi = locelements.Get(i).PNum(j);

		    int oldts = trigsonnode.Size();
		    if (gpi >= oldts+PointIndex::BASE)
		      {
			trigsonnode.SetSize (gpi+1-PointIndex::BASE);
			illegalpoint.SetSize (gpi+1-PointIndex::BASE);
			for (int k = oldts+PointIndex::BASE; 
			     k <= gpi; k++)
			  {
			    trigsonnode[k] = 0;
			    illegalpoint[k] = 0;
			  }
		      }

		    trigsonnode[gpi]++;
		  
		    if (trigsonnode[gpi] > 20)
		      {
			illegalpoint[gpi] = 1;
			//		      cout << "illegal point: " << gpi << endl;
			(*testout) << "illegal point: " << gpi << endl;
		      }

		    static int mtonnode = 0;
		    if (trigsonnode[gpi] > mtonnode)
		      mtonnode = trigsonnode[gpi];
		  }
		//	      cout << "els = " << cntelem << " trials = " << trials << endl;
		//	      if (trials > 100)		return;
	      }
	      
	    for (int i = 1; i <= dellines.Size(); i++)
	      adfront.DeleteLine (lindex.Get(dellines.Get(i)));
	      
	    //	  rname = rules.Get(rulenr)->Name();
#ifdef MYGRAPH
	    if (silentflag<3) 
	      {
		plotsurf.DrawPnL(locpoints, loclines);
		plotsurf.Plot(testmode, testmode);
	      }
#endif

	    if (morerisc)
	      {
		cout << "generated due to morerisc" << endl;
		//	      multithread.drawing = 1;
		//	      glrender(1);
	      }



	  
	    if ( debugparam.haltsuccess || debugflag )
	      {
		// adfront.PrintOpenSegments (*testout);
		cout << "success of rule" << rules[rulenr-1]->Name() << endl;
		multithread.drawing = 1;
		multithread.testmode = 1;
		multithread.pause = 1;


		/*
		  extern STLGeometry * stlgeometry;
		  stlgeometry->ClearMarkedSegs();
		  for (i = 1; i <= loclines.Size(); i++)
		  {
		  stlgeometry->AddMarkedSeg(locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()));
		  }
		*/

		(*testout) << "success of rule" << rules[rulenr-1]->Name() << endl;
		(*testout) << "trials = " << trials << endl;

		(*testout) << "locpoints " << endl;
		for (int i = 1; i <= pindex.Size(); i++)
		  (*testout) << adfront.GetGlobalIndex (pindex.Get(i)) << endl;

		(*testout) << "old number of lines = " << oldnl << endl;

                UpdateVisSurfaceMeshData(oldnl);

		for (int i = 1; i <= loclines.Size(); i++)
		  {
		    (*testout) << "line ";
		    for (int j = 1; j <= 2; j++)
		      {
			int hi = 0;
			if (loclines.Get(i).I(j) >= 1 &&
			    loclines.Get(i).I(j) <= pindex.Size())
			  hi = adfront.GetGlobalIndex (pindex.Get(loclines.Get(i).I(j)));

			(*testout) << hi << " ";
		      }
		    (*testout) << " : " 
			       << plainpoints.Get(loclines.Get(i).I1()) << " - "
			       << plainpoints.Get(loclines.Get(i).I2()) << " 3d: "
			       << locpoints.Get(loclines.Get(i).I1()) << " - "
			       << locpoints.Get(loclines.Get(i).I2()) 
			       << endl;
		  }



		glrender(1);
	      }
	  }
	else
	  {
	    adfront.IncrementClass (lindex.Get(1));

	    if ( debugparam.haltnosuccess || debugflag )
	      {
		cout << "Problem with seg " << gpi1 << " - " << gpi2
		     << ", class = " << qualclass << endl;

		(*testout) << "Problem with seg " << gpi1 << " - " << gpi2
			   << ", class = " << qualclass << endl;

		multithread.drawing = 1;
		multithread.testmode = 1;
		multithread.pause = 1;


		/*
		  extern STLGeometry * stlgeometry;
		  stlgeometry->ClearMarkedSegs();
		  for (i = 1; i <= loclines.Size(); i++)
		  {
		  stlgeometry->AddMarkedSeg(locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()));
		  }
		*/

		for (int i = 1; i <= loclines.Size(); i++)
		  {
		    (*testout) << "line ";
		    for (int j = 1; j <= 2; j++)
		      {
			int hi = 0;
			if (loclines.Get(i).I(j) >= 1 &&
			    loclines.Get(i).I(j) <= pindex.Size())
			  hi = adfront.GetGlobalIndex (pindex.Get(loclines.Get(i).I(j)));

			(*testout) << hi << " ";
		      }
		    (*testout) << " : " 
			       << plainpoints.Get(loclines.Get(i).I1()) << " - "
			       << plainpoints.Get(loclines.Get(i).I2()) << " 3d: "
			       << locpoints.Get(loclines.Get(i).I1()) << " - "
			       << locpoints.Get(loclines.Get(i).I2()) 
			       << endl;
		  }


		/*
		  cout << "p1gi = " << blgeominfo[0].trignum 
		  << ", p2gi = " << blgeominfo[1].trignum << endl;
		*/

		glrender(1);
	      }

	  
#ifdef MYGRAPH      
	    if (silentflag<3)
	      {
		if (testmode || trials%2 == 0)
		  {
		    plotsurf.DrawPnL(locpoints, loclines);
		    plotsurf.Plot(testmode, testmode);
		  }
	      }
#endif
	  }

      }
    }
    PrintMessage (3, "Surface meshing done");


    adfront.PrintOpenSegments (*testout);

    multithread.task = savetask;


    EndMesh ();


    if (!adfront.Empty())
      return MESHING2_GIVEUP;
    
    return MESHING2_OK;
  }


}


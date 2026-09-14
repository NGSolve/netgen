#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{


  /// quad badness of a local element, via the mesh Element2d implementation
  static double CalcJacobianBadness (const MiniElement2d & elem,
                                     const Array<Point<2>, LocalPointIndex> & points)
  {
    Array<Point<2>, PointIndex> hpoints(points.Size());
    for (LocalPointIndex pi : points.Range())
      hpoints[pi.Nr0()+IndexBASE<PointIndex>()] = points[pi];

    Element2d hel(elem.GetNP());
    for (int j = 1; j <= elem.GetNP(); j++)
      hel.PNum(j) = elem.PNum(j).Nr0()+IndexBASE<PointIndex>();

    return hel.CalcJacobianBadness (hpoints);
  }

  static double CalcElementBadness (const Array<Point<2>, LocalPointIndex> & points,
				    const MiniElement2d & elem)
  {
    // badness = sqrt(3) /36 * circumference^2 / area - 1 +
    //           h / li + li / h - 2

    Vec<2> v12, v13, v23;
    double l12, l13, l23, cir, area;
    static const double c = sqrt(3.0) / 36;

    v12 = points[elem.PNum(2)] - points[elem.PNum(1)];
    v13 = points[elem.PNum(3)] - points[elem.PNum(1)];
    v23 = points[elem.PNum(3)] - points[elem.PNum(2)];

    l12 = v12.Length();
    l13 = v13.Length();
    l23 = v23.Length();

    cir = l12 + l13 + l23;
    area = 0.5 * (v12[0] * v13[1] - v12[1] * v13[0]);
    if (area < 1e-6)
      {
	return 1e8;
      }

    if (testmode)
      {
	(*testout) << "l = " << l12 << " + " << l13 << " + " << l23 << " = " 
		   << cir << ", area = " << area << endl;
	(*testout) << "shapeerr = " << 10 * (c * cir * cir / area - 1) << endl
		   << "sizeerr = " << 1/l12 + l12 + 1/l13 + l13 + 1/l23 + l23 - 6
		   << endl;
      }

    return 10 * (c * cir * cir / area - 1)
      + 1/l12 + l12 + 1/l13 + l13 + 1/l23 + l23 - 6;
  }



  int Meshing2 ::ApplyRules (Array<Point<2>, LocalPointIndex> & lpoints, 
			     Array<int, LocalPointIndex> & legalpoints,
			     int maxlegalpoint,
			     Array<IVec<2,LocalPointIndex>> & llines1,
			     int maxlegalline,
			     Array<MiniElement2d> & elements,
			     Array<INDEX> & dellines, int tolerance,
			     const MeshingParameters & mp)
  {
    // static Timer timer ("meshing2::ApplyRules"); RegionTimer reg (timer);


    double maxerr = 0.5 + 0.3 * tolerance;
    double minelerr = 2 + 0.5 * tolerance * tolerance;

    int noldlp = lpoints.Size();
    int noldll = llines1.Size();


    Array<int,LocalPointIndex> pused(maxlegalpoint);
    ArrayMem<int,100> lused(maxlegalline);
    Array<int,LocalPointIndex> pnearness(noldlp);
    ArrayMem<int,100> lnearness(llines1.Size());

    ArrayMem<LocalPointIndex, 20, RulePointIndex> pmap;   // rule point -> local point
    ArrayMem<bool, 20, RulePointIndex> pfixed;
    ArrayMem<int, 20> lmap;
  
    ArrayMem<Point<2>,100> tempnewpoints;
    ArrayMem<IVec<2,LocalPointIndex>,100> tempnewlines;
    ArrayMem<int,100> tempdellines;
    ArrayMem<MiniElement2d,100> tempelements;

    // a least 2 * maximal number of old points in rules,
    // what is actually 4 now
    double oldumem[20];  

    elements.SetSize (0);
    dellines.SetSize (0);

    testmode = debugparam.debugoutput;

#ifdef LOCDEBUG
    int loctestmode = testmode;

    if (loctestmode)
      {
	(*testout) << endl << endl << "Check new environment" << endl;
	(*testout) << "tolerance = " << tolerance << endl;
	for (int i = 1; i <= lpoints.Size(); i++)
	  (*testout) << "P" << i << " = " << lpoints[i] << endl;
	(*testout) << endl;
	for (int i = 1; i <= llines1.Size(); i++)
	  (*testout) << "(" << llines1.Get(i).I1() << "-" << llines1.Get(i).I2() << ")" << endl;
      }
#endif

    // check every rule

    int found = 0;   // rule number

    pnearness = 1000;
  
    for (int j = 0; j < 2; j++)
      pnearness[llines1[0][j]] = 0;



    enum { MAX_NEARNESS = 3 };

    for (int cnt = 0; cnt < MAX_NEARNESS; cnt++)
      {
	bool ok = true;
	for (int i = 0; i < maxlegalline; i++)
	  {
	    const IVec<2,LocalPointIndex> & hline = llines1[i];

	    int minn = min2 (pnearness[hline[0]],  pnearness[hline[1]]);

	    for (int j = 0; j < 2; j++)
	      if (pnearness[hline[j]] > minn+1)
		{
		  ok = false;
		  pnearness[hline[j]] = minn+1;
		}
	  }
	if (!ok) break;
      }


    for (int i = 0; i < maxlegalline; i++)
      lnearness[i] = pnearness[llines1[i][0]] + pnearness[llines1[i][1]];


    // resort lines after lnearness
    Array<IVec<2,LocalPointIndex>> llines(llines1.Size());
    Array<int> sortlines(llines1.Size());
    int lnearness_class[MAX_NEARNESS];

    for (int j = 0; j < MAX_NEARNESS; j++)
      lnearness_class[j] = 0;
    for (int i = 0; i < maxlegalline; i++)
      if (lnearness[i] < MAX_NEARNESS)
	lnearness_class[lnearness[i]]++;
    
    int cumm = 0;
    for (int j = 0; j < MAX_NEARNESS; j++)
      {
	int hcnt = lnearness_class[j];
	lnearness_class[j] = cumm;
	cumm += hcnt;
      }

    for (int i = 0; i < maxlegalline; i++)
      if (lnearness[i] < MAX_NEARNESS)
	{
	  llines[lnearness_class[lnearness[i]]] = llines1[i];
	  sortlines[lnearness_class[lnearness[i]]] = i+1;
	  lnearness_class[lnearness[i]]++;
	}
      else
	{
	  llines[cumm] = llines1[i];
	  sortlines[cumm] = i+1;
	  cumm++;
	}

    for (int i = maxlegalline; i < llines1.Size(); i++)
      {
	llines[cumm] = llines1[i];
	sortlines[cumm] = i+1;
	cumm++;
      }

    for (int i = 0; i < maxlegalline; i++)
      lnearness[i] = pnearness[llines[i][0]] + pnearness[llines[i][1]];




    static bool firsttime = true;
    // static int timers[100];
    // static int timers2[100];
    // static int timers3[100];
    if (firsttime)
      {
	/*
	for (int ri = 0; ri < rules.Size(); ri++)
	  timers[ri] = NgProfiler::CreateTimer (string("netrule ")+rules[ri]->Name());
	for (int ri = 0; ri < rules.Size(); ri++)
	  timers2[ri] = NgProfiler::CreateTimer (string("netrule,mapped ")+rules[ri]->Name());
	for (int ri = 0; ri < rules.Size(); ri++)
	  timers3[ri] = NgProfiler::CreateTimer (string("netrule,lines mapped ")+rules[ri]->Name());
	*/
	firsttime = false;
      }

    lused = 0;
    pused = 0;


    static Timer timer1("meshing2::ApplyRules 1");
    RegionTimer reg1 (timer1);


    for (int ri = 1; ri <= rules.Size(); ri++)
      {
	// NgProfiler::RegionTimer reg(timers[ri-1]);
	netrule * rule = rules[ri-1].get();

#ifdef LOCDEBUG
	if (loctestmode)
	  (*testout) << "Rule " << rule->Name() << endl;
#endif

	if (rule->GetQuality() > tolerance) continue;

	pmap.SetSize (rule->GetNP());
	lmap.SetSize (rule->GetNL());
      
	for (auto & p : pmap) p.Invalidate();
	lmap = 0;

	lused[0] = 1; 
	lmap[0] = 1;  

	for (int j = 0; j < 2; j++)
	  {
	    pmap[rule->GetLine(1)[j]] = llines[0][j];
	    pused[llines[0][j]]++;
	  }



	int nlok = 2;


	bool ok = false;

	while (nlok >= 2)
	  {

	    if (nlok <= rule->GetNOldL())

	      {
		ok = 0;
		
		int maxline = (rule->GetLNearness(nlok) < MAX_NEARNESS) ? lnearness_class[rule->GetLNearness(nlok)] : maxlegalline;
		// int maxline = maxlegalline;

		while (!ok && lmap[nlok-1] < maxline)
		  {
		    lmap[nlok-1]++;
		    int locli = lmap[nlok-1];

		    if (lnearness[locli-1] > rule->GetLNearness (nlok) ) continue;
		    if (lused[locli-1]) continue;


		    ok = 1;

		    IVec<2,LocalPointIndex> loclin = llines[locli-1];
		    auto linevec = lpoints[loclin[1]] - lpoints[loclin[0]];

		    if (rule->CalcLineError (nlok, linevec) > maxerr)
		      {
			ok = 0;
#ifdef LOCDEBUG
			if(loctestmode)
			  (*testout) << "not ok pos1" << endl;
#endif
			continue;
		      }

		    for (int j = 0; j < 2; j++)
		      {
			RulePointIndex refpi = rule->GetLine(nlok)[j];

			if (pmap[refpi].IsValid())
			  {
			    if (pmap[refpi] != loclin[j])
			      {
				ok = 0;
#ifdef LOCDEBUG
				if(loctestmode)
				  (*testout) << "not ok pos2" << endl;
#endif
				break;
			      }
			  }
			else
			  {
			    if (rule->CalcPointDist (refpi, lpoints[loclin[j]]) > maxerr
				|| !legalpoints[loclin[j]]
				|| pused[loclin[j]])
			      {
				ok = 0;
#ifdef LOCDEBUG
				if(loctestmode)
				  {
				    (*testout) << "nok pos3" << endl;
				    //if(rule->CalcPointDist (refpi, lpoints[loclin[j]]) > maxerr)
				    //(*testout) << "r1" << endl;
				    //if(!legalpoints[loclin[j]])
				    //(*testout) << "r2 legalpoints " << legalpoints << " loclin " << loclin << " j " << j << endl;
				    //if(pused[loclin[j]])
				    //(*testout) << "r3" << endl;
				  }
#endif
				break;
			      }
			  }
		      }
		  }

		if (ok)
		  {
		    int locli = lmap[nlok-1];
		    IVec<2,LocalPointIndex> loclin = llines[locli-1];

		    lused[locli-1] = 1;
		    for (int j = 0; j < 2; j++)
		      {
			pmap[rule->GetLine (nlok)[j]] = loclin[j];
			pused[loclin[j]]++;
		      }

		    nlok++;
		  }
		else
		  {
		    lmap[nlok-1] = 0;
		    nlok--;

		    lused[lmap[nlok-1]-1] = 0;
		    for (int j = 0; j < 2; j++)
		      {
			pused[llines[lmap[nlok-1]-1][j]] --;
			if (! pused[llines[lmap[nlok-1]-1][j]])
			  pmap[rule->GetLine (nlok)[j]].Invalidate();
		      }
		  }
	      }

	    else

	      {
		// NgProfiler::RegionTimer reg(timers3[ri-1]);

		// all lines are mapped !!

		// map also all points:

		RulePointIndex npok = IndexBASE<RulePointIndex>();
		int incnpok = 1;

		pfixed.SetSize (pmap.Size());
		for (auto i : pmap.Range())
		  pfixed[i] = pmap[i].IsValid();
 
		while (npok >= IndexBASE<RulePointIndex>())
		  {

		    if (npok <= RuleP(rule->GetNOldP()))

		      {
			if (pfixed[npok])

			  {
			    if (incnpok)
			      npok++;
			    else
			      npok--;
			  }

			else

			  {
			    ok = 0;

			    if (pmap[npok].IsValid())
			      pused[pmap[npok]]--;

			    while (!ok && pmap[npok] < maxlegalpoint+IndexBASE<LocalPointIndex>()-1)
			      {
				ok = 1;

				pmap[npok]++;

				if (pused[pmap[npok]])
				  {
				    ok = 0;
				  }
				else
				  {
				    if (rule->CalcPointDist (npok, lpoints[pmap[npok]]) > maxerr 
					|| !legalpoints[pmap[npok]]) 
                                    
				      ok = 0;
				  }
			      }

			    if (ok)
			      {
				pused[pmap[npok]]++;
				npok++;
				incnpok = 1;
			      }

			    else

			      {
				pmap[npok].Invalidate();
				npok--;
				incnpok = 0;
			      }
			  }
		      }

		    else

		      {
			// NgProfiler::RegionTimer reg(timers2[ri-1]);

			npok = RuleP(rule->GetNOldP());
			incnpok = 0;

			if (ok)
			  foundmap[ri-1]++; 

#ifdef LOCDEBUG
			if (loctestmode)
			  (*testout) << "lines and points mapped" << endl;
#endif

			ok = 1;

			// check orientations

			for (int i = 1; i <= rule->GetNOrientations(); i++)
			  {
			    if (CW (lpoints[pmap[rule->GetOrientation(i).i1]],
				    lpoints[pmap[rule->GetOrientation(i).i2]],
				    lpoints[pmap[rule->GetOrientation(i).i3]]) )
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "Orientation " << i << " not ok" << endl;
#endif
				break;
			      }
			  }


			if (!ok) continue;

			// Vector oldu (2 * rule->GetNOldP());
                        Vector oldu (2 * rule->GetNOldP(), &oldumem[0]);
		      
			for (auto pi : pmap.Range().Modify(0, rule->GetNOldP()-pmap.Size()))
			  {
			    Vec<2> ui(rule->GetPoint(pi), lpoints[pmap[pi]]);
			    int i = pi.Nr1();
			    oldu (2*i-2) = ui(0);
			    oldu (2*i-1) = ui(1);
			  }
		      
			rule -> SetFreeZoneTransformation (oldu, tolerance);

		      
			if (!ok) continue;
			if (!rule->ConvexFreeZone())
			  {
			    ok = 0;
#ifdef LOCDEBUG
			    if (loctestmode) 
			      (*testout) << "freezone not convex" << endl;
#endif
			    /*
			      static int cnt = 0;
			      cnt++;
			      if (cnt % 100 == 0)
			      {
			      cout << "freezone not convex, cnt = " << cnt << "; rule = " << rule->Name() << endl;
			      (*testout) << "freezone not convex, cnt = " << cnt << "; rule = " << rule->Name() << endl;
			      (*testout) << "tol = " << tolerance << endl;
			      (*testout) << "maxerr = " << maxerr << "; minerr = " << minelerr << endl;
			      (*testout) << "freezone = " << rule->GetTransFreeZone() << endl;
			      }
			    */
			  }

			// check freezone:
			if (!ok) continue;
			for (auto i : lpoints.Range().Modify(0, maxlegalpoint-lpoints.Size()))
			  {
			    if (!ok) break;
			    if ( !pused[i] &&
				 rule->IsInFreeZone (lpoints[i]) )
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "Point " << i << " in freezone" << endl;
#endif
				break;
			      }
			  }

			if (!ok) continue;
			for (auto i : lpoints.Range().Modify(maxlegalpoint, 0))
			  {
			    if ( rule->IsInFreeZone (lpoints[i]) )
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "Point " << i << " in freezone" << endl;
#endif
				break;
			      }
			  }


			if (!ok) continue;
			for (int i = 1; i <= maxlegalline; i++)
			  {
			    if (!lused[i-1] && 
				rule->IsLineInFreeZone (lpoints[llines[i-1][0]],
							lpoints[llines[i-1][1]]))
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "line " << llines.Get(i)[0] << "-"
					     << llines.Get(i)[1] << " in freezone" << endl;
#endif
				break;
			      }
			  }

			if (!ok) continue;

			for (int i = maxlegalline+1; i <= llines.Size(); i++)
			  {
			    if (rule->IsLineInFreeZone (lpoints[llines[i-1][0]],
							lpoints[llines[i-1][1]]))
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "line " << llines.Get(i)[0] << "-"
					     << llines.Get(i)[1] << " in freezone" << endl;
#endif
				break;
			      }
			  }


			/*
			// check orientations

			for (i = 1; i <= rule->GetNOrientations() && ok; i++)
			{
			if (CW (lpoints[pmap[rule->GetOrientation(i).i1]],
			lpoints[pmap[rule->GetOrientation(i).i2]],
			lpoints[pmap[rule->GetOrientation(i).i3]]) )
			{
			ok = 0;
			if (loctestmode)
			(*testout) << "Orientation " << i << " not ok" << endl;
			}
			}
			*/


			if (!ok) continue;

#ifdef LOCDEBUG
			if (loctestmode)
			  (*testout) << "rule ok" << endl;
#endif

			// Setze neue Punkte:
			if (rule->GetNOldP() < rule->GetNP())
			  {
			    Vector newu(rule->GetOldUToNewU().Height());
			    rule->GetOldUToNewU().Mult (oldu, newu);
			    
			    int oldnp = rule->GetNOldP();
			    for (auto pi : pmap.Range().Modify(oldnp, 0))
			      {
				auto np = rule->GetPoint(pi);
				int i = pi.Nr1();
				np[0] += newu (2 * (i-oldnp) - 2);
				np[1] += newu (2 * (i-oldnp) - 1);

                                lpoints.Append (np);
				pmap[pi] = lpoints.Range().Next()-1;
			      }
			  }

			// Setze neue Linien:

			for (int i = rule->GetNOldL() + 1; i <= rule->GetNL(); i++)
			  {
			    llines.Append (IVec<2,LocalPointIndex> (pmap[rule->GetLine (i)[0]],
                                                                   pmap[rule->GetLine (i)[1]]));
			  }


			// delete old lines:
			for (int i = 1; i <= rule->GetNDelL(); i++)
			  dellines.Append (sortlines[lmap[(rule->GetDelLine(i))-1]-1]);
			// dellines.Append (lmap.Get(rule->GetDelLine(i))));

			// dellines.Append (lmap.Elem(rule->GetDelLines()));
			// lmap[rule->GetDelLines()];


			// insert new elements:

			for (int i = 1; i <= rule->GetNE(); i++)
			  {
			    const RuleElement2d & rel = rule->GetElement(i);
			    MiniElement2d el(rel.GetNP());
			    for (int j = 1; j <= rel.GetNP(); j++)
			      el.PNum(j) = pmap[rel.PNum(j)];   // rule nr -> local nr
			    elements.Append (el);
			  }


			double elerr = 0;
			for (int i = 1; i <= elements.Size(); i++)
			  {
			    double hf;
			    if (!mp.quad)
			      hf = CalcElementBadness (lpoints, elements[i-1]);
			    else
			      hf = CalcJacobianBadness (elements[i-1], lpoints) * 5;
#ifdef LOCDEBUG
			    if (loctestmode)
			      (*testout) << "r " << rule->Name() << "bad = " << hf << endl;
#endif
			    if (hf > elerr) elerr = hf;
			  }

#ifdef LOCDEBUG
			if (loctestmode)
			  (*testout) << "error = " << elerr;
#endif

			canuse[ri-1] ++;

			if (elerr < 0.99*minelerr)
			  {
#ifdef LOCDEBUG
			    if (loctestmode)
			      {
				(*testout) << "rule = " << rule->Name() << endl;
				(*testout) << "class = " << tolerance << endl;
				(*testout) << "lpoints: " << endl;
				for (int i = 1; i <= lpoints.Size(); i++)
				  (*testout) << lpoints[i] << endl;
				(*testout) << "llines: " << endl;
				for (int i = 1; i <= llines.Size(); i++)
				  (*testout) << llines.Get(i)[0] << " " << llines.Get(i)[1] << endl;

				(*testout) << "Freezone: ";
				for (int i = 1; i <= rule -> GetTransFreeZone().Size(); i++)
				  (*testout) << rule->GetTransFreeZone().Get(i) << endl;
			      }
#endif

			    minelerr = elerr;
			    found = ri;

			    tempnewpoints = lpoints.Range (noldlp, lpoints.Size());
			    tempnewlines = llines.Range (noldll, llines.Size());
			    tempdellines = dellines;
			    tempelements = elements;
			  }

			lpoints.SetSize (noldlp);
			llines.SetSize (noldll);
			dellines.SetSize (0);
			elements.SetSize (0);
			ok = 0;
		      }
		  }

		nlok = rule->GetNOldL();

		lused[lmap[nlok-1]-1] = 0;

		for (int j = 1; j <= 2; j++)
		  {
		    RulePointIndex refpi = rule->GetPointNr (nlok, j);
		    if (!pmap[refpi].IsValid()) continue;   // point not mapped
		    pused[pmap[refpi]]--;

		    if (pused[pmap[refpi]] == 0)
		      pmap[refpi] = LocalPointIndex::INVALID;
		  }
	      }
	  }
      }


    if (found)
      {
	lpoints.Append (tempnewpoints);
	llines1.Append (tempnewlines);
	dellines.Append (tempdellines);
	elements.Append (tempelements);
      }


    return found;
  }





}

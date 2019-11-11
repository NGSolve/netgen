/*
  Advancing front class for surfaces
*/

#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{
  FrontPoint2 :: FrontPoint2 (const Point<3> & ap, PointIndex agi,
			      MultiPointGeomInfo * amgi, bool aonsurface)
  {
    p = ap;
    globalindex = agi;
    nlinetopoint = 0;
    frontnr = INT_MAX-10;
    onsurface = aonsurface;

    if (amgi)
      {
	mgi = new MultiPointGeomInfo (*amgi);
	for (int i = 1; i <= mgi->GetNPGI(); i++)
	  if (mgi->GetPGI(i).trignum <= 0)
	    cout << "WARNING: Add FrontPoint2, illegal geominfo = " << mgi->GetPGI(i).trignum << endl;
      }
    else
      mgi = NULL;
  }


  AdFront2 :: AdFront2 (const Box3d & aboundingbox)
    : boundingbox(aboundingbox), 
      linesearchtree(boundingbox.PMin(), boundingbox.PMax()),
      pointsearchtree(boundingbox.PMin(), boundingbox.PMax()),
      cpointsearchtree(boundingbox.PMin(), boundingbox.PMax())
  {
    nfl = 0;
    allflines = 0;

    minval = 0;
    // starti = lines.Begin();
    starti = *lines.Range().begin();
  }

  AdFront2 :: ~AdFront2 ()
  {
    delete allflines;
  }


  void AdFront2 :: PrintOpenSegments (ostream & ost) const
  {
    if (nfl > 0)
      {
	ost << nfl << " open front segments left:" << endl;
	// for (int i = lines.Begin(); i < lines.End(); i++)
        for (int i : lines.Range())
	  if (lines[i].Valid())
	    ost << i << ": " 
                << GetGlobalIndex (lines[i].L().I1()) << "-"
		<< GetGlobalIndex (lines[i].L().I2()) << endl;
      }
  }

  /*
  void AdFront2 :: GetPoints (NgArray<Point<3> > & apoints) const
  {
    apoints.Append (points);
    // for (int i = 0; i < points.Size(); i++)
    // apoints.Append (points[i].P());
  }
  */



  int AdFront2 :: AddPoint (const Point<3> & p, PointIndex globind, 
                            MultiPointGeomInfo * mgi,
                            bool pointonsurface)
  {
    // inserts at empty position or resizes array
    int pi;

    if (delpointl.Size() != 0)
      {
	pi = delpointl.Last();
	delpointl.DeleteLast ();

	points[pi] = FrontPoint2 (p, globind, mgi, pointonsurface);
      }
    else
      {
	points.Append (FrontPoint2 (p, globind, mgi, pointonsurface));
        pi = points.Size()-1;
      }

    if (mgi)
      cpointsearchtree.Insert (p, pi);

    if (pointonsurface)
      pointsearchtree.Insert (p, pi);
    
    return pi;
  }


  int AdFront2 :: AddLine (int pi1, int pi2,
                           const PointGeomInfo & gi1, const PointGeomInfo & gi2)
  {
    int minfn;
    int li;

    FrontPoint2 & p1 = points[pi1];
    FrontPoint2 & p2 = points[pi2];


    nfl++;

    p1.AddLine();
    p2.AddLine();

    minfn = min2 (p1.FrontNr(), p2.FrontNr());
    p1.DecFrontNr (minfn+1);
    p2.DecFrontNr (minfn+1);

    if (dellinel.Size() != 0)
      {
	li = dellinel.Last();
	dellinel.DeleteLast ();
	lines[li] = FrontLine (INDEX_2(pi1, pi2));
      }
    else
      {
	lines.Append(FrontLine (INDEX_2(pi1, pi2)));
        li = lines.Size()-1;
      }

  
    if (!gi1.trignum || !gi2.trignum)
      {
	cout << "WARNING: in AdFront::AddLine, illegal geominfo" << endl;
      }
  
    lines[li].SetGeomInfo (gi1, gi2);

    Box3d lbox;
    lbox.SetPoint(p1.P());
    lbox.AddPoint(p2.P());

    linesearchtree.Insert (lbox.PMin(), lbox.PMax(), li);

    if (allflines)
      {
	if (allflines->Used (INDEX_2 (GetGlobalIndex (pi1), 
				      GetGlobalIndex (pi2))))
	  {
	    cerr << "ERROR Adfront2::AddLine: line exists" << endl;
	    (*testout) << "ERROR Adfront2::AddLine: line exists" << endl;
	  }

	allflines->Set (INDEX_2 (GetGlobalIndex (pi1), 
				 GetGlobalIndex (pi2)), 1);
      }

    return li;
  }


  void AdFront2 :: DeleteLine (int li)
  {
    int pi;

    nfl--;

    for (int i = 1; i <= 2; i++)
      {
	pi = lines[li].L().I(i);
	points[pi].RemoveLine();

	if (!points[pi].Valid())
	  {
	    delpointl.Append (pi);
	    if (points[pi].mgi)
	      {
		cpointsearchtree.DeleteElement (pi);
		delete points[pi].mgi;
		points[pi].mgi = NULL;
	      }

            pointsearchtree.DeleteElement (pi);
	  }
      }

    if (allflines)
      {
	allflines->Set (INDEX_2 (GetGlobalIndex (lines[li].L().I1()),
				 GetGlobalIndex (lines[li].L().I2())), 2);
      }

    lines[li].Invalidate();
    linesearchtree.DeleteElement (li);

    dellinel.Append (li);
  }


  int AdFront2 :: ExistsLine (int pi1, int pi2)
  {
    if (!allflines)
      return 0;
    if (allflines->Used (INDEX_2(pi1, pi2)))
      return allflines->Get (INDEX_2 (pi1, pi2));
    else
      return 0;
  }


  int AdFront2 :: SelectBaseLine (Point<3>  & p1, Point<3>  & p2, 
				  const PointGeomInfo *& geominfo1,
				  const PointGeomInfo *& geominfo2,
				  int & qualclass)
  {
    int baselineindex = -1; 
    
    // for (int i = starti; i < lines.End(); i++)
    for (int i = starti; i < *lines.Range().end(); i++)
      {
	if (lines[i].Valid())
	  {
	    int hi = lines[i].LineClass() +
	      points[lines[i].L().I1()].FrontNr() +
	      points[lines[i].L().I2()].FrontNr();
	  
	    if (hi <= minval)
	      {
		minval = hi;
		baselineindex = i;
		break;
	      }
	  }
      }
  
    if (baselineindex == -1)
      {
	minval = INT_MAX;
	// for (int i = lines.Begin(); i < lines.End(); i++)
        for (int i : lines.Range())
	  if (lines[i].Valid())
	    {
	      int hi = lines[i].LineClass() +
		points[lines[i].L().I1()].FrontNr() +
		points[lines[i].L().I2()].FrontNr();
	    
	      if (hi < minval)
		{
		  minval = hi;
		  baselineindex = i;
		}
	    }
      }
    starti = baselineindex+1;

    p1 = points[lines[baselineindex].L().I1()].P();
    p2 = points[lines[baselineindex].L().I2()].P();
    geominfo1 = &lines[baselineindex].GetGeomInfo(1);
    geominfo2 = &lines[baselineindex].GetGeomInfo(2);

    qualclass = lines[baselineindex].LineClass();

    return baselineindex;
  }




  int AdFront2 :: GetLocals (int baselineindex,
			     NgArray<Point<3>> & locpoints,
			     NgArray<MultiPointGeomInfo> & pgeominfo,
			     NgArray<INDEX_2> & loclines,   // local index
			     NgArray<INDEX> & pindex,
			     NgArray<INDEX> & lindex,
			     double xh)
  {
    static Timer timer("adfront2::GetLocals"); RegionTimer reg (timer);
    
    int pstind;
    Point<3>  midp, p0;

    pstind = lines[baselineindex].L().I1();
    p0 = points[pstind].P();

    loclines.Append(lines[baselineindex].L());
    lindex.Append(baselineindex);  

    ArrayMem<int, 1000> nearlines(0);
    NgArrayMem<int, 1000> nearpoints(0);

    // dominating costs !!
    linesearchtree.GetIntersecting (p0 - Vec3d(xh, xh, xh),
				    p0 + Vec3d(xh, xh, xh),
				    nearlines);

    // only special points that are not in adfront,
    // other points are from linesearchtree
    cpointsearchtree.GetIntersecting(p0 - Vec3d(xh, xh, xh),
                                     p0 + Vec3d(xh, xh, xh),
                                     nearpoints);

    for(auto i : nearlines)
      {
	if (lines[i].Valid() && i != baselineindex) 
	  {
            loclines.Append(lines[i].L());
            lindex.Append(i);
	  }
      }

    // static NgArray<int> invpindex;
    invpindex.SetSize (points.Size()); 
    // invpindex = -1;
    for(auto pi : nearpoints)
      invpindex[pi] = -1;

    for(const auto& li : loclines)
      {
	invpindex[li.I1()] = 0;
	invpindex[li.I2()] = 0;
      }


    for(auto& line : loclines)
      {
        for(auto i : Range(2))
          {
            auto& pi = line[i];
	    if (invpindex[pi] == 0)
	      {
		pindex.Append (pi);
		invpindex[pi] = pindex.Size();
                locpoints.Append (points[pi].P());
		pi = locpoints.Size();
	      }
	    else
	      pi = invpindex[pi];
	  }
      }


    // double xh2 = xh*xh;
    for(auto i : nearpoints)
      {
	if (points[i].Valid() && 
	    points[i].OnSurface() &&
	    // Dist2 (points.Get(i).P(), p0) <= xh2 &&
	    invpindex[i] <= 0)
	  {
            locpoints.Append (points[i].P());
	    invpindex[i] = locpoints.Size();
	    pindex.Append(i);
	  }
      }
    /*
    double xh2 = xh*xh;
    for (i = 1; i <= points.Size(); i++)
      {
	if (points.Get(i).Valid() && 
	    points.Get(i).OnSurface() &&
	    Dist2 (points.Get(i).P(), p0) <= xh2 &&
	    invpindex.Get(i) <= 0)
	  {
	    invpindex.Elem(i) =
	      locpoints.Append (points.Get(i).P());
	    pindex.Append(i);
	  }
      }
    */

    pgeominfo.SetSize (locpoints.Size());
    for (int i = 0; i < pgeominfo.Size(); i++)
      pgeominfo[i].Init();


    for (int i = 0; i < loclines.Size(); i++)
      for (int j = 0; j < 2; j++)
	{
	  int lpi = loclines[i][j];
	
	  const PointGeomInfo & gi = 
	    lines[lindex[i]].GetGeomInfo (j+1);
	  pgeominfo.Elem(lpi).AddPointGeomInfo (gi);
	
	  /*
	    if (pgeominfo.Elem(lpi).cnt == MULTIPOINTGEOMINFO_MAX)
	    break;

	    const PointGeomInfo & gi = 
	    lines.Get(lindex.Get(i)).GetGeomInfo (j);
	
	    PointGeomInfo * pgi = pgeominfo.Elem(lpi).mgi;

	    int found = 0;
	    for (k = 0; k < pgeominfo.Elem(lpi).cnt; k++)
	    if (pgi[k].trignum == gi.trignum)
	    found = 1;

	    if (!found)
	    {
	    pgi[pgeominfo.Elem(lpi).cnt] = gi;
	    pgeominfo.Elem(lpi).cnt++;
	    }
	  */
	}

    for (int i = 0; i < locpoints.Size(); i++)
      {
	int pi = pindex[i];
      
	if (points[pi].mgi)
	  for (int j = 1; j <= points[pi].mgi->GetNPGI(); j++)
	    pgeominfo[i].AddPointGeomInfo (points[pi].mgi->GetPGI(j));
      }
   
    if (loclines.Size() == 1)
      {
	cout << "loclines.Size = 1" << endl;
	(*testout) << "loclines.size = 1" << endl
		   << " h = " << xh << endl
		   << " nearline.size = " << nearlines.Size() << endl
		   << " p0 = " << p0 << endl;
      }

    return lines[baselineindex].LineClass();
  }



  void AdFront2 :: SetStartFront ()
  {
    // for (int i = lines.Begin(); i < lines.End(); i++)
    for (int i : lines.Range())
      if (lines[i].Valid())
	for (int j = 1; j <= 2; j++)
	  points[lines[i].L().I(j)].DecFrontNr(0);
  }


  void AdFront2 :: Print (ostream & ost) const
  {
    ost << points.Size() << " Points: " << endl;
    // for (int i = points.Begin(); i < points.End(); i++)
    for (int i : points.Range())
      if (points[i].Valid())
	ost << i << "  " << points[i].P() << endl;

    ost << nfl << " Lines: " << endl;
    // for (int i = lines.Begin(); i < lines.End(); i++)
    for (int i : lines.Range())
      if (lines[i].Valid())
	ost << lines[i].L().I1() << " - " << lines[i].L().I2() << endl;

    ost << flush;
  }


  bool AdFront2 :: Inside (const Point<2> & p) const
  {
    int cnt;
    Vec<2> n;
    Vec<3> v1;
    DenseMatrix a(2), ainv(2);
    Vector b(2), u(2);
    
    // quasi-random numbers:
    n(0) = 0.123871;
    n(1) = 0.15432;
    
    cnt = 0;
    for (int i = 0; i < lines.Size(); i++)
      if (lines[i].Valid())
	{
	  const Point<3> & p1 = points[lines[i].L().I1()].P();
	  const Point<3> & p2 = points[lines[i].L().I2()].P();
	  
	  v1 = p2 - p1;
	  
	  a(0, 0) = v1(0);
	  a(1, 0) = v1(1);
	  
	  a(0, 1) = -n(0);
	  a(1, 1) = -n(1);

	  b(0) = p(0) - p1(0);
	  b(1) = p(1) - p1(1);
	  
	  CalcInverse (a, ainv);
	  ainv.Mult (b, u);
	  
	  if (u(0) >= 0 && u(0) <= 1 && u(1) > 0)
	    cnt++;
	}
    
    return ((cnt % 2) != 0);
  }

  bool AdFront2 :: SameSide (const Point<2> & lp1, const Point<2> & lp2, 
                             const NgArray<int> * testfaces) const
  {
    int cnt = 0;

    if (testfaces)
      {
        for (int ii = 0; ii < testfaces->Size(); ii++)
          if (lines[(*testfaces)[ii]].Valid())
            {
              int i = (*testfaces)[ii];
              const Point<3> & p13d = points[lines[i].L().I1()].P();
              const Point<3> & p23d = points[lines[i].L().I2()].P();
              
              Point<2> p1(p13d(0), p13d(1));
              Point<2> p2(p23d(0), p23d(1));
              
              // p1 + alpha v = lp1 + beta vl
              Vec<2> v = p2-p1;
              Vec<2> vl = lp2 - lp1;
              Mat<2,2> mat, inv;
              Vec<2> rhs, sol;
              mat(0,0) = v(0);
              mat(1,0) = v(1);
              mat(0,1) = -vl(0);
              mat(1,1) = -vl(1);
              rhs = lp1-p1;
              
              if (Det(mat) == 0) continue;
              CalcInverse (mat, inv);
              sol = inv * rhs;
              if ( (sol(0) >= 0) && (sol(0) <= 1) && (sol(1) >= 0) && (sol(1) <= 1))
                { cnt++; }
            }

      }
    else
      {
        for (int i = 0; i < lines.Size(); i++)
          if (lines[i].Valid())
            {
              const Point<3> & p13d = points[lines[i].L().I1()].P();
              const Point<3> & p23d = points[lines[i].L().I2()].P();
              
              Point<2> p1(p13d(0), p13d(1));
              Point<2> p2(p23d(0), p23d(1));
              
              // p1 + alpha v = lp1 + beta vl
              Vec<2> v = p2-p1;
              Vec<2> vl = lp2 - lp1;
              Mat<2,2> mat, inv;
              Vec<2> rhs, sol;
              mat(0,0) = v(0);
              mat(1,0) = v(1);
              mat(0,1) = -vl(0);
              mat(1,1) = -vl(1);
              rhs = lp1-p1;
              
              if (Det(mat) == 0) continue;
              CalcInverse (mat, inv);
              sol = inv * rhs;
              if ((sol(0) >= 0) && (sol(0) <= 1) && (sol(1) >= 0) && (sol(1) <= 1))
                { cnt++; }
            }
      }
    return ((cnt % 2) == 0);
  }
}

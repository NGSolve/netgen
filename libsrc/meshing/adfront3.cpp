#include <mystdlib.h>

#include <gprim/geomtest3d.hpp>
#include "adfront3.hpp"

/* ********************** FrontPoint ********************** */

namespace netgen
{

FrontPoint3 :: FrontPoint3 () 
{ 
  globalindex.Invalidate(); //  = -1;
  nfacetopoint = 0; 
  frontnr = 1000; 
  cluster = Front3PointIndex::INVALID;
}


FrontPoint3 :: FrontPoint3 (const Point<3> & ap, PointIndex agi)
{ 
  p = ap; 
  globalindex = agi;
  nfacetopoint = 0; 
  frontnr = 1000; 
  cluster = Front3PointIndex::INVALID;
}



/* ********************** FrontFace ********************** */

FrontFace :: FrontFace () 
{ 
  qualclass = 1; 
  oldfront = 0; 
  hashvalue = 0;
  cluster = Front3PointIndex::INVALID;
}

FrontFace :: FrontFace (const FrontElement2d & af)
{ 
  f = af; 
  oldfront = 0; 
  qualclass = 1; 
  hashvalue = 0;
}

void FrontFace :: Invalidate ()
{ 
  f.Delete(); 
  oldfront = 0; 
  qualclass = 1000; 
}




/* ********************** AddFront ********************** */
 

AdFront3 :: AdFront3 ()
{
  nff = 0;
  nff4 = 0;
  vol = 0;

  hashon = 1;
  hashcreated = 0;
  if (hashon) 
    hashtable.Init(&points, &faces);

  facetree = NULL;
  // connectedpairs = NULL;

  rebuildcounter = -1;
  lasti = 0;
  minval = -1;
}


AdFront3 :: ~AdFront3 ()
{
  delete facetree;
  // delete connectedpairs;
}

void AdFront3 :: GetPoints (Array<Point<3> > & apoints) const
{
  /*
  for (Front3PointIndex pi = points.Begin(); pi < points.End(); pi++)
    
    apoints.Append (points[pi].P());
  */
  for (auto & p : points)
    apoints.Append(p.P());
}


Front3PointIndex AdFront3 :: AddPoint (const Point<3> & p, PointIndex globind)
{
  if (delpointl.Size())
    {
      Front3PointIndex pi = delpointl.Last();
      delpointl.DeleteLast ();
      
      points[pi] = FrontPoint3 (p, globind);
      return pi;
    }
  else
    {
      points.Append (FrontPoint3 (p, globind));
      // return --points.End();
      return *points.Range().end()-1;
      // return points.Size()-1+PointIndex::BASE;
    }
}


INDEX AdFront3 :: AddFace (const FrontElement2d & aface)
{
  int i, minfn;

  nff++;

  for (i = 0; i < aface.GetNP(); i++)
    points[aface[i]].AddFace();

  const Point<3> & p1 = points[aface[0]].P();
  const Point<3> & p2 = points[aface[1]].P();
  const Point<3> & p3 = points[aface[2]].P();

  vol += 1.0/6.0 * (p1(0) + p2(0) + p3(0)) *
    ( (p2(1) - p1(1)) * (p3(2) - p1(2)) -
      (p2(2) - p1(2)) * (p3(1) - p1(1)) );

  if (aface.GetNP() == 4)
    {
      nff4++;
      const Point<3> & p4 = points[aface[3]].P();      
      vol += 1.0/6.0 * (p1(0) + p3(0) + p4(0)) *
	( (p3(1) - p1(1)) * (p4(2) - p1(2)) -
	  (p3(2) - p1(2)) * (p4(1) - p1(1)) );
    }


  minfn = 1000;
  for (i = 0; i < aface.GetNP(); i++)
    {
      int fpn = points[aface[i]].FrontNr();
      if (i == 0 || fpn < minfn)
	minfn = fpn;
    }


  Front3PointIndex cluster = Front3PointIndex::INVALID;
  for (i = 1; i <= aface.GetNP(); i++)
    {
      if (points[aface.PNum(i)].cluster.IsValid())
	cluster = points[aface.PNum(i)].cluster;
    }
  for (i = 1; i <= aface.GetNP(); i++)
    points[aface.PNum(i)].cluster = cluster;


  for (i = 1; i <= aface.GetNP(); i++)
    points[aface.PNum(i)].DecFrontNr (minfn+1);
  
  faces.Append(FrontFace (aface));
  int nfn = faces.Size();
  faces[nfn-1].cluster = cluster;

  if (hashon && hashcreated) 
    hashtable.AddElem(aface, nfn);

  return nfn;
}



void AdFront3 :: DeleteFace (INDEX fi)
{
  nff--;

  /*
  for (int i = 1; i <= faces.Get(fi).Face().GetNP(); i++)
    {
      Front3PointIndex pi = faces.Get(fi).Face().PNum(i);
  */
  for (Front3PointIndex pi : faces[fi-1].Face().PNums())
    {
      points[pi].RemoveFace();
      if (!points[pi].Valid())
	delpointl.Append (pi);
    }

  const FrontElement2d & face = faces[fi-1].Face();
  const Point<3> & p1 = points[face.PNum(1)].P();
  const Point<3> & p2 = points[face.PNum(2)].P();
  const Point<3> & p3 = points[face.PNum(3)].P();

  vol -= 1.0/6.0 * (p1(0) + p2(0) + p3(0)) *
    ( (p2(1) - p1(1)) * (p3(2) - p1(2)) -
      (p2(2) - p1(2)) * (p3(1) - p1(1)) );

  if (face.GetNP() == 4)
    {
      const Point<3> & p4 = points[face.PNum(4)].P();      
      vol -= 1.0/6.0 * (p1(0) + p3(0) + p4(0)) *
	( (p3(1) - p1(1)) * (p4(2) - p1(2)) -
	  (p3(2) - p1(2)) * (p4(1) - p1(1)) );

      nff4--;
    }

  faces[fi-1].Invalidate();
}


INDEX AdFront3 :: AddConnectedPair (IVec<2,Front3PointIndex> apair)
{
  if (!connectedpairs)
    connectedpairs = make_unique<DynamicTable<Front3PointIndex, Front3PointIndex>> (GetNP());

  connectedpairs->Add (apair[0], apair[1]);
  connectedpairs->Add (apair[1], apair[0]);

  return 0;
}


void AdFront3 :: CreateTrees ()
{
  int i, j;
  Front3PointIndex pi;
  Point<3> pmin, pmax;

  for (pi = IndexBASE<Front3PointIndex>(); 
       pi < GetNP()+IndexBASE<Front3PointIndex>(); pi++)
    {
      const Point<3> & p = GetPoint(pi);
      if (pi == IndexBASE<Front3PointIndex>())
	{
	  pmin = p;
	  pmax = p;
	}
      else
	{
	  SetToMin (pmin, p);
	  SetToMax (pmax, p);
	}
    }

  pmax = pmax + 0.5 * (pmax - pmin);
  pmin = pmin + 0.5 * (pmin - pmax);

  delete facetree;
  facetree = new BoxTree<3> (pmin, pmax);
  
  for (i = 1; i <= GetNF(); i++)
    {
      const FrontElement2d & el = GetFace(i);
      pmin = GetPoint (el[0]);
      pmax = pmin;
      for (j = 1; j < 3; j++)
	{
	  const Point<3> & p = GetPoint (el[j]);
	  SetToMin (pmin, p);
	  SetToMax (pmax, p);
	}
      pmax = pmax + 0.01 * (pmax - pmin);
      pmin = pmin + 0.01 * (pmin - pmax);
      //      (*testout) << "insert " << i << ": " << pmin << " - " << pmax << "\n";
      facetree -> Insert (pmin, pmax, i);
    }
}


void AdFront3 :: GetIntersectingFaces (const Point<3> & pmin, const Point<3> & pmax, 
				       Array<int> & ifaces) const
{
  facetree -> GetIntersecting (pmin, pmax, ifaces);
}

void AdFront3 :: GetFaceBoundingBox (int i, Box3d & box) const
{
  const FrontFace & face = faces[i-1];
  box.SetPoint (points[face.f[0]].p);
  box.AddPoint (points[face.f[1]].p);
  box.AddPoint (points[face.f[2]].p);
}

void AdFront3 :: RebuildInternalTables ()
{
  static Timer timer_a("Adfront3::RebuildInternal A");
  static Timer timer_b("Adfront3::RebuildInternal B");
  static Timer timer_c("Adfront3::RebuildInternal C");
  static Timer timer_d("Adfront3::RebuildInternal D");


  timer_a.Start();	  
  int hi = 0;
  for (int i = 1; i <= faces.Size(); i++)
    if (faces[i-1].Valid())
      {
	hi++;
	if (hi < i)
	  faces[hi-1] = faces[i-1];
      }
  
  faces.SetSize (nff);

  int np = points.Size();

  // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
  for (Front3PointIndex pi : points.Range())
    points[pi].cluster = pi;
  
  timer_a.Stop();	  
  timer_b.Start();

  int change;
  do
    {
      change = 0;
      for (int i = 1; i <= faces.Size(); i++)
	{
	  const FrontElement2d & el = faces[i-1].Face();

	  Front3PointIndex mini = points[el.PNum(1)].cluster;
	  Front3PointIndex maxi = mini;
	  
	  for (int j = 2; j <= 3; j++)
	    {
	      Front3PointIndex ci = points[el.PNum(j)].cluster;
	      if (ci < mini) mini = ci;
	      if (ci > maxi) maxi = ci;
	    }

	  if (mini < maxi)
	    {
	      change = 1;
	      for (int j = 1; j <= 3; j++)
		points[el.PNum(j)].cluster = mini;
	    }
	}
    }
  while (change);


  timer_b.Stop();
  timer_c.Start();




  Array<bool, Front3PointIndex> usecl(np);
  usecl = false;
  for (int i = 1; i <= faces.Size(); i++)
    {
      usecl[points[faces[i-1].Face().PNum(1)].cluster] = true;
      faces[i-1].cluster =
	points[faces[i-1].Face().PNum(1)].cluster;
    }
  /*
  int cntcl = 0;
  for (int i = PointIndex::BASE; 
       i < np+PointIndex::BASE; i++)
    if (usecl[i])
      cntcl++;
  */
  
  Array<double, Front3PointIndex> clvol (np);
  clvol = 0.0;

  for (int i = 1; i <= faces.Size(); i++)
    {
      const FrontElement2d & face = faces[i-1].Face();

      const Point<3> p1 = points[face.PNum(1)].P();      
      const Point<3> p2 = points[face.PNum(2)].P();      
      const Point<3> p3 = points[face.PNum(3)].P();      
      
      double vi = 1.0/6.0 * (p1(0) + p2(0) + p3(0)) *
	( (p2(1) - p1(1)) * (p3(2) - p1(2)) -
	  (p2(2) - p1(2)) * (p3(1) - p1(1)) );
      
      if (face.GetNP() == 4)
	{
	  const Point<3> p4 = points[face.PNum(4)].P();      
	  vi += 1.0/6.0 * (p1(0) + p3(0) + p4(0)) *
	    ( (p3(1) - p1(1)) * (p4(2) - p1(2)) -
	      (p3(2) - p1(2)) * (p4(1) - p1(1)) );
	}
     
      clvol[faces[i-1].cluster] += vi;
    }

  timer_c.Stop();	  
  timer_d.Start();



  bool negvol = false;
  for (auto i : clvol.Range())
    if (clvol[i] < 0)
      negvol = true;
  
  if (negvol)
    {
      for (int i = 1; i <= faces.Size(); i++)
	faces[i-1].cluster = IndexBASE<Front3PointIndex>();
      for (Front3PointIndex pi : points.Range())
	points[pi].cluster = IndexBASE<Front3PointIndex>();
    }

  if (hashon) 
    hashtable.Create();

  timer_d.Stop();
}



int AdFront3 :: SelectBaseElement ()
{
  /*
  static int minval = -1;
  static int lasti = 0;
  static int counter = 0;
  */
  if (rebuildcounter <= 0)
    {
      RebuildInternalTables();
      rebuildcounter = nff / 10 + 1;
      
      lasti = 0;
    }
  rebuildcounter--;

  /*
  if (faces.Size() > 2 * nff)
    {
      // compress facelist

      RebuildInternalTables ();
      lasti = 0;
    }
    */
  
  int fstind = 0;

  for (int i = lasti+1; i <= faces.Size() && !fstind; i++)
    if (faces[i-1].Valid())
      {
	int hi = faces[i-1].QualClass() +
	  points[faces[i-1].Face().PNum(1)].FrontNr() +
	  points[faces[i-1].Face().PNum(2)].FrontNr() +
	  points[faces[i-1].Face().PNum(3)].FrontNr();
	
	if (hi <= minval)
	  {
	    minval = hi;
	    fstind = i;
	    lasti = fstind;
	  }
      }
  
  if (!fstind)
    {
      minval = INT_MAX;
      for (int i = 1; i <= faces.Size(); i++)
	if (faces[i-1].Valid())
	  {
	    int hi = faces[i-1].QualClass() +
	      points[faces[i-1].Face().PNum(1)].FrontNr() +
	      points[faces[i-1].Face().PNum(2)].FrontNr() +
	      points[faces[i-1].Face().PNum(3)].FrontNr();
	    
	    if (hi <= minval)
	      {
		minval = hi;
		fstind = i;
		lasti = 0;
	      }
	  }
    }


  return fstind;
}



int AdFront3 :: GetLocals (int fstind,
			   Array<Point<3>, LocalPointIndex> & locpoints,
			   Array<MiniElement2d> & locfaces,   // local index
			   Array<Front3PointIndex, LocalPointIndex> & pindex,
			   Array<INDEX> & findex,
			   INDEX_2_HASHTABLE<int> & getconnectedpairs,
			   float xh,
			   float relh,
			   INDEX& facesplit)
{
  // static Timer timer("AdFront3::GetLocals");
  // RegionTimer reg (timer);


  if (hashon && faces.Size() < 500) { hashon=0; }
  if (hashon && !hashcreated) 
    {
      hashtable.Create(); 
      hashcreated=1;
    }

  INDEX i;
  Front3PointIndex pstind;
  Point<3> midp, p0;

  //  static Array<int, PointIndex::BASE> invpindex;
  
  Array<FrontElement2d> locfaces2;          // all front faces in radius xh
  Array<int> locfaces3;           // all faces in outer radius relh
  Array<INDEX> findex2;

  locfaces2.SetSize(0);
  locfaces3.SetSize(0);
  findex2.SetSize(0);

  Front3PointIndex cluster = faces[fstind-1].cluster;

  pstind = faces[fstind-1].Face().PNum(1);
  p0 = points[pstind].P();
  
  locfaces2.Append(faces[fstind-1].Face());
  findex2.Append(fstind);


  Box3d b1 (p0 - Vec<3>(xh, xh, xh), p0 + Vec<3> (xh, xh, xh));

  if (hashon)
    {
      hashtable.GetLocals(locfaces2, findex2, fstind, p0, xh);
    }
  else
    {
      for (i = 1; i <= faces.Size(); i++)
	{
	  const FrontElement2d & face = faces[i-1].Face();
	  if (faces[i-1].cluster == cluster && faces[i-1].Valid() && i != fstind)
	    {
	      Box3d b2;
	      b2.SetPoint (points[face[0]].P());
	      b2.AddPoint (points[face[1]].P());
	      b2.AddPoint (points[face[2]].P());

	      if (b1.Intersect (b2))
		{
		  locfaces2.Append(faces[i-1].Face());
		  findex2.Append(i);
		}
	    }
	}
    }

  Array<FrontElement2d> frontfaces;         // the selected faces, front numbering

  //local faces for inner radius:
  for (i = 1; i <= locfaces2.Size(); i++)
    {
      const FrontElement2d & face = locfaces2[i-1];
      const Point<3> & p1 = points[face[0]].P();
      const Point<3> & p2 = points[face[1]].P();
      const Point<3> & p3 = points[face[2]].P();

      midp = Center (p1, p2, p3);

      if (Dist2 (midp, p0) <= relh * relh || i == 1)
	{
          frontfaces.Append(locfaces2[i-1]);
	  findex.Append(findex2[i-1]);
	}
      else
	locfaces3.Append (i);
    }
  
  facesplit=frontfaces.Size();
  
  
  //local faces for outer radius:
  for (i = 1; i <= locfaces3.Size(); i++)
    {
      frontfaces.Append (locfaces2[locfaces3[i-1]-1]);
      findex.Append (findex2[locfaces3[i-1]-1]);
    }


  invpindex.SetSize (points.Size());
  /*
  for (i = 1; i <= locfaces.Size(); i++)
    for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
      {
	PointIndex pi = locfaces.Get(i).PNum(j);
        invpindex[pi] = PointIndex::INVALID;
      }
  */
  for (auto & f : frontfaces)
    for (int j = 1; j <= f.GetNP(); j++)
      invpindex[f.PNum(j)] = LocalPointIndex::INVALID;

  for (auto & f : frontfaces)
    {
      MiniElement2d locface(f.GetNP());
      for (int j = 1; j <= f.GetNP(); j++)
	{
          Front3PointIndex pi = f.PNum(j);
	  if (!invpindex[pi].IsValid())
	    {
	      pindex.Append (pi);
              locpoints.Append (points[pi].P());
	      invpindex[pi] = pindex.Size()-1+IndexBASE<LocalPointIndex>();
            }
          locface.PNum(j) = invpindex[pi];
	}
      locfaces.Append (locface);
    }



  if (connectedpairs)
    {
      // for (i = 1; i <= locpoints.Size(); i++)
      for (auto i : locpoints.Range())
	{
	  Front3PointIndex pind = pindex[i]; // .Get(i);
	  // if (pind.IsValid() && pind <= connectedpairs->Size ())
          if (connectedpairs->Range().Contains(pind))
	    {
	      // for (int j = 1; j <= connectedpairs->EntrySize(pind); j++)
              for (auto j : (*connectedpairs)[pind].Range())
		{
		  //PointIndex oi = connectedpairs->Get(pind, j);
                  Front3PointIndex oi = (*connectedpairs)[pind][j];
		  LocalPointIndex other = invpindex[oi];
		  // if (other >= 1 && other <= pindex.Size() &&
                  if (pindex.Range().Contains(other) &&
		      pindex[other] == oi)
		    {
		      // INDEX_2 coned(i, other);
		      // coned.Sort();
		      // (*testout) << "connected: " << locpoints.Get(i) << "-" << locpoints.Get(other) << endl;
		      getconnectedpairs.Set (INDEX_2::Sort (i.Nr0(),
							    other.Nr0()), 1);
		    }
		}
	    }
	}
    }
  

  /*
    // add isolated points
  for (i = 1; i <= points.Size(); i++)
    if (points.Elem(i).Valid() && Dist (points.Elem(i).P(), p0) <= xh)
      {
	if (!invpindex.Get(i))
	  {
	    locpoints.Append (points.Get(i).P());
	    pindex.Append (i);
	    invpindex.Elem(i) = pindex.Size();
	  }
      }
      */
  return faces[fstind-1].QualClass();
}


// returns all points connected with fi
void AdFront3 :: GetGroup (int fi,
			   Array<MeshPoint, LocalPointIndex> & grouppoints,
			   Array<MiniElement2d> & groupelements,
			   Array<Front3PointIndex, LocalPointIndex> & pindex,
			   Array<INDEX> & findex) 
{
  // static Array<char> pingroup;
  int changed;

  pingroup.SetSize(points.Size());

  pingroup = 0;
  for (int j = 1; j <= 3; j++)
    pingroup[faces[fi-1].Face().PNum(j)] = 1;

  do
    {
      changed = 0;

      /*
      for (i = 1; i <= faces.Size(); i++)
	if (faces.Get(i).Valid())
	  {
	    const MiniElement2d & face = faces.Get(i).Face();

	    int fused = 0;
	    for (j = 1; j <= 3; j++)
	      if (pingroup.Elem(face.PNum(j))) 
		fused++;
            
	    if (fused >= 2)
	      for (j = 1; j <= 3; j++)
		if (!pingroup.Elem(face.PNum(j)))
		  {
		    pingroup.Elem(face.PNum(j)) = 1;
		    changed = 1;
		  }
	  }
      */
      for (auto & f : faces)
	if (f.Valid())
	  {
	    const FrontElement2d & face = f.Face();

	    int fused = 0;
	    for (int j = 1; j <= 3; j++)
	      if (pingroup[face.PNum(j)]) 
		fused++;
            
	    if (fused >= 2)
	      for (int j = 1; j <= 3; j++)
		if (!pingroup[face.PNum(j)])
		  {
		    pingroup[face.PNum(j)] = 1;
		    changed = 1;
		  }
	  }

    }
  while (changed);

  invpindex.SetSize (points.Size());

  // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
  for (Front3PointIndex pi : points.Range())
    if (points[pi].Valid())
      {
	grouppoints.Append (points[pi].P());
        pindex.Append (pi);
	invpindex[pi] = pindex.Size()-1 + IndexBASE<LocalPointIndex>();
      }

  for (int i = 1; i <= faces.Size(); i++)
    if (faces[i-1].Valid())
      {
	int fused = 0;
	for (int j = 1; j <= 3; j++)
	  if (pingroup[faces[i-1].Face().PNum(j)])
	    fused++;

	if (fused >= 2)
	  {
	    const FrontElement2d & f = faces[i-1].Face();
	    MiniElement2d ge(f.GetNP());
	    for (int j = 1; j <= f.GetNP(); j++)
	      ge.PNum(j) = invpindex[f.PNum(j)];
	    groupelements.Append (ge);
	    findex.Append (i);
	  }
      }

}


void AdFront3 :: SetStartFront (int /* baseelnp */)
{
  for (INDEX i = 1; i <= faces.Size(); i++)
    if (faces[i-1].Valid())
      {
	const FrontElement2d & face = faces[i-1].Face();
	for (int j = 1; j <= 3; j++)
	  points[face.PNum(j)].DecFrontNr(0);
      }

  /*
  if (baseelnp)
    {
      for (i = 1; i <= faces.Size(); i++)
	if (faces.Get(i).Valid() && faces.Get(i).Face().GetNP() != baseelnp)
	  faces.Elem(i).qualclass = 1000;
    }
    */
}

bool AdFront3 :: PointInsideGroup(const Array<Front3PointIndex, LocalPointIndex> &grouppindex,
                                  const Array<MiniElement2d> &groupfaces) const
{
  for(auto pi : Range(points))
    {
      const auto& p = points[pi].P();
      bool found = false;
      for(const auto& f : groupfaces)
        {
          for(auto i : Range(3))
            if(grouppindex[f.PNum(i+1)] == pi)
              {
                found = true;
                break;
              }
        }
      if(found)
        continue;

      // "random" direction
      Vec<3> dir = { 0.123871, 0.15432,-0.43989 };
      DenseMatrix a(3), ainv(3);
      Vector b(3), u(3);

      int count = 0;
      for(const auto& f : groupfaces)
        {
          const auto& p1 = points[grouppindex[f.PNum(1)]].P();
          auto v1 = points[grouppindex[f.PNum(2)]].P() - p1;
          auto v2 = points[grouppindex[f.PNum(3)]].P() - p1;
          for(auto i : Range(3))
            {
              a(i,0) = v1[i];
              a(i,1) = v2[i];
              a(i,2) = -dir[i];
              b(i) = p[i] - p1[i];
            }
          CalcInverse (a, ainv);
          ainv.Mult (b, u);
          if (u(0) >= 0 && u(1) >= 0 && u(0)+u(1) <= 1 &&
	    u(2) > 0)
	    count++;
        }
        if (count % 2 == 1)
          return true;
    }
  return false;
}

bool AdFront3 :: Inside (const Point<3> & p) const
{
  static Timer timer("AdFront3::Inside"); RegionTimer rt(timer);
  int cnt;
  Vec<3> n, v1, v2;
  DenseMatrix a(3), ainv(3);
  Vector b(3), u(3);

  // random numbers:
  n(0) = 0.123871;
  n(1) = 0.15432;
  n(2) = -0.43989;

  cnt = 0;
  for (int i = 1; i <= faces.Size(); i++)
    if (faces[i-1].Valid())
      {
	const Point<3> & p1 = points[faces[i-1].Face().PNum(1)].P();
	const Point<3> & p2 = points[faces[i-1].Face().PNum(2)].P();
	const Point<3> & p3 = points[faces[i-1].Face().PNum(3)].P();

	v1 = p2 - p1;
	v2 = p3 - p1;

	a(0, 0) = v1(0);
	a(1, 0) = v1(1);
	a(2, 0) = v1(2);
	a(0, 1) = v2(0);
	a(1, 1) = v2(1);
	a(2, 1) = v2(2);
	a(0, 2) = -n(0);
	a(1, 2) = -n(1);
	a(2, 2) = -n(2);

	b(0) = p(0) - p1(0);
	b(1) = p(1) - p1(1);
	b(2) = p(2) - p1(2);

	CalcInverse (a, ainv);
	ainv.Mult (b, u);

	if (u(0) >= 0 && u(1) >= 0 && u(0)+u(1) <= 1 &&
	    u(2) > 0)
	  {
	    cnt++;
	  }
      }

  return ((cnt % 2) != 0);
}





int AdFront3 :: SameSide (const Point<3> & lp1, const Point<3> & lp2,
			  const Array<int> * testfaces) const
{
  const Point<3> *line[2];
  line[0] = &lp1;
  line[1] = &lp2;


  Point<3> pmin(lp1);
  Point<3> pmax(lp1);
  SetToMin (pmin, lp2);
  SetToMax (pmax, lp2);
  
  ArrayMem<int, 100> aprif;
  aprif.SetSize(0);
  
  if (!testfaces)
    facetree->GetIntersecting (pmin, pmax, aprif);
  else
    for (int i = 1; i <= testfaces->Size(); i++)
      aprif.Append (testfaces->operator[](i-1));

  int cnt = 0;
  for (int ii = 1; ii <= aprif.Size(); ii++)
    {
      int i = aprif[ii-1];
      
      if (faces[i-1].Valid())
	{
	  const Point<3> *tri[3];
	  tri[0] = &points[faces[i-1].Face().PNum(1)].P();
	  tri[1] = &points[faces[i-1].Face().PNum(2)].P();
	  tri[2] = &points[faces[i-1].Face().PNum(3)].P();
	  	  
	  if (IntersectTriangleLine (&tri[0], &line[0]))
	    cnt++;
	}
    }

  return ((cnt+1) % 2);
}
}

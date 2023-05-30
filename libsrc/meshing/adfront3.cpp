#include <mystdlib.h>
#include "meshing.hpp"


/* ********************** FrontPoint ********************** */

namespace netgen
{

FrontPoint3 :: FrontPoint3 () 
{ 
  globalindex.Invalidate(); //  = -1;
  nfacetopoint = 0; 
  frontnr = 1000; 
  cluster = 0;
}


FrontPoint3 :: FrontPoint3 (const Point<3> & ap, PointIndex agi)
{ 
  p = ap; 
  globalindex = agi;
  nfacetopoint = 0; 
  frontnr = 1000; 
  cluster = 0;
}



/* ********************** FrontFace ********************** */

FrontFace :: FrontFace () 
{ 
  qualclass = 1; 
  oldfront = 0; 
  hashvalue = 0;
  cluster = 0;
}

FrontFace :: FrontFace (const MiniElement2d & af)
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
  connectedpairs = NULL;

  rebuildcounter = -1;
  lasti = 0;
  minval = -1;
}


AdFront3 :: ~AdFront3 ()
{
  delete facetree;
  delete connectedpairs;
}

void AdFront3 :: GetPoints (NgArray<Point<3> > & apoints) const
{
  /*
  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
    
    apoints.Append (points[pi].P());
  */
  for (auto & p : points)
    apoints.Append(p.P());
}


PointIndex AdFront3 :: AddPoint (const Point<3> & p, PointIndex globind)
{
  if (delpointl.Size())
    {
      PointIndex pi = delpointl.Last();
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


INDEX AdFront3 :: AddFace (const MiniElement2d & aface)
{
  int i, minfn;

  nff++;

  for (i = 0; i < aface.GetNP(); i++)
    points[aface[i]].AddFace();

  const Point3d & p1 = points[aface[0]].P();
  const Point3d & p2 = points[aface[1]].P();
  const Point3d & p3 = points[aface[2]].P();

  vol += 1.0/6.0 * (p1.X() + p2.X() + p3.X()) *
    ( (p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) -
      (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()) );

  if (aface.GetNP() == 4)
    {
      nff4++;
      const Point3d & p4 = points[aface[3]].P();      
      vol += 1.0/6.0 * (p1.X() + p3.X() + p4.X()) *
	( (p3.Y() - p1.Y()) * (p4.Z() - p1.Z()) -
	  (p3.Z() - p1.Z()) * (p4.Y() - p1.Y()) );
    }


  minfn = 1000;
  for (i = 0; i < aface.GetNP(); i++)
    {
      int fpn = points[aface[i]].FrontNr();
      if (i == 0 || fpn < minfn)
	minfn = fpn;
    }


  int cluster = 0;
  for (i = 1; i <= aface.GetNP(); i++)
    {
      if (points[aface.PNum(i)].cluster)
	cluster = points[aface.PNum(i)].cluster;
    }
  for (i = 1; i <= aface.GetNP(); i++)
    points[aface.PNum(i)].cluster = cluster;


  for (i = 1; i <= aface.GetNP(); i++)
    points[aface.PNum(i)].DecFrontNr (minfn+1);
  
  faces.Append(FrontFace (aface));
  int nfn = faces.Size();
  faces.Elem(nfn).cluster = cluster;

  if (hashon && hashcreated) 
    hashtable.AddElem(aface, nfn);

  return nfn;
}



void AdFront3 :: DeleteFace (INDEX fi)
{
  nff--;

  for (int i = 1; i <= faces.Get(fi).Face().GetNP(); i++)
    {
      PointIndex pi = faces.Get(fi).Face().PNum(i);
      points[pi].RemoveFace();
      if (!points[pi].Valid())
	delpointl.Append (pi);
    }

  const MiniElement2d & face = faces.Get(fi).Face();
  const Point3d & p1 = points[face.PNum(1)].P();
  const Point3d & p2 = points[face.PNum(2)].P();
  const Point3d & p3 = points[face.PNum(3)].P();

  vol -= 1.0/6.0 * (p1.X() + p2.X() + p3.X()) *
    ( (p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) -
      (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()) );

  if (face.GetNP() == 4)
    {
      const Point3d & p4 = points[face.PNum(4)].P();      
      vol -= 1.0/6.0 * (p1.X() + p3.X() + p4.X()) *
	( (p3.Y() - p1.Y()) * (p4.Z() - p1.Z()) -
	  (p3.Z() - p1.Z()) * (p4.Y() - p1.Y()) );

      nff4--;
    }

  faces.Elem(fi).Invalidate();
}


INDEX AdFront3 :: AddConnectedPair (const INDEX_2 & apair)
{
  if (!connectedpairs)
    connectedpairs = new TABLE<int, PointIndex::BASE> (GetNP());

  connectedpairs->Add (apair.I1(), apair.I2());
  connectedpairs->Add (apair.I2(), apair.I1());

  return 0;
}


void AdFront3 :: CreateTrees ()
{
  int i, j;
  PointIndex pi;
  Point3d pmin, pmax;

  for (pi = PointIndex::BASE; 
       pi < GetNP()+PointIndex::BASE; pi++)
    {
      const Point<3> & p = GetPoint(pi);
      if (pi == PointIndex::BASE)
	{
	  pmin = p;
	  pmax = p;
	}
      else
	{
	  pmin.SetToMin (p);
	  pmax.SetToMax (p);
	}
    }

  pmax = pmax + 0.5 * (pmax - pmin);
  pmin = pmin + 0.5 * (pmin - pmax);

  delete facetree;
  facetree = new BoxTree<3> (pmin, pmax);
  
  for (i = 1; i <= GetNF(); i++)
    {
      const MiniElement2d & el = GetFace(i);
      pmin = GetPoint (el[0]);
      pmax = pmin;
      for (j = 1; j < 3; j++)
	{
	  const Point<3> & p = GetPoint (el[j]);
	  pmin.SetToMin (p);
	  pmax.SetToMax (p);
	}
      pmax = pmax + 0.01 * (pmax - pmin);
      pmin = pmin + 0.01 * (pmin - pmax);
      //      (*testout) << "insert " << i << ": " << pmin << " - " << pmax << "\n";
      facetree -> Insert (pmin, pmax, i);
    }
}


void AdFront3 :: GetIntersectingFaces (const Point<3> & pmin, const Point<3> & pmax, 
				       NgArray<int> & ifaces) const
{
  facetree -> GetIntersecting (pmin, pmax, ifaces);
}

void AdFront3 :: GetFaceBoundingBox (int i, Box3d & box) const
{
  const FrontFace & face = faces.Get(i);
  box.SetPoint (points[face.f[0]].p);
  box.AddPoint (points[face.f[1]].p);
  box.AddPoint (points[face.f[2]].p);
}

void AdFront3 :: RebuildInternalTables ()
{
  static int timer_a = NgProfiler::CreateTimer ("Adfront3::RebuildInternal A");
  static int timer_b = NgProfiler::CreateTimer ("Adfront3::RebuildInternal B");
  static int timer_c = NgProfiler::CreateTimer ("Adfront3::RebuildInternal C");
  static int timer_d = NgProfiler::CreateTimer ("Adfront3::RebuildInternal D");


  NgProfiler::StartTimer (timer_a);	  
  int hi = 0;
  for (int i = 1; i <= faces.Size(); i++)
    if (faces.Get(i).Valid())
      {
	hi++;
	if (hi < i)
	  faces.Elem(hi) = faces.Get(i);
      }
  
  faces.SetSize (nff);

  int np = points.Size();

  // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
  for (PointIndex pi : points.Range())
    points[pi].cluster = pi;
  
  NgProfiler::StopTimer (timer_a);	  
  NgProfiler::StartTimer (timer_b);	  

  int change;
  do
    {
      change = 0;
      for (int i = 1; i <= faces.Size(); i++)
	{
	  const MiniElement2d & el = faces.Get(i).Face();

	  int mini = points[el.PNum(1)].cluster;
	  int maxi = mini;
	  
	  for (int j = 2; j <= 3; j++)
	    {
	      int ci = points[el.PNum(j)].cluster;
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


  NgProfiler::StopTimer (timer_b);	  
  NgProfiler::StartTimer (timer_c);	  




  Array<bool, PointIndex> usecl(np);
  usecl = false;
  for (int i = 1; i <= faces.Size(); i++)
    {
      usecl[points[faces.Get(i).Face().PNum(1)].cluster] = true;
      faces.Elem(i).cluster =
	points[faces.Get(i).Face().PNum(1)].cluster;
    }
  int cntcl = 0;
  for (int i = PointIndex::BASE; 
       i < np+PointIndex::BASE; i++)
    if (usecl[i])
      cntcl++;

  NgArray<double, PointIndex::BASE> clvol (np);
  clvol = 0.0;

  for (int i = 1; i <= faces.Size(); i++)
    {
      const MiniElement2d & face = faces.Get(i).Face();

      const Point3d p1 = points[face.PNum(1)].P();      
      const Point3d p2 = points[face.PNum(2)].P();      
      const Point3d p3 = points[face.PNum(3)].P();      
      
      double vi = 1.0/6.0 * (p1.X() + p2.X() + p3.X()) *
	( (p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) -
	  (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()) );
      
      if (face.GetNP() == 4)
	{
	  const Point3d p4 = points[face.PNum(4)].P();      
	  vi += 1.0/6.0 * (p1.X() + p3.X() + p4.X()) *
	    ( (p3.Y() - p1.Y()) * (p4.Z() - p1.Z()) -
	      (p3.Z() - p1.Z()) * (p4.Y() - p1.Y()) );
	}
     
      clvol[faces.Get(i).cluster] += vi;
    }

  NgProfiler::StopTimer (timer_c);	  
  NgProfiler::StartTimer (timer_d);	  



  int negvol = 0;
  for (int i = PointIndex::BASE; 
       i < clvol.Size()+PointIndex::BASE; i++)
    {
      if (clvol[i] < 0)
	negvol = 1;
    }
  
  if (negvol)
    {
      for (int i = 1; i <= faces.Size(); i++)
	faces.Elem(i).cluster = 1;
      // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
      for (PointIndex pi : points.Range())
	points[pi].cluster = 1;
    }

  if (hashon) 
    hashtable.Create();

  NgProfiler::StopTimer (timer_d);	  
}



int AdFront3 :: SelectBaseElement ()
{
  int i, hi, fstind;

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
  
  fstind = 0;

  for (i = lasti+1; i <= faces.Size() && !fstind; i++)
    if (faces.Elem(i).Valid())
      {
	hi = faces.Get(i).QualClass() +
	  points[faces.Get(i).Face().PNum(1)].FrontNr() +
	  points[faces.Get(i).Face().PNum(2)].FrontNr() +
	  points[faces.Get(i).Face().PNum(3)].FrontNr();
	
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
      for (i = 1; i <= faces.Size(); i++)
	if (faces.Elem(i).Valid())
	  {
	    hi = faces.Get(i).QualClass() +
	      points[faces.Get(i).Face().PNum(1)].FrontNr() +
	      points[faces.Get(i).Face().PNum(2)].FrontNr() +
	      points[faces.Get(i).Face().PNum(3)].FrontNr();
	    
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
			   NgArray<Point3d, PointIndex::BASE> & locpoints,
			   NgArray<MiniElement2d> & locfaces,   // local index
			   NgArray<PointIndex, PointIndex::BASE> & pindex,
			   NgArray<INDEX> & findex,
			   INDEX_2_HASHTABLE<int> & getconnectedpairs,
			   float xh,
			   float relh,
			   INDEX& facesplit)
{
  // static int timer = NgProfiler::CreateTimer ("AdFront3::GetLocals");
  // NgProfiler::RegionTimer reg (timer);


  if (hashon && faces.Size() < 500) { hashon=0; }
  if (hashon && !hashcreated) 
    {
      hashtable.Create(); 
      hashcreated=1;
    }

  INDEX i, j;
  PointIndex pstind;
  Point3d midp, p0;

  //  static NgArray<int, PointIndex::BASE> invpindex;
  
  NgArray<MiniElement2d> locfaces2;           //all local faces in radius xh
  NgArray<int> locfaces3;           // all faces in outer radius relh
  NgArray<INDEX> findex2;

  locfaces2.SetSize(0);
  locfaces3.SetSize(0);
  findex2.SetSize(0);

  int cluster = faces.Get(fstind).cluster;

  pstind = faces.Get(fstind).Face().PNum(1);
  p0 = points[pstind].P();
  
  locfaces2.Append(faces.Get(fstind).Face());
  findex2.Append(fstind);


  Box3d b1 (p0 - Vec3d(xh, xh, xh), p0 + Vec3d (xh, xh, xh));

  if (hashon)
    {
      hashtable.GetLocals(locfaces2, findex2, fstind, p0, xh);
    }
  else
    {
      for (i = 1; i <= faces.Size(); i++)
	{
	  const MiniElement2d & face = faces.Get(i).Face();
	  if (faces.Get(i).cluster == cluster && faces.Get(i).Valid() && i != fstind)
	    {
	      Box3d b2;
	      b2.SetPoint (points[face[0]].P());
	      b2.AddPoint (points[face[1]].P());
	      b2.AddPoint (points[face[2]].P());

	      if (b1.Intersect (b2))
		{
		  locfaces2.Append(faces.Get(i).Face());
		  findex2.Append(i);
		}
	    }
	}
    }

  //local faces for inner radius:
  for (i = 1; i <= locfaces2.Size(); i++)
    {
      const MiniElement2d & face = locfaces2.Get(i);
      const Point3d & p1 = points[face[0]].P();
      const Point3d & p2 = points[face[1]].P();
      const Point3d & p3 = points[face[2]].P();

      midp = Center (p1, p2, p3);

      if (Dist2 (midp, p0) <= relh * relh || i == 1)
	{
          locfaces.Append(locfaces2.Get(i));
	  findex.Append(findex2.Get(i));
	}
      else
	locfaces3.Append (i);
    }
  
  facesplit=locfaces.Size();
  
  
  //local faces for outer radius:
  for (i = 1; i <= locfaces3.Size(); i++)
    {
      locfaces.Append (locfaces2.Get(locfaces3.Get(i)));
      findex.Append (findex2.Get(locfaces3.Get(i)));
    }


  invpindex.SetSize (points.Size());
  for (i = 1; i <= locfaces.Size(); i++)
    for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
      {
	PointIndex pi = locfaces.Get(i).PNum(j);
        invpindex[pi] = PointIndex::INVALID;
      }

  for (i = 1; i <= locfaces.Size(); i++)
    {
      for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
	{
	  PointIndex pi = locfaces.Get(i).PNum(j);
	  if (!invpindex[pi].IsValid())
	    {
	      pindex.Append (pi);
              locpoints.Append (points[pi].P());
	      invpindex[pi] = pindex.Size()-1+PointIndex::BASE;
            }
          // locfaces.Elem(i).PNum(j) = locpoints.Append (points[pi].P());
          // }
	  // else
          locfaces.Elem(i).PNum(j) = invpindex[pi];
	}
    }



  if (connectedpairs)
    {
      for (i = 1; i <= locpoints.Size(); i++)
	{
	  int pind = pindex.Get(i);
	  if (pind >= 1 && pind <= connectedpairs->Size ())
	    {
	      for (j = 1; j <= connectedpairs->EntrySize(pind); j++)
		{
		  int oi = connectedpairs->Get(pind, j);
		  int other = invpindex.Get(oi);
		  if (other >= 1 && other <= pindex.Size() && 
		      pindex.Get(other) == oi)
		    {
		      // INDEX_2 coned(i, other);
		      // coned.Sort();
		      // (*testout) << "connected: " << locpoints.Get(i) << "-" << locpoints.Get(other) << endl;
		      getconnectedpairs.Set (INDEX_2::Sort (i, other), 1);
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
  return faces.Get(fstind).QualClass();
}


// returns all points connected with fi
void AdFront3 :: GetGroup (int fi,
			   NgArray<MeshPoint, PointIndex::BASE> & grouppoints,
			   NgArray<MiniElement2d> & groupelements,
			   NgArray<PointIndex, PointIndex::BASE> & pindex,
			   NgArray<INDEX> & findex) 
{
  // static NgArray<char> pingroup;
  int changed;

  pingroup.SetSize(points.Size());

  pingroup = 0;
  for (int j = 1; j <= 3; j++)
    pingroup[faces.Get(fi).Face().PNum(j)] = 1;

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
	    const MiniElement2d & face = f.Face();

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
  for (PointIndex pi : points.Range())
    if (points[pi].Valid())
      {
	grouppoints.Append (points[pi].P());
	pindex.Append (pi);
	invpindex[pi] = pindex.Size();
      }

  for (int i = 1; i <= faces.Size(); i++)
    if (faces.Get(i).Valid())
      {
	int fused = 0;
	for (int j = 1; j <= 3; j++)
	  if (pingroup[faces.Get(i).Face().PNum(j)])
	    fused++;

	if (fused >= 2)
	  {
	    groupelements.Append (faces.Get(i).Face());
	    findex.Append (i);
	  }
      }

  /*
  for (int i = 1; i <= groupelements.Size(); i++)
    for (int j = 1; j <= 3; j++)
      {
	groupelements.Elem(i).PNum(j) =
	  invpindex.Get(groupelements.Elem(i).PNum(j));
      }
  */
  for (auto & e : groupelements)
    for (int j = 1; j <= 3; j++)
      e.PNum(j) = invpindex[e.PNum(j)];
}


void AdFront3 :: SetStartFront (int /* baseelnp */)
{
  for (INDEX i = 1; i <= faces.Size(); i++)
    if (faces.Get(i).Valid())
      {
	const MiniElement2d & face = faces.Get(i).Face();
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


bool AdFront3 :: Inside (const Point<3> & p) const
{
  static Timer timer("AdFront3::Inside"); RegionTimer rt(timer);
  int cnt;
  Vec3d n, v1, v2;
  DenseMatrix a(3), ainv(3);
  Vector b(3), u(3);

  // random numbers:
  n.X() = 0.123871;
  n.Y() = 0.15432;
  n.Z() = -0.43989;

  cnt = 0;
  for (int i = 1; i <= faces.Size(); i++)
    if (faces.Get(i).Valid())
      {
	const Point<3> & p1 = points[faces.Get(i).Face().PNum(1)].P();
	const Point<3> & p2 = points[faces.Get(i).Face().PNum(2)].P();
	const Point<3> & p3 = points[faces.Get(i).Face().PNum(3)].P();

	v1 = p2 - p1;
	v2 = p3 - p1;

	a(0, 0) = v1.X();
	a(1, 0) = v1.Y();
	a(2, 0) = v1.Z();
	a(0, 1) = v2.X();
	a(1, 1) = v2.Y();
	a(2, 1) = v2.Z();
	a(0, 2) = -n.X();
	a(1, 2) = -n.Y();
	a(2, 2) = -n.Z();

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
			  const NgArray<int> * testfaces) const
{
  const Point<3> *line[2];
  line[0] = &lp1;
  line[1] = &lp2;


  Point3d pmin(lp1);
  Point3d pmax(lp1);
  pmin.SetToMin (lp2);
  pmax.SetToMax (lp2);
  
  NgArrayMem<int, 100> aprif;
  aprif.SetSize(0);
  
  if (!testfaces)
    facetree->GetIntersecting (pmin, pmax, aprif);
  else
    for (int i = 1; i <= testfaces->Size(); i++)
      aprif.Append (testfaces->Get(i));

  int cnt = 0;
  for (int ii = 1; ii <= aprif.Size(); ii++)
    {
      int i = aprif.Get(ii);
      
      if (faces.Get(i).Valid())
	{
	  const Point<3> *tri[3];
	  tri[0] = &points[faces.Get(i).Face().PNum(1)].P();
	  tri[1] = &points[faces.Get(i).Face().PNum(2)].P();
	  tri[2] = &points[faces.Get(i).Face().PNum(3)].P();
	  	  
	  if (IntersectTriangleLine (&tri[0], &line[0]))
	    cnt++;
	}
    }

  return ((cnt+1) % 2);
}
}

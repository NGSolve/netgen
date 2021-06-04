#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{
  // typedef BoxTree<3> DTREE;
  typedef DelaunayTree<3> DTREE;

  static const int deltetfaces[][3] = 
    { { 1, 2, 3 },
      { 2, 0, 3 },
      { 0, 1, 3 },
      { 1, 0, 2 } };


  class DelaunayTet
  {
    PointIndex pnums[4];
    int nb[4];

  public:
    DelaunayTet () { ; }

    DelaunayTet (const DelaunayTet & el)
    {
      for (int i = 0; i < 4; i++)
	pnums[i] = el[i];
    }

    DelaunayTet (const Element & el)
    {
      for (int i = 0; i < 4; i++)
	pnums[i] = el[i];
    }
    
    PointIndex & operator[] (int i) { return pnums[i]; }
    PointIndex operator[] (int i) const { return pnums[i]; }

    int & NB(int i) { return nb[i]; }
    int NB(int i) const { return nb[i]; }


    int FaceNr (INDEX_3 & face) const  // which face nr is it ?
    {
      for (int i = 0; i < 3; i++)
	if (pnums[i] != face.I1() && 
	    pnums[i] != face.I2() && 
	    pnums[i] != face.I3())
	  return i;
      return 3;
    }

    INDEX_3 GetFace (int i) const
    {
      return INDEX_3 (pnums[deltetfaces[i][0]],
		      pnums[deltetfaces[i][1]],
		      pnums[deltetfaces[i][2]]);
    }

    void GetFace (int i, Element2d & face) const
    {
      // face.SetType(TRIG);
      face[0] = pnums[deltetfaces[i][0]];
      face[1] = pnums[deltetfaces[i][1]];
      face[2] = pnums[deltetfaces[i][2]];
    }
  };







  /*
    Table to maintain neighbour elements
  */
  class MeshNB
  {
    // face nodes -> one element
    INDEX_3_CLOSED_HASHTABLE<int> faces;

    // 
    NgArray<DelaunayTet> & tets;

  public:

    // estimated number of points
    MeshNB (NgArray<DelaunayTet> & atets, int np)
      : faces(200), tets(atets)
    { ; }

    // add element with 4 nodes
    void Add (int elnr);

    // delete element with 4 nodes
    void Delete (int elnr)
    {
      DelaunayTet & el = tets.Elem(elnr);
      for (int i = 0; i < 4; i++)
	faces.Set (el.GetFace(i).Sort(), el.NB(i));
    }

    // get neighbour of element elnr in direction fnr 
    int GetNB (int elnr, int fnr)
    { 
      return tets.Get(elnr).NB(fnr); 
    }

    //
    void ResetFaceHT (int size)
    {
      faces.SetSize (size);
    }
  };



  void MeshNB :: Add (int elnr)
  {
    DelaunayTet & el = tets.Elem(elnr);

    for (int i = 0; i < 4; i++)
      {
	INDEX_3 i3 = INDEX_3::Sort (el.GetFace(i));

	int posnr;
	
	if (!faces.PositionCreate (i3, posnr))
	  {
	    // face already in use
	    int othertet = faces.GetData (posnr);

	    el.NB(i) = othertet;
	    if (othertet)
	      {
		int fnr = tets.Get(othertet).FaceNr (i3);
		tets.Elem(othertet).NB(fnr) = elnr;
	      }
	  }
	else
	  {
	    faces.SetData (posnr, elnr);	
	    el.NB(i) = 0;
	  }
      }
  }





  /*
    connected lists of cosphereical elements
  */
  class SphereList 
  {
    NgArray<int> links;
  public:
    SphereList () 
    { ; }

    void AddElement (int elnr)
    {
      if (elnr > links.Size())
	links.Append (1);
      links.Elem(elnr) = elnr;
    }

    void DeleteElement (int elnr)
    {
      links.Elem(elnr) = 0;
    }    
  
    void ConnectElement (int eli, int toi)
    {
      links.Elem (eli) = links.Get (toi);
      links.Elem (toi) = eli;
    }
      
    void GetList (int eli, NgArray<int> & linked) const;

    template <typename TFUNC>
    void IterateList (int eli, TFUNC func)
    {
      int pi = eli;
      do
        {
          func(pi);
          pi = links.Get(pi);
        }
      while (pi != eli);
    }
    
  };


  void SphereList :: GetList (int eli, NgArray<int> & linked) const
  {
    linked.SetSize (0);
    int pi = eli;

    do
      {
#ifdef DELAUNAY_DEBUG
	if (pi <= 0 || pi > links.Size())
	  {
	    cerr << "link, error " << endl;
	    cerr << "pi = " << pi << " linked.s = " << linked.Size() << endl;
	    exit(1);
	  }
	if (linked.Size() > links.Size())
	  {
	    cerr << "links have loop" << endl;
	    exit(1);
	  }
#endif
	linked.Append (pi);
	pi = links.Get(pi);
      }
    while (pi != eli);
  }





  void AddDelaunayPoint (PointIndex newpi, const Point3d & newp, 
			 NgArray<DelaunayTet> & tempels, 
			 Mesh & mesh,
			 DTREE & tettree,
			 MeshNB & meshnb,
			 NgArray<Point<3> > & centers, NgArray<double> & radi2,
			 NgArray<int> & connected, NgArray<int> & treesearch, 
			 NgArray<int> & freelist, SphereList & list,
			 IndexSet & insphere, IndexSet & closesphere)
  {
    static Timer t("Meshing3::AddDelaunayPoint");// RegionTimer reg(t);
    // static Timer tsearch("addpoint, search");
    // static Timer tinsert("addpoint, insert");

    /*
      find any sphere, such that newp is contained in
    */
    // DelaunayTet el;
    int cfelind = -1;

    const Point<3> * pp[4];
    Point<3> pc;
    Point3d tpmin, tpmax;



    /*
    // stop search if intersecting point is close enough to center
    tettree.GetFirstIntersecting (newp, newp, treesearch, [&](const auto pi)
            {
                double quot = Dist2 (centers.Get(pi), newp);
                return (quot < 0.917632 * radi2.Get(pi));
            } );
    // tettree.GetIntersecting (newp, newp, treesearch);

    double quot,minquot(1e20);

    for (auto jjj : treesearch)
      {
	quot = Dist2 (centers.Get(jjj), newp) / radi2.Get(jjj);
	
	if((cfelind == -1 || quot < 0.99*minquot) && quot < 1)
	  {
	    minquot = quot;
	    el = tempels.Get(jjj);
	    cfelind = jjj;
	    if(minquot < 0.917632)
	      break;
	  }
      }
    */

    // tsearch.Start();
    double minquot{1e20};
    tettree.GetFirstIntersecting
      (newp, newp, [&](const auto pi)
       {
         double rad2 = radi2.Get(pi);
         double d2 = Dist2 (centers.Get(pi), newp); //  / radi2.Get(pi);
         if (d2 >= rad2) return false;
         
         // if (d2 < 0.917632 * rad2)
         if (d2 < 0.99 * rad2)
           {
             cfelind = pi;
             return true;
           }
         
	if (cfelind == -1 || d2 < 0.99*minquot*rad2) 
	  {
	    minquot = d2/rad2;
	    cfelind = pi;
          }
        return false;
       } );
    // tsearch.Stop();

    if (cfelind == -1)
      {
	PrintWarning ("Delaunay, point not in any sphere");
	return;
      }
	
    /*
      insphere:     point is in sphere -> delete element
      closesphere:  point is close to sphere -> considered for same center
    */

    // save overestimate
    insphere.SetMaxIndex (2 * tempels.Size() + 5 * mesh.GetNP());
    closesphere.SetMaxIndex (2 * tempels.Size() + 5 * mesh.GetNP());

    insphere.Clear();
    closesphere.Clear();


    insphere.Add (cfelind);
      
    bool changed = true;
    int nstarti = 1, starti;

    while (changed)
      {
	changed = false;
	starti = nstarti;
	nstarti = insphere.GetArray().Size()+1;


	// if point in sphere, then it is also closesphere
	for (int j = starti; j < nstarti; j++)
	  {
	    int helind = insphere.GetArray().Get(j);
	    if (!closesphere.IsIn (helind))
	      closesphere.Add (helind);
	  }

	// add connected spheres to insphere - list
	for (int j = starti; j < nstarti; j++)
	  {
	    list.IterateList (insphere.GetArray().Get(j),
                              [&](int celind)
                              {
                                if (tempels.Get(celind)[0] != PointIndex(PointIndex::INVALID) && 
                                    !insphere.IsIn (celind))
                                  {
                                    changed = true;
                                    insphere.Add (celind);
                                  }
                              });
	  }
	
	// check neighbour-tets
	for (int j = starti; j < nstarti; j++)
	  for (int k = 0; k < 4; k++)
	    {
	      int helind = insphere.GetArray().Get(j);
	      int nbind = meshnb.GetNB (helind, k);

	      if (nbind && !insphere.IsIn (nbind) )
		{
                  double d2 = Dist2 (centers.Get(nbind), newp);
		  if (d2 < radi2.Get(nbind) * (1+1e-8) )
		    closesphere.Add (nbind);
		    
		  if (d2 < radi2.Get(nbind) * (1 + 1e-12))
		    {
		      // point is in sphere -> remove tet
		      insphere.Add (nbind);
		      changed = true;
		    }
		  else
		    {
		      INDEX_3 i3 = tempels.Get(helind).GetFace (k);
		      const Point<3> & p1 = mesh[PointIndex (i3.I1())];
		      const Point<3> & p2 = mesh[PointIndex (i3.I2())];
		      const Point<3> & p3 = mesh[PointIndex (i3.I3())];

		      Vec<3> n = Cross (p2-p1, p3-p1);
                      n /= n.Length();

		      double dist = n * (newp-p1);
                      double scal = n * (mesh.Point (tempels.Get(helind)[k])-p1);
		      if (scal > 0) dist *= -1;

		      if (dist > -1e-10)  // 1e-10
			{
			  insphere.Add (nbind);
			  changed = true;
			}
		    }
		}
	    }
      } // while (changed)

    // NgArray<Element> newels;
    static NgArray<DelaunayTet> newels;
    newels.SetSize(0);
    
    Element2d face(TRIG);

    for (int celind : insphere.GetArray())
      for (int k = 0; k < 4; k++)
	{
	  int nbind = meshnb.GetNB (celind, k);

	  if (!nbind || !insphere.IsIn (nbind))
	    {
	      tempels.Get (celind).GetFace (k, face);
		
	      // Element newel(TET);
              DelaunayTet newel;
	      for (int l = 0; l < 3; l++)
                newel[l] = face[l];
              newel[3] = newpi;

	      newels.Append (newel);

#ifdef DEBUG_DELAUNAY
              Vec<3> v1 = mesh[face[1]] - mesh[face[0]];
              Vec<3> v2 = mesh[face[2]] - mesh[face[0]];
	      Vec<3> n = Cross (v1, v2);

              n.Normalize();
	      if (n * Vec3d(mesh.Point (face[0]), 
			    mesh.Point (tempels.Get(celind)[k]))
		  > 0)
		n *= -1;

              double hval = n *  ( newp - mesh[face[0]]);
		
	      if (hval > -1e-12)
		{
		  cerr << "vec to outer" << endl;
		  (*testout) << "vec to outer, hval = " << hval << endl;
		  (*testout) << "v1 x v2 = " << Cross (v1, v2) << endl;
		  (*testout) << "facep: "
			     << mesh.Point (face[0]) << " "
			     << mesh.Point (face[1]) << " "
			     << mesh.Point (face[2]) << endl;
		}
#endif
	    }
	}

    meshnb.ResetFaceHT (10*insphere.GetArray().Size()+1);

    for (auto celind : insphere.GetArray())
      {
	meshnb.Delete (celind); 
	list.DeleteElement (celind);
	  
	for (int k = 0; k < 4; k++)
        tempels.Elem(celind)[k] = PointIndex::INVALID;
        
        tettree.DeleteElement (celind);
	freelist.Append (celind);
      }

    bool hasclose = false;
    for (int ind : closesphere.GetArray())
      {
	if (!insphere.IsIn(ind) &&
	    fabs (Dist2 (centers.Get (ind), newp) - radi2.Get(ind)) < 1e-8 )
	  hasclose = true;
      }

    /*
    for (int j = 1; j <= newels.Size(); j++)
      {
        const auto & newel = newels.Get(j);
    */
    for (const auto & newel : newels)
      {
	int nelind;

	if (!freelist.Size())
	  {
	    tempels.Append (newel);
	    nelind = tempels.Size();
	  }
	else
	  {
	    nelind = freelist.Last();
	    freelist.DeleteLast();

	    tempels.Elem(nelind) = newel;
	  }

	meshnb.Add (nelind);
	list.AddElement (nelind);

	for (int k = 0; k < 4; k++)
	  pp[k] = &mesh.Point (newel[k]);

	if (CalcSphereCenter (&pp[0], pc) )
	  {
#ifdef DEBUG_DELAUNAY
	    PrintSysError ("Delaunay: New tet is flat");

	    (*testout) << "new tet is flat" << endl;
	    for (int k = 0; k < 4; k++)
	      (*testout) << newel[k] << " ";
	    (*testout) << endl;
	    for (int k = 0; k < 4; k++)
	      (*testout) << *pp[k-1] << " ";
	    (*testout) << endl;
#endif
            ;
	  }

	double r2 = Dist2 (*pp[0], pc);
	if (hasclose)
          /*
	  for (int k = 1; k <= closesphere.GetArray().Size(); k++)
	    {
	      int csameind = closesphere.GetArray().Get(k); 
          */
          for (int csameind : closesphere.GetArray())
            {
	      if (!insphere.IsIn(csameind) &&
		  fabs (r2 - radi2.Get(csameind)) < 1e-10 && 
		  Dist2 (pc, centers.Get(csameind)) < 1e-20)
		{
		  pc = centers.Get(csameind);
		  r2 = radi2.Get(csameind);
		  list.ConnectElement (nelind, csameind);
		  break;
		}
	    }
      
	if (centers.Size() < nelind)
	  {
	    centers.Append (pc);
	    radi2.Append (r2);
	  }
	else
	  {
	    centers.Elem(nelind) = pc;
	    radi2.Elem(nelind) = r2;
	  }

	closesphere.Add (nelind);
	  
	tpmax = tpmin = *pp[0];
	for (int k = 1; k <= 3; k++)
	  {
	    tpmin.SetToMin (*pp[k]);
	    tpmax.SetToMax (*pp[k]);
	  }
	tpmax = tpmax + 0.01 * (tpmax - tpmin);
        // tinsert.Start();
	tettree.Insert (tpmin, tpmax, nelind);
        // tinsert.Stop();
      }
  }






  void Delaunay1 (Mesh & mesh, const MeshingParameters & mp, AdFront3 * adfront,
		  NgArray<DelaunayTet> & tempels,
		  int oldnp, DelaunayTet & startel, Point3d & pmin, Point3d & pmax)
  {
    static Timer t("Meshing3::Delaunay1"); RegionTimer reg(t);
    
    NgArray<Point<3>> centers;
    NgArray<double> radi2;
  
    Box<3> bbox(Box<3>::EMPTY_BOX);

    for (auto & face : adfront->Faces())
      for (PointIndex pi : face.Face().PNums())      
        bbox.Add (mesh.Point(pi));

    for (PointIndex pi : mesh.LockedPoints())
      bbox.Add (mesh.Point (pi));
    
    pmin = bbox.PMin();
    pmax = bbox.PMax();


    Vec<3> vdiag = pmax-pmin;
    // double r1 = vdiag.Length();
    double r1 = sqrt (3.0) * max3(vdiag(0), vdiag(1), vdiag(2));
    vdiag = Vec<3> (r1, r1, r1);
    //double r2;

    Point<3> pmin2 = pmin - 8 * vdiag;
    Point<3> pmax2 = pmax + 8 * vdiag;

    Point<3> cp1(pmin2), cp2(pmax2), cp3(pmax2), cp4(pmax2);
    cp2(0) = pmin2(0);
    cp3(1) = pmin2(1);
    cp4(2) = pmin2(2);


    size_t np = mesh.GetNP();

    startel[0] = mesh.AddPoint (cp1);
    startel[1] = mesh.AddPoint (cp2);
    startel[2] = mesh.AddPoint (cp3);
    startel[3] = mesh.AddPoint (cp4);

    // flag points to use for Delaunay:
    Array<bool, PointIndex> usep(np);
    usep = false;

    for (auto & face : adfront->Faces())
      for (PointIndex pi : face.Face().PNums())      
        usep[pi] = true;
    
    for (size_t i = oldnp + PointIndex::BASE; 
	 i < np + PointIndex::BASE; i++)
      usep[i] = true;

    for (PointIndex pi : mesh.LockedPoints())
      usep[pi] = true;
    

    NgArray<int> freelist;

    int cntp = 0;

    MeshNB meshnb (tempels, mesh.GetNP() + 5);
    SphereList list;

    pmin2 = pmin2 + 0.1 * (pmin2 - pmax2);
    pmax2 = pmax2 + 0.1 * (pmax2 - pmin2);

    DTREE tettree(pmin2, pmax2);


    tempels.Append (startel);
    meshnb.Add (1);
    list.AddElement (1);
    NgArray<int> connected, treesearch;

    Box<3> tbox(Box<3>::EMPTY_BOX);
    for (size_t k = 0; k < 4; k++)
      tbox.Add (mesh.Point(startel[k]));
    Point<3> tpmin = tbox.PMin();
    Point<3> tpmax = tbox.PMax();
    
    tpmax = tpmax + 0.01 * (tpmax - tpmin);
    tettree.Insert (tpmin, tpmax, 1);

    Point<3> pc;
	  
    const Point<3> * pp[4];
    for (int k = 0; k < 4; k++)
      pp[k] = &mesh.Point (startel[k]);
    CalcSphereCenter (&pp[0], pc);
    
    centers.Append (pc);
    radi2.Append (Dist2 (*pp[0], pc));


    IndexSet insphere(mesh.GetNP());
    IndexSet closesphere(mesh.GetNP());

    // "random" reordering of points  (speeds a factor 3 - 5 !!!)
    NgArray<PointIndex, PointIndex::BASE, PointIndex> mixed(np);
    // int prims[] = { 11, 13, 17, 19, 23, 29, 31, 37 };
    // int prims[] = { 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
    int prims[] = { 211, 223, 227, 229, 233, 239, 241, 251, 257, 263 }; 
    int prim;
  
    {
      int i = 0;
      while (np % prims[i] == 0) i++;
      prim = prims[i];
    }

    // for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End()-4; pi++)
    for (PointIndex pi : mesh.Points().Range().Modify(0, -4))
      mixed[pi] = PointIndex ( (prim * pi) % np + PointIndex::BASE );

    // for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End()-4; pi++)
    for (PointIndex pi : mesh.Points().Range().Modify(0, -4))      
      {
	if (pi % 1000 == 0)
	  {
	    if (pi % 10000 == 0)
	      PrintDot ('+');
	    else
	      PrintDot ('.');
	  }

	multithread.percent = 100.0 * pi / np;
	if (multithread.terminate)
	  break;

	PointIndex newpi = mixed[pi];

	if (!usep[newpi])
	  continue;

	cntp++;

	const MeshPoint & newp = mesh[newpi];
      
	AddDelaunayPoint (newpi, newp, tempels, mesh,
			  tettree, meshnb, centers, radi2, 
			  connected, treesearch, freelist, list, insphere, closesphere);

      }
    
    for (int i = tempels.Size(); i >= 1; i--)
      if (tempels.Get(i)[0] <= 0)
	tempels.DeleteElement (i);

    PrintDot ('\n');

    PrintMessage (3, "Points: ", cntp);
    PrintMessage (3, "Elements: ", tempels.Size());
    PrintMessage (3, "Tree data entries per element: ", 1.0*tettree.N*tettree.GetNLeaves() / tempels.Size());
    PrintMessage (3, "Tree nodes per element: ", 1.0*tettree.GetNNodes() / tempels.Size());
    //   (*mycout) << cntp << " / " << tempels.Size() << " points/elements" << endl;

    /*
      cout << "tempels: ";
      tempels.PrintMemInfo(cout);
      cout << "Searchtree: ";
      tettree.Tree().PrintMemInfo(cout);
      cout << "MeshNB: ";
      meshnb.PrintMemInfo(cout);
    */
  }


  void DelaunayRemoveDegenerated( const Mesh::T_POINTS & points, NgArray<DelaunayTet> & tempels, int np )
  {
    static Timer tdegenerated("Delaunay - remove degenerated"); RegionTimer rt(tdegenerated);

    NgBitArray badnode(points.Size());
    badnode.Clear();

    int ndeg = 0;
    for (int i = 1; i <= tempels.Size(); i++)
      {
	Element el(4);
	for (int j = 0; j < 4; j++)
	  el[j] = tempels.Elem(i)[j];
	//      Element & el = tempels.Elem(i);
	const Point3d & lp1 = points[el[0]];
	const Point3d & lp2 = points[el[1]];
	const Point3d & lp3 = points[el[2]];
	const Point3d & lp4 = points[el[3]];
	Vec3d v1(lp1, lp2);
	Vec3d v2(lp1, lp3);
	Vec3d v3(lp1, lp4);
	Vec3d n = Cross (v1, v2);
	double vol = n * v3;

	double h = v1.Length() + v2.Length() + v3.Length();
	if (fabs (vol) < 1e-8 * (h * h * h) &&
	    (el[0] <= np && el[1] <= np &&
	     el[2] <= np && el[3] <= np) )   // old: 1e-12
	  {
	    badnode.Set(el[0]);
	    badnode.Set(el[1]);
	    badnode.Set(el[2]);
	    badnode.Set(el[3]);
	    ndeg++;
	    (*testout) << "vol = " << vol << " h = " << h << endl;
	  }

	if (vol > 0)
	  Swap (el[2], el[3]);
      }

    auto ne = tempels.Size();
    for (int i = ne; i >= 1; i--)
      {
	const DelaunayTet & el = tempels.Get(i);
	if (badnode.Test(el[0]) ||
	    badnode.Test(el[1]) ||
	    badnode.Test(el[2]) ||
	    badnode.Test(el[3]) )
	  tempels.DeleteElement(i);
      }

  
    PrintMessage (3, ndeg, " degenerated elements removed");
  }

  // Remove flat tets containing two adjacent surface trigs
  void DelaunayRemoveTwoTriaTets( const Mesh & mesh, NgArray<DelaunayTet> & tempels, NgArray<int> & openels )
  {
    static Timer topenel("Delaunay - find openel"); RegionTimer rt(topenel);

    // find surface triangles which are no face of any tet
    BitArray bnd_points( mesh.GetNP() + PointIndex::BASE );
    bnd_points.Clear();

    for (int i = 1; i <= mesh.GetNOpenElements(); i++)
      {
        const Element2d & tri = mesh.OpenElement(i);
        bnd_points.SetBit(tri[0]);
        bnd_points.SetBit(tri[1]);
        bnd_points.SetBit(tri[2]);
      }

    auto ne = tempels.Size();
    Array<int> tets_with_3_bnd_points(ne);
    atomic<int> cnt = 0;

    // table of tets with >= 2 boundary points, store in extra array tets with >=3 boundary points
    auto p2el = ngcore::CreateSortedTable<int, PointIndex>( ne,
               [&](auto & table, int ei)
               {
                 const auto & el = tempels[ei];
                 int num_bnd_points = 0;

                 for( auto i : Range(4) )
                   if(bnd_points[el[i]])
                     num_bnd_points++;

                 if(num_bnd_points>1)
                   {
                     table.Add (el[0], ei);
                     table.Add (el[1], ei);
                     table.Add (el[2], ei);
                     table.Add (el[3], ei);
                   }

                 // table creator is running this code 2 times, only store tets on last run
                 if(table.GetMode()==3 && num_bnd_points>2)
                   tets_with_3_bnd_points[cnt++] = ei;
               }, mesh.GetNP());

    tets_with_3_bnd_points.SetSize(cnt);

    static Timer t1("Build face table"); t1.Start();
    ngcore::ClosedHashTable< ngcore::INT<3>, int > face_table( 4*cnt + 3 );
    for(auto ei : tets_with_3_bnd_points)
        for(auto j : Range(4))
        {
            auto i3_ = tempels[ei].GetFace (j);
            ngcore::INT<3> i3 = {i3_[0], i3_[1], i3_[2]};
            if(bnd_points[i3[0]] && bnd_points[i3[1]] && bnd_points[i3[2]])
            {
                i3.Sort();
                face_table.Set( i3, true );
            }
        }
    t1.Stop();

    static Timer t2("check faces"); t2.Start();

    openels.SetSize(0);
    for (int i = 1; i <= mesh.GetNOpenElements(); i++)
      {
	const Element2d & tri = mesh.OpenElement(i);
        ngcore::INT<3> i3(tri[0], tri[1], tri[2]);
	i3.Sort();
        if(!face_table.Used(i3))
            openels.Append(i);
      }

    t2.Stop();

    auto p2sel = ngcore::CreateSortedTable<int, PointIndex>( Range(openels.Size()),
           [&](auto & table, int i)
           {
             auto openel_i = openels[i];
             const Element2d & tri = mesh.OpenElement(openel_i);
             table.Add(tri[0], openel_i);
             table.Add(tri[1], openel_i);
             table.Add(tri[2], openel_i);
           }, mesh.GetNP());

    ngcore::BitArray badnode(mesh.GetNP()+PointIndex::BASE);
    badnode.Clear();

    ngcore::ParallelForRange(openels.Size(), [&] (auto myrange) {
        for (auto i_ : myrange)
          {
            auto i = openels[i_];
            const Element2d & tri = mesh.OpenElement(i);

            for( auto edge : Range(3) )
            {
                auto pi0 = tri[edge];
                auto pi1 = tri[(edge+1)%3];
                if(pi0>pi1)
                    Swap(pi0, pi1);

                // find other trig with edge pi0, pi1
                int i_other = -1;
                for(auto ii : p2sel[pi0])
                {
                    if(ii==i)
                        continue;
                    auto & tri_other = mesh.OpenElement(ii);
                    if(tri_other[0]==pi1 || tri_other[1]==pi1 || tri_other[2]==pi1)
                    {
                        i_other = ii;
                        break;
                    }
                }

                if(i_other>i)
                {
                    auto & tri_other = mesh.OpenElement(i_other);
                    PointIndex pi2 = tri[(edge+2)%3];
                    PointIndex pi3 = tri_other[0]+tri_other[1]+tri_other[2] - pi0 - pi1;
                    if(pi2>pi3)
                        Swap(pi2, pi3);

                    // search for tet with edge pi2-pi3
                    for(auto ei : p2el[pi2])
                    {
                        auto & el = tempels[ei];

                        if(el[0]==pi3 || el[1]==pi3 || el[2]==pi3 || el[3]==pi3)
                        {
                            const Point3d & p1 = mesh[pi0];
                            const Point3d & p2 = mesh[pi1];
                            const Point3d & p3 = mesh[pi2];
                            const Point3d & p4 = mesh[pi3];
                            Vec3d v1(p1, p2);
                            Vec3d v2(p1, p3);
                            Vec3d v3(p1, p4);
                            Vec3d n = Cross (v1, v2);
                            double vol = n * v3;

                            double h = v1.Length() + v2.Length() + v3.Length();
                            if (fabs (vol) < 1e-4 * (h * h * h))   // old: 1e-12
                            {
                                badnode.SetBitAtomic(pi2);
                                badnode.SetBitAtomic(pi3);
                            }
                            break;
                        }
                    }
                }
            }
          }
    });

    for (int i = ne; i >= 1; i--)
      {
	const DelaunayTet & el = tempels.Get(i);
	if (badnode[el[0]] ||
	    badnode[el[1]] ||
	    badnode[el[2]] ||
	    badnode[el[3]] )
	  tempels.DeleteElement(i);
      }
  }

  void DelaunayRemoveIntersecting( const Mesh & mesh, NgArray<DelaunayTet> & tempels, NgArray<int> & openels, Point3d pmin, Point3d pmax )
  {
    static Timer trem_intersect("Delaunay - remove intersecting"); RegionTimer rt(trem_intersect);

    // find intersecting:
    PrintMessage (3, "Remove intersecting");
    if (openels.Size())
      {
	BoxTree<3> setree(pmin, pmax);

	/*      
		cout << "open elements in search tree: " << openels.Size() << endl;
		cout << "pmin, pmax = " << pmin << " - " << pmax << endl;
	*/

	for (int i = 1; i <= openels.Size(); i++)
	  {
	    int fnr;
	    fnr = openels.Get(i);
	    if (fnr)
	      {
		const Element2d & tri = mesh.OpenElement(fnr);
	      
		Point3d ltpmin (mesh.Point(tri[0]));
		Point3d ltpmax (ltpmin);
	      
		for (int k = 2; k <= 3; k++)
		  {
		    ltpmin.SetToMin (mesh.Point (tri.PNum(k)));
		    ltpmax.SetToMax (mesh.Point (tri.PNum(k)));
		  }
		setree.Insert (ltpmin, ltpmax, fnr);
	      }
	  }
      
	NgArray<int> neartrias;
	for (int i = 1; i <= tempels.Size(); i++)
	  {
	    const Point<3> *pp[4];
	    int tetpi[4];
	    DelaunayTet & el = tempels.Elem(i);
	  
	    int intersect = 0;
	  
	    for (int j = 0; j < 4; j++)
	      {
		pp[j] = &mesh.Point(el[j]);
		tetpi[j] = el[j];
	      }
	  
	    Point3d tetpmin(*pp[0]);
	    Point3d tetpmax(tetpmin);
	    for (int j = 1; j < 4; j++)
	      {
		tetpmin.SetToMin (*pp[j]);
		tetpmax.SetToMax (*pp[j]);
	      }
	    tetpmin = tetpmin + 0.01 * (tetpmin - tetpmax);
	    tetpmax = tetpmax + 0.01 * (tetpmax - tetpmin);
	  
	    setree.GetIntersecting (tetpmin, tetpmax, neartrias);
	  
	  
	    //      for (j = 1; j <= mesh.GetNSE(); j++)
	    //	{
	    for (int jj = 1; jj <= neartrias.Size(); jj++)
	      {
		int j = neartrias.Get(jj);
	      
		const Element2d & tri = mesh.OpenElement(j);
		const Point<3> *tripp[3];
		int tripi[3];
	      
		for (int k = 1; k <= 3; k++)
		  {
		    tripp[k-1] = &mesh.Point (tri.PNum(k));
		    tripi[k-1] = tri.PNum(k);
		  }
	      
		if (IntersectTetTriangle (&pp[0], &tripp[0], tetpi, tripi))
		  {
		    /*
		    int il1, il2;
		    (*testout) << "intersect !" << endl;
		    (*testout) << "triind: ";
		    for (il1 = 0; il1 < 3; il1++)
		      (*testout) << " " << tripi[il1];
		    (*testout) << endl;
		    (*testout) << "tetind: ";
		    for (il2 = 0; il2 < 4; il2++)
		      (*testout) << " " << tetpi[il2];
		    (*testout) << endl;
		  
		    (*testout) << "trip: ";
		    for (il1 = 0; il1 < 3; il1++)
		      (*testout) << " " << *tripp[il1];
		    (*testout) << endl;
		    (*testout) << "tetp: ";
		    for (il2 = 0; il2 < 4; il2++)
		      (*testout) << " " << *pp[il2];
		    (*testout) << endl;
		    */
		  
		  
		    intersect = 1;
		    break;
		  }
	      }
	  
	  
	    if (intersect)
	      {
		tempels.DeleteElement(i);
		i--;
	      }
	  }
      }
  }

  void DelaunayRemoveOuter( const Mesh & mesh, NgArray<DelaunayTet> & tempels, AdFront3 * adfront )
  {
    static Timer trem_outer("Delaunay - remove outer"); RegionTimer rt(trem_outer);


    PrintMessage (3, "Remove outer");

    // find connected tets (with no face between, and no hole due
    // to removed intersecting tets.
    //  INDEX_3_HASHTABLE<INDEX_2> innerfaces(np);

  
    INDEX_3_HASHTABLE<int> boundaryfaces(mesh.GetNOpenElements()/3+1);
    /*
    for (int i = 1; i <= mesh.GetNOpenElements(); i++)
      {
	const Element2d & tri = mesh.OpenElement(i);
	INDEX_3 i3 (tri[0], tri[1], tri[2]);
	i3.Sort();
	boundaryfaces.PrepareSet (i3);
      }
    */
    for (const Element2d & tri : mesh.OpenElements())
      {
	INDEX_3 i3 (tri[0], tri[1], tri[2]);
	i3.Sort();
	boundaryfaces.PrepareSet (i3);
      }
    boundaryfaces.AllocateElements();
    for (int i = 1; i <= mesh.GetNOpenElements(); i++)
      {
	const Element2d & tri = mesh.OpenElement(i);
	INDEX_3 i3 (tri[0], tri[1], tri[2]);
	i3.Sort();
	boundaryfaces.Set (i3, 1);
      }

    /*
    for (int i = 0; i < tempels.Size(); i++)
      for (int j = 0; j < 4; j++)
	tempels[i].NB(j) = 0;
    */
    for (auto & el : tempels)
      for (int j = 0; j < 4; j++)
	el.NB(j) = 0;
      
    TABLE<int,PointIndex::BASE> elsonpoint(mesh.GetNP());
    /*
    for (int i = 0; i < tempels.Size(); i++)
      {
	const DelaunayTet & el = tempels[i];
    */
    for (const DelaunayTet & el : tempels)
      {
	INDEX_4 i4(el[0], el[1], el[2], el[3]);
	i4.Sort();
	elsonpoint.IncSizePrepare (i4.I1());
	elsonpoint.IncSizePrepare (i4.I2());
      }

    elsonpoint.AllocateElementsOneBlock();

    for (int i = 0; i < tempels.Size(); i++)
      {
	const DelaunayTet & el = tempels[i];
	INDEX_4 i4(el[0], el[1], el[2], el[3]);
	i4.Sort();
	elsonpoint.Add (i4.I1(), i+1);
	elsonpoint.Add (i4.I2(), i+1);
      }

    //  cout << "elsonpoint mem: ";
    //  elsonpoint.PrintMemInfo(cout);

    INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht(100);   
  
    Element2d hel(TRIG);
    // for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
    for (PointIndex pi : mesh.Points().Range())
      {
	faceht.SetSize (4 * elsonpoint[pi].Size());
	for (int ii = 0; ii < elsonpoint[pi].Size(); ii++)
	  {
	    int i = elsonpoint[pi][ii];
	    const DelaunayTet & el = tempels.Get(i);

	    for (int j = 1; j <= 4; j++)
	      {
		el.GetFace (j-1, hel);
		hel.Invert();
		hel.NormalizeNumbering();
	      
		if (hel[0] == pi)
		  {
		    INDEX_3 i3(hel[0], hel[1], hel[2]);
		  
		    if (!boundaryfaces.Used (i3))
		      {
			if (faceht.Used (i3))
			  {
			    INDEX_2 i2 = faceht.Get(i3);
			  
			    tempels.Elem(i).NB(j-1) = i2.I1();
			    tempels.Elem(i2.I1()).NB(i2.I2()-1) = i;
			  }
			else
			  {
			    hel.Invert();
			    hel.NormalizeNumbering();
			    INDEX_3 i3i(hel[0], hel[1], hel[2]);
			    INDEX_2 i2(i, j);
			    faceht.Set (i3i, i2);
			  }
		      }
		  }
	      }
	  }
      }
  
    /*
      for (i = 1; i <= tempels.Size(); i++)
      {
      const DelaunayTet & el = tempels.Get(i);
      for (j = 1; j <= 4; j++)
      {
      INDEX_3 i3;
      Element2d face;
      el.GetFace1 (j, face);
      for (int kk = 1; kk <= 3; kk++)
      i3.I(kk) = face.PNum(kk);

      i3.Sort();
      if (!boundaryfaces.Used (i3))
      {
      if (innerfaces.Used(i3))
      {
      INDEX_2 i2;
      i2 = innerfaces.Get(i3);
      i2.I2() = i;
      innerfaces.Set (i3, i2);
      }
      else
      {
      INDEX_2 i2;
      i2.I1() = i;
      i2.I2() = 0;
      innerfaces.Set (i3, i2);
      }
      }
      }
      }
    */

    /*
      (*testout) << "nb elements:" << endl;
      for (i = 1; i <= tempels.Size(); i++)
      {
      (*testout) << i << " ";
      for (j = 1; j <= 4; j++)
      (*testout) << tempels.Get(i).NB1(j) << " ";
      (*testout) << endl;
      }
  
      (*testout) << "pairs:" << endl;
      for (i = 1; i <= innerfaces.GetNBags(); i++)
      for (j = 1; j <= innerfaces.GetBagSize(i); j++)
      {
      INDEX_3 i3;
      INDEX_2 i2;
      innerfaces.GetData (i, j, i3, i2);
      (*testout) << i2 << endl;
      }
    */







    /*
      cout << "innerfaces: ";
      innerfaces.PrintMemInfo (cout);
    */

    //  cout << "boundaryfaces: ";
    //  boundaryfaces.PrintMemInfo (cout);


    PrintMessage (5, "tables filled");
 

    auto ne = tempels.Size();
    NgBitArray inner(ne), outer(ne);
    inner.Clear();
    outer.Clear();
    NgArray<int> elstack;

    /*
      int starti = 0;
      for (i = 1; i <= ne; i++)
      {
      const Element & el = tempels.Get(i);
      for (j = 1; j <= 4; j++)
      for (k = 1; k <= 4; k++)
      if (el.PNum(j) == startel.PNum(k))
      {
      outer.Set(i);
      starti = i;
      }
      }
    */

    int lowest_undefined_el = 1;
    while (1)
      {
	int inside;
	bool done = 1;

	int i;
	for (i = lowest_undefined_el; i <= ne; i++)
	  if (!inner.Test(i) && !outer.Test(i))
	    {
              lowest_undefined_el = i+1;
	      done = 0;
	      break;
	    }

	if (done) break;
      
	const DelaunayTet & el = tempels.Get(i);
	const Point3d & p1 = mesh.Point (el[0]);
	const Point3d & p2 = mesh.Point (el[1]);
	const Point3d & p3 = mesh.Point (el[2]);
	const Point3d & p4 = mesh.Point (el[3]);
      
	Point3d ci = Center (p1, p2, p3, p4);

	inside = adfront->Inside (ci);

	/*
	  cout << "startel: " << i << endl;
	  cout << "inside = " << inside << endl;
	  cout << "ins2 = " << adfront->Inside (Center (ci, p1)) << endl;
	  cout << "ins3 = " << adfront->Inside (Center (ci, p2)) << endl;
	*/
      
	elstack.SetSize(0);
	elstack.Append (i);
  
	while (elstack.Size())
	  {
	    int ei = elstack.Last();
	    elstack.DeleteLast();
	  
	    if (!inner.Test(ei) && !outer.Test(ei))
	      {
		if (inside)
		  inner.Set(ei);
		else
		  outer.Set(ei);


		for (int j = 1; j <= 4; j++)
		  {
		    INDEX_3 i3 = tempels.Get(ei).GetFace(j-1);
		    /*
		    Element2d face;
		    tempels.Get(ei).GetFace(j, face);
		    for (int kk = 1; kk <= 3; kk++)
		      i3.I(kk) = face.PNum(kk);
		    */
		    i3.Sort();
		  

		    if (tempels.Get(ei).NB(j-1))
		      elstack.Append (tempels.Get(ei).NB(j-1));

		    /*
		      if (innerfaces.Used(i3))
		      {
		      INDEX_2 i2 = innerfaces.Get(i3);
		      int other = i2.I1() + i2.I2() - ei;

		      if (other != tempels.Get(ei).NB1(j))
		      cerr << "different1 !!" << endl;

		      if (other)
		      {
		      elstack.Append (other);
		      }
		      }
		      else
		      if (tempels.Get(ei).NB1(j))
		      cerr << "different2 !!" << endl;
		    */

		  }
	      }
	  }
      }



    // check outer elements
    if (debugparam.slowchecks)
      {
	for (int i = 1; i <= ne; i++)
	  {
	    const DelaunayTet & el = tempels.Get(i);
	    const Point3d & p1 = mesh.Point (el[0]);
	    const Point3d & p2 = mesh.Point (el[1]);
	    const Point3d & p3 = mesh.Point (el[2]);
	    const Point3d & p4 = mesh.Point (el[3]);
	  
	    Point3d ci = Center (p1, p2, p3, p4);
	  
	    //       if (adfront->Inside (ci) != adfront->Inside (Center (ci, p1)))
	    // 	cout << "ERROR: outer test unclear !!!" << endl;	
	  
	    if (inner.Test(i) != adfront->Inside (ci))
	      {
		/*
		  cout << "ERROR: outer test wrong !!!" 
		  << "inner = " << int(inner.Test(i))
		  << "outer = " << int(outer.Test(i))
		  << endl;
	      
		  cout << "Vol = " << Determinant(Vec3d(p1, p2),
		  Vec3d(p1, p3),
		  Vec3d(p1, p4)) << endl;
	      
		*/	      
		for (int j = 1; j <= 4; j++)
		  {
		    Point3d hp;
		    switch (j)
		      {
		      case 1: hp = Center (ci, p1); break;
		      case 2: hp = Center (ci, p2); break;
		      case 3: hp = Center (ci, p3); break;
		      case 4: hp = Center (ci, p4); break;
		      }
		    //		  cout << "inside(" << hp << ") = " << adfront->Inside(hp) << endl;
		  }
	      
	      }
	  
	    if (adfront->Inside(ci))
	      outer.Clear(i);
	    else
	      outer.Set(i);
	  }
      }


    /*

    // find bug in innerfaces

    tempmesh.DeleteVolumeElements();

    for (i = 1; i <= innerfaces.GetNBags(); i++)
    for (j = 1; j <= innerfaces.GetBagSize(i); j++)
    {
    INDEX_3 i3;
    INDEX_2 i2;
    innerfaces.GetData (i, j, i3, i2);
    if (i2.I2())
    {
    if (outer.Test(i2.I1()) != outer.Test(i2.I2()))
    {
    tempmesh.AddVolumeElement (tempels.Get(i2.I1()));
    tempmesh.AddVolumeElement (tempels.Get(i2.I2()));
    cerr << "outer flag different for connected els" << endl;
    }
    }
    }


    cout << "Check intersectiong once more" << endl;

    for (i = 1; i <= openels.Size(); i++)
    {
    tempmesh.SurfaceElement(2*openels.Get(i)).SetIndex(2);
    tempmesh.SurfaceElement(2*openels.Get(i)-1).SetIndex(2);
    }

    //  for (i = 1; i <= tempmesh.GetNE(); i++)
    //    for (j = 1; j <= tempmesh.GetNSE(); j++)
    i = 6; j = 403;
    if (i <= tempmesh.GetNE() && j <= tempmesh.GetNSE())
    if (tempmesh.SurfaceElement(j).GetIndex()==2)
    {
    const Element & el = tempmesh.VolumeElement(i);
    const Element2d & sel = tempmesh.SurfaceElement(j);

    const Point3d *tripp[3];
    const Point3d *pp[4];
    int tetpi[4], tripi[3];

    for (k = 1; k <= 4; k++)
    {
    pp[k-1] = &tempmesh.Point(el.PNum(k));
    tetpi[k-1] = el.PNum(k);
    }

    for (k = 1; k <= 3; k++)
    {
    tripp[k-1] = &tempmesh.Point (sel.PNum(k));
    tripi[k-1] = sel.PNum(k);
    }

    (*testout) << "Check Triangle " << j << ":";
    for (k = 1; k <= 3; k++)
    (*testout) << " " << sel.PNum(k);
    for (k = 1; k <= 3; k++)
    (*testout) << " " << tempmesh.Point(sel.PNum(k));
    (*testout) << endl;

    (*testout) << "Check Tet " << i << ":";
    for (k = 1; k <= 4; k++)
    (*testout) << " " << el.PNum(k);
    for (k = 1; k <= 4; k++)
    (*testout) << " " << tempmesh.Point(el.PNum(k));
    (*testout) << endl;

    if (IntersectTetTriangle (&pp[0], &tripp[0], tetpi, tripi))
    {
    cout << "Intesection detected !!" << endl;
    }
    }

    tempmesh.Save ("temp.vol");

    // end bug search
    */


    for (int i = ne; i >= 1; i--)
      {
	if (outer.Test(i))
	  tempels.DeleteElement(i);
      }


    // mesh.points.SetSize(mesh.points.Size()-4);

    PrintMessage (5, "outer removed");

  }



  void Meshing3 :: Delaunay (Mesh & mesh, int domainnr, const MeshingParameters & mp)
  {
    static Timer t("Meshing3::Delaunay"); RegionTimer reg(t);
    
    int np, ne;

    PrintMessage (1, "Delaunay meshing");
    PrintMessage (3, "number of points: ", mesh.GetNP());
    PushStatus ("Delaunay meshing");


    NgArray<DelaunayTet> tempels;
    Point3d pmin, pmax;

    DelaunayTet startel;

    int oldnp = mesh.GetNP();
    if (mp.blockfill)
      {
	BlockFillLocalH (mesh, mp);
	PrintMessage (3, "number of points: ", mesh.GetNP());
      }

    np = mesh.GetNP();

    Delaunay1 (mesh, mp, adfront, tempels, oldnp, startel, pmin, pmax);

    {
      // improve delaunay - mesh by swapping !!!!

      Mesh tempmesh;
      tempmesh.GetMemoryTracer().SetName("delaunay-tempmesh");

      for (auto & meshpoint : mesh.Points())
        tempmesh.AddPoint (meshpoint);
      
      for (auto & tempel : tempels)
	{   
	  Element el(4);
	  for (int j = 0; j < 4; j++)
            el[j] = tempel[j];

	  el.SetIndex (1);

	  const Point<3> & lp1 = mesh.Point (el[0]);
	  const Point<3> & lp2 = mesh.Point (el[1]);
	  const Point<3> & lp3 = mesh.Point (el[2]);
	  const Point<3> & lp4 = mesh.Point (el[3]);
	  Vec<3> v1 = lp2-lp1;
	  Vec<3> v2 = lp3-lp1;
	  Vec<3> v3 = lp4-lp1;

	  Vec<3> n = Cross (v1, v2);
	  double vol = n * v3;
	  if (vol > 0) swap (el[2], el[3]);

	  tempmesh.AddVolumeElement (el);
	}

      tempels.DeleteAll();

      MeshQuality3d (tempmesh);

      tempmesh.AddFaceDescriptor (FaceDescriptor (1, 1, 0, 0));
      tempmesh.AddFaceDescriptor (FaceDescriptor (2, 1, 0, 0));


    
      for (int i = 1; i <= mesh.GetNOpenElements(); i++)
	{
	  Element2d sel = mesh.OpenElement(i);
	  sel.SetIndex(1);
	  tempmesh.AddSurfaceElement (sel);
	  swap (sel[1], sel[2]);
	  tempmesh.AddSurfaceElement (sel);
	}


      for (int i = 1; i <= 4; i++)
	{
	  Element2d self(TRIG);
	  self.SetIndex (1);
	  startel.GetFace (i-1, self);
	  tempmesh.AddSurfaceElement (self);
	}

      
      //  for (i = mesh.GetNP() - 3; i <= mesh.GetNP(); i++)
      //    tempmesh.AddLockedPoint (i);
      for (auto pi : tempmesh.Points().Range())
        tempmesh.AddLockedPoint (pi);
      
      //  tempmesh.PrintMemInfo(cout);
      // tempmesh.Save ("tempmesh.vol");

      {
        RegionTaskManager rtm(mp.parallel_meshing ? mp.nthreads : 0);
        for (int i = 1; i <= 4; i++)
          {
            tempmesh.FindOpenElements ();

            PrintMessage (5, "Num open: ", tempmesh.GetNOpenElements());
            tempmesh.CalcSurfacesOfNode ();

            tempmesh.FreeOpenElementsEnvironment (1);

            MeshOptimize3d meshopt(mp);
            // tempmesh.CalcSurfacesOfNode();
            meshopt.SwapImprove(tempmesh, OPT_CONFORM);
          }
      }
    
      MeshQuality3d (tempmesh);
    
      tempels.SetSize(tempmesh.GetNE());
      tempels.SetSize(0);
      for (auto & el : tempmesh.VolumeElements())
        tempels.Append (el);
    }

    DelaunayRemoveDegenerated(mesh.Points(), tempels, np);

    NgArray<int> openels;
    DelaunayRemoveTwoTriaTets(mesh, tempels, openels);
    DelaunayRemoveIntersecting(mesh, tempels, openels, pmin, pmax);
    DelaunayRemoveOuter(mesh, tempels, adfront);

    for (int i = 0; i < tempels.Size(); i++)
      {
	Element el(4);
	for (int j = 0; j < 4; j++)
	  el[j] = tempels[i][j];
	mesh.AddVolumeElement (el);
      }

    mesh.FindOpenElements(domainnr);

    mesh.Compress();
    PopStatus ();
  }
}

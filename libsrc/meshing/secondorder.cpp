 #include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{

  void Refinement :: MakeSecondOrder (Mesh & mesh) const
  {
    const_cast<Refinement&> (*this).MakeSecondOrder(mesh);
  }

  
  void Refinement :: MakeSecondOrder (Mesh & mesh)
  {
    /*
      Berlin, 2014: if we have curved surface elements, keep them !
    */

    mesh.ComputeNVertices();
    // mesh.SetNP(mesh.GetNV());
    mesh.SetNP(mesh.GetNP());  // setup multilevel-table

    INDEX_2_HASHTABLE<PointIndex> between(mesh.GetNP() + 5);

    for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
      {
        auto & seg = mesh[si];
        if (seg.GetType() == SEGMENT3)
          between.Set(INDEX_2::Sort(seg[0],seg[1]), seg[2]);
      }
    
    
    for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
      {
	const Element2d & el = mesh[sei];

	static int betw_trig[3][3] =
	  { { 1, 2, 3 }, { 0, 2, 4 }, { 0, 1, 5 } };
	static int betw_quad6[2][3] =
	  { { 0, 1, 4 }, { 3, 2, 5 } };
	static int betw_quad8[4][3] =
	  { { 0, 1, 4 }, { 3, 2, 5 },
	    { 0, 3, 6 }, { 1, 2, 7 } };
      
	int onp = 0;
	int (*betw)[3] = NULL;
	switch (el.GetType())
	  {
	  case TRIG6:
	    {
	      betw = betw_trig;
	      onp = 3;
	      break;
	    }
	  case QUAD6: 
            {
              betw = betw_quad6;
              onp = 4;
              break;
            }
	  case QUAD8: 
            {
              betw = betw_quad8;
              onp = 4;
              break;
            }
          default:
            ;
          }
        
        if (betw)
          for (int j = 0; j < el.GetNP()-onp; j++)
            {
              int pi1 = el[betw[j][0]];
              int pi2 = el[betw[j][1]];
              INDEX_2 i2 = INDEX_2::Sort (pi1, pi2);
              between.Set (i2, el[onp+j]);
            }
      }
      

    bool thinlayers = 0;
    for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      if (mesh[ei].GetType() == PRISM ||
	  mesh[ei].GetType() == PRISM12)
	thinlayers = 1;
    

    int nseg = mesh.GetNSeg();
    for (SegmentIndex si = 0; si < nseg; si++)
      {
	Segment & el = mesh.LineSegment(si);

	INDEX_2 i2 = INDEX_2::Sort (el[0], el[1]);

	if (between.Used(i2))
	  el[2] = between.Get(i2);
	else
	  {
	    Point<3> pb;
	    EdgePointGeomInfo ngi;
            geo.PointBetweenEdge(mesh.Point (el[0]),
                                 mesh.Point (el[1]), 0.5,
                                 el.surfnr1, el.surfnr2,
                                 el.epgeominfo[0], el.epgeominfo[1],
                                 pb, ngi);
	  
	    el[2] = mesh.AddPoint (pb, mesh.Point(el[0]).GetLayer(), 
				   EDGEPOINT);
	    between.Set (i2, el[2]);
	  }
      }

    // refine surface elements
    for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
      {
	const Element2d & el = mesh[sei];

	int onp = 0;
      
	Element2d newel(TRIG);
	newel.SetIndex (el.GetIndex());

	static int betw_trig[3][3] =
	  { { 1, 2, 3 },
	    { 0, 2, 4 },
	    { 0, 1, 5 } };
	static int betw_quad6[2][3] =
	  { { 0, 1, 4 },
	    { 3, 2, 5 } };
	static int betw_quad8[4][3] =
	  { { 0, 1, 4 },
	    { 3, 2, 5 },
	    { 0, 3, 6 },
	    { 1, 2, 7 } };
	int (*betw)[3] = NULL;
      
	switch (el.GetType())
	  {
	  case TRIG:
	  case TRIG6:
	    {
	      betw = betw_trig;
	      newel.SetType (TRIG6);
	      onp = 3;
	      break;
	    }
	  case QUAD:
	  case QUAD6: 
	  case QUAD8:
	    {
	      if (thinlayers)
		{
		  betw = betw_quad6;
		  newel.SetType (QUAD6);
		}
	      else
		{
		  betw = betw_quad8;
		  newel.SetType (QUAD8);
		}
	      onp = 4;
	      break;
	    }
	  default:
	    PrintSysError ("Unhandled element in secondorder:", int(el.GetType()));
	  }

	for (int j = 0; j < onp; j++)
	  newel[j] = el[j];
      
	int nnp = newel.GetNP();
	for (int j = 0; j < nnp-onp; j++)
	  {
	    int pi1 = newel[betw[j][0]];
	    int pi2 = newel[betw[j][1]];
	  
	    INDEX_2 i2 = INDEX_2::Sort (pi1, pi2);
	  
	    if (between.Used(i2))
	      newel[onp+j] = between.Get(i2);
	    else
	      {
		Point<3> pb;
		PointGeomInfo newgi;
		geo.PointBetween(mesh.Point (pi1),
                                 mesh.Point (pi2), 0.5, 
                                 mesh.GetFaceDescriptor(el.GetIndex ()).SurfNr(),
                                 el.GeomInfoPi (betw[j][0]+1),
                                 el.GeomInfoPi (betw[j][1]+1),
                                 pb, newgi);

		newel[onp+j] = mesh.AddPoint (pb, mesh.Point(pi1).GetLayer(), 
					      SURFACEPOINT);
		between.Set (i2, newel[onp+j]);
	      }
	  }
      
	mesh[sei] = newel;
      }

 
    //    int i, j;



    // refine volume elements
    for (int i = 1; i <= mesh.GetNE(); i++)
      {
	const Element & el = mesh.VolumeElement(i);
	int onp = 0;

	Element newel(TET);
	newel.SetIndex (el.GetIndex());

	static int betw_tet[6][3] =
	  { { 0, 1, 4 },
	    { 0, 2, 5 },
	    { 0, 3, 6 },
	    { 1, 2, 7 },
	    { 1, 3, 8 },
	    { 2, 3, 9 } };
	static int betw_prism[6][3] =
	  {
	    { 0, 2, 6 },
	    { 0, 1, 7 },
	    { 1, 2, 8 },
	    { 3, 5, 9 },
	    { 3, 4, 10 },
	    { 4, 5, 11 },
	  };
        static int betw_prism15[9][3] =
          {
           { 0, 1, 6 },
           { 0, 2, 7 },
           { 1, 2, 8 },
           { 0, 3, 9 },
           { 1, 4, 10 },
           { 2, 5, 11 },
           { 3, 4, 12 },
           { 3, 5, 13 },
           { 4, 5, 14 }
          };
        static int betw_pyramid[8][3] =
          {
           { 0, 1, 5 },
           { 3, 2, 6 },
           { 3, 0, 7 },
           { 1, 2, 8 },
           { 0, 4, 9 },
           { 1, 4, 10 },
           { 2, 4, 11 },
           { 3, 4, 12 }
          };
	static int betw_hex[12][3] =
	  {
            { 0, 1, 8 },
            { 2, 3, 9 },
            { 3, 0, 10 },
            { 1, 2, 11 },
            { 4, 5, 12 },
            { 6, 7, 13 },
            { 7, 4, 14 },
            { 5, 6, 15 },
            { 0, 4, 16 },
            { 1, 5, 17 },
            { 2, 6, 18 },
            { 3, 7, 19 },
          };
        
	int (*betw)[3] = NULL;

	switch (el.GetType())
	  {
	  case TET:
	  case TET10:
	    {
	      betw = betw_tet;
	      newel.SetType (TET10);
	      onp = 4;
	      break;
	    }
	  case PRISM:
	  case PRISM12:
	    {
	      betw = betw_prism;
	      newel.SetType (PRISM12);
	      onp = 6;
	      break;
	    }
          case PRISM15:
            {
              betw = betw_prism15;
              newel.SetType(PRISM15);
              onp = 6;
              break;
            }
          case PYRAMID:
          case PYRAMID13:
            {
              betw = betw_pyramid;
              newel.SetType(PYRAMID13);
              onp = 5;
              break;
            }
          case HEX:
          case HEX20:
	    {
	      betw = betw_hex;
	      newel.SetType (HEX20);
	      onp = 8;
	      break;
	    }
	  default:
	    PrintSysError ("MakeSecondOrder, illegal vol type ", int(el.GetType()));
	  }


	for (int j = 1; j <= onp; j++)
	  newel.PNum(j) = el.PNum(j);
	int nnp = newel.GetNP();

	for (int j = 0; j < nnp-onp; j++)
	  {
	    INDEX_2 i2(newel[betw[j][0]],
		       newel[betw[j][1]]);
	    i2.Sort();
	  
	    if (between.Used(i2))
	      newel.PNum(onp+1+j) = between.Get(i2);
	    else
	      {
		newel.PNum(onp+1+j) = mesh.AddPoint
		  (Center (mesh.Point(i2.I1()),
			   mesh.Point(i2.I2())),
		   mesh.Point(i2.I1()).GetLayer(), 
		   INNERPOINT);

		between.Set (i2, newel.PNum(onp+1+j));
	      }
	  }

	mesh.VolumeElement (i) = newel;
      }


    // makes problems after linear mesh refinement, since
    // 2nd order identifications are not removed
    // update identification tables
    for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
      {
	NgArray<int,PointIndex::BASE> identmap;
	mesh.GetIdentifications().GetMap (i, identmap);

	for (INDEX_2_HASHTABLE<PointIndex>::Iterator it = between.Begin();
	     it != between.End(); it++)
	  {
	      INDEX_2 i2;
	      PointIndex newpi;
	      between.GetData (it, i2, newpi);
	      INDEX_2 oi2(identmap.Get(i2.I1()),
			  identmap.Get(i2.I2()));
	      oi2.Sort();
	      if (between.Used (oi2))
		{
		  PointIndex onewpi = between.Get(oi2);
		  mesh.GetIdentifications().Add (newpi, onewpi, i);
		}
	  }

	/*
	for (int j = 1; j <= between.GetNBags(); j++)
	  for (int k = 1; k <= between.GetBagSize(j); k++)
	    {
	      INDEX_2 i2;
	      int newpi;
	      between.GetData (j, k, i2, newpi);
	      INDEX_2 oi2(identmap.Get(i2.I1()),
			  identmap.Get(i2.I2()));
	      oi2.Sort();
	      if (between.Used (oi2))
		{
		  int onewpi = between.Get(oi2);
		  mesh.GetIdentifications().Add (newpi, onewpi, i);
		}
	    }
	*/
      }


    //  mesh.mglevels++;
    int oldsize = mesh.mlbetweennodes.Size();
    mesh.mlbetweennodes.SetSize(mesh.GetNP());
    for (int i = oldsize; i < mesh.GetNP(); i++)
      mesh.mlbetweennodes[i] = INDEX_2(0,0);

    /*
    for (i = 1; i <= between.GetNBags(); i++)
      for (j = 1; j <= between.GetBagSize(i); j++)
	{
	  INDEX_2 oldp;
	  int newp;
	  between.GetData (i, j, oldp, newp);
	  mesh.mlbetweennodes.Elem(newp) = oldp;
	}
    */

    for (INDEX_2_HASHTABLE<PointIndex>::Iterator it = between.Begin();
	 it != between.End(); it++)
      {
	mesh.mlbetweennodes[between.GetData (it)] = between.GetHash(it);
      }

    mesh.ComputeNVertices();
    mesh.RebuildSurfaceElementLists();
    //  ValidateSecondOrder (mesh);
  }


  void Refinement :: ValidateSecondOrder (Mesh & mesh)
  {
    PrintMessage (3, "Validate mesh");
    int np = mesh.GetNP();
    int ne = mesh.GetNE();
    // int i, j;
    NgArray<INDEX_2> parents(np);
  
    for (int i = 1; i <= np; i++)
      parents.Elem(i) = INDEX_2(0,0);

    for (int i = 1; i <= ne; i++)
      {
	const Element & el = mesh.VolumeElement(i);
	if (el.GetType() == TET10)
	  {
	    static int betweentab[6][3] =
	      { { 1, 2, 5 },
		{ 1, 3, 6 },
		{ 1, 4, 7 },
		{ 2, 3, 8 },
		{ 2, 4, 9 },
		{ 3, 4, 10 } };
	    for (int j = 0; j < 6; j++)
	      {
		int f1 = el.PNum (betweentab[j][0]);
		int f2 = el.PNum (betweentab[j][1]);
		int son = el.PNum (betweentab[j][2]);
		parents.Elem(son).I1() = f1;
		parents.Elem(son).I2() = f2;
	      }
	  }
      }

    ValidateRefinedMesh (mesh, parents);
  }


  void Refinement ::
  ValidateRefinedMesh (Mesh & mesh, 
		       NgArray<INDEX_2> & parents)
  {
    // int i, j, k;
  
    // homotopy method

    int ne = mesh.GetNE();

    int cnttrials = 100;
    int wrongels = 0;
    for (int i = 1; i <= ne; i++)
      if (mesh.VolumeElement(i).CalcJacobianBadness (mesh.Points()) > 1e10)
	{
	  wrongels++;
	  mesh.VolumeElement(i).flags.badel = 1;
	}
      else
	mesh.VolumeElement(i).flags.badel = 0;

    double facok = 0;
    double factry;

    NgBitArray illegalels(ne);
    illegalels.Clear();

      
    if (wrongels)
      {
	cout << "WARNING: " << wrongels << " illegal element(s) found" << endl;

	int np = mesh.GetNP();
	NgArray<Point<3> > should(np);
	NgArray<Point<3> > can(np);

	for (int i = 1; i <= np; i++)
	  {
	    should.Elem(i) = can.Elem(i) = mesh.Point(i);
	  }

	for (int i = 1; i <= parents.Size(); i++)
	  {
	    if (parents.Get(i).I1())
	      can.Elem(i) = Center (can.Elem(parents.Get(i).I1()),
				    can.Elem(parents.Get(i).I2()));
	  }

	NgBitArray boundp(np);
	boundp.Clear();
	for (int i = 1; i <= mesh.GetNSE(); i++)
	  {
	    const Element2d & sel = mesh.SurfaceElement(i);
	    for (int j = 1; j <= sel.GetNP(); j++)
	      boundp.Set(sel.PNum(j));
	  }


	(*testout) << "bpoints:" << endl;
	for (int i = 1; i <= np; i++)
	  if (boundp.Test(i))
	    (*testout) << i << endl;

	double lam = 0.5;

	while (facok < 1-1e-8 && cnttrials > 0)
	  {
	    lam *= 4;
	    if (lam > 2) lam = 2;

	    do
	      {
		//	      cout << "trials: " << cnttrials << endl;
		lam *= 0.5;
		cnttrials--;

		cout << "lam = " << lam << endl;

		factry = lam + (1-lam) * facok;
		cout << "trying: " << factry << endl;

		for (int i = 1; i <= np; i++)
		  if (boundp.Test(i))
		    {
		      for (int j = 0; j < 3; j++)
			mesh.Point(i)(j) = 
			  lam * should.Get(i)(j) +
			  (1-lam) * can.Get(i)(j);
		    }
		  else
		    mesh.Point(i) = Point<3> (can.Get(i));
	      
		//	      (*testout) << "bad els: " << endl;
		wrongels = 0;
		for (int i = 1; i <= ne; i++)
		  {
		    if (!illegalels.Test(i) && 
			mesh.VolumeElement(i).
			CalcJacobianBadness(mesh.Points()) > 1e10)
		      {
			wrongels++;
			Element & el = mesh.VolumeElement(i);
			el.flags.badel = 1;
		     
		      
			if (lam < 1e-4)
			  illegalels.Set(i);
 

			/*
			  (*testout) << i << ": ";
			  for (j = 1; j <= el.GetNP(); j++)
			  (*testout) << el.PNum(j) << " ";
			  (*testout) << endl;
			*/
		      }
		    else
		      mesh.VolumeElement(i).flags.badel = 0;
		  }
		cout << "wrongels = " << wrongels << endl;
	      }
	    while (wrongels && cnttrials > 0);

	    mesh.CalcSurfacesOfNode();
	    MeshingParameters dummymp;
	    mesh.ImproveMeshJacobian (dummymp, OPT_WORSTCASE);	      
	  
	    facok = factry;
	    for (int i = 1; i <= np; i++)
	      can.Elem(i) = mesh.Point(i);
	  }
      }


      
    for (int i = 1; i <= ne; i++)
      {
	if (illegalels.Test(i))
	  {
	    cout << "illegal element: " << i << endl;
	    mesh.VolumeElement(i).flags.badel = 1;	
	  }
	else
	  mesh.VolumeElement(i).flags.badel = 0;	
      }
  
    /*
      if (cnttrials <= 0)
      {
      cerr << "ERROR: Sorry, illegal elements:" << endl;
      }
    */
  }

}

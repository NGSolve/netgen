#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{

  void Refinement :: Refine (Mesh & mesh) const
  {
    const_cast<Refinement&> (*this).Refine(mesh);
  }


  void Refinement :: Refine (Mesh & mesh)
  {
    PrintMessage (3, "Refine mesh");

    mesh.SetNextMajorTimeStamp();
    
    if (ntasks > 1 && id == 0)
      return;


    // reduce 2nd order
    mesh.ComputeNVertices();
    mesh.SetNP(mesh.GetNV());

    if (mesh.mlbetweennodes.Size() < mesh.GetNV())
      {
        mesh.mlbetweennodes.SetSize(mesh.GetNV());
        mesh.mlbetweennodes = INDEX_2(PointIndex::BASE-1,PointIndex::BASE-1);
      }

    
    INDEX_2_HASHTABLE<PointIndex> between(mesh.GetNP() + 5);


    // new version with consistent ordering across sub-domains

    NgArray<INDEX_2> parents;
    for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
      {
	const Segment & el = mesh[si];
	INDEX_2 i2 = INDEX_2::Sort(el[0], el[1]);
        if (!between.Used(i2))
          {
            between.Set (i2, 0);          
            parents.Append(i2);
          }
      }
    for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
      {
	const Element2d & el = mesh[sei];
	switch (el.GetType())
	  {
	  case TRIG:
	  case TRIG6:
	    {
              static int betw[3][3] =
		{ { 1, 2, 3 },
		  { 0, 2, 4 },
		  { 0, 1, 5 } };
              for (int j = 0; j < 3; j++)
                {
                  auto i2 = PointIndices<2>::Sort(el[betw[j][0]],el[betw[j][1]]);
                  if (!between.Used(i2))
                    {
                      between.Set (i2, 0);          
                      parents.Append(i2);
                    }
                }
              break;
            }
          default:
            throw NgException ("currently refinement for quad-elements is not supported");
          }
      }
    for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      {
	const Element & el = mesh[ei];
	switch (el.GetType())
	  {
	  case TET:
	  case TET10:
	    {
	      static int betw[6][3] =
                { { 1, 2, 5 },
                  { 1, 3, 6 },
                  { 1, 4, 7 },
                  { 2, 3, 8 },
                  { 2, 4, 9 },
                  { 3, 4, 10 } };

              for (int j = 0; j < 6; j++)
                {
                  INDEX_2 i2 = INDEX_2::Sort(el.PNum(betw[j][0]),el.PNum(betw[j][1]));
                  if (!between.Used(i2))
                    {
                      between.Set (i2, 0);          
                      parents.Append(i2);
                    }
                }
              break;
            }
          default:
            throw NgException ("currently refinement for non-tet elements is not supported");
          }
      }

    PrintMessage (5, "have points");
    
    NgArray<int> par_nr(parents.Size());
    for (int i = 0; i < par_nr.Size(); i++)
      par_nr[i] = i;
    QuickSort (parents, par_nr);
    mesh.mlbetweennodes.SetSize(mesh.GetNV()+parents.Size());
    for (int i = 0; i < parents.Size(); i++)
      {
        between.Set (parents[i], mesh.GetNV()+i+PointIndex::BASE);
        mesh.mlbetweennodes[mesh.GetNV()+i+PointIndex::BASE] = parents[i];
      }
    
    mesh.SetNP(mesh.GetNV() + parents.Size());
    NgArray<bool, PointIndex::BASE> pointset(mesh.GetNP());
    pointset = false;
    
    PrintMessage (5, "sorting complete");
    
    // refine edges
    NgArray<EdgePointGeomInfo,PointIndex::BASE> epgi;

    int oldns = mesh.GetNSeg();
    for (SegmentIndex si = 0; si < oldns; si++)
      {
	const Segment & el = mesh.LineSegment(si);

	INDEX_2 i2 = INDEX_2::Sort(el[0], el[1]);
	PointIndex pinew = between.Get(i2);
	EdgePointGeomInfo ngi;

	if (pointset[pinew])
	  {
	    // pinew = between.Get(i2);
	    ngi = epgi[pinew]; 
	  }
	else
	  {
            pointset[pinew] = true;
	    Point<3> pnew;
	    geo.PointBetweenEdge(mesh.Point (el[0]),
                                 mesh.Point (el[1]), 0.5,
                                 el.surfnr1, el.surfnr2,
                                 el.epgeominfo[0], el.epgeominfo[1],
                                 pnew, ngi);

	    // pinew = mesh.AddPoint (pnew);
            mesh.Point(pinew) = pnew;
	    // between.Set (i2, pinew);

	    if (pinew >= epgi.Size()+PointIndex::BASE)
	      epgi.SetSize (pinew+1-PointIndex::BASE);
	    epgi[pinew] = ngi;
	  }

	Segment ns1 = el;
	Segment ns2 = el;
	ns1[1] = pinew;
	ns1.epgeominfo[1] = ngi;
	ns2[0] = pinew;
	ns2.epgeominfo[0] = ngi;

	mesh.LineSegment(si) = ns1;
	mesh.AddSegment (ns2);
      }

    PrintMessage (5, "have 1d elements");
    
    // refine surface elements
    NgArray<PointGeomInfo,PointIndex::BASE> surfgi (8*mesh.GetNP());
    for (int i = PointIndex::BASE;
	 i < surfgi.Size()+PointIndex::BASE; i++)
      surfgi[i].trignum = -1;


    int oldnf = mesh.GetNSE();
    for (SurfaceElementIndex sei = 0; sei < oldnf; sei++)
      {
	const Element2d & el = mesh[sei];

	switch (el.GetType())
	  {
	  case TRIG:
	  case TRIG6:
	    {
	      NgArrayMem<PointIndex,6> pnums(6);
	      NgArrayMem<PointGeomInfo,6> pgis(6);

	      static int betw[3][3] =
		{ { 2, 3, 4 },
		  { 1, 3, 5 },
		  { 1, 2, 6 } };

	      for (int j = 1; j <= 3; j++)
		{
		  pnums.Elem(j) = el.PNum(j);
		  pgis.Elem(j) = el.GeomInfoPi(j);
		}

	      for (int j = 0; j < 3; j++)
		{
		  PointIndex pi1 = pnums.Elem(betw[j][0]);
		  PointIndex pi2 = pnums.Elem(betw[j][1]);

		  INDEX_2 i2 (pi1, pi2);
		  i2.Sort();

		  Point<3> pb;
		  PointGeomInfo pgi;
		  geo.PointBetween(mesh.Point (pi1),
                                   mesh.Point (pi2), 0.5,
                                   mesh.GetFaceDescriptor(el.GetIndex ()).SurfNr(),
                                   el.GeomInfoPi (betw[j][0]),
                                   el.GeomInfoPi (betw[j][1]),
                                   pb, pgi);


		  pgis.Elem(4+j) = pgi;
                  PointIndex pinew = between.Get(i2); 
                  pnums.Elem(4+j) = pinew;
                  if (!pointset[pinew])
                    {
                      pointset[pinew] = true;
                      mesh.Point(pinew) = pb;                      
                    }
                  /*
		  if (between.Used(i2))
		    pnums.Elem(4+j) = between.Get(i2);
		  else
		    {
		      pnums.Elem(4+j) = mesh.AddPoint (pb);
		      between.Set (i2, pnums.Get(4+j));
		    }
                  */
		  if (surfgi.Size() < pnums.Elem(4+j))
		    surfgi.SetSize (pnums.Elem(4+j));
		  surfgi.Elem(pnums.Elem(4+j)) = pgis.Elem(4+j);
		}


	      static int reftab[4][3] =
		{ { 1, 6, 5 },
		  { 2, 4, 6 },
		  { 3, 5, 4 },
		  { 6, 4, 5 } };

	      int ind = el.GetIndex();
	      for (int j = 0; j < 4; j++)
		{
		  Element2d nel(TRIG);
		  for (int k = 1; k <= 3; k++)
		    {
		      nel.PNum(k) = pnums.Get(reftab[j][k-1]);
		      nel.GeomInfoPi(k) = pgis.Get(reftab[j][k-1]);
		    }
		  nel.SetIndex(ind);

		  if (j == 0)
		    mesh[sei] = nel;
		  else
		    mesh.AddSurfaceElement(nel);
		}
	      break;
	    }
	  case QUAD:
	  case QUAD6:
	  case QUAD8:
	    {
	      NgArrayMem<PointIndex,9> pnums(9);
	      NgArrayMem<PointGeomInfo,9> pgis(9);

	      static int betw[5][3] =
		{ { 1, 2, 5 },
		  { 2, 3, 6 },
		  { 3, 4, 7 },
		  { 1, 4, 8 },
		  { 5, 7, 9 } };

	      for (int j = 1; j <= 4; j++)
		{
		  pnums.Elem(j) = el.PNum(j);
		  pgis.Elem(j) = el.GeomInfoPi(j);
		}

	      for (int j = 0; j < 5; j++)
		{
		  int pi1 = pnums.Elem(betw[j][0]);
		  int pi2 = pnums.Elem(betw[j][1]);

		  INDEX_2 i2 (pi1, pi2);
		  i2.Sort();

		  if (between.Used(i2))
		    {
		      pnums.Elem(5+j) = between.Get(i2);
		      pgis.Elem(5+j) = surfgi.Get(pnums.Elem(4+j));
		    }
		  else
		    {
		      Point<3> pb;
		      geo.PointBetween(mesh.Point (pi1),
                                       mesh.Point (pi2), 0.5,
                                       mesh.GetFaceDescriptor(el.GetIndex ()).SurfNr(),
                                       el.GeomInfoPi (betw[j][0]),
                                       el.GeomInfoPi (betw[j][1]),
                                       pb, pgis.Elem(5+j));

		      pnums.Elem(5+j) = mesh.AddPoint (pb);

		      between.Set (i2, pnums.Get(5+j));
		      
		      if (surfgi.Size() < pnums.Elem(5+j))
			surfgi.SetSize (pnums.Elem(5+j));
		      surfgi.Elem(pnums.Elem(5+j)) = pgis.Elem(5+j);
		    }
		}

	      static int reftab[4][4] =
		{
		  { 1, 5, 9, 8 },
		  { 5, 2, 6, 9 },
		  { 8, 9, 7, 4 },
		  { 9, 6, 3, 7 } };

	      int ind = el.GetIndex();
	      for (int j = 0; j < 4; j++)
		{
		  Element2d nel(QUAD);
		  for (int k = 1; k <= 4; k++)
		    {
		      nel.PNum(k) = pnums.Get(reftab[j][k-1]);
		      nel.GeomInfoPi(k) = pgis.Get(reftab[j][k-1]);
		    }
		  nel.SetIndex(ind);

		  if (j == 0)
		    mesh[sei] = nel;
		  else
		    mesh.AddSurfaceElement(nel);
		}
	      break;
	    }
	  default:
	    PrintSysError ("Refine: undefined surface element type ", int(el.GetType()));
	  }
      }

    PrintMessage (5, "have 2d elements");
    // cout << "id = " << id << ", ne = " << mesh.GetNE() << endl;
    // refine volume elements
    int oldne = mesh.GetNE();
    mesh.VolumeElements().SetAllocSize(8*oldne);
    for (ElementIndex ei = 0; ei < oldne; ei++)
      {
	const Element & el = mesh[ei];
	switch (el.GetType())
	  {
	  case TET:
	  case TET10:
	    {
	     NgArrayMem<PointIndex,10> pnums(10);
	     static int betw[6][3] =
	     { { 1, 2, 5 },
	       { 1, 3, 6 },
	       { 1, 4, 7 },
	       { 2, 3, 8 },
	       { 2, 4, 9 },
	       { 3, 4, 10 } };

	     int elrev = el.flags.reverse;

	     for (int j = 1; j <= 4; j++)
               pnums.Elem(j) = el.PNum(j);
	     if (elrev)
	     swap (pnums.Elem(3), pnums.Elem(4));

	     for (int j = 0; j < 6; j++)
               {
                 PointIndex pi1 = pnums.Get(betw[j][0]);
                 PointIndex pi2 = pnums.Get(betw[j][1]);
                 INDEX_2 i2 (pi1, pi2);
                 i2.Sort();

	       /*
	       if (between.Used(i2))
	          pnums.Elem(5+j) = between.Get(i2);
	       else
	       {
		  pnums.Elem(5+j) = mesh.AddPoint
		  (Center (mesh.Point(i2.I1()),
			   mesh.Point(i2.I2())));
		  between.Set (i2, pnums.Elem(5+j));
	       }
	       */
	       PointIndex pinew = between.Get(i2);
	       pnums.Elem(j+5) = pinew;
	       if (!pointset[pinew])
		 {
		   pointset[pinew] = true;
		   mesh.Point(pinew) = Center(mesh.Point(pi1),
					      mesh.Point(pi2));
		 }
	    }

	    static int reftab[8][4] =
	    { { 1, 5, 6, 7 },
	      { 5, 2, 8, 9 },
	      { 6, 8, 3, 10 },
	      { 7, 9, 10, 4 },
	      { 5, 6, 7, 9 },
	      { 5, 6, 9, 8 },
	      { 6, 7, 9, 10 },
	      { 6, 8, 10, 9 } };
	/*
	  { { 1, 5, 6, 7 },
	  { 5, 2, 8, 9 },
	  { 6, 8, 3, 10 },
	  { 7, 9, 10, 4 },
	  { 5, 6, 7, 9 },
	  { 5, 6, 8, 9 },
	  { 6, 7, 9, 10 },
	  { 6, 8, 9, 10 } };
	*/
	   static bool reverse[8] =
	   {
	      false, false, false, false, false, true, false, true
	   };

	   int ind = el.GetIndex();
	   for (int j = 0; j < 8; j++)
	   {
             Element nel(TET);
	      for (int k = 1; k <= 4; k++)
	        nel.PNum(k) = pnums.Get(reftab[j][k-1]);
	      nel.SetIndex(ind);
	      nel.flags.reverse = reverse[j];
	      if (elrev)
	      {
		nel.flags.reverse = !nel.flags.reverse;
		swap (nel.PNum(3), nel.PNum(4));
	      }

	      if (j == 0)
	        mesh.VolumeElement(ei) = nel;
	      else
	        mesh.AddVolumeElement (nel);
	    }
	    break;
          }
          case HEX:
          {
	     NgArrayMem<PointIndex,27> pnums(27);
	     static int betw[13][3] =
	     { { 1, 2, 9 },
	       { 3, 4, 10 },
	       { 4, 1, 11 },
               { 2, 3, 12 },
	       { 5, 6, 13 },
	       { 7, 8, 14 },
	       { 8, 5, 15 },
	       { 6, 7, 16 },
	       { 1, 5, 17 },
	       { 2, 6, 18 },
	       { 3, 7, 19 },
	       { 4, 8, 20 },
	       { 2, 8, 21 },
	       };

             /*
	     static int fbetw[12][3] =
	     { { 1, 3, 22 },
	       { 2, 4, 22 },
	       { 5, 7, 23 },
               { 6, 8, 23 },
	       { 1, 6, 24 },
	       { 2, 5, 24 },
	       { 2, 7, 25 },
	       { 3, 6, 25 },
	       { 3, 8, 26 },
	       { 4, 7, 26 },
	       { 1, 8, 27 },
	       { 4, 5, 27 },
	       };
             */
             
             // updated by anonymous supporter, donations please to Karo W.
             static int fbetw[12][3] =
               { { 11, 12, 22 },
                 { 9, 10, 22 },
                 { 13, 14, 23 },
                 { 15, 16, 23 },
                 { 9, 13, 24 },
                 { 17, 18, 24 },
                 { 12, 16, 25 },
                 { 18, 19, 25 },
                 { 19, 20, 26 },
                 { 10, 14, 26 },
                 { 11, 15, 27 },
                 { 17, 20, 27 },
               };

	     pnums = PointIndex(-1);

	     for (int j = 1; j <= 8; j++)
               pnums.Elem(j) = el.PNum(j);


	     for (int j = 0; j < 13; j++)
	     {
	       INDEX_2 i2;
	       i2.I1() = pnums.Get(betw[j][0]);
	       i2.I2() = pnums.Get(betw[j][1]);
	       i2.Sort();

	       if (between.Used(i2))
	          pnums.Elem(9+j) = between.Get(i2);
	       else
	       {
		  pnums.Elem(9+j) = mesh.AddPoint
		  (Center (mesh.Point(i2.I1()),
			   mesh.Point(i2.I2())));
		  between.Set (i2, pnums.Elem(9+j));
	       }
	    }

	    for (int j = 0; j < 6; j++)
	    {
	       INDEX_2 i2a, i2b;
	       i2a.I1() = pnums.Get(fbetw[2*j][0]);
	       i2a.I2() = pnums.Get(fbetw[2*j][1]);
	       i2a.Sort();
	       i2b.I1() = pnums.Get(fbetw[2*j+1][0]);
	       i2b.I2() = pnums.Get(fbetw[2*j+1][1]);
	       i2b.Sort();

	       if (between.Used(i2a))
		 pnums.Elem(22+j) = between.Get(i2a);
	       else if (between.Used(i2b))
		 pnums.Elem(22+j) = between.Get(i2b);
	       else
		 {
		   pnums.Elem(22+j) = mesh.AddPoint
		     (Center (mesh.Point(i2a.I1()),
			      mesh.Point(i2a.I2())));

		   between.Set (i2a, pnums.Elem(22+j));
		 }
	    }

	    static int reftab[8][8] =
	    { { 1, 9, 22, 11, 17, 24, 21, 27 },
	      { 9, 2, 12, 22, 24, 18, 25, 21 },
	      { 11, 22, 10, 4, 27, 21, 26, 20},
	      { 22, 12, 3, 10, 21, 25, 19, 26},
	      { 17, 24, 21, 27, 5, 13, 23, 15},
	      { 24, 18, 25, 21, 13, 6, 16, 23},
	      { 27, 21, 26, 20, 15, 23, 14, 8},
	      { 21, 25, 19, 26, 23, 16, 7, 14} };


	   int ind = el.GetIndex();
	   for (int j = 0; j < 8; j++)
	   {
	      Element nel(HEX);
	      for (int k = 1; k <= 8; k++)
	        nel.PNum(k) = pnums.Get(reftab[j][k-1]);
	      nel.SetIndex(ind);

              if (j == 0)
	        mesh.VolumeElement(ei) = nel;
	      else
	        mesh.AddVolumeElement (nel);
           }
           break;
	  }
	  case PRISM:
          {
	     NgArrayMem<PointIndex,18> pnums(18);
	     static int betw[9][3] =
	     { { 3, 1, 7 },
	       { 1, 2, 8 },
	       { 3, 2, 9 },
               { 6, 4, 10 },
	       { 4, 5, 11 },
	       { 6, 5, 12 },
	       { 1, 4, 13 },
	       { 3, 6, 14 },
	       { 2, 5, 15 },
	       };

// he: 15.jul 08, old version is wrong
//                produces double points ad quad faces and inconsistent mesh
// 	     static int fbetw[6][3] =
// 	     { { 1, 6, 16 },
// 	       { 3, 4, 16 },
// 	       { 1, 5, 17 },
//                { 2, 4, 17 },
// 	       { 2, 6, 18 },
// 	       { 3, 5, 18 },
// 	       };
           
           static int fbetw[6][3] =
           { { 7, 10, 16 },
           { 14, 13, 16 },
           { 11, 8, 17 },
           { 13, 15, 17 },
           { 12, 9, 18 },
           { 14, 15, 18 },
           };

	     //int elrev = el.flags.reverse;
           pnums = PointIndex(-1);
           
           for (int j = 1; j <= 6; j++)
	     pnums.Elem(j) = el.PNum(j);
	    // if (elrev)
	    // swap (pnums.Elem(3), pnums.Elem(4));

	     for (int j = 0; j < 9; j++)
	     {
	       INDEX_2 i2;
	       i2.I1() = pnums.Get(betw[j][0]);
	       i2.I2() = pnums.Get(betw[j][1]);
	       i2.Sort();

	       if (between.Used(i2))
	          pnums.Elem(7+j) = between.Get(i2);
	       else
	       {
		  pnums.Elem(7+j) = mesh.AddPoint
		  (Center (mesh.Point(i2.I1()),
			   mesh.Point(i2.I2())));
		  between.Set (i2, pnums.Elem(7+j));
	       }
	    }

	    for (int j = 0; j < 3; j++)
	    {
	       INDEX_2 i2a, i2b;
	       i2a.I1() = pnums.Get(fbetw[2*j][0]);
	       i2a.I2() = pnums.Get(fbetw[2*j][1]);
	       i2a.Sort();
	       i2b.I1() = pnums.Get(fbetw[2*j+1][0]);
	       i2b.I2() = pnums.Get(fbetw[2*j+1][1]);
	       i2b.Sort();

	       if (between.Used(i2a))
		 pnums.Elem(16+j) = between.Get(i2a);
	       else if (between.Used(i2b))
		 pnums.Elem(16+j) = between.Get(i2b);
	       else
		 {
		   pnums.Elem(16+j) = mesh.AddPoint
		     (Center (mesh.Point(i2a.I1()),
			      mesh.Point(i2a.I2())));

		   between.Set (i2a, pnums.Elem(16+j));
		 }
	    }


	    static int reftab[8][6] =
	    { { 1, 8, 7, 13, 17, 16 },
	      { 7, 8, 9, 16, 17, 18 },
	      { 7, 9, 3, 16, 18, 14 },
	      { 8, 2, 9, 17, 15, 18 },
	      { 13, 17, 16, 4, 11, 10 },
	      { 16, 17, 18, 10, 11, 12 },
	      { 16, 18, 14, 10, 12, 6 },
	      { 17, 15, 18, 11, 5, 12 } };


	   int ind = el.GetIndex();
	   for (int j = 0; j < 8; j++)
	   {
	      Element nel(PRISM);
	      for (int k = 1; k <= 6; k++)
	        nel.PNum(k) = pnums.Get(reftab[j][k-1]);
	      nel.SetIndex(ind);


	      //nel.flags.reverse = reverse[j];
	      //if (elrev)
	     // {
		//nel.flags.reverse = 1 - nel.flags.reverse;
		//swap (nel.PNum(3), nel.PNum(4));


	      if (j == 0)
	        mesh.VolumeElement(ei) = nel;
	      else
	        mesh.AddVolumeElement (nel);
           }
           break;
	  }
	  default:
	    PrintSysError ("Refine: undefined volume element type ", int(el.GetType()));
        }
      }


    // update identification tables
    for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
      {
	NgArray<int,PointIndex::BASE> identmap;
	mesh.GetIdentifications().GetMap (i, identmap);

	for (int j = 1; j <= between.GetNBags(); j++)
	  for (int k = 1; k <= between.GetBagSize(j); k++)
	    {
	      INDEX_2 i2;
	      PointIndex newpi;
	      between.GetData (j, k, i2, newpi);
	      INDEX_2 oi2(identmap.Get(i2.I1()),
			  identmap.Get(i2.I2()));
	      oi2.Sort();
	      if (between.Used (oi2))
		{
		  PointIndex onewpi = between.Get(oi2);
		  mesh.GetIdentifications().Add (newpi, onewpi, i);
		}
	    }

      }

    PrintMessage (5, "have 3d elements");
    mesh.ComputeNVertices();
    mesh.RebuildSurfaceElementLists();
    PrintMessage (5, "mesh updates complete");
    return;

    int cnttrials = 10;
    int wrongels = 0;
    for (int i = 1; i <= mesh.GetNE(); i++)
      if (mesh.VolumeElement(i).Volume(mesh.Points()) < 0)
	{
	  wrongels++;
	  mesh.VolumeElement(i).flags.badel = 1;
	}
      else
	mesh.VolumeElement(i).flags.badel = 0;

    if (wrongels)
      {
	cout << "WARNING: " << wrongels << " with wrong orientation found" << endl;

	int np = mesh.GetNP();
	NgArray<Point<3> > should(np);
	NgArray<Point<3> > can(np);
	for (int i = 1; i <= np; i++)
	  {
	    should.Elem(i) = can.Elem(i) = mesh.Point(i);
	  }
	for (int i = 1; i <= between.GetNBags(); i++)
	  for (int j = 1; j <= between.GetBagSize(i); j++)
	    {
	      INDEX_2 parent;
	      PointIndex child;
	      between.GetData (i, j, parent, child);
	      can.Elem(child) = Center (can.Elem(parent.I1()),
					can.Elem(parent.I2()));
	    }

	NgBitArray boundp(np);
	boundp.Clear();
	for (auto & sel : mesh.SurfaceElements())
          for (auto pi : sel.PNums())
            boundp.Set(pi);


	double lam = 0.5;

	while (lam < 0.9 && cnttrials > 0)
	  {
	    lam = 2;
	    do
	      {
		lam *= 0.5;
		cnttrials--;

		cout << "lam = " << lam << endl;

		for (int i = 1; i <= np; i++)
		  if (boundp.Test(i))
		    {
		      for (int j = 0; j < 3; j++)
			mesh.Point(i)(j) = 
			  lam * should.Get(i)(j) +
			  (1-lam) * can.Get(i)(j);
		    }
		  else
		    mesh.Point(i) = can.Get(i);
	      

		NgBitArray free (mesh.GetNP()), fhelp(mesh.GetNP());
		free.Clear();
		for (int i = 1; i <= mesh.GetNE(); i++)
		  {
		    const Element & el = mesh.VolumeElement(i);
		    if (el.Volume(mesh.Points()) < 0)
		      for (int j = 1; j <= el.GetNP(); j++)
			free.Set (el.PNum(j));
		  }
		for (int k = 1; k <= 3; k++)
		  {
		    fhelp.Clear();
		    for (int i = 1; i <= mesh.GetNE(); i++)
		      {
			const Element & el = mesh.VolumeElement(i);
			int freeel = 0;
			for (int j = 1; j <= el.GetNP(); j++)
			  if (free.Test(el.PNum(j)))
			    freeel = 1;
			if (freeel)
			  for (int j = 1; j <= el.GetNP(); j++)
			    fhelp.Set (el.PNum(j));
		      }
		    free.Or (fhelp);
		  }

		(*testout) << "smooth points: " << endl;
		for (int i = 1; i <= free.Size(); i++)
		  if (free.Test(i))
		    (*testout) << "p " << i << endl;

		(*testout) << "surf points: " << endl;
		for (auto & sel : mesh.SurfaceElements())
		  for (auto pi : sel.PNums())
		    (*testout) << pi << endl;

		mesh.CalcSurfacesOfNode();
		free.Invert();
		mesh.FixPoints (free);
		MeshingParameters dummymp;
		mesh.ImproveMesh (dummymp, OPT_REST);


		wrongels = 0;
		for (int i = 1; i <= mesh.GetNE(); i++)
		  {
		    if (mesh.VolumeElement(i).Volume(mesh.Points()) < 0)
		      {
			wrongels++;
			mesh.VolumeElement(i).flags.badel = 1;
			(*testout) << "wrong el: ";
			for (int j = 1; j <= 4; j++)
			  (*testout) << mesh.VolumeElement(i).PNum(j) << " ";
			(*testout) << endl;
		      }
		    else
		      mesh.VolumeElement(i).flags.badel = 0;
		  }
		cout << "wrongels = " << wrongels << endl;
	      }
	    while (wrongels && cnttrials > 0);
	  
	    for (int i = 1; i <= np; i++)
	      can.Elem(i) = mesh.Point(i);
	  }
      }

    if (cnttrials <= 0)
      {
	cerr << "ERROR: Sorry, reverted elements" << endl;
      }
 
    mesh.ComputeNVertices();
  }
}

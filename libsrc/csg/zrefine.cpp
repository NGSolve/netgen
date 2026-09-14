#include <mystdlib.h>
#include "meshing.hpp"

#include <csg.hpp>

namespace netgen
{

  // find singular edges
  void SelectSingularEdges (const Mesh & mesh, const CSGeometry & geom, 
			    ClosedHashTable<SortedPointIndices<2>, int> & singedges,
			    ZRefinementOptions & opt)
  {
    // edges selected in csg input file
    for (int i = 1; i <= geom.singedges.Size(); i++)
      {
	//if(geom.singedges.Get(i)->maxhinit > 0)
	//  continue; //!!!!

	const SingularEdge & se = *geom.singedges[i-1];
	for (int j = 1; j <= se.segms.Size(); j++)
	  {
	    PointIndices<2> i2 = se.segms[j-1];
	    singedges.Set (i2, 1);
	  }
      }

    // edges interactively selected
    for (int i = 1; i <= mesh.GetNSeg(); i++)
      {
	const Segment & seg = mesh.LineSegment(i);
	auto & ed = mesh.GetEdgeDescriptor(seg.GetIndex());
	if (ed.SingEdgeLeft() || ed.SingEdgeRight())
	  {
	    PointIndices<2> i2(seg[0], seg[1]);
	    singedges.Set (i2, 1);
	  }
      }
  }


  /**
     Convert elements (vol-tets, surf-trigs) into prisms/quads
  */
  void MakePrismsSingEdge (Mesh & mesh, ClosedHashTable<SortedPointIndices<2>, int> & singedges)
  {
    // volume elements
    // for (int i = 1; i <= mesh.GetNE(); i++)
    for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      {
	Element & el = mesh.VolumeElement(ei);
	if (el.GetType() != TET) continue;

	for (int j = 1; j <= 3; j++)
	  for (int k = j+1; k <= 4; k++)
	    {
	      SortedPointIndices<2> edge(el.PNum(j), el.PNum(k));
	      if (singedges.Used (edge))
		{
		  int pi3 = 1, pi4 = 1;
		  while (pi3 == j || pi3 == k) pi3++;
		  pi4 = 10 - j - k - pi3;
		
		  PointIndex p3 = el.PNum(pi3);
		  PointIndex p4 = el.PNum(pi4);

		  el.SetType(PRISM);
		  el.PNum(1) = edge.I1();
		  el.PNum(2) = p3;
		  el.PNum(3) = p4;
		  el.PNum(4) = edge.I2();
		  el.PNum(5) = p3;
		  el.PNum(6) = p4;
		}
	    }
      }

    // surface elements
    for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
      {
	Element2d & el = mesh.SurfaceElement(sei);
	if (el.GetType() != TRIG) continue;

	for (int j = 1; j <= 3; j++)
	  {
	    int k = (j % 3) + 1;
	    SortedPointIndices<2> edge(el.PNum(j), el.PNum(k));

	    if (singedges.Used (edge))
	      {
		int pi3 = 6-j-k;
		PointIndex p3 = el.PNum(pi3);
		PointIndex p1 = el.PNum(j);
		PointIndex p2 = el.PNum(k);

		el.SetType(QUAD);
		el.PNum(1) = p2;
		el.PNum(2) = p3;
		el.PNum(3) = p3;
		el.PNum(4) = p1;
	      }
	  }
      }
  }


  /*
    Convert tets and pyramids next to close (identified) points into prisms
  */
  void MakePrismsClosePoints (Mesh & mesh)
  {
    // int i, j, k;
    for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      {
	Element & el = mesh.VolumeElement(ei);
	if (el.GetType() == TET)
	  {
	    for (int j = 1; j <= 3; j++)
	      for (int k = j+1; k <= 4; k++)
		{
		  SortedPointIndices<2> edge(el.PNum(j), el.PNum(k));
		  if (mesh.GetIdentifications().UsedSymmetric (el.PNum(j), el.PNum(k)))
		    {
		      int pi3 = 1, pi4 = 1;
		      while (pi3 == j || pi3 == k) pi3++;
		      pi4 = 10 - j - k - pi3;
		    
		      PointIndex p3 = el.PNum(pi3);
		      PointIndex p4 = el.PNum(pi4);
		    
		      el.SetType(PRISM);
		      el.PNum(1) = edge.I1();
		      el.PNum(2) = p3;
		      el.PNum(3) = p4;
		      el.PNum(4) = edge.I2();
		      el.PNum(5) = p3;
		      el.PNum(6) = p4;
		    }
		}
	  }

	if (el.GetType() == PYRAMID)
	  {
	    // pyramid, base face = 1,2,3,4
	  
	    for (int j = 0; j <= 1; j++)
	      {
		PointIndex pi1 = el.PNum( (j+0) % 4 + 1);
		PointIndex pi2 = el.PNum( (j+1) % 4 + 1);
		PointIndex pi3 = el.PNum( (j+2) % 4 + 1);
		PointIndex pi4 = el.PNum( (j+3) % 4 + 1);
		PointIndex pi5 = el.PNum(5);

		if (mesh.GetIdentifications().UsedSymmetric (pi1, pi4) &&
		    mesh.GetIdentifications().UsedSymmetric (pi2, pi3))
		  {
		    //int p3 = el.PNum(pi3);
		    //int p4 = el.PNum(pi4);
		  
		    el.SetType(PRISM);
		    el.PNum(1) = pi1;
		    el.PNum(2) = pi2;
		    el.PNum(3) = pi5;
		    el.PNum(4) = pi4;
		    el.PNum(5) = pi3;
		    el.PNum(6) = pi5;
		  }
	      }
	  }
      }
  
    for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
      {
	Element2d & el = mesh.SurfaceElement(sei);
	if (el.GetType() != TRIG) continue;

	for (int j = 1; j <= 3; j++)
	  {
	    int k = (j % 3) + 1;
	    if (mesh.GetIdentifications().UsedSymmetric (el.PNum(j), el.PNum(k)))
	      {
		int pi3 = 6-j-k;
		PointIndex p3 = el.PNum(pi3);
		PointIndex p1 = el.PNum(j);
		PointIndex p2 = el.PNum(k);

		el.SetType(QUAD);
		el.PNum(1) = p2;
		el.PNum(2) = p3;
		el.PNum(3) = p3;
		el.PNum(4) = p1;
	      }
	  }
      }
  }



#ifdef OLD
  void MakeCornerNodes (Mesh & mesh,
			INDEX_HASHTABLE<int> & cornernodes)
  {
    int i, j;
    int nseg = mesh.GetNSeg();
    Array<int> edgesonpoint(mesh.GetNP());
    for (i = 1; i <= mesh.GetNP(); i++)
      edgesonpoint.Elem(i) = 0;

    for (i = 1; i <= nseg; i++)
      {
	for (j = 1; j <= 2; j++)
	  {
	    int pi = (j == 1) ? 
	      mesh.LineSegment(i)[0] :
	      mesh.LineSegment(i)[1];
	    edgesonpoint.Elem(pi)++;
	  }
      }

    /*
      cout << "cornernodes: ";
      for (i = 1; i <= edgesonpoint.Size(); i++)
      if (edgesonpoint.Get(i) >= 6)
      {
      cornernodes.Set (i, 1);
      cout << i << " ";
      }
      cout << endl;
    */
    //  cornernodes.Set (5, 1);
  }
#endif


  void RefinePrisms (Mesh & mesh, const CSGeometry * geom, 
		     ZRefinementOptions & opt)
  {
    // int i, j;
    bool found, change;
    int cnt = 0;


    // markers for z-refinement: the edge pts[0]-pts[1] is to be refined
    struct RefEdge { PointIndices<2> pts; int levels; };
    struct RefSliceEdge { PointIndices<2> pts; int idnr, slicenr; };
    Array<RefEdge> ref_uniform;
    Array<RefEdge> ref_singular;
    Array<RefSliceEdge> ref_slices;

    BitArray first_id(geom->identifications.Size());
    first_id.Set();

  
    // if (mesh.GetIdentifications().HasIdentifiedPoints())
      {
        auto & identpts =
          mesh.GetIdentifications().GetIdentifiedPoints ();

        /*
	for (int i = 1; i <= identpts.GetNBags(); i++)
	  for (int j = 1; j <= identpts.GetBagSize(i); j++)
	    {
	      INDEX_3 pair;
	      int dummy;
	      identpts.GetData(i, j, pair, dummy);
        */
        for (auto [hash, val] : identpts)\
          {
            auto [hash_pts, idnr] = hash;
            auto [pi1, pi2] = hash_pts;
            // auto idnr = pair[2];
            
	      const CloseSurfaceIdentification * csid = 
		dynamic_cast<const CloseSurfaceIdentification*> 
		(geom->identifications[idnr-1]);
	      if (csid)
		{
		  if (!csid->GetSlices().Size())
		    {
		      if (first_id.Test (idnr))
			{
			  first_id.Clear(idnr);
                          /*
			  ref_uniform.Append (INDEX_3 (pair.I1(), pair.I2(), csid->RefLevels()));
			  ref_singular.Append (INDEX_3 (pair.I1(), pair.I2(), csid->RefLevels1()));
			  ref_singular.Append (INDEX_3 (pair.I2(), pair.I1(), csid->RefLevels2()));
                          */
			  ref_uniform.Append ( { { pi1, pi2 }, csid->RefLevels() } );
			  ref_singular.Append ( { { pi1, pi2 }, csid->RefLevels1() } );
			  ref_singular.Append ( { { pi2, pi1 }, csid->RefLevels2() } );
                          
			}
		    }
		  else
		    {   
		      //const Array<double> & slices = csid->GetSlices();
		      ref_slices.Append ( { { pi1, pi2 }, idnr,
                                            int(csid->GetSlices().Size()) } );
		    }
		}
	    }
      }

  
  
    Array<EdgePointGeomInfo> epgi;

    while (1)
      {
	cnt++;
	PrintMessage (3, "Z-Refinement, level = ", cnt);
	ClosedHashTable<SortedPointIndices<2>, PointIndex> refedges(mesh.GetNSE()+1);


	found = 0;
	// mark prisms due to close surface flags:
	size_t oldsize = ref_uniform.Size();
	for (size_t i = 0; i < oldsize; i++)
	  {
	    PointIndex pi1 = ref_uniform[i].pts[0];
	    PointIndex pi2 = ref_uniform[i].pts[1];
	    int levels = ref_uniform[i].levels;

	    if (levels > 0)
	      {
		const Point3d & p1 = mesh[pi1];
		const Point3d & p2 = mesh[pi2];
		PointIndex npi = PointIndex::INVALID;
	      
		SortedPointIndices<2> edge(pi1, pi2);
		if (!refedges.Used(edge))
		  {
		    Point3d np = Center (p1, p2);
		    npi = mesh.AddPoint (np);
		    refedges.Set (edge, npi);
		    found = 1;
		  }

		ref_uniform[i] = { { pi1, npi }, levels-1 };
		ref_uniform.Append ( { { pi2, npi }, levels-1 } );
	      }
	  }
	for (size_t i = 0; i < ref_singular.Size(); i++)
	  {
	    PointIndex pi1 = ref_singular[i].pts[0];
	    PointIndex pi2 = ref_singular[i].pts[1];
	    int levels = ref_singular[i].levels;

	    if (levels > 0)
	      {
		const Point3d & p1 = mesh[pi1];
		const Point3d & p2 = mesh[pi2];
		PointIndex npi;
	      
		SortedPointIndices<2> edge(pi1, pi2);
		if (!refedges.Used(edge))
		  {
		    Point3d np = Center (p1, p2);
		    npi = mesh.AddPoint (np);
		    refedges.Set (edge, npi);
		    found = 1;
		  }
		else
		  npi = refedges.Get (edge);

		ref_singular[i] = { { pi1, npi }, levels-1 };
	      }
	  }

	for (size_t i = 0; i < ref_slices.Size(); i++)
	  {
	    PointIndex pi1 = ref_slices[i].pts[0];
	    PointIndex pi2 = ref_slices[i].pts[1];
	    int idnr = ref_slices[i].idnr;
	    int slicenr = ref_slices[i].slicenr;

	    if (slicenr > 0)
	      {
		const Point3d & p1 = mesh[pi1];
		const Point3d & p2 = mesh[pi2];
		PointIndex npi;

		const CloseSurfaceIdentification * csid = 
		  dynamic_cast<const CloseSurfaceIdentification*> 
		  (geom->identifications[idnr-1]);

	      
		SortedPointIndices<2> edge(pi1, pi2);
		if (!refedges.Used(edge))
		  {
		    const auto& slices = csid->GetSlices();
		    //(*testout) << "idnr " << idnr << " i " << i << endl;
		    //(*testout) << "slices " << slices << endl;
		    double slicefac = slices[slicenr-1];
		    double slicefaclast = 
		      (slicenr == slices.Size()) ? 1 : slices[slicenr];
		    
		    Point3d np = p1 + (slicefac / slicefaclast) * (p2-p1);
		    //(*testout) << "slicenr " << slicenr << " slicefac " << slicefac << " quot " << (slicefac / slicefaclast) << " np " << np << endl;
		    npi = mesh.AddPoint (np);
		    refedges.Set (edge, npi);
		    found = 1;
		  }
		else
		  npi = refedges.Get (edge);
		
		ref_slices[i].pts[1] = npi;
		ref_slices[i].slicenr--;
	      }
	  }




	for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	  {
	    Element & el = mesh.VolumeElement (ei);
	    if (el.GetType() != PRISM)
	      continue;

	    for (int j = 1; j <= 3; j++)
	      {
		PointIndex pi1 = el.PNum(j);
		PointIndex pi2 = el.PNum(j+3);
		const Point3d & p1 = mesh[pi1];
		const Point3d & p2 = mesh[pi2];

		bool ref = 0;

		/*
		  if (Dist (p1, p2) > mesh.GetH (Center (p1, p2)))
		  ref = 1;
		*/

		/*
		  if (cnt <= opt.minref)
		  ref = 1;
		*/

		/*
		  if ((pi1 == 460 || pi2 == 460 ||
		  pi1 == 461 || pi2 == 461) && cnt <= 8) ref = 1;
		*/
		if (ref == 1)
		  {
		    SortedPointIndices<2> edge(pi1, pi2);
		    if (!refedges.Used(edge))
		      {
			Point3d np = Center (p1, p2);
			PointIndex npi = mesh.AddPoint (np);
			refedges.Set (edge, npi);
			found = 1;
		      }
		  }
	      }
	  }
      
	if (!found) break;

	// build closure:
	PrintMessage (5, "start closure");
	do
	  {
	    PrintMessage (5, "start loop");
	    change = 0;
	    for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	      {
		Element & el = mesh.VolumeElement (ei);
		if (el.GetType() != PRISM)
		  continue;
	      
		bool hasref = 0, hasnonref = 0;
		for (int j = 1; j <= 3; j++)
		  {
		    PointIndex pi1 = el.PNum(j);
		    PointIndex pi2 = el.PNum(j+3);
		    if (pi1 != pi2)
		      {
			SortedPointIndices<2> edge(pi1, pi2);
			if (refedges.Used(edge))
			  hasref = 1;
			else 
			  hasnonref = 1;
		      }
		  }

		if (hasref && hasnonref)
		  {
		    //		  cout << "el " << i << " in closure" << endl;
		    change = 1;
		    for (int j = 1; j <= 3; j++)
		      {
			PointIndex pi1 = el.PNum(j);
			PointIndex pi2 = el.PNum(j+3);
			const Point3d & p1 = mesh[pi1];
			const Point3d & p2 = mesh[pi2];
		      
			SortedPointIndices<2> edge(pi1, pi2);
			if (!refedges.Used(edge))
			  {
			    Point3d np = Center (p1, p2);
			    PointIndex npi = mesh.AddPoint (np);
			    refedges.Set (edge, npi);
			  }
		      }
		  }
	      }
	  }
	while (change);

	PrintMessage (5, "Do segments");

	//      (*testout) << "closure formed, np = " << mesh.GetNP() << endl;

	int oldns = mesh.GetNSeg();

	for (int i = 1; i <= oldns; i++)
	  {
	    const Segment & el = mesh.LineSegment(i);

	    SortedPointIndices<2> i2(el[0], el[1]);
	  
	    PointIndex pnew;
	    EdgePointGeomInfo ngi;
      
	    if (refedges.Used(i2))
	      {
		pnew = refedges.Get(i2);
		//	      ngi = epgi.Get(pnew);
	      }
	    else
	      {
		continue;

		// 	      Point3d pb;

		// 	      /*
		// 	      geom->PointBetween (mesh.Point (el[0]),
		// 				  mesh.Point (el[1]),
		// 				  el.surfnr1, el.surfnr2,
		// 				  el.epgeominfo[0], el.epgeominfo[1],
		// 				  pb, ngi);
		// 	      */
		// 	      pb = Center (mesh.Point (el[0]), mesh.Point (el[1]));

		// 	      pnew = mesh.AddPoint (pb);
	      
		// 	      refedges.Set (i2, pnew);
	      
		// 	      if (pnew > epgi.Size())
		// 		epgi.SetSize (pnew);
		// 	      epgi.Elem(pnew) = ngi;
	      }
	  
	    Segment ns1 = el;
	    Segment ns2 = el;
	    ns1[1] = pnew;
	    ns1.EPGeomInfo(1) = ngi;
	    ns2[0] = pnew;
	    ns2.EPGeomInfo(0) = ngi;

	    mesh.LineSegment(i) = ns1;
	    mesh.AddSegment (ns2);
	  }
      
	PrintMessage (5, "Segments done, NSeg = ", mesh.GetNSeg());

	// do refinement
	int oldne = mesh.GetNE();
	for (ElementIndex ei = 0; ei < oldne; ei++)
	  {
	    Element & el = mesh.VolumeElement (ei);
	    if (el.GetNP() != 6)
	      continue;

	    PointIndex npi[3];
	    for (int j = 1; j <= 3; j++)
	      {
		PointIndex pi1 = el.PNum(j);
		PointIndex pi2 = el.PNum(j+3);

		if (pi1 == pi2)
		  npi[j-1] = pi1;
		else
		  {
		    SortedPointIndices<2> edge(pi1, pi2);
		    if (refedges.Used (edge))
		      npi[j-1] = refedges.Get(edge);
		    else
		      {
			/*
			  (*testout) << "ERROR: prism " << i << " has hanging node !!" 
			  << ", edge = " << edge << endl;
			  cerr << "ERROR: prism " << i << " has hanging node !!" << endl;
			*/
			npi[j-1] = PointIndex::INVALID;
		      }
		  }
	      }

	    if (npi[0].IsValid())
	      {
		Element nel1(6), nel2(6);
		for (int j = 1; j <= 3; j++)
		  {
		    nel1.PNum(j) = el.PNum(j);
		    nel1.PNum(j+3) = npi[j-1];
		    nel2.PNum(j) = npi[j-1];
		    nel2.PNum(j+3) = el.PNum(j+3);
		  }
		nel1.SetIndex (el.GetIndex());
		nel2.SetIndex (el.GetIndex());
		mesh.VolumeElement (ei) = nel1;
		mesh.AddVolumeElement (nel2);
	      }
	  }

      
	PrintMessage (5, "Elements done, NE = ", mesh.GetNE());


	// do surface elements
	int oldnse = mesh.GetNSE();
	//      cout << "oldnse = " << oldnse << endl;
	for (SurfaceElementIndex sei = 0; sei < oldnse; sei++)
	  {
	    Element2d & el = mesh.SurfaceElement (sei);
	    if (el.GetType() != QUAD)
	      continue;

	    int index = el.GetIndex();
	    PointIndex npi[2];
	    for (int j = 1; j <= 2; j++)
	      {
		PointIndex pi1, pi2;

		if (j == 1)
		  {
		    pi1 = el.PNum(1);
		    pi2 = el.PNum(4);
		  }
		else
		  {
		    pi1 = el.PNum(2);
		    pi2 = el.PNum(3);
		  }

		if (pi1 == pi2)
		  npi[j-1] = pi1;
		else
		  {
		    SortedPointIndices<2> edge(pi1, pi2);
		    if (refedges.Used (edge))
		      npi[j-1] = refedges.Get(edge);
		    else
		      {
			npi[j-1] = PointIndex::INVALID;
		      }
		  }
	      }

	    if (npi[0].IsValid())
	      {
		Element2d nel1(QUAD), nel2(QUAD);
		for (int j = 1; j <= 4; j++)
		  {
		    nel1.PNum(j) = el.PNum(j);
		    nel2.PNum(j) = el.PNum(j);
		  }
		nel1.PNum(3) = npi[1];
		nel1.PNum(4) = npi[0];
		nel2.PNum(1) = npi[0];
		nel2.PNum(2) = npi[1];
		/*
		  for (j = 1; j <= 2; j++)
		  {
		  nel1.PNum(j) = el.PNum(j);
		  nel1.PNum(j+2) = npi[j-1];
		  nel2.PNum(j) = npi[j-1];
		  nel2.PNum(j+2) = el.PNum(j+2);
		  }
		*/
		nel1.SetIndex (el.GetIndex());
		nel2.SetIndex (el.GetIndex());

		mesh.SurfaceElement (sei) = nel1;
		mesh.AddSurfaceElement (nel2);

		int si = mesh.GetFaceDescriptor (index).SurfNr();

		Point<3> hp = mesh[npi[0]];
		geom->GetSurface(si)->Project (hp);
		mesh[npi[0]].SetPoint (hp);

		hp = mesh[npi[1]];
		geom->GetSurface(si)->Project (hp);
		mesh[npi[1]].SetPoint (hp);

		//	      geom->GetSurface(si)->Project (mesh[PointIndex(npi[0])]);
		//	      geom->GetSurface(si)->Project (mesh[PointIndex(npi[1])]);
	      }
	  }

        mesh.RebuildSurfaceElementLists();
	PrintMessage (5, "Surface elements done, NSE = ", mesh.GetNSE());
      }
    
  }

  void CombineSingularPrisms(Mesh& mesh)
  {
    for(ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
      {
        Element& el = mesh.VolumeElement(ei);
        if(el.GetType() != PRISM)
          continue;
        if(el.PNum(3) == el.PNum(6))
          {
            if(el.PNum(2) == el.PNum(5))
              {
                el.SetType(TET);
              }
            else
              {
                el.SetType(PYRAMID);
                PointIndex pnr5 = el.PNum(3);
                el.PNum(3) = el.PNum(5);
                el.PNum(5) = pnr5;
              }
          }
      }
  }

  void ZRefinement (Mesh & mesh, const NetgenGeometry * hgeom,
		    ZRefinementOptions & opt)
  {
    const CSGeometry * geom = dynamic_cast<const CSGeometry*> (hgeom);
    if (!geom) return;

    ClosedHashTable<SortedPointIndices<2>, int> singedges(mesh.GetNSeg());

    SelectSingularEdges (mesh, *geom, singedges, opt);
    //MakePrismsSingEdge (mesh, singedges);
    MakePrismsClosePoints (mesh);

    RefinePrisms (mesh, geom, opt);

    CombineSingularPrisms(mesh);
  }



  ZRefinementOptions :: ZRefinementOptions()
  {
    minref = 0;
  }

}

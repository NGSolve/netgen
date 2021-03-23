#include <meshing.hpp>
#include <geometry2d.hpp>

namespace netgen
{
  // extern DLL_HEADER MeshingParameters mparam;

  extern void Optimize2d (Mesh & mesh, MeshingParameters & mp);




  void CalcPartition (const SplineSegExt & spline, 
		      // double l, 
		      MeshingParameters & mp, Mesh & mesh, 
		      // double h, double h1, double h2, double hcurve, 
		      double elto0, NgArray<double> & points)
  {
    double fperel, oldf, f;

    int n = 1;
    NgArray<Point<2> > xi;
    NgArray<double> hi;

    // do one extra step
    int not_fine_enough = 2;

    while(not_fine_enough && n < 10000)
    {
      not_fine_enough--;
      n*=4;

      xi.SetSize(n);
      hi.SetSize(n);

      for (int i = 0; i < n; i++)
      {
        xi[i] = spline.GetPoint ( (i+0.5) / n );
        hi[i] = mesh.GetH (Point<3> (xi[i](0), xi[i](1), 0));
      }

      // limit slope
      double gradh = min(1/elto0,mp.grading);
      for (int i = 0; i < n-1; i++)
      {
        double hnext = hi[i] + gradh * (xi[i+1]-xi[i]).Length();
        if(hnext > 2*hi[i])
          not_fine_enough = 2;
        hi[i+1] = min(hi[i+1], hnext);
      } 
      for (int i = n-1; i > 1; i--)
      {
        double hnext = hi[i] + gradh * (xi[i-1]-xi[i]).Length();
        if(hnext > 2*hi[i])
          not_fine_enough = 2;
        hi[i-1] = min(hi[i-1], hnext);
      }
    }

    points.SetSize (0);

    double len = spline.Length();
    double dt =  len / n;

    double sum = 0;
    for (int i = 1; i <= n; i++)
      {
	// double t = (i-0.5)*dt;
	double fun = hi[i-1];
	sum += dt / fun;
      }

    int nel = int (sum+0.5);
    if (nel == 0) nel = 1;
    fperel = sum / nel;

    points.Append (0);

    int i = 1;
    oldf = 0;

    for (int j = 1; j <= n && i < nel; j++)
      {
	double t = (j-0.5)*dt;
	double fun = hi[j-1];

	f = oldf + dt / fun;

	while (i * fperel < f && i < nel)
	  {
	    points.Append ( dt * (j-1) +  (i * fperel - oldf) * fun);
	    i++;
	  }
	oldf = f;
	t += dt;
      }
    points.Append (len);
  }









  // partitionizes spline curve
  void Partition (const SplineSegExt & spline,
		  MeshingParameters & mp, double hxxx, double elto0,
		  Mesh & mesh, Point3dTree & searchtree, int segnr) 
  {
    int n = 100;

    Point<2> mark, oldmark;
    NgArray<double> curvepoints;
    double edgelength, edgelengthold;

    CalcPartition (spline, mp, mesh, elto0, curvepoints);

    double dt = 1.0 / n;

    int j = 1;

    Point<2> pold = spline.GetPoint (0);
    double lold = 0;
    oldmark = pold;
    edgelengthold = 0;
    NgArray<int> locsearch;
    
    for (int i = 1; i <= n; i++)
      {
	Point<2> p = spline.GetPoint (i*dt);
	double l = lold + Dist (p, pold);
	while (j < curvepoints.Size() && (l >= curvepoints[j] || i == n))
	  {
	    double frac = (curvepoints[j]-lold) / (l-lold);
	    edgelength = i*dt + (frac-1)*dt;
	    mark = spline.GetPoint (edgelength);
	  
	    {
	      PointIndex pi1 = -1, pi2 = -1;
	  
	      Point3d mark3(mark(0), mark(1), 0);
	      Point3d oldmark3(oldmark(0), oldmark(1), 0);

	      double h = mesh.GetH (Point<3> (oldmark(0), oldmark(1), 0));
	      Vec<3> v (1e-4*h, 1e-4*h, 1e-4*h);
	      searchtree.GetIntersecting (oldmark3 - v, oldmark3 + v, locsearch);

	      for (int k = 0; k < locsearch.Size(); k++)
		if ( mesh[PointIndex(locsearch[k])].GetLayer() == spline.layer)
		  pi1 = locsearch[k];
	      
	      searchtree.GetIntersecting (mark3 - v, mark3 + v, locsearch);
	      for (int k = 0; k < locsearch.Size(); k++)
		if ( mesh[PointIndex(locsearch[k])].GetLayer() == spline.layer)
		  pi2 = locsearch[k];

	      if (pi1 == -1)
		{
		  pi1 = mesh.AddPoint(oldmark3, spline.layer);
		  searchtree.Insert (oldmark3, pi1);
		}
	      if (pi2 == -1)
		{
		  pi2 = mesh.AddPoint(mark3, spline.layer);
		  searchtree.Insert (mark3, pi2);
		}

	      Segment seg;
	      seg.edgenr = segnr;
	      seg.si = spline.bc; // segnr;
	      seg[0] = pi1;
	      seg[1] = pi2;
	      seg.domin = spline.leftdom;
	      seg.domout = spline.rightdom;
	      seg.epgeominfo[0].edgenr = segnr;
	      seg.epgeominfo[0].dist = edgelengthold;
	      seg.epgeominfo[1].edgenr = segnr;
	      seg.epgeominfo[1].dist = edgelength;
	      seg.singedge_left = spline.hpref_left;
	      seg.singedge_right = spline.hpref_right;
	      mesh.AddSegment (seg);
	    }
	
	    oldmark = mark;
	    edgelengthold = edgelength;
	    j++;
	  }
    
	pold = p;
	lold = l;
      }
  }



  void SplineGeometry2d :: PartitionBoundary (MeshingParameters & mp, double h, Mesh & mesh2d)
  {
    enum { D = 2 };
    Box<D> bbox;
    GetBoundingBox (bbox);
    double dist = Dist (bbox.PMin(), bbox.PMax());
    Point<3> pmin;
    Point<3> pmax;
  
    pmin(2) = -dist; pmax(2) = dist;
    for(int j=0;j<D;j++)
      {
	pmin(j) = bbox.PMin()(j);
	pmax(j) = bbox.PMax()(j);
      }

    Point3dTree searchtree (pmin, pmax);

    for (int i = 0; i < splines.Size(); i++)
      for (int side = 0; side <= 1; side++)
	{
	  int dom = (side == 0) ? GetSpline(i).leftdom : GetSpline(i).rightdom;
	  if (dom != 0) GetSpline(i).layer = GetDomainLayer (dom);
	}


    // mesh size restrictions ...

    for (auto & point : geompoints)
      mesh2d.RestrictLocalH (Point<3> (point(0), point(1), 0), point.hmax);
    
    for (int i = 0; i < splines.Size(); i++)
      {
	const SplineSegExt & spline = GetSpline(i);
	const GeomPoint<2> & p1 = spline.StartPI();
	const GeomPoint<2> & p2 = spline.EndPI();

	double h1 = min (p1.hmax, h/p1.refatpoint);
	mesh2d.RestrictLocalH (Point<3>(p1(0),p1(1),0), h1);
	double h2 = min (p2.hmax, h/p2.refatpoint);
	mesh2d.RestrictLocalH (Point<3>(p2(0),p2(1),0), h2);

	double len = spline.Length();
	mesh2d.RestrictLocalHLine (Point<3>(p1(0),p1(1),0), 
				   Point<3>(p2(0),p2(1),0), len/mp.segmentsperedge);

	double hcurve = min (spline.hmax, h/spline.reffak);
	double hl = GetDomainMaxh (spline.leftdom);
	if (hl > 0) hcurve = min2 (hcurve, hl);
	double hr = GetDomainMaxh (spline.rightdom);
	if (hr > 0) hcurve = min2 (hcurve, hr);

        // skip curvature restrictions for straight lines
        if(spline.MaxCurvature()==0)
        {
          mesh2d.RestrictLocalHLine (Point<3>(p1(0),p1(1),0), 
                                     Point<3>(p2(0),p2(1),0), hcurve);
        }
        else
        {
          int np = 1000;
          for (double t = 0.5/np; t < 1; t += 1.0/np)
            {
              Point<2> x = spline.GetPoint(t);
              double hc = 1.0/mp.curvaturesafety / (1e-99+spline.CalcCurvature (t));
              mesh2d.RestrictLocalH (Point<3> (x(0), x(1), 0), min2(hc, hcurve));
            }
        }
      }
    
    for (auto mspnt : mp.meshsize_points)
      mesh2d.RestrictLocalH (mspnt.pnt, mspnt.h);

    // add point elements
    for (auto & point : geompoints)
      if (point.name.length())
        {
          Point<3> newp(point(0), point(1), 0);
          PointIndex npi = mesh2d.AddPoint (newp, 1, FIXEDPOINT);
          mesh2d.AddLockedPoint(npi);
          Element0d el(npi, npi);
          el.name = point.name;
          mesh2d.SetCD2Name(npi, point.name);
          mesh2d.pointelements.Append (el);
          searchtree.Insert (newp, npi);          
        }

    // first add all vertices (for compatible orientation on periodic bnds)
    {
      double diam2 = Dist2(pmin, pmax);
      for (int i = 0; i < splines.Size(); i++)
        for (int j : { 1, 2 })
          {
            Point<2> hnewp = (j == 1) ? splines[i]->StartPI() : splines[i]->EndPI();
            Point<3> newp(hnewp(0), hnewp(1), 0);
            int layer = GetSpline(i).layer;
            int npi = -1;
            for (PointIndex pi = PointIndex::BASE; 
                 pi < mesh2d.GetNP()+PointIndex::BASE; pi++)
              if (Dist2 (mesh2d.Point(pi), newp) < 1e-12 * diam2 && mesh2d.Point(pi).GetLayer() == layer)
                npi = pi;
            
	    if (npi == -1)
	      {
		npi = mesh2d.AddPoint (newp, layer);
		searchtree.Insert (newp, npi);
                mesh2d.AddLockedPoint(npi);
                Element0d el(npi, npi);
                el.name = "";
                mesh2d.SetCD2Name(npi, "");
                mesh2d.pointelements.Append (el);
	      }
          }
    }
    
    for (int i = 0; i < splines.Size(); i++)
      if (GetSpline(i).copyfrom == -1)
	{
	  // astrid - set boundary meshsize to  domain meshsize h
	  // if no domain mesh size is given, the max h value from the bounding box is used
	  double hl = GetDomainMaxh ( GetSpline(i).leftdom );
	  double hr = GetDomainMaxh ( GetSpline(i).rightdom );

	  double useh = h;
	  if (hl > 0) useh = min2 (h, hl);
	  if (hr > 0) useh = min2 (h, hr);
	  Partition(GetSpline(i), mp, useh, elto0, mesh2d, searchtree, i+1);	    
	}
      else
	{
	  CopyEdgeMesh (GetSpline(i).copyfrom, i+1, mesh2d, searchtree);
	}
  }






  void SplineGeometry2d :: CopyEdgeMesh (int from, int to, Mesh & mesh, Point3dTree & searchtree)
  {
    Array<PointIndex, PointIndex> mappoints (mesh.GetNP());
    Array<double, PointIndex> param (mesh.GetNP());
    mappoints = PointIndex::INVALID;
    param = 0;

    Point3d pmin, pmax;
    mesh.GetBox (pmin, pmax);
    double diam2 = Dist2(pmin, pmax);

    PrintMessage(3, string("Copy edge, from ") + ToString(from) + " to " + ToString(to));
  
    for (const auto& seg : mesh.LineSegments())
      {
	if (seg.edgenr == from)
	  {
	    mappoints[seg[0]] = 1;
	    param[seg[0]] = seg.epgeominfo[0].dist;

	    mappoints[seg[1]] = 1;
	    param[seg[1]] = seg.epgeominfo[1].dist;
	  }
      }

    bool mapped = false;
    for (auto i : Range(mappoints))
      {
	if (mappoints[i].IsValid())
	  {
	    Point<2> newp = splines.Get(to)->GetPoint (param[i]);
	    Point<3> newp3 (newp(0), newp(1), 0);
	  
	    PointIndex npi = PointIndex::INVALID;
	  
	    for (auto pi : Range(mesh.Points()))
	      if (Dist2 (mesh.Point(pi), newp3) < 1e-12 * diam2)
		npi = pi;
	  
	    if (!npi.IsValid())
	      {
		npi = mesh.AddPoint (newp3);
		searchtree.Insert (newp3, npi);
	      }

	    mappoints[i] = npi;

	    mesh.GetIdentifications().Add (i, npi, to);
	    mapped = true;
	  }
      }
    if(mapped)
      mesh.GetIdentifications().SetType(to,Identifications::PERIODIC);

    // copy segments
    int oldnseg = mesh.GetNSeg();
    for (int i = 1; i <= oldnseg; i++)
      {
	const Segment & seg = mesh.LineSegment(i);
	if (seg.edgenr == from)
	  {
	    Segment nseg;
	    nseg.edgenr = to;
	    nseg.si = GetSpline(to-1).bc;      // splines.Get(to)->bc;
	    nseg[0] = mappoints[seg[0]];
	    nseg[1] = mappoints[seg[1]];
	    nseg.domin = GetSpline(to-1).leftdom;
	    nseg.domout = GetSpline(to-1).rightdom;
	  
	    nseg.epgeominfo[0].edgenr = to;
	    nseg.epgeominfo[0].dist = param[seg[0]];
	    nseg.epgeominfo[1].edgenr = to;
	    nseg.epgeominfo[1].dist = param[seg[1]];
	    mesh.AddSegment (nseg);
	  }
      }
  }





  void MeshFromSpline2D (SplineGeometry2d & geometry,
			 shared_ptr<Mesh> & mesh, 
			 MeshingParameters & mp)
  {
    static Timer tall("MeshFromSpline2D"); RegionTimer rtall(tall);
    static Timer t_h("SetH");
    static Timer t_tensor("tensor domain meshing");
    static Timer t_part_boundary("PartitionBoundary");
    static Timer t_hpref("mark hpref points");
    PrintMessage (1, "Generate Mesh from spline geometry");

    Box<2> bbox = geometry.GetBoundingBox ();

    t_h.Start();
    if (bbox.Diam() < mp.maxh) 
      mp.maxh = bbox.Diam();

    // double h = mp.maxh;

    // mesh = make_shared<Mesh>();
    mesh->SetDimension (2);

    Point3d pmin(bbox.PMin()(0), bbox.PMin()(1), -bbox.Diam());
    Point3d pmax(bbox.PMax()(0), bbox.PMax()(1), bbox.Diam());

    mesh->SetLocalH (pmin, pmax, mp.grading);
    mesh->SetGlobalH (mp.maxh);

    

    t_part_boundary.Start();
    geometry.PartitionBoundary (mp, mp.maxh, *mesh);
    t_part_boundary.Stop();
    
    PrintMessage (3, "Boundary mesh done, np = ", mesh->GetNP());


    t_hpref.Start();
    // marks mesh points for hp-refinement
    for (int i = 0; i < geometry.GetNP(); i++)
      if (geometry.GetPoint(i).hpref)
	{
	  double mindist = 1e99;
	  PointIndex mpi(0);
	  Point<2> gp = geometry.GetPoint(i);
	  Point<3> gp3(gp(0), gp(1), 0);
	  for (PointIndex pi = PointIndex::BASE; 
	       pi < mesh->GetNP()+PointIndex::BASE; pi++)
	    if (Dist2(gp3, (*mesh)[pi]) < mindist)
	      {
		mpi = pi;
		mindist = Dist2(gp3, (*mesh)[pi]);
	      }
	  (*mesh)[mpi].Singularity(geometry.GetPoint(i).hpref);
	}
    t_hpref.Stop();


    int maxdomnr = 0;
    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      {
	if ( (*mesh)[si].domin > maxdomnr) maxdomnr = (*mesh)[si].domin;
	if ( (*mesh)[si].domout > maxdomnr) maxdomnr = (*mesh)[si].domout;
      }

    TableCreator<const Segment*> dom2seg_creator(maxdomnr+1);
    for ( ; !dom2seg_creator.Done(); dom2seg_creator++)
      for (const Segment & seg : mesh->LineSegments())
        {
          dom2seg_creator.Add (seg.domin, &seg);
          if (seg.domin != seg.domout)
            dom2seg_creator.Add (seg.domout, &seg);
        }
    auto dom2seg = dom2seg_creator.MoveTable();
    
    mesh->ClearFaceDescriptors();
    for (int i = 1; i <= maxdomnr; i++)
      mesh->AddFaceDescriptor (FaceDescriptor (i, 0, 0, i));

    // set NgArray<string*> bcnames... 
    // number of bcnames
    int maxsegmentindex = 0;
    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      if ( (*mesh)[si].si > maxsegmentindex) maxsegmentindex = (*mesh)[si].si;

    mesh->SetNBCNames(maxsegmentindex);

    for ( int sindex = 0; sindex < maxsegmentindex; sindex++ )
      mesh->SetBCName ( sindex, geometry.GetBCName( sindex+1 ) );

    mesh->CalcLocalH(mp.grading);
    t_h.Stop();

    int bnp = mesh->GetNP(); // boundary points
    auto BndPntRange = mesh->Points().Range();

    int hquad = mp.quad;


    for (int domnr = 1; domnr <= maxdomnr; domnr++)
      if (geometry.GetDomainTensorMeshing (domnr))
        { // tensor product mesh
          RegionTimer rt(t_tensor);
          
          NgArray<PointIndex, PointIndex::BASE> nextpi(bnp);
          NgArray<int, PointIndex::BASE> si1(bnp), si2(bnp);
          // PointIndex firstpi;
          
          nextpi = -1;
          si1 = -1;
          si2 = -1;
          for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
            {
              int p1 = -1, p2 = -2;

              if ( (*mesh)[si].domin == domnr)
                { p1 = (*mesh)[si][0]; p2 = (*mesh)[si][1]; }
              if ( (*mesh)[si].domout == domnr)
                { p1 = (*mesh)[si][1]; p2 = (*mesh)[si][0]; }
              
              if (p1 == -1) continue;

              nextpi[p1] = p2;       // counter-clockwise
              
              int index = (*mesh)[si].si;
              if (si1[p1] != index && si2[p1] != index)
                { si2[p1] = si1[p1]; si1[p1] = index; }
              if (si1[p2] != index && si2[p2] != index)
                { si2[p2] = si1[p2]; si1[p2] = index; }
            }

          PointIndex c1(0), c2, c3, c4;  // 4 corner points
          int nex = 1, ney = 1;

          for (PointIndex pi = 1; pi <= si2.Size(); pi++)
            if (si2[pi] != -1)
              { c1 = pi; break; }      

          for (c2 = nextpi[c1]; si2[c2] == -1; c2 = nextpi[c2], nex++);
          for (c3 = nextpi[c2]; si2[c3] == -1; c3 = nextpi[c3], ney++);
          for (c4 = nextpi[c3]; si2[c4] == -1; c4 = nextpi[c4]);



          NgArray<PointIndex> pts ( (nex+1) * (ney+1) );   // x ... inner loop
          pts = -1;

          for (PointIndex pi = c1, i = 0; pi != c2; pi = nextpi[pi], i++)
            pts[i] = pi;
          for (PointIndex pi = c2, i = 0; pi != c3; pi = nextpi[pi], i++)
            pts[(nex+1)*i+nex] = pi;
          for (PointIndex pi = c3, i = 0; pi != c4; pi = nextpi[pi], i++)
            pts[(nex+1)*(ney+1)-i-1] = pi;
          for (PointIndex pi = c4, i = 0; pi != c1; pi = nextpi[pi], i++)
            pts[(nex+1)*(ney-i)] = pi;


          for (PointIndex pix = nextpi[c1], ix = 0; pix != c2; pix = nextpi[pix], ix++)
          {
            Point<3> px = (*mesh)[pix];
            for (PointIndex piy = nextpi[c2], iy = 0; piy != c3; piy = nextpi[piy], iy++)
              {
                double lam = Dist((*mesh)[piy],(*mesh)[c2]) / Dist((*mesh)[c3],(*mesh)[c2]);
                auto pix1 = pts[(nex+1)*ney+ix+1];
                auto pnew = px + lam*((*mesh)[pix1]-px);
                pts[(nex+1)*(iy+1) + ix+1] = mesh -> AddPoint (pnew, 1, FIXEDPOINT);
              }
          }

          for (int i = 0; i < ney; i++)
            for (int j = 0; j < nex; j++)
              {
                Element2d el(QUAD);
                el[0] = pts[i*(nex+1)+j];
                el[1] = pts[i*(nex+1)+j+1];
                el[2] = pts[(i+1)*(nex+1)+j+1];
                el[3] = pts[(i+1)*(nex+1)+j];
                el.SetIndex (domnr);

                mesh -> AddSurfaceElement (el);
              }
          char* material;
          geometry.GetMaterial(domnr, material);
          if(material)
            mesh->SetMaterial(domnr, material);
        }




    static Timer t_domain("Mesh domain");
    static Timer t_points("Mesh domain - find points");
    for (int domnr = 1; domnr <= maxdomnr; domnr++)
      {
        RegionTimer rt(t_domain);
        if (geometry.GetDomainTensorMeshing (domnr)) continue;
        
        double h = mp.maxh;
	if ( geometry.GetDomainMaxh ( domnr ) > 0 )
	  h = geometry.GetDomainMaxh(domnr);


	PrintMessage (3, "Meshing domain ", domnr, " / ", maxdomnr);

	int oldnf = mesh->GetNSE();

        mp.quad = hquad || geometry.GetDomainQuadMeshing (domnr);

	Meshing2 meshing (geometry, mp, Box<3> (pmin, pmax));

	NgArray<int, PointIndex::BASE> compress(mesh->GetNP());
	compress = -1;
	int cnt = 0;

        t_points.Start();
        /*
        for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
        {
          const auto & s = (*mesh)[si];
          if ( s.domin==domnr || s.domout==domnr )
          {
            for (auto pi : {s[0], s[1]})
            {
              if(compress[pi]==-1)
              {
                meshing.AddPoint((*mesh)[pi], pi);
                cnt++;
                compress[pi] = cnt;
              }
            }
          }
        }
        */
        for (const Segment * seg : dom2seg[domnr])
          if (seg->domin==domnr || seg->domout==domnr )
            for (auto pi : {(*seg)[0], (*seg)[1]})
              if (compress[pi]==-1)
                {
                  meshing.AddPoint((*mesh)[pi], pi);
                  cnt++;
                  compress[pi] = cnt;
                }


	PointGeomInfo gi;
	gi.trignum = 1;
        /*
	for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
	  {
	    if ( (*mesh)[si].domin == domnr)
	      {
		meshing.AddBoundaryElement (compress[(*mesh)[si][0]], 
                                            compress[(*mesh)[si][1]], gi, gi);
	      }
	    if ( (*mesh)[si].domout == domnr)
	      {
		meshing.AddBoundaryElement (compress[(*mesh)[si][1]],
                                            compress[(*mesh)[si][0]], gi, gi);
	      }
	  }
        */

	for (const Segment * seg : dom2seg[domnr])
	  {
	    if (seg->domin == domnr)
              meshing.AddBoundaryElement (compress[(*seg)[0]], 
                                          compress[(*seg)[1]], gi, gi);
            
	    if (seg->domout == domnr)
              meshing.AddBoundaryElement (compress[(*seg)[1]],
                                          compress[(*seg)[0]], gi, gi);
	  }

        
        t_points.Stop();

        if(mp.delaunay2d && cnt>1)
          meshing.Delaunay(*mesh, domnr, mp);
        else
        {
          // mp.checkoverlap = 0;
          auto res = meshing.GenerateMesh (*mesh, mp, h, domnr);
          if (res != 0)
            throw NgException("meshing failed");
        }
        
	for (SurfaceElementIndex sei = oldnf; sei < mesh->GetNSE(); sei++)
	  (*mesh)[sei].SetIndex (domnr);

	// astrid
	char * material;
	geometry.GetMaterial (domnr, material);
	if (material)
          mesh->SetMaterial (domnr, material);
      }

    mesh->Compress();

    mp.quad = hquad;



    int hsteps = mp.optsteps2d;

    mp.optimize2d = "smcm"; 
    mp.optsteps2d = hsteps/2;
    Optimize2d (*mesh, mp);

    mp.optimize2d = "Smcm"; 
    mp.optsteps2d = (hsteps+1)/2;
    Optimize2d (*mesh, mp);

    mp.optsteps2d = hsteps;

    mesh->Compress();
    mesh->OrderElements();
    mesh -> SetNextMajorTimeStamp();

    mp.Render();
  }








}

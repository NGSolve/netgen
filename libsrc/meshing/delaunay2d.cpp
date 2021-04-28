#include <mystdlib.h>
#include "meshing.hpp"

#include <geom2d/csg2d.hpp>


// not yet working ....

namespace netgen
{
  using ngcore::INT;
  extern void Optimize2d (Mesh & mesh, MeshingParameters & mp);

  static inline Point<2> P2( Point<3> p )
  {
    return {p[0], p[1]};
  }

  static inline Point<3> P3( Point<2> p )
  {
    return {p[0], p[1], 0};
  }


  class DelaunayTrig
  {
    PointIndex pnums[3];
    Point<2> c;
  public:    
    double r;
    double rad2;
    DelaunayTrig () = default;
    DelaunayTrig (int p1, int p2, int p3)
    { 
      pnums[0] = p1;
      pnums[1] = p2;
      pnums[2] = p3;
    }

    PointIndex & operator[] (int j) { return pnums[j]; }
    const PointIndex & operator[] (int j) const { return pnums[j]; }

    void CalcCenter (Mesh & mesh)
    {
      Point<2> p1 = P2(mesh[pnums[0]]);
      Point<2> p2 = P2(mesh[pnums[1]]);
      Point<2> p3 = P2(mesh[pnums[2]]);

      Vec<2> v1 = p2-p1;
      Vec<2> v2 = p3-p1;
      /*
      Mat<2,2> mat, inv;
      mat(0,0) = v1*v1;
      mat(0,1) = v1*v2;
      mat(1,0) = v2*v1;
      mat(1,1) = v2*v2;
      Vec<2> rhs, sol;
      rhs(0) = 0.5 * v1*v1;
      rhs(1) = 0.5 * v2*v2;
      CalcInverse (mat, inv);
      sol = inv * rhs;
      
      c = p1 + sol(0) * v1 + sol(1) * v2;
      rad2 = Dist2(c, p1);
      r = sqrt(rad2);
      */

      // without normal equation ...
      Mat<2,2> mat, inv;
      mat(0,0) = v1(0); mat(0,1) = v1(1);
      mat(1,0) = v2(0); mat(1,1) = v2(1);
      CalcInverse (mat, inv);
      Vec<2> rhs, sol;
      rhs(0) = 0.5 * v1*v1;
      rhs(1) = 0.5 * v2*v2;
      sol = inv * rhs;
      c = p1 + sol;

      rad2 = Dist2(c, p1);
      r = sqrt(rad2);
    }

    Point<2> Center() const { return c; }
    double Radius2() const { return rad2; }
    Box<2> BoundingBox() const { return Box<2> (c-Vec<2>(r,r), c+Vec<2>(r,r)); }

    mutable PointIndex visited_pi = -1;
  };

  class DelaunayMesh
  {
    ngcore::ClosedHashTable<INT<2>, INT<2>> edge_to_trig;
    Array<DelaunayTrig> trigs;
    unique_ptr<DelaunayTree<2>> tree;
    Mesh & mesh;

    Array<int> closeels;
    Array<int> intersecting;
    Array<INT<2>> edges;

    int GetNeighbour( int eli, int edge )
    {
      auto p0 = trigs[eli][(edge+1)%3];
      auto p1 = trigs[eli][(edge+2)%3];
      if(p1<p0)
        Swap(p0,p1);

      INT<2> hash = {p0,p1};
      /*
      if(!edge_to_trig.Used(hash))
        return -1;

      auto i2 = edge_to_trig.Get({p0,p1});

      return i2[0] == eli ? i2[1] : i2[0];
      */
      auto pos = edge_to_trig.Position(hash);
      if (pos == -1) return -1;
      auto i2 = edge_to_trig.GetData(pos);
      return i2[0] == eli ? i2[1] : i2[0];
    }

    void SetNeighbour( int eli, int edge )
    {
      auto p0 = trigs[eli][(edge+1)%3];
      auto p1 = trigs[eli][(edge+2)%3];
      if(p1<p0)
        Swap(p0,p1);

      INT<2> hash = {p0,p1};
      auto pos = edge_to_trig.Position(hash);
      if (pos == -1)
        edge_to_trig[hash] = {eli, -1};
      else
        {
          auto i2 = edge_to_trig.GetData(pos);
          if(i2[0]==-1)
            i2[0] = eli;
          else
            {
              if(i2[1]==-1)
                i2[1] = eli;
            }
          edge_to_trig.SetData (pos, i2);
        }
      /*
      if(!edge_to_trig.Used(hash))
        edge_to_trig[hash] = {eli, -1};
      else
      {

        auto i2 = edge_to_trig.Get({p0,p1});

        if(i2[0]==-1)
          i2[0] = eli;
        else
        {
          if(i2[1]==-1)
            i2[1] = eli;
        }

        edge_to_trig[hash] = i2;
      }
      */
    }

    void UnsetNeighbours( int eli )
    {
      for(int edge : Range(3))
      {
        auto p0 = trigs[eli][(edge+1)%3];
        auto p1 = trigs[eli][(edge+2)%3];
        if(p1<p0)
          Swap(p0,p1);

        INT<2> hash = {p0,p1};
        // auto i2 = edge_to_trig.Get({p0,p1});
        auto pos = edge_to_trig.Position(hash);
        auto i2 = edge_to_trig.GetData(pos);
          
        if(i2[0]==eli)
          i2[0] = i2[1];
        i2[1] = -1;

        // edge_to_trig[hash] = i2;
        edge_to_trig.SetData (pos, i2);
      }
    }


    void AppendTrig( int pi0, int pi1, int pi2 )
    {
      DelaunayTrig el;
      el[0] = pi0;
      el[1] = pi1;
      el[2] = pi2;

      el.CalcCenter(mesh);

      trigs.Append(el);
      int ti = trigs.Size()-1;
      tree->Insert(el.BoundingBox(), ti);

      for(int i : Range(3))
        SetNeighbour(ti, i);
    }

    public:
    DelaunayMesh( Mesh & mesh_, Box<2> box  ) 
      : mesh(mesh_)
    {
      Vec<2> vdiag = box.PMax()-box.PMin();

      double w = vdiag(0);
      double h = vdiag(1);

      Point<2> p0 = box.PMin() + Vec<2> ( -3*h, -h);
      Point<2> p1 = box.PMin() + Vec<2> (w+3*h, -h);
      Point<2> p2 = box.Center() + Vec<2> (0, 1.5*h+0.5*w);

      box.Add( p0 );
      box.Add( p1 );
      box.Add( p2 );

      tree = make_unique<DelaunayTree<2>>(box);

      auto pi0 = mesh.AddPoint (P3(p0));
      auto pi1 = mesh.AddPoint (P3(p1));
      auto pi2 = mesh.AddPoint (P3(p2));
      AppendTrig(pi0, pi1, pi2);
    }

    void AddPoint( PointIndex pi_new )
    {
      static Timer t("AddPoint"); RegionTimer reg(t);
      Point<2> newp = P2(mesh[pi_new]);
      intersecting.SetSize(0);
      edges.SetSize(0);

      int definitive_overlapping_trig = -1;

      double minquot{1e20};
      tree->GetFirstIntersecting (newp, newp, [&] (const auto i_trig)
          {
             const auto trig = trigs[i_trig];
             double rad2 = trig.Radius2();
             double d2 = Dist2 (trig.Center(), newp);
             if (d2 >= rad2) return false;
             
             if (d2 < 0.999 * rad2)
               {
                 definitive_overlapping_trig = i_trig;
                 return true;
               }
             
            if (definitive_overlapping_trig == -1 || d2 < 0.99*minquot*rad2) 
              {
                minquot = d2/rad2;
                definitive_overlapping_trig = i_trig;
              }
            return false;
          });

    if(definitive_overlapping_trig==-1)
    {
      static Timer t("slow check"); RegionTimer reg(t);
      PrintMessage (5, "Warning in delaunay tree - didn't find overlapping circle, check all trigs again");
      for(auto i_trig : trigs.Range())
      {
        const auto trig = trigs[i_trig];

        if(trig[0]==-1)
          continue;

        double rad2 = trig.Radius2();
        double d2 = Dist2 (trig.Center(), newp);

        // if (d2 < 0.999 * rad2)
        if (d2 < (1-1e-10)*rad2)
        {
          definitive_overlapping_trig = i_trig;
          break;
        }
      }
    }

    // static Timer tvis("trig visited");
    // tvis.Start();
    // BitArray trig_visited(trigs.Size());
    // trig_visited.Clear();
    if(definitive_overlapping_trig==-1)
      throw Exception("point not in any circle "+ ToString(pi_new));
    // tvis.Stop();
    // static Timer t2("addpoint - rest"); RegionTimer r2(t2);
    Array<int> trigs_to_visit;
    trigs_to_visit.Append(definitive_overlapping_trig);
    intersecting.Append(definitive_overlapping_trig);
    // trig_visited.SetBit(definitive_overlapping_trig);
    trigs[definitive_overlapping_trig].visited_pi = pi_new;
    
    while(trigs_to_visit.Size())
    {
      int ti = trigs_to_visit.Last();
      trigs_to_visit.DeleteLast();

      // trig_visited.SetBit(ti);

      auto & trig = trigs[ti];
      trig.visited_pi = pi_new;
      
      for(auto ei : Range(3))
      {
        auto nb = GetNeighbour(ti, ei);
        if(nb==-1)
          continue;
        // if(trig_visited.Test(nb))
        // continue;

        const auto & trig_nb = trigs[nb];
        if (trig_nb.visited_pi == pi_new)
          continue;

        // trig_visited.SetBit(nb);
        trig_nb.visited_pi = pi_new;
        
        bool is_intersecting = Dist2(newp, trig_nb.Center()) < trig_nb.Radius2()*(1+1e-12);

        if(!is_intersecting)
        {
          const Point<2> p0 = P2(mesh[PointIndex (trig[(ei+1)%3])]);
          const Point<2> p1 = P2(mesh[PointIndex (trig[(ei+2)%3])]);
          const Point<2> p2 = P2(mesh[PointIndex (trig[ei])]);
          auto v = p1-p0;

          Vec<2> n = {-v[1], v[0]};
          n /= n.Length();

          double dist = n * (newp-p1);
          double scal = n * (p2 - p1);
          if (scal > 0) dist *= -1;

          if (dist > -1e-10)
            is_intersecting = true;
        }

        if(is_intersecting)
        {
          trigs_to_visit.Append(nb);
          intersecting.Append(nb);
        }
      }
    }

    // find outer edges
    for (auto j : intersecting)
    {
      const DelaunayTrig & trig = trigs[j];
      for (int k = 0; k < 3; k++)
      {
        int p1 = trig[k];
        int p2 = trig[(k+1)%3];
        INT<2> edge{p1,p2};
        edge.Sort();
        bool found = false;
        for (int l = 0; l < edges.Size(); l++)
          if (edges[l] == edge)
          {
            edges.RemoveElement(l); 
            found = true;
            break;
          }
        if (!found)
          edges.Append (edge);
      }
    }

    for (int j : intersecting)
    {
      UnsetNeighbours(j);
      trigs[j][0] = -1;
      trigs[j][1] = -1;
      trigs[j][2] = -1;
    }

    for (auto edge : edges)
      AppendTrig( edge[0], edge[1], pi_new );

    for (int j : intersecting)
      tree->DeleteElement (j);

    static int counter=0;
    if(0)
    {
      Mesh m;
      m = mesh;
      m.ClearSurfaceElements();
      for (DelaunayTrig & trig : trigs)
      {
        if (trig[0] < 0) continue;

        Vec<3> n = Cross (mesh[trig[1]]-mesh[trig[0]], 
            mesh[trig[2]]-mesh[trig[0]]);
        if (n(2) < 0) Swap (trig[1], trig[2]);

        Element2d el(trig[0], trig[1], trig[2]);
        el.SetIndex (1);
        m.AddSurfaceElement (el);
      }
      m.Save("meshes/mesh_"+ToString(counter++)+".vol.gz");
    }

    }

    Array<DelaunayTrig> & GetElements() { return trigs; }

  };

  ostream & operator<< (ostream & ost, DelaunayTrig trig)
  {
    ost << trig[0] << "-" << trig[1] << "-" << trig[2] << endl;
    return ost;
  }

  
  void Meshing2 :: BlockFillLocalH (Mesh & mesh, const MeshingParameters & mp)
  {
    static Timer timer("Meshing2::BlockFill");
    static Timer timer1("Meshing2::BlockFill 1");
    static Timer timer2("Meshing2::BlockFill 2");
    static Timer timer3("Meshing2::BlockFill 3");
    static Timer timer4("Meshing2::BlockFill 4");
    RegionTimer reg (timer);

    timer1.Start();

    double filldist = mp.filldist;
    
    PrintMessage (6, "blockfill local h");

    NgArray<Point<3> > npoints;
    
    // adfront -> CreateTrees();

    Box<3> bbox ( Box<3>::EMPTY_BOX );
    double maxh = 0;
    
    for (int i = 0; i < adfront.GetNFL(); i++)
      {
	const FrontLine & line = adfront.GetLine (i);

	const Point<3> & p1 = adfront.GetPoint(line.L().I1());
	const Point<3> & p2 = adfront.GetPoint(line.L().I2());
	
        maxh = max (maxh, Dist (p1, p2));
	
	bbox.Add (p1);
	bbox.Add (p2);
      }

    
    // Point<3> mpc = bbox.Center();
    bbox.Increase (bbox.Diam()/2);
    Box<3> meshbox = bbox;

    timer1.Stop();
    timer2.Start();


    LocalH loch2 (bbox, 1, 2);
    
    if (mp.maxh < maxh) maxh = mp.maxh;
    
    bool changed;
    do 
      {
        static Timer tcf("clear flags");
        tcf.Start();
        // mesh.LocalHFunction().ClearFlags();
        mesh.LocalHFunction().ClearRootFlags();
	tcf.Stop();
        
        static Timer tcut("tcut");
        tcut.Start();
	for (int i = 0; i < adfront.GetNFL(); i++)
	  {
	    const FrontLine & line = adfront.GetLine(i);
	    
	    Box<3> bbox (adfront.GetPoint (line.L().I1()));
	    bbox.Add (adfront.GetPoint (line.L().I2()));

	    
	    double filld = filldist * bbox.Diam();
	    bbox.Increase (filld);
	    
	    mesh.LocalHFunction().CutBoundary (bbox); 
	  }
	tcut.Stop();

	mesh.LocalHFunction().FindInnerBoxes (&adfront, NULL);
	
	npoints.SetSize(0);
	mesh.LocalHFunction().GetInnerPoints (npoints);

	changed = false;
	for (int i = 0; i < npoints.Size(); i++)
	  {
	    if (mesh.LocalHFunction().GetH(npoints[i]) > 1.2 * maxh)
	      {
		mesh.LocalHFunction().SetH (npoints[i], maxh);
		changed = true;
	      }
	  }
      }
    while (changed);

    timer2.Stop();
    timer3.Start();


    if (debugparam.slowchecks)
      {
        (*testout) << "Blockfill with points: " << endl;
        *testout << "loch = " << mesh.LocalHFunction() << endl;

        *testout << "npoints = " << endl << npoints << endl;
      }


    int prims[] = { 211, 223, 227, 229, 233, 239, 241, 251, 257, 263 }; 
    int prim;
  
    {
      int i = 0;
      if (npoints.Size())
        while (npoints.Size() % prims[i] == 0) i++;
      prim = prims[i];
    }

    for (int i = 0; i < npoints.Size(); i++)
      {
        size_t hi = (size_t(prim) * size_t(i)) % npoints.Size();
        
	if (meshbox.IsIn (npoints[hi]))
	  {
	    PointIndex gpnum = mesh.AddPoint (npoints[hi]);
	    adfront.AddPoint (npoints[hi], gpnum);
	    
	    if (debugparam.slowchecks)
	      {
		(*testout) << npoints[hi] << endl;

		Point<2> p2d (npoints[hi](0), npoints[hi](1));
		if (!adfront.Inside(p2d))
		  {
		    cout << "add outside point" << endl;
		    (*testout) << "outside" << endl;
		  }
	      }
	    
	  }
      }
    
    timer3.Stop();
    timer4.Start();
  

  // find outer points
  
    loch2.ClearFlags();

    for (int i = 0; i < adfront.GetNFL(); i++)
      {
	const FrontLine & line = adfront.GetLine(i);
	
	Box<3> bbox (adfront.GetPoint (line.L().I1()));
	bbox.Add (adfront.GetPoint (line.L().I2()));
	
	loch2.SetH (bbox.Center(), bbox.Diam());
      }


    for (int i = 0; i < adfront.GetNFL(); i++)
      {
	const FrontLine & line = adfront.GetLine(i);
	
	Box<3> bbox (adfront.GetPoint (line.L().I1()));
	bbox.Add (adfront.GetPoint (line.L().I2()));

	bbox.Increase (filldist * bbox.Diam());
	loch2.CutBoundary (bbox);
      }
    
    loch2.FindInnerBoxes (&adfront, NULL);

      // outer points : smooth mesh-grading
    npoints.SetSize(0);
    loch2.GetOuterPoints (npoints);

    /*
    for (int i = 1; i <= npoints.Size(); i++)
      {
	if (meshbox.IsIn (npoints.Get(i)))
	  {
	    PointIndex gpnum = mesh.AddPoint (npoints.Get(i));
	    adfront.AddPoint (npoints.Get(i), gpnum);
	  }
      }  
    */

    for (const Point<3> p : npoints)
      if (meshbox.IsIn(p))
        {
          PointIndex gpnum = mesh.AddPoint (p);
          adfront.AddPoint (p, gpnum);
        }
    timer4.Stop();
  }



  
  void Meshing2 :: Delaunay (Mesh & mesh, int domainnr, const MeshingParameters & mp)
  {
    static Timer timer("Meshing2::Delaunay");
    static Timer t1("Meshing2::Delaunay1");
    static Timer t2("Meshing2::Delaunay2");
    static Timer t3("Meshing2::Delaunay3");
    static Timer timer_addpoints("add points");
    RegionTimer reg (timer);

    PrintMessage (4, "2D Delaunay meshing");

    auto first_point_blockfill = mesh.Points().Range().Next();

    BlockFillLocalH (mesh, mp);

    auto last_point_blockfill = mesh.Points().Range().Next();

    t1.Start();
    // Bounding box for starting trig in delaunay
    Box<2> bbox (Box<2>::EMPTY_BOX);

    for (int i = 0; i < adfront.GetNFL(); i++)
      {
	const FrontLine & line = adfront.GetLine(i);
        bbox.Add (P2(Point<3> (adfront.GetPoint (line.L()[0]))));
        bbox.Add (P2(Point<3> (adfront.GetPoint (line.L()[1]))));
      }

    for (PointIndex pi : Range(first_point_blockfill, last_point_blockfill))
      bbox.Add(P2(mesh[pi]));

    for (int i = 0; i < mesh.LockedPoints().Size(); i++)
      bbox.Add (P2(mesh.Point (mesh.LockedPoints()[i])));
    t1.Stop();

    t2.Start();
    Array<PointIndex> old_points;
    BitArray add_point(mesh.Points().Size()+1);
    Array<PointIndex> addpoints;
    add_point.Clear();
    /*
    for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
    {
      const auto & s = mesh[si];
      if ( s.domin==domainnr || s.domout==domainnr )
      {
        add_point.SetBit(s[0]);
        add_point.SetBit(s[1]);
      }
    }
    */
    /*
    for (int i = 0; i < adfront.GetNFL(); i++)
      {
	const FrontLine & line = adfront.GetLine(i);
        for (int j = 0; j < 2; j++)
          add_point.SetBit (adfront.GetGlobalIndex (line.L()[j]))adfront.GetGlobalIndex (line.L()[j]));
      }
    */
    for (const auto & line : adfront.GetLines())
      for (int j = 0; j < 2; j++)
        {
          PointIndex pnum = adfront.GetGlobalIndex (line.L()[j]);
          if (!add_point.Test(pnum))
            addpoints.Append(pnum);
          add_point.SetBit (pnum);
        }
      
    
    t2.Stop();

    t3.Start();
    Mesh tempmesh;
    tempmesh.AddFaceDescriptor (FaceDescriptor (1, 1, 0, 0));
    tempmesh.AddFaceDescriptor (FaceDescriptor (2, 1, 0, 0));
    tempmesh.AddFaceDescriptor (FaceDescriptor (3, 1, 0, 0));

    Array<PointIndex, PointIndex> compress;
    Array<PointIndex, PointIndex> icompress(mesh.Points().Size());

    /*
    for(auto pi : mesh.Points().Range())
      if(add_point.Test(pi))
    */
    for (PointIndex pi : addpoints)
      {
        icompress[pi] = tempmesh.AddPoint(mesh[pi]);
        compress.Append(pi);
      }

    for (PointIndex pi : Range(first_point_blockfill, last_point_blockfill))
      {
        icompress[pi] = tempmesh.AddPoint(mesh[pi]);
        compress.Append(pi);
      }
    t3.Stop();
    // DelaunayMesh adds surrounding trig (don't add the last 3 points to delaunay AGAIN!
    auto tempmesh_points = tempmesh.Points().Range();

    DelaunayMesh dmesh(tempmesh, bbox);

    timer_addpoints.Start();

//     // reorder points
//     NgArray<PointIndex, PointIndex::BASE, PointIndex> mixed(old_points.Size());
//     int prims[] = { 11, 13, 17, 19, 23, 29, 31, 37 };
//     int prim;
//   
//     {
//       int i = 0;
//       while (old_points.Size() % prims[i] == 0) i++;
//       prim = prims[i];
//     }
// 
//     for (PointIndex pi : old_points)
//       mixed[pi] = PointIndex ( (prim * pi) % old_points.Size() + PointIndex::BASE );
    
    for (auto pi : tempmesh_points)
      dmesh.AddPoint(pi);

    timer_addpoints.Stop();

    static Timer taddseg("addseg");
    taddseg.Start();

    /*
    for (auto seg : mesh.LineSegments())
    {
      if ( seg.domin == domainnr || seg.domout == domainnr )
      {
        if(seg.domin==domainnr)
          seg.domout = 0;
        if(seg.domout==domainnr)
          seg.domin = 0;
        seg[0] = icompress[seg[0]];
        seg[1] = icompress[seg[1]];
        tempmesh.AddSegment(seg);
      }
    }
    */
    for (const auto & line : adfront.GetLines())
      {
        Segment seg;
        for (int j = 0; j < 2; j++)
          seg[j] = icompress [adfront.GetGlobalIndex (line.L()[j])];
        seg.domin = domainnr;
        seg.domout = 0;
        tempmesh.AddSegment(seg);
      }
           
    taddseg.Stop();
    
    for (auto & trig : dmesh.GetElements())
    {
      if (trig[0] < 0) continue;

      Element2d el(trig[0], trig[1], trig[2]);
      el.SetIndex (1);
      tempmesh.AddSurfaceElement (el);
    }

    bool conforming = false;
    while(!conforming)
    {
      conforming = true;
      BitArray marked_points(tempmesh.Points().Size()+1);
      marked_points = false;
      // Check for trigs cutting a boundary edge (non-conforming mesh)
      auto point_to_trigs = tempmesh.CreatePoint2SurfaceElementTable( 0 );
      for (auto & seg : tempmesh.LineSegments())
      {
        int count_adjacent = 0;;
        PointIndex pi0 = seg[0];
        PointIndex pi1 = seg[1];
        if(marked_points.Test(pi0)) continue;
        if(marked_points.Test(pi1)) continue;

        for(auto sei : point_to_trigs[pi0])
          for( auto i : Range(3))
            if(tempmesh[sei][i] == pi1)
              count_adjacent++;

        if(count_adjacent==2)
          continue;

        PointIndex pi2;
        PointIndex pi3;
        ArrayMem<SurfaceElementIndex, 2> cutting_trigs;
        for(auto sei : point_to_trigs[pi0])
        {
          auto & el = tempmesh[sei];
          pi2 = el[0] == pi0 ? el[1] : el[0];
          pi3 = el[2] == pi0 ? el[1] : el[2];
          double alpha, beta;
          auto itype = intersect( P2(tempmesh[pi0]), P2(tempmesh[pi1]), P2(tempmesh[pi2]), P2(tempmesh[pi3]), alpha, beta );
          if(itype == X_INTERSECTION)
          {
            cutting_trigs.Append(sei);
            break;
          }
        }
        if(cutting_trigs.Size()==0)
          continue;
        for(auto sei : point_to_trigs[pi2])
        {
          if(sei==cutting_trigs[0])
            continue;
          for(auto i : IntRange(3))
            if(tempmesh[sei][i]==pi3)
              cutting_trigs.Append(sei);
        }

        // Found two trigs cutting a boundary edge -> perform swap
        if(cutting_trigs.Size()==2)
        {
          conforming = false;
          if(marked_points.Test(pi2)) continue;
          if(marked_points.Test(pi3)) continue;

          auto & el0 = tempmesh[cutting_trigs[0]];
          auto & el1 = tempmesh[cutting_trigs[1]];

          pi1 = el1[0]+el1[1]+el1[2] - pi2-pi3;

          if(marked_points.Test(pi1)) continue;

          marked_points.SetBit(pi0);
          marked_points.SetBit(pi1);
          marked_points.SetBit(pi2);
          marked_points.SetBit(pi3);

          el0[0] = pi2;
          el0[1] = pi1;
          el0[2] = pi0;

          el1[0] = pi3;
          el1[1] = pi0;
          el1[2] = pi1;
        }
      }
    }

    auto point_to_trigs = tempmesh.CreatePoint2SurfaceElementTable( 0 );

    // Mark edges and trigs as inside or outside, starting with boundary edges
    enum POSITION { UNKNOWN, BOUNDARY, INSIDE, OUTSIDE };
    Array<POSITION, SurfaceElementIndex> trig_pos(tempmesh.SurfaceElements().Size());
    ngcore::ClosedHashTable<INT<2>, POSITION> edge_pos(3*tempmesh.SurfaceElements().Size());
    trig_pos = UNKNOWN;

    for (auto & seg : tempmesh.LineSegments())
    {
      ArrayMem<SurfaceElementIndex, 2> els;
      INT<2> edge{seg[0], seg[1]};
      edge.Sort();
      edge_pos[edge] = BOUNDARY;

      for(auto sei : point_to_trigs[seg[0]])
        for( auto i : Range(3))
          if(tempmesh[sei][i] == seg[1])
            els.Append(sei);

      for(auto sei : els)
      {
        auto & el = tempmesh[sei];
        PointIndex pi2 = el[0]+el[1]+el[2] - seg[0] - seg[1];
        bool is_left = ::netgen::Area(P2(tempmesh[seg[0]]), P2(tempmesh[seg[1]]), P2(tempmesh[pi2]))>0.0;
        POSITION pos;

        if(is_left == (seg.domin==domainnr))
          pos = INSIDE;
        else
          pos = OUTSIDE;

        INT<2> e1{seg[0], pi2};
        INT<2> e2{seg[1], pi2};
        e1.Sort();
        e2.Sort();
        if(!edge_pos.Used(e1))
          edge_pos[e1] = pos;
        if(!edge_pos.Used(e2))
          edge_pos[e2] = pos;
        trig_pos[sei] = pos;
      }
    }

    // Advance from boundary edges/trigs to all others
    bool have_unknown_trigs = true;
    while(have_unknown_trigs)
    {
      have_unknown_trigs = false;

      for (auto sei : Range(tempmesh.SurfaceElements()))
      {
        auto & el = tempmesh[sei];

        if(trig_pos[sei] == UNKNOWN)
        {
          have_unknown_trigs = true;

          // any edge of unkown trig already marked?
          for(auto i : IntRange(3))
          {
            INT<2> edge{el[(i+1)%3], el[(i+2)%3]};
            edge.Sort();
            if(edge_pos.Used(edge) && edge_pos[edge]!=BOUNDARY)
            {
              trig_pos[sei] = edge_pos[edge];
              break;
            }
          }
        }

        // if we could mark the trig -> also mark all edges
        if(trig_pos[sei] != UNKNOWN)
          for(auto i : IntRange(3))
          {
            INT<2> edge{el[(i+1)%3], el[(i+2)%3]};
            edge.Sort();
            if(!edge_pos.Used(edge) || edge_pos[edge]==BOUNDARY)
              edge_pos[edge] = trig_pos[sei];
          }
      }
    }

    // add inside trigs to actual mesh
    for (auto sei : Range(tempmesh.SurfaceElements()))
    {
      if(trig_pos[sei] == INSIDE)
      {
        auto el = tempmesh[sei];

        Vec<3> n = Cross (tempmesh[el[1]]-tempmesh[el[0]],
            tempmesh[el[2]]-tempmesh[el[0]]);
        if (n(2) < 0) Swap (el[1], el[2]);

        el[0] = compress[el[0]];
        el[1] = compress[el[1]];
        el[2] = compress[el[2]];
        el.SetIndex(domainnr);
        mesh.AddSurfaceElement(el);
      }
    }

    // mesh.Compress();  // don't compress whole mesh after every sub-domain
  }

}

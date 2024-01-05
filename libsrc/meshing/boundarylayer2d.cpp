#include <mystdlib.h>
#include "meshing.hpp"
#include "meshing2.hpp"
#include "../geom2d/csg2d.hpp"

namespace netgen
{
   void InsertVirtualBoundaryLayer (Mesh & mesh)
   {
      cout << "Insert virt. b.l." << endl;

      int surfid;

      cout << "Boundary Nr:";
      cin >> surfid;

      int i;
      int np = mesh.GetNP();

      cout << "Old NP: " << mesh.GetNP() << endl;
      cout << "Trigs: " << mesh.GetNSE() << endl;

      NgBitArray bndnodes(np);
      NgArray<int> mapto(np);

      bndnodes.Clear();
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         cout << "snr = " << snr << endl;
         if (snr == surfid)
         {
            bndnodes.Set (mesh.LineSegment(i)[0]);
            bndnodes.Set (mesh.LineSegment(i)[1]);
         }
      }
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         if (snr != surfid)
         {
            bndnodes.Clear (mesh.LineSegment(i)[0]);
            bndnodes.Clear (mesh.LineSegment(i)[1]);
         }
      }

      for (i = 1; i <= np; i++)
        {
          if (bndnodes.Test(i))
            mapto.Elem(i) = mesh.AddPoint (mesh.Point (i));
          else
            mapto.Elem(i) = 0;
        }

      for (i = 1; i <= mesh.GetNSE(); i++)
      {
         Element2d & el = mesh.SurfaceElement(i);
         for (int j = 1; j <= el.GetNP(); j++)
            if (mapto.Get(el.PNum(j)))
               el.PNum(j) = mapto.Get(el.PNum(j));
      }


      int nq = 0;
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         if (snr == surfid)
         {
            int p1 = mesh.LineSegment(i)[0];
            int p2 = mesh.LineSegment(i)[1];
            int p3 = mapto.Get (p1);
            if (!p3) p3 = p1;
            int p4 = mapto.Get (p2);
            if (!p4) p4 = p2;

            Element2d el(QUAD);
            el.PNum(1) = p1;
            el.PNum(2) = p2;
            el.PNum(3) = p3;
            el.PNum(4) = p4;
            el.SetIndex (2);
            mesh.AddSurfaceElement (el);
            nq++;
         }
      }

      cout << "New NP: " << mesh.GetNP() << endl;
      cout << "Quads: " << nq << endl;
   }


  void AddDirection( Vec<3> & a, Vec<3> b )
  {
     if(a.Length2()==0.)
     {
        a = b;
        return;
     }

     if(b.Length2()==0.)
        return;

     auto ab = a * b;
     if(fabs(ab)>1-1e-8)
        return;

     Mat<2> m;
     m(0,0) = a[0];
     m(0,1) = a[1];
     m(1,0) = b[0];
     m(1,1) = b[1];
     Vec<2> lam;
     Vec<2> rhs;
     rhs[0] = a[0]-b[0];
     rhs[1] = a[1]-b[1];

     const auto Dot = [](Vec<3> a, Vec<3> b)
     { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; };

     rhs[0] = Dot(a,a);
     rhs[1] = Dot(b,b);

     m.Solve(rhs, lam);
     a[0] = lam[0];
     a[1] = lam[1];
     a[2] = 0.0;
     return;
  }

  static void Generate2dMesh( Mesh & mesh, int domain )
  {
     Box<3> box{Box<3>::EMPTY_BOX};
     for(const auto & seg : mesh.LineSegments())
        if (seg.domin == domain || seg.domout == domain)
           for (auto pi : {seg[0], seg[1]})
              box.Add(mesh[pi]);

     MeshingParameters mp;
     Meshing2 meshing (*mesh.GetGeometry(), mp, box);

     Array<PointIndex, PointIndex> compress(mesh.GetNP());
     compress = PointIndex::INVALID;

     PointIndex cnt = PointIndex::BASE;

     auto p2sel = mesh.CreatePoint2SurfaceElementTable();
     PointGeomInfo gi;
     gi.u = 0.0;
     gi.v = 0.0;
     gi.trignum = domain;
     for(auto seg : mesh.LineSegments())
     {
         if(seg.domin == domain || seg.domout == domain)
         for (auto pi : {seg[0], seg[1]})
            if (compress[pi]==PointIndex{PointIndex::INVALID})
            {
               meshing.AddPoint(mesh[pi], pi);
               compress[pi] = cnt++;
            }
         if(seg.domin == domain)
             meshing.AddBoundaryElement (compress[seg[0]], compress[seg[1]], gi, gi);
         if(seg.domout == domain)
             meshing.AddBoundaryElement (compress[seg[1]], compress[seg[0]], gi, gi);
     }

     auto oldnf = mesh.GetNSE();
     // auto res =
     meshing.GenerateMesh (mesh, mp, mp.maxh, domain);
     for (SurfaceElementIndex sei : Range(oldnf, mesh.GetNSE()))
        mesh[sei].SetIndex (domain);
     
     // int hsteps = mp.optsteps2d;

     const char * optstr = mp.optimize2d.c_str();
     MeshOptimize2d meshopt(mesh);
     meshopt.SetFaceIndex(domain);
     meshopt.SetMetricWeight (mp.elsizeweight);
     for (size_t j = 1; j <= strlen(optstr); j++)
     {
        switch (optstr[j-1])
        {
           case 's':
              {  // topological swap
                 meshopt.EdgeSwapping (0);
                 break;
              }
           case 'S':
              {  // metric swap
                 meshopt.EdgeSwapping (1);
                 break;
              }
           case 'm':
              {
                 meshopt.ImproveMesh(mp);
                 break;
              }
           case 'c':
              {
                 meshopt.CombineImprove();
                 break;
              }
           default:
              cerr << "Optimization code " << optstr[j-1] << " not defined" << endl;
        }
     }

     mesh.Compress();
     mesh.CalcSurfacesOfNode();
     mesh.OrderElements();
     mesh.SetNextMajorTimeStamp();
  }

  int GenerateBoundaryLayer2 (Mesh & mesh, int domain, const Array<double> & thicknesses, bool should_make_new_domain, const Array<int> & boundaries)
  {
     mesh.GetTopology().SetBuildVertex2Element(true);
     mesh.UpdateTopology();
     const auto & line_segments = mesh.LineSegments();
     SegmentIndex first_new_seg = mesh.LineSegments().Range().Next();

     int np = mesh.GetNP();
     int nseg = line_segments.Size();
     // int ne = mesh.GetNSE();
     mesh.UpdateTopology();

     double total_thickness = 0.0;
     for(auto thickness : thicknesses)
        total_thickness += thickness;

     Array<Array<PointIndex>, PointIndex> mapto(np);

     // Bit array to keep track of segments already processed
     BitArray segs_done(nseg);
     segs_done.Clear();

     // moved segments
     Array<SegmentIndex> moved_segs;

     Array<Vec<3>, PointIndex> growthvectors(np);
     growthvectors = 0.;

     auto & meshtopo = mesh.GetTopology();

     Array<SegmentIndex> segments;

    // surface index map
    Array<int> si_map(mesh.GetNFD()+2);
    si_map = -1;

    // int fd_old = mesh.GetNFD();

    int max_edge_nr = -1;
    int max_domain = -1;

    for(const auto& seg : line_segments)
    {
      if(seg.epgeominfo[0].edgenr > max_edge_nr)
        max_edge_nr = seg.epgeominfo[0].edgenr;
      if(seg.domin > max_domain)
         max_domain = seg.domin;
      if(seg.domout > max_domain)
         max_domain = seg.domout;
    }

    int new_domain = max_domain+1;

    BitArray active_boundaries(max_edge_nr+1);
    BitArray active_segments(nseg);
    active_boundaries.Clear();
    active_segments.Clear();

    if(boundaries.Size() == 0)
       active_boundaries.Set();
    else
       for(auto edgenr : boundaries)
          active_boundaries.SetBit(edgenr);

    for(auto segi : Range(line_segments))
    {
       const auto seg = line_segments[segi];
       if(active_boundaries.Test(seg.epgeominfo[0].edgenr) && (seg.domin==domain || seg.domout==domain))
          active_segments.SetBit(segi);
    }

    {
        FaceDescriptor new_fd(0, 0, 0, -1);
        new_fd.SetBCProperty(new_domain);
        // int new_fd_index =
        mesh.AddFaceDescriptor(new_fd);
        if(should_make_new_domain)
           mesh.SetBCName(new_domain-1, "mapped_" + mesh.GetBCName(domain-1));
    }

    for(auto segi : Range(line_segments))
      {
        if(segs_done[segi]) continue;
        segs_done.SetBit(segi);
        const auto& seg = line_segments[segi];
        if(seg.domin != domain && seg.domout != domain) continue;
        if(!active_boundaries.Test(seg.epgeominfo[0].edgenr))
           continue;
        moved_segs.Append(segi);
      }

     // calculate growth vectors (average normal vectors of adjacent segments at each point)
     for (auto si : moved_segs)
     {
       auto & seg = line_segments[si];

       auto n = mesh[seg[1]] - mesh[seg[0]];
       n = {-n[1], n[0], 0};
       n.Normalize();

       if(seg.domout == domain)
           n = -n;

       AddDirection(growthvectors[seg[0]], n);
       AddDirection(growthvectors[seg[1]], n);
     }

     //////////////////////////////////////////////////////////////////////////
     // average growthvectors along straight lines to avoid overlaps in corners
     BitArray points_done(np+1);
     points_done.Clear();

     for(auto si : moved_segs)
     {
        auto current_seg = line_segments[si];
        auto current_si = si;

        auto first = current_seg[0];
        auto current = -1;
        auto next =  current_seg[1];

        if(points_done.Test(first))
           continue;

        Array<PointIndex> chain;
        chain.Append(first);

        // first find closed loops of segments
        while(next != current && next != first)
        {
           current = next;
           points_done.SetBit(current);
           chain.Append(current);
           for(auto sj : meshtopo.GetVertexSegments( current ))
           {
              if(!active_segments.Test(sj))
                 continue;

              if(sj!=current_si)
              {
                 current_si = sj;
                 current_seg = mesh[sj];

                 next = current_seg[0] + current_seg[1] - current;
                 break;
              }
           }
        }

        auto ifirst = 0;
        auto n = chain.Size();

        // angle of adjacent segments at points a[i-1], a[i], a[i+1]
        auto getAngle = [&mesh, &growthvectors] (FlatArray<PointIndex> a, size_t i)
        {
           auto n = a.Size();
           auto v0 = growthvectors[a[(i+n-1)%n]];
           auto v1 = growthvectors[a[i]];
           auto v2 = growthvectors[a[(i+1)%n]];

           auto p0 = mesh[a[(i+n-1)%n]];
           auto p1 = mesh[a[i]];
           auto p2 = mesh[a[(i+1)%n]];

           v0 = p1-p0;
           v1 = p2-p1;

           auto angle = abs(atan2(v1[0], v1[1]) - atan2(v0[0], v0[1]));
           if(angle>M_PI)
              angle = 2*M_PI-angle;

           return angle;
        };

        // find first corner point
        while(getAngle(chain, ifirst) < 1e-5 )
           ifirst = (ifirst+1)%n;

        // Copy points of closed loop in correct order, starting with a corner
        Array<PointIndex> pis(n+1);
        pis.Range(0, n-ifirst) = chain.Range(ifirst, n);
        pis.Range(n-ifirst, n) = chain.Range(0, n-ifirst);
        pis[n] = pis[0];

        Array<double> lengths(n);

        for(auto i : Range(n))
           lengths[i] = (mesh[pis[(i+1)%n]] - mesh[pis[i]]).Length();

        auto averageGrowthVectors = [&] (size_t first, size_t last)
        {
           if(first+1 >= last)
              return;

           double total_len = 0.0;
           for(auto l : lengths.Range(first, last))
              total_len += l;

           double len = lengths[first];
           auto v0 = growthvectors[pis[first]];
           auto v1 = growthvectors[pis[last]];

           for(auto i : Range(first+1, last))
           {
              auto pi = pis[i];
              growthvectors[pi] = (len/total_len)*v1 + (1.0-len/total_len)*v0;
              len += lengths[i];
           }
        };

        auto icurrent = 0;

        while(icurrent<n)
        {
           auto ilast = icurrent+1;

           while(getAngle(pis, ilast) < 1e-5 && ilast < n)
              ilast++;

           // found straight line -> average growth vectors between end points
           if(icurrent!=ilast)
              averageGrowthVectors(icurrent, ilast);

           icurrent = ilast;
        }
     }

     //////////////////////////////////////////////////////////////////////
     // reduce growthvectors where necessary to avoid overlaps/slim regions
     const auto getSegmentBox = [&] (SegmentIndex segi)
     {
        PointIndex pi0=mesh[segi][0], pi1=mesh[segi][1];
        Box<3> box( mesh[pi0], mesh[pi1] );
        box.Add( mesh[pi0]+growthvectors[pi0] );
        box.Add( mesh[pi1]+growthvectors[pi1] );
        return box;
     };

     Array<double, PointIndex> growth(np);
     growth = 1.0;

     const auto Dot = [](auto a, auto b)
     { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; };

     const auto restrictGrowthVectors = [&] (SegmentIndex segi0, SegmentIndex segi1)
     {
        if(!active_segments.Test(segi0))
           return;

        const auto & seg0 = mesh[segi0];
        const auto & seg1 = mesh[segi1];

        if( (seg0.domin != domain && seg0.domout != domain) ||
            (seg1.domin != domain && seg1.domout != domain) )
            return;

        if(segi0 == segi1)
           return;

        if(seg0[0]==seg1[0] || seg0[0]==seg1[1] || seg0[1]==seg1[0] || seg0[1] == seg1[1])
           return;

        auto n = mesh[seg0[0]] - mesh[seg0[1]];
        n = {-n[1], n[0], 0};
        n.Normalize();
        if(Dot(n, growthvectors[seg0[0]])<0) n = -n;
        if(Dot(n, growthvectors[seg0[1]])<0) n = -n;

        auto n1 = mesh[seg1[0]] - mesh[seg1[1]];
        n1 = {-n1[1], n1[0], 0};
        n1.Normalize();
        if(Dot(n1, growthvectors[seg1[0]])<0) n1 = -n;
        if(Dot(n1, growthvectors[seg1[1]])<0) n1 = -n;

        auto p10 = mesh[seg1[0]];
        auto p11 = mesh[seg1[1]];

        for ( auto pi : {seg0[0], seg0[1]} )
        {
           if(growthvectors[pi].Length2() == 0.0)
              continue;

           PointIndex pi1 = seg0[0] + seg0[1] - pi;
           auto p1 = mesh[pi1];
           auto p = mesh[pi];

           Point<3> points[] = { p10, p11, p10+total_thickness*growthvectors[seg1[0]], p11+total_thickness*growthvectors[seg1[1]], p1+total_thickness*growthvectors[pi1] };

           Vec<3> gn{ growthvectors[pi][1], -growthvectors[pi][0], 0.0 };
           if(Dot(gn, p1-p) < 0)
              gn = -gn;

           double d0 = Dot(gn, p);
           double d1 = Dot(gn, p1);
           if(d0>d1)
              Swap(d0,d1);

           bool all_left=true, all_right=true;

           for (auto i: Range(4))
           {
              auto p_other = points[i];
              auto dot = Dot(gn,p_other);
              if(dot>d0) all_left = false;
              if(dot<d1) all_right = false;
           }

           if(all_left || all_right)
              return;

           //for ( auto pi : {seg0[0], seg0[1]} )
           {
              double safety = 1.3;
              double t = safety*total_thickness;
              if(growthvectors[pi].Length2() == 0.0)
                 continue;

              Point<3> points[] = { p10, p10+t*growthvectors[seg1[0]], p11, p11+t*growthvectors[seg1[1]] };
              auto p0 = mesh[pi];
              auto p1 = p0 + t*growthvectors[pi];
              auto P2 = [](Point<3> p) { return Point<2>{p[0], p[1]}; };
              ArrayMem<pair<double, double>, 4> intersections;

              double alpha, beta;

              auto checkIntersection = [] (Point<2> p0, Point<2> p1, Point<2> q0, Point<2> q1, double & alpha, double & beta) {
                auto intersection_type = intersect( p0, p1, q0, q1, alpha, beta );
                return intersection_type == X_INTERSECTION || intersection_type == T_INTERSECTION_P || intersection_type == T_INTERSECTION_Q;
              };

              if(checkIntersection( P2(p0), P2(p1), P2(points[0]), P2(points[2]), alpha, beta ))
                 intersections.Append({alpha, 0.0});

              if(checkIntersection( P2(p0), P2(p1), P2(points[1]), P2(points[3]), alpha, beta ))
                 intersections.Append({alpha, 1.0});

              if(checkIntersection( P2(p0), P2(p1), P2(points[0]), P2(points[1]), alpha, beta ))
                 intersections.Append({alpha, beta});

              if(checkIntersection( P2(p0), P2(p1), P2(points[2]), P2(points[3]), alpha, beta ))
                 intersections.Append({alpha, beta});

              QuickSort(intersections);
              for(auto [alpha,beta] : intersections)
              {
                 if(!active_segments.Test(segi1))
                    growth[pi] = min(growth[pi], alpha);
                 else
                 {
                    double mean = 0.5*(alpha+beta);
                    growth[pi] = min(growth[pi], mean);
                    growth[seg1[0]] = min(growth[seg1[0]], mean);
                    growth[seg1[1]] = min(growth[seg1[1]], mean);
                 }
              }
           }
        }
     };

     Box<3> box(Box<3>::EMPTY_BOX);
     for (auto segi : Range(mesh.LineSegments()))
     {
        auto segbox = getSegmentBox( segi );
        box.Add(segbox.PMin());
        box.Add(segbox.PMax());
     }
     BoxTree<3> segtree(box);

     for (auto segi : Range(mesh.LineSegments()))
     {
        auto p2 = [](Point<3> p) { return Point<2>{p[0], p[1]}; };

        auto seg = line_segments[segi];
        double alpha,beta;
        intersect( p2(mesh[seg[0]]), p2(mesh[seg[0]]+total_thickness*growthvectors[seg[0]]), p2(mesh[seg[1]]), p2(mesh[seg[1]]+total_thickness*growthvectors[seg[1]]), alpha, beta );

        if(beta>0 && alpha>0 && alpha<1.1)
           growth[seg[0]] = min(growth[seg[0]], 0.8*alpha);
        if(alpha>0 && beta>0 && beta<1.1)
           growth[seg[1]] = min(growth[seg[1]], 0.8*beta);

        for (auto segj : Range(mesh.LineSegments()))
           if(segi!=segj)
              restrictGrowthVectors(segi, segj);
     }

     for( auto pi : Range(growthvectors))
        growthvectors[pi] *= growth[pi];


     // insert new points
     for(PointIndex pi : Range(mesh.Points()))
        if(growthvectors[pi].Length2()!=0)
        {

           auto & pnew = mapto[pi];
           auto dist = 0.0;
           for(auto t : thicknesses)
           {
              dist+=t;
              pnew.Append( mesh.AddPoint( mesh[pi] + dist*growthvectors[pi] ) );
              mesh[pnew.Last()].SetType(FIXEDPOINT);
           }
        }

     map<pair<PointIndex, PointIndex>, int> seg2edge;

     // insert new elements ( and move old ones )
     for(auto si : moved_segs)
     {
        auto seg = line_segments[si];

        bool swap = false;
        auto & pm0 = mapto[seg[0]];
        auto & pm1 = mapto[seg[1]];

        // auto newindex = si_map[domain];

        Segment s = seg;
        s.geominfo[0] = {};
        s.geominfo[1] = {};
        s[0] = pm0.Last();
        s[1] = pm1.Last();
        s[2] = PointIndex::INVALID;
        auto pair = s[0] < s[1] ? make_pair(s[0], s[1]) : make_pair(s[1], s[0]);
        if(seg2edge.find(pair) == seg2edge.end())
           seg2edge[pair] = ++max_edge_nr;
        s.edgenr = seg2edge[pair];
        s.si = seg.si;
        mesh.AddSegment(s);

        for ( auto i : Range(thicknesses))
        {
           PointIndex pi0, pi1, pi2, pi3;

           if(i==0)
           {
              pi0 = seg[0];
              pi1 = seg[1];
           }
           else
           {
              pi0 = pm0[i-1];
              pi1 = pm1[i-1];
           }

           pi2 = pm1[i];
           pi3 = pm0[i];

           if(i==0)
           {
              auto p0 = mesh[pi0];
              auto p1 = mesh[pi1];
              auto q0 = mesh[pi2];
              // auto q1 = mesh[pi3];

              Vec<2> n = {-p1[1]+p0[1], p1[0]-p0[0]};
              Vec<2> v = { q0[0]-p0[0], q0[1]-p0[1]};
              if(n[0]*v[0]+n[1]*v[1]<0)
                 swap = true;
           }

           Element2d newel;
           newel.SetType(QUAD);
           newel[0] = pi0;
           newel[1] = pi1;
           newel[2] = pi2;
           newel[3] = pi3;
           newel.SetIndex(new_domain);
           newel.GeomInfo() = PointGeomInfo{};

            if(swap)
            {
               Swap(newel[0], newel[1]);
               Swap(newel[2], newel[3]);
            }

           for(auto i : Range(4))
           {
               newel.GeomInfo()[i].u = 0.0;
               newel.GeomInfo()[i].v = 0.0;
           }
           mesh.AddSurfaceElement(newel);

        }
        // segment now adjacent to new 2d-domain!
        if(line_segments[si].domin == domain)
            line_segments[si].domin = new_domain;
        if(line_segments[si].domout == domain)
            line_segments[si].domout = new_domain;
     }

     for(auto pi : Range(mapto))
     {
        if(mapto[pi].Size() == 0)
           continue;
        auto pnew = mapto[pi].Last();
        for(auto old_sei : meshtopo.GetVertexSurfaceElements( pi ))
        {
           if(mesh[old_sei].GetIndex() == domain)
           {
              auto & old_el = mesh[old_sei];
              for(auto i : IntRange(old_el.GetNP()))
                 if(old_el[i]==pi)
                    old_el[i] = pnew;
           }
        }
     }

     for(auto & sel : mesh.SurfaceElements())
        if(sel.GetIndex() == domain)
           sel.Delete();

     mesh.Compress();
     mesh.CalcSurfacesOfNode();

     Generate2dMesh(mesh, domain);

     // even without new domain, we need temporarily a new domain to mesh the remaining area, without confusing the meshes with quads -> add segments temporarily and reset domain number and segments afterwards
     if(!should_make_new_domain)
     {
        // map new domain back to old one
        for(auto & sel : mesh.SurfaceElements())
           if(sel.GetIndex()==new_domain)
              sel.SetIndex(domain);

        // remove (temporary) inner segments
        for(auto segi : Range(first_new_seg, mesh.LineSegments().Range().Next()))
        {
           mesh[segi][0].Invalidate();
           mesh[segi][1].Invalidate();
        }

        for(auto segi : moved_segs)
        {
            if(mesh[segi].domin == new_domain)
                mesh[segi].domin = domain;
            if(mesh[segi].domout == new_domain)
                mesh[segi].domout = domain;
        }

        mesh.Compress();
        mesh.CalcSurfacesOfNode();
     }

     return new_domain;
   }

} // namespace netgen

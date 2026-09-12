#include <mystdlib.h>
#include <regex>

#include "boundarylayer.hpp"
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
      Array<PointIndex, PointIndex> mapto(np);

      bndnodes.Clear();
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.GetEdgeDescriptor(mesh.LineSegment(i).GetIndex()).EdgeNr();
         cout << "snr = " << snr << endl;
         if (snr == surfid)
         {
            bndnodes.Set (mesh.LineSegment(i)[0]);
            bndnodes.Set (mesh.LineSegment(i)[1]);
         }
      }
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.GetEdgeDescriptor(mesh.LineSegment(i).GetIndex()).EdgeNr();
         if (snr != surfid)
         {
            bndnodes.Clear (mesh.LineSegment(i)[0]);
            bndnodes.Clear (mesh.LineSegment(i)[1]);
         }
      }

      for (PointIndex pi : mesh.Points().Range())
        mapto[pi] = bndnodes.Test(pi) ? mesh.AddPoint (mesh[pi]) : PointIndex(PointIndex::INVALID);

      for (i = 1; i <= mesh.GetNSE(); i++)
      {
         Element2d & el = mesh.SurfaceElement(i);
         for (int j = 1; j <= el.GetNP(); j++)
            if (mapto[el.PNum(j)].IsValid())
               el.PNum(j) = mapto[el.PNum(j)];
      }


      int nq = 0;
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.GetEdgeDescriptor(mesh.LineSegment(i).GetIndex()).EdgeNr();
         if (snr == surfid)
         {
            PointIndex p1 = mesh.LineSegment(i)[0];
            PointIndex p2 = mesh.LineSegment(i)[1];
            PointIndex p3 = mapto[p1];
            if (!p3.IsValid()) p3 = p1;
            PointIndex p4 = mapto[p2];
            if (!p4.IsValid()) p4 = p2;

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


  static void SetEdgeDomains (EdgeDescriptor & ed, int dom0, int dom1)
  {
     ed.SetSurfNr(0, dom0);
     ed.SetSurfNr(1, dom1);
     ed.SetDomainIn(dom0);
     ed.SetDomainOut(dom1);
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

  BoundaryLayer2dInfo InsertBoundaryLayer2d (Mesh & mesh, int domain, const Array<double> & thicknesses, bool should_make_new_domain, const Array<int> & boundaries)
  {
     mesh.GetTopology().SetBuildVertex2Element(true);
     mesh.UpdateTopology();
     const auto & line_segments = mesh.LineSegments();

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
     std::map<int, int> edge_to_new_edge;

     Array<Vec<3>, PointIndex> growthvectors(np);
     growthvectors = 0.;

     auto & meshtopo = mesh.GetTopology();

     Array<SegmentIndex> segments;

    // int fd_old = mesh.GetNFD();

    int max_edge_nr = -1;
    int max_domain = -1;

    for(const auto& seg : line_segments)
    {
      if(mesh.GetEdgeDescriptor(seg.GetIndex()).EdgeNr() > max_edge_nr)
        max_edge_nr = mesh.GetEdgeDescriptor(seg.GetIndex()).EdgeNr();
      const auto & ed2 = mesh.GetEdgeDescriptor(seg.GetIndex());
      if(ed2.SurfNr(0) > max_domain)
         max_domain = ed2.SurfNr(0);
      if(ed2.SurfNr(1) > max_domain)
         max_domain = ed2.SurfNr(1);
    }

    int new_domain = max_domain+1;
    int next_edge_nr = max_edge_nr+1;

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
       const auto & ed3 = mesh.GetEdgeDescriptor(seg.GetIndex());
       if(active_boundaries.Test(mesh.GetEdgeDescriptor(seg.GetIndex()).EdgeNr()) && (ed3.SurfNr(0)==domain || ed3.SurfNr(1)==domain))
          active_segments.SetBit(segi);
    }

    {
        while(mesh.GetNFD() < new_domain)
        {
           FaceDescriptor fd(0, 0, 0, -1);
           fd.SetBCProperty(mesh.GetNFD()+1);
           mesh.AddFaceDescriptor(fd);
        }
        if(should_make_new_domain)
           mesh.SetMaterial(new_domain, "layer_" + mesh.GetMaterial(domain));
    }

    for(auto segi : Range(line_segments))
      {
        if(segs_done[segi]) continue;
        segs_done.SetBit(segi);
        const auto& seg = line_segments[segi];
        const auto & ed4 = mesh.GetEdgeDescriptor(seg.GetIndex());
        if(ed4.SurfNr(0) != domain && ed4.SurfNr(1) != domain) continue;
        if(!active_boundaries.Test(mesh.GetEdgeDescriptor(seg.GetIndex()).EdgeNr()))
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

       if(mesh.GetEdgeDescriptor(seg.GetIndex()).SurfNr(1) == domain)
           n = -n;

       AddDirection(growthvectors[seg[0]], n);
       AddDirection(growthvectors[seg[1]], n);
     }

     BitArray is_junction(np+1);
     is_junction.Clear();

     for(auto segi : Range(line_segments))
     {
        if(active_segments.Test(segi)) continue;
        const auto & seg = line_segments[segi];
        const auto & ed = mesh.GetEdgeDescriptor(seg.GetIndex());
        if(ed.SurfNr(0) != domain && ed.SurfNr(1) != domain) continue;

        for(auto pi : {seg[0], seg[1]})
        {
           bool on_layer = false;
           for(auto sj : meshtopo.GetVertexSegments(pi))
              if(active_segments.Test(sj))
                 on_layer = true;
           if(!on_layer) continue;

           auto n = growthvectors[pi];
           if(n.Length2() == 0.0) continue;
           n.Normalize();

           auto other = seg[0] == pi ? seg[1] : seg[0];
           auto t = mesh[other] - mesh[pi];
           if(t.Length2() == 0.0) continue;
           t.Normalize();

           auto tn = t[0]*n[0] + t[1]*n[1] + t[2]*n[2];
           if(tn < 1e-8)
              throw Exception("Boundary layer in 2d: boundary " + ToString(ed.EdgeNr())
                              + " has no layer and does not lead into domain "
                              + ToString(domain) + " at the end of the layer");

           growthvectors[pi] = (1.0/tn) * t;
           is_junction.SetBit(pi);
        }
     }

     //////////////////////////////////////////////////////////////////////////
     // average growthvectors along straight lines to avoid overlaps in corners
     BitArray points_done(np+1);
     points_done.Clear();

     auto nextActive = [&] (PointIndex pi, SegmentIndex from, SegmentIndex & to)
     {
        for(auto sj : meshtopo.GetVertexSegments(pi))
        {
           if(!active_segments.Test(sj) || sj == from) continue;
           to = sj;
           const auto & sg = mesh[sj];
           return PointIndex(sg[0] - pi + sg[1]);
        }
        return PointIndex(PointIndex::INVALID);
     };

     for(auto si : moved_segs)
     {
        if(points_done.Test(line_segments[si][0]))
           continue;

        Array<PointIndex> chain;
        bool closed = false;
        {
           chain.Append(line_segments[si][0]);
           chain.Append(line_segments[si][1]);

           auto current_si = si;
           auto current = line_segments[si][1];
           while(true)
           {
              SegmentIndex next_si;
              auto next = nextActive(current, current_si, next_si);
              if(!next.IsValid()) break;
              if(next == chain[0]) { closed = true; break; }
              chain.Append(next);
              current = next;
              current_si = next_si;
           }

           if(!closed)
           {
              current_si = si;
              current = line_segments[si][0];
              while(true)
              {
                 SegmentIndex next_si;
                 auto next = nextActive(current, current_si, next_si);
                 if(!next.IsValid()) break;
                 chain.Insert(0, next);
                 current = next;
                 current_si = next_si;
              }
           }
        }

        for(auto pi : chain)
           points_done.SetBit(pi);

        auto n = chain.Size();

        // angle between the segments meeting at a[i], the ends of an open chain
        // count as corners
        auto getAngle = [&mesh, closed] (FlatArray<PointIndex> a, size_t i)
        {
           auto n = a.Size();
           if(!closed && (i == 0 || i+1 >= n))
              return M_PI;

           auto p0 = mesh[a[(i+n-1)%n]];
           auto p1 = mesh[a[i]];
           auto p2 = mesh[a[(i+1)%n]];

           auto v0 = p1-p0;
           auto v1 = p2-p1;

           auto angle = abs(atan2(v1[0], v1[1]) - atan2(v0[0], v0[1]));
           if(angle>M_PI)
              angle = 2*M_PI-angle;

           return angle;
        };

        Array<PointIndex> pis;
        if(closed)
        {
           // start the loop at a corner
           size_t ifirst = 0;
           while(ifirst < n && getAngle(chain, ifirst) < 1e-5)
              ifirst++;
           if(ifirst == n)
              ifirst = 0;  // no corner at all, e.g. a circle

           pis.SetSize(n+1);
           pis.Range(0, n-ifirst) = chain.Range(ifirst, n);
           pis.Range(n-ifirst, n) = chain.Range(0, n-ifirst);
           pis[n] = pis[0];
        }
        else
           pis = chain;

        auto nseg_chain = pis.Size()-1;
        Array<double> lengths(nseg_chain);
        for(auto i : Range(nseg_chain))
           lengths[i] = (mesh[pis[i+1]] - mesh[pis[i]]).Length();

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
              if(is_junction.Test(pi))
                 continue;
              growthvectors[pi] = (len/total_len)*v1 + (1.0-len/total_len)*v0;
              len += lengths[i];
           }
        };

        size_t icurrent = 0;
        while(icurrent < nseg_chain)
        {
           auto ilast = icurrent+1;

           while(ilast < nseg_chain && getAngle(pis, ilast) < 1e-5)
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

        const auto & ed0 = mesh.GetEdgeDescriptor(seg0.GetIndex());
        const auto & ed1 = mesh.GetEdgeDescriptor(seg1.GetIndex());
        if( (ed0.SurfNr(0) != domain && ed0.SurfNr(1) != domain) ||
            (ed1.SurfNr(0) != domain && ed1.SurfNr(1) != domain) )
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

           PointIndex pi1 = seg0[0] - pi + seg0[1];
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
        double alpha=0.0;
        double beta=0.0;
        if (intersect(p2(mesh[seg[0]]), p2(mesh[seg[0]] + total_thickness * growthvectors[seg[0]]), p2(mesh[seg[1]]), p2(mesh[seg[1]] + total_thickness * growthvectors[seg[1]]), alpha, beta))
          {
            if (beta > 0 && alpha > 0 && alpha < 1.1)
              growth[seg[0]] = min(growth[seg[0]], 0.8 * alpha);
            if (alpha > 0 && beta > 0 && beta < 1.1)
              growth[seg[1]] = min(growth[seg[1]], 0.8 * beta);
          }

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

     std::map<int, int> ed_to_bl_ed;
     auto getBLEdgeDescriptor = [&] (int old_ed_idx)
     {
        if(ed_to_bl_ed.find(old_ed_idx) == ed_to_bl_ed.end())
        {
           const auto & ed = mesh.GetEdgeDescriptor(old_ed_idx);
           EdgeDescriptor new_ed;
           new_ed.SetEdgeNr(ed.EdgeNr());
           SetEdgeDomains(new_ed,
                          ed.SurfNr(0) == domain ? new_domain : ed.SurfNr(0),
                          ed.SurfNr(1) == domain ? new_domain : ed.SurfNr(1));
           new_ed.SetName(ed.GetName());
           ed_to_bl_ed[old_ed_idx] = mesh.AddEdgeDescriptor(new_ed);
        }
        return ed_to_bl_ed[old_ed_idx];
     };

     for(PointIndex pi = PointIndex::BASE; pi < np + PointIndex::BASE; pi++)
     {
        if(!is_junction.Test(pi)) continue;
        if(mapto[pi].Size() == 0) continue;

        auto gv = growthvectors[pi];
        double gv_len = gv.Length();
        if(gv_len < 1e-12) continue;
        auto gv_dir = (1.0/gv_len) * gv;
        double bl_len = total_thickness * gv_len;

        // walk along the unlayered boundary, away from the layer
        Array<SegmentIndex> wsegs;
        Array<PointIndex> wpts;
        Array<double> wdist, wparam;
        wpts.Append(pi);
        wdist.Append(0.0);
        wparam.Append(0.0);

        PointIndex current = pi, previous(PointIndex::INVALID);
        int ed_idx = -1;
        for(int iter = 0; iter < nseg; iter++)
        {
           SegmentIndex found(0);
           PointIndex next(PointIndex::INVALID);
           for(SegmentIndex segi(0); segi < nseg; segi++)
           {
              if(active_segments.Test(segi)) continue;
              const auto & sg = line_segments[segi];
              if(sg[0] != current && sg[1] != current) continue;
              const auto & ed = mesh.GetEdgeDescriptor(sg.GetIndex());
              if(ed.SurfNr(0) != domain && ed.SurfNr(1) != domain) continue;
              auto other = sg[0] == current ? sg[1] : sg[0];
              if(other == previous) continue;
              found = segi;
              next = other;
              break;
           }
           if(!next.IsValid()) break;

           const auto & sg = line_segments[found];
           if(ed_idx == -1)
           {
              ed_idx = sg.GetIndex();
              wparam[0] = sg.EPGeomInfo(sg[0] == pi ? 0 : 1).dist;
           }
           else if(sg.GetIndex() != ed_idx)
              break;  // reached the next boundary of the geometry

           auto v = mesh[next] - mesh[pi];
           double d = Dot(v, gv_dir);
           if((v - d*gv_dir).Length() > 1e-8 * max2(d, 1.0))
              throw Exception("Boundary layer in 2d: boundary " + ToString(mesh.GetEdgeDescriptor(ed_idx).EdgeNr())
                              + " is not straight where the layer of domain " + ToString(domain)
                              + " ends on it");

           wsegs.Append(found);
           wpts.Append(next);
           wdist.Append(d);
           wparam.Append(sg.EPGeomInfo(sg[0] == next ? 0 : 1).dist);

           previous = current;
           current = next;
           if(d > 2*bl_len) break;
        }

        // first point beyond the layer end
        size_t j = 1;
        while(j < wdist.Size() && wdist[j] <= bl_len + 1e-12)
           j++;
        if(j >= wdist.Size())
           throw Exception("Boundary layer in 2d: the layer of domain " + ToString(domain)
                           + " is thicker than the boundary it ends on");

        // do not leave a sliver of boundary between layer end and next point
        if(j+1 < wdist.Size() && wdist[j] - bl_len < 0.3*(wdist[j] - wdist[j-1]))
           j++;

        auto paramAt = [&] (double d)
        {
           size_t k = 0;
           while(k+1 < wdist.Size()-1 && wdist[k+1] < d) k++;
           double dd = wdist[k+1] - wdist[k];
           double f = dd > 0 ? (d - wdist[k])/dd : 0.0;
           return (1.0-f)*wparam[k] + f*wparam[k+1];
        };

        bool forward = line_segments[wsegs[0]][0] == pi;
        int bl_ed_idx = getBLEdgeDescriptor(ed_idx);

        // pi -> layer points -> first point beyond the layer
        Array<PointIndex> newpts;
        Array<double> newparam;
        newpts.Append(pi);
        newparam.Append(wparam[0]);
        double dist = 0.0;
        for(auto i : Range(thicknesses))
        {
           dist += thicknesses[i];
           newpts.Append(mapto[pi][i]);
           newparam.Append(paramAt(dist*gv_len));
        }
        newpts.Append(wpts[j]);
        newparam.Append(wparam[j]);

        for(auto i : Range(newpts.Size()-1))
        {
           Segment s = line_segments[wsegs[0]];
           // all but the last piece are under the layer
           s.SetIndex(i+2 < newpts.Size() ? bl_ed_idx : ed_idx);
           int a = forward ? 0 : 1;
           s[a] = newpts[i];
           s[1-a] = newpts[i+1];
           s.EPGeomInfo(a).dist = newparam[i];
           s.EPGeomInfo(1-a).dist = newparam[i+1];
           s.GeomInfo(0) = {};
           s.GeomInfo(1) = {};

           if(i < j)
              mesh[wsegs[i]] = s;   // reuse the segments we replace
           else
              mesh.AddSegment(s);
        }
        // drop the ones we did not reuse
        for(size_t i = newpts.Size()-1; i < size_t(j); i++)
        {
           mesh[wsegs[i]][0].Invalidate();
           mesh[wsegs[i]][1].Invalidate();
        }
     }

     // insert new elements ( and move old ones )
     for(auto si : moved_segs)
     {
        auto seg = line_segments[si];

        bool swap = false;
        auto & pm0 = mapto[seg[0]];
        auto & pm1 = mapto[seg[1]];

        Segment s = seg;
        s.GeomInfo(0) = {};
        s.GeomInfo(1) = {};
        s[0] = pm0.Last();
        s[1] = pm1.Last();
        s[2] = PointIndex::INVALID;


        if(edge_to_new_edge.find(seg.GetIndex()) == edge_to_new_edge.end())
        {
          EdgeDescriptor ed;
          ed.SetEdgeNr(next_edge_nr++);
          auto & orig_ed = mesh.GetEdgeDescriptor(seg.GetIndex());
          // the moved segment separates the layer from the rest of the domain
          SetEdgeDomains(ed,
                         orig_ed.SurfNr(0) == domain ? domain : new_domain,
                         orig_ed.SurfNr(1) == domain ? domain : new_domain);
          ed.SetName("moved_" + orig_ed.GetName());
          edge_to_new_edge[seg.GetIndex()] = mesh.AddEdgeDescriptor(ed);
        }

        s.SetIndex(edge_to_new_edge[seg.GetIndex()]);
        // auto pair = s[0] < s[1] ? make_pair(s[0], s[1]) : make_pair(s[1], s[0]);
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
        auto & old_ed = mesh.GetEdgeDescriptor(line_segments[si].GetIndex());
        SetEdgeDomains(old_ed,
                       old_ed.SurfNr(0) == domain ? new_domain : old_ed.SurfNr(0),
                       old_ed.SurfNr(1) == domain ? new_domain : old_ed.SurfNr(1));
     }

     BoundaryLayer2dInfo info;
     info.domain = domain;
     info.new_domain = new_domain;
     info.make_new_domain = should_make_new_domain;
     info.n_edge_descriptors = max_edge_nr;
     for(auto si : moved_segs)
        if(!info.moved_edge_descriptors.Contains(line_segments[si].GetIndex()))
           info.moved_edge_descriptors.Append(line_segments[si].GetIndex());
     for(auto & [orig_idx, new_idx] : edge_to_new_edge)
        info.front_edge_descriptors.Append(new_idx);
     for(auto & [orig_idx, bl_idx] : ed_to_bl_ed)
     {
        info.bl_edge_descriptors.Append(bl_idx);
        info.bl_edge_descriptors_orig.Append(orig_idx);
     }

     mesh.Compress();
     mesh.CalcSurfacesOfNode();
     mesh.ComputeNVertices();

     return info;
   }

  Array<BoundaryLayer2dInfo> InsertBoundaryLayers2d (Mesh & mesh, const MeshingParameters & mp)
  {
     Array<BoundaryLayer2dInfo> infos;

     for(const auto & blp : mp.boundary_layers)
     {
        Array<double> thicknesses;
        if(auto t = get_if<double>(&blp.thickness); t)
           thicknesses.Append(*t);
        else
           for(auto t : *get_if<std::vector<double>>(&blp.thickness))
              thicknesses.Append(t);

        Array<int> domains;
        if(auto d = get_if<int>(&blp.domain); d)
           domains.Append(*d);
        else if(auto s = get_if<string>(&blp.domain); s)
        {
           std::regex pattern(*s);
           for(int i = 1; i <= mesh.GetNDomains(); i++)
              if(std::regex_match(mesh.GetMaterial(i), pattern))
                 domains.Append(i);
           if(auto geo = mesh.GetGeometry(); geo)
              for(auto i : Range(geo->GetNFaces()))
                 if(std::regex_match(geo->GetFace(i).properties.GetName(), pattern)
                    && !domains.Contains(int(i)+1))
                    domains.Append(int(i)+1);
           if(domains.Size() == 0)
              throw Exception("Boundary layer in 2d: no domain matches '" + *s + "'");
        }
        else
           for(auto d : *get_if<std::vector<int>>(&blp.domain))
              domains.Append(d);

        Array<int> boundaries;
        if(auto b = get_if<int>(&blp.boundary); b)
           boundaries.Append(*b);
        else if(auto s = get_if<string>(&blp.boundary); s)
        {
           std::regex pattern(*s);
           for(int i = 1; i <= mesh.GetNED(); i++)
           {
              const auto & ed = mesh.GetEdgeDescriptor(i);
              if(std::regex_match(ed.GetName(), pattern) && !boundaries.Contains(ed.EdgeNr()))
                 boundaries.Append(ed.EdgeNr());
           }
        }
        else
           for(auto b : *get_if<std::vector<int>>(&blp.boundary))
              boundaries.Append(b);

        if(thicknesses.Size() == 0)
           throw Exception("Boundary layer in 2d: no thickness given");

        bool make_new_domain = blp.new_material.has_value();

        for(auto domain : domains)
        {
           auto info = InsertBoundaryLayer2d(mesh, domain, thicknesses, make_new_domain, boundaries);
           if(make_new_domain)
              if(auto mat = get_if<string>(&*blp.new_material); mat)
                 mesh.SetMaterial(info.new_domain, *mat);
           infos.Append(std::move(info));
        }
     }

     return infos;
  }

  void FinalizeBoundaryLayers2d (Mesh & mesh, FlatArray<BoundaryLayer2dInfo> infos)
  {
     bool any = false;
     int n_edge_descriptors = mesh.EdgeDescriptors().Size();

     for(const auto & info : infos)
     {
        if(info.make_new_domain)
           continue;
        any = true;
        n_edge_descriptors = min2(n_edge_descriptors, info.n_edge_descriptors);

        for(auto & sel : mesh.SurfaceElements())
           if(sel.GetIndex() == info.new_domain)
              sel.SetIndex(info.domain);

        // the segments in front of the layer are interior now
        for(auto segi : Range(mesh.LineSegments()))
           if(info.front_edge_descriptors.Contains(mesh[segi].GetIndex()))
           {
              mesh[segi][0].Invalidate();
              mesh[segi][1].Invalidate();
           }

        // the boundary under the layer belongs to the domain again
        for(auto i : Range(info.bl_edge_descriptors))
           for(auto segi : Range(mesh.LineSegments()))
              if(mesh[segi].GetIndex() == info.bl_edge_descriptors[i])
                 mesh[segi].SetIndex(info.bl_edge_descriptors_orig[i]);

        for(auto edi : info.moved_edge_descriptors)
        {
           auto & ed = mesh.GetEdgeDescriptor(edi);
           SetEdgeDomains(ed,
                          ed.SurfNr(0) == info.new_domain ? info.domain : ed.SurfNr(0),
                          ed.SurfNr(1) == info.new_domain ? info.domain : ed.SurfNr(1));
        }
     }

     if(!any)
        return;

     mesh.Compress();
     mesh.EdgeDescriptors().SetSize(n_edge_descriptors);
     mesh.CalcSurfacesOfNode();
  }

} // namespace netgen

#include <mystdlib.h>
#include "meshing.hpp"

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

  void GenerateBoundaryLayer(Mesh& mesh, const BoundaryLayerParameters& blp)
  {
    static Timer timer("Create Boundarylayers");
    RegionTimer regt(timer);

    int max_edge_nr = -1;
    for(const auto& seg : mesh.LineSegments())
      if(seg.edgenr > max_edge_nr)
        max_edge_nr = seg.edgenr;

    int new_mat_nr = mesh.GetNDomains() +1;
    mesh.SetMaterial(new_mat_nr, blp.new_mat);

    auto domains = blp.domains;
    if(!blp.outside)
      domains.Invert();

    mesh.UpdateTopology();
    auto& meshtopo = mesh.GetTopology();

    int np = mesh.GetNP();
    int ne = mesh.GetNE();
    int nse = mesh.GetNSE();
    int nseg = mesh.GetNSeg();

    Array<Array<PointIndex>, PointIndex> mapto(np);

    Array<Vec<3>, PointIndex> growthvectors(np);
    growthvectors = 0.;

    Array<double> surfacefacs(mesh.GetNFD()+1);
    surfacefacs = 0.;

    auto getSurfaceNormal = [&mesh] (const Element2d& el)
    {
      auto v0 = mesh[el[0]];
      return Cross(mesh[el[1]]-v0, mesh[el[2]]-v0).Normalize();
    };

    // surface index map
    Array<int> si_map(mesh.GetNFD()+1);
    si_map = -1;

    int fd_old = mesh.GetNFD();

    // create new FaceDescriptors
    for(auto i : Range(1, fd_old+1))
      {
        const auto& fd = mesh.GetFaceDescriptor(i);
        string name = fd.GetBCName();
        if(blp.surfid.Contains(i))
          {
            if(auto isIn = domains.Test(fd.DomainIn()); isIn != domains.Test(fd.DomainOut()))
              {
                int new_si = mesh.GetNFD()+1;
                surfacefacs[i] = isIn ? 1. : -1.;
                // -1 surf nr is so that curving does not do anything
                FaceDescriptor new_fd(-1, isIn ? new_mat_nr : fd.DomainIn(),
                                      isIn ? fd.DomainOut() : new_mat_nr, -1);
                new_fd.SetBCProperty(new_si);
                mesh.AddFaceDescriptor(new_fd);
                si_map[i] = new_si;
                mesh.SetBCName(new_si-1, "mapped_" + name);
              }
          }
      }

    // mark points for remapping
    for(const auto& sel : mesh.SurfaceElements())
      {
        auto n = surfacefacs[sel.GetIndex()] * getSurfaceNormal(sel);
        if(n.Length2() != 0.)
          {
            for(auto pi : sel.PNums())
              {
                auto & np = growthvectors[pi];
                if(np.Length() == 0) { np = n; continue; }
                auto npn = np * n;
                auto npnp = np * np;
                auto nn = n * n;
                if(nn-npn*npn/npnp == 0) { np = n; continue; }
                np += (nn - npn)/(nn - npn*npn/npnp) * (n - npn/npnp * np);
              }
          }
      }

    // Bit array to keep track of segments already processed
    BitArray segs_done(nseg);
    segs_done.Clear();

    // map for all segments with same points
    // points to pair of SegmentIndex, int
    // int is type of other segment, either:
    // 0 == adjacent surface grows layer
    // 1 == adjacent surface doesn't grow layer, but layer ends on it
    // 2 == adjacent surface is interior surface that ends on layer
    // 3 == adjacent surface is exterior surface that ends on layer (not allowed yet)
    Array<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap(mesh.GetNSeg());

    // moved segments
    Array<SegmentIndex> moved_segs;

    // boundaries to project endings to
    BitArray project_boundaries(fd_old+1);
    BitArray move_boundaries(fd_old+1);
    project_boundaries.Clear();
    move_boundaries.Clear();

    Array<SurfaceElementIndex, SegmentIndex> seg2surfel(mesh.GetNSeg());
    for(auto si : Range(mesh.SurfaceElements()))
      {
        NgArray<int> surfeledges;
        meshtopo.GetSurfaceElementEdges(si+1, surfeledges);
        for(auto edgenr : surfeledges)
          for(auto sei : Range(mesh.LineSegments()))
            if(meshtopo.GetEdge(sei)+1 == edgenr &&
               mesh[sei].si == mesh[si].GetIndex())
              seg2surfel[sei] = si;
      }

    for(auto si : Range(mesh.LineSegments()))
      {
        if(segs_done[si]) continue;
        const auto& segi = mesh[si];
        if(si_map[segi.si] == -1) continue;
        segs_done.SetBit(si);
        segmap[si].Append(make_pair(si, 0));
        moved_segs.Append(si);
        for(auto sj : Range(mesh.LineSegments()))
          {
            if(segs_done.Test(sj)) continue;
            const auto& segj = mesh[sj];
            if((segi[0] == segj[0] && segi[1] == segj[1]) ||
               (segi[0] == segj[1] && segi[1] == segj[0]))
              {
                segs_done.SetBit(sj);
                int type;
                if(si_map[segj.si] != -1)
                  type = 0;
                else if(const auto& fd = mesh.GetFaceDescriptor(segj.si); domains.Test(fd.DomainIn()) && domains.Test(fd.DomainOut()))
                  {
                    type = 2;
                    if(fd.DomainIn() == 0 || fd.DomainOut() == 0)
                      project_boundaries.SetBit(segj.si);
                  }
                else if(const auto& fd = mesh.GetFaceDescriptor(segj.si); !domains.Test(fd.DomainIn()) && !domains.Test(fd.DomainOut()))
                  {
                    type = 3;
                    if(fd.DomainIn() == 0 || fd.DomainOut() == 0)
                      project_boundaries.SetBit(segj.si);
                    move_boundaries.SetBit(segj.si);
                  }
                else
                  {
                    type = 1;
                    // in case 1 we project the growthvector onto the surface
                    project_boundaries.SetBit(segj.si);
                  }
                segmap[si].Append(make_pair(sj, type));
              }
          }
      }

    BitArray in_surface_direction(fd_old+1);
    in_surface_direction.Clear();

    // project growthvector on surface for inner angles
    if(blp.grow_edges)
      {
        for(const auto& sel : mesh.SurfaceElements())
          if(project_boundaries.Test(sel.GetIndex()))
            {
              auto n = getSurfaceNormal(sel);
              for(auto i : Range(sel.PNums()))
                {
                  auto pi = sel.PNums()[i];
                  if(growthvectors[pi].Length2() == 0.)
                    continue;
                  auto next = sel.PNums()[(i+1)%sel.GetNV()];
                  auto prev = sel.PNums()[i == 0 ? sel.GetNV()-1 : i-1];
                  auto v1 = (mesh[next] - mesh[pi]).Normalize();
                  auto v2 = (mesh[prev] - mesh[pi]).Normalize();
                  auto v3 = growthvectors[pi];
                  v3.Normalize();
                  if((v1 * v3 > 1e-12) || (v2 * v3 > 1e-12))
                    in_surface_direction.SetBit(sel.GetIndex());

                  auto& g = growthvectors[pi];
                  auto ng = n * g;
                  auto gg = g * g;
                  auto nn = n * n;
                  // if(fabs(ng*ng-nn*gg) < 1e-12 || fabs(ng) < 1e-12) continue;
                  auto a = -ng*ng/(ng*ng-nn * gg);
                  auto b = ng*gg/(ng*ng-nn*gg);
                  g += a*g + b*n;
                }
            }
      }
    else
      {
        for(const auto& seg : mesh.LineSegments())
          {
            int count = 0;
            for(const auto& seg2 : mesh.LineSegments())
              if(((seg[0] == seg2[0] && seg[1] == seg2[1]) || (seg[0] == seg2[1] && seg[1] == seg2[0])) && blp.surfid.Contains(seg2.si))
                count++;
            if(count == 1)
              {
                growthvectors[seg[0]] = {0., 0., 0.};
                growthvectors[seg[1]] = {0., 0., 0.};
              }
          }
      }

    // insert new points
    for (PointIndex pi = 1; pi <= np; pi++)
      if (growthvectors[pi].Length2() != 0)
        {
          Point<3> p = mesh[pi];
          for(auto i : Range(blp.heights))
            {
              p += blp.heights[i] * growthvectors[pi];
              mapto[pi].Append(mesh.AddPoint(p));
            }
        }

    // add 2d quads on required surfaces
    map<pair<PointIndex, PointIndex>, int> seg2edge;
    if(blp.grow_edges)
      {
        for(auto sei : moved_segs)
          {
            // copy here since we will add segments and this would
            // invalidate a reference!
            auto segi = mesh[sei];
            for(auto [sej, type] : segmap[sei])
              {
                auto segj = mesh[sej];
                if(type == 0)
                  {
                    Segment s;
                    s[0] = mapto[segj[0]].Last();
                    s[1] = mapto[segj[1]].Last();
                    s[2] = PointIndex::INVALID;
                    auto pair = s[0] < s[1] ? make_pair(s[0], s[1]) : make_pair(s[1], s[0]);
                    if(seg2edge.find(pair) == seg2edge.end())
                      seg2edge[pair] = ++max_edge_nr;
                    s.edgenr = seg2edge[pair];
                    s.si = si_map[segj.si];
                    mesh.AddSegment(s);
                  }
                // here we need to grow the quad elements
                else if(type == 1)
                  {
                    PointIndex pp1 = segj[1];
                    PointIndex pp2 = segj[0];
                    if(in_surface_direction.Test(segj.si))
                      {
                        Swap(pp1, pp2);
                        move_boundaries.SetBit(segj.si);
                      }
                    PointIndex p1 = pp1;
                    PointIndex p2 = pp2;
                    PointIndex p3, p4;
                    Segment s0;
                    s0[0] = p1;
                    s0[1] = p2;
                    s0[2] = PointIndex::INVALID;
                    s0.edgenr = segj.edgenr;
                    s0.si = segj.si;
                    mesh.AddSegment(s0);

                    for(auto i : Range(blp.heights))
                      {
                        Element2d sel(QUAD);
                        p3 = mapto[pp2][i];
                        p4 = mapto[pp1][i];
                        sel[0] = p1;
                        sel[1] = p2;
                        sel[2] = p3;
                        sel[3] = p4;
                        sel.SetIndex(segj.si);
                        mesh.AddSurfaceElement(sel);

                        // TODO: Too many, would be enough to only add outermost ones
                        Segment s1;
                        s1[0] = p2;
                        s1[1] = p3;
                        s1[2] = PointIndex::INVALID;
                        auto pair = make_pair(p2, p3);
                        if(seg2edge.find(pair) == seg2edge.end())
                          seg2edge[pair] = ++max_edge_nr;
                        s1.edgenr = seg2edge[pair];
                        s1.si = segj.si;
                        mesh.AddSegment(s1);
                        Segment s2;
                        s2[0] = p4;
                        s2[1] = p1;
                        s2[2] = PointIndex::INVALID;
                        pair = make_pair(p1, p4);
                        if(seg2edge.find(pair) == seg2edge.end())
                          seg2edge[pair] = ++max_edge_nr;
                        s2.edgenr = seg2edge[pair];
                        s2.si = segj.si;
                        mesh.AddSegment(s2);
                        p1 = p4;
                        p2 = p3;
                      }
                    Segment s3;
                    s3[0] = p3;
                    s3[1] = p4;
                    s3[2] = PointIndex::INVALID;
                    auto pair = p3 < p4 ? make_pair(p3, p4) : make_pair(p4, p3);
                    if(seg2edge.find(pair) == seg2edge.end())
                      seg2edge[pair] = ++max_edge_nr;
                    s3.edgenr = seg2edge[pair];
                    s3.si = segj.si;
                    mesh.AddSegment(s3);
                  }
              }
          }
      }

    BitArray fixed_points(np+1);
    fixed_points.Clear();
    BitArray moveboundarypoint(np+1);
    moveboundarypoint.Clear();
    for(SurfaceElementIndex si = 0; si < nse; si++)
      {
        // copy because surfaceels array will be resized!
        auto sel = mesh[si];
        if(si_map[sel.GetIndex()] != -1)
          {
            Array<PointIndex> points(sel.PNums());
            if(surfacefacs[sel.GetIndex()] > 0) Swap(points[0], points[2]);
            for(auto j : Range(blp.heights))
              {
                auto eltype = points.Size() == 3 ? PRISM : HEX;
                Element el(eltype);
                for(auto i : Range(points))
                    el[i] = points[i];
                for(auto i : Range(points))
                  points[i] = mapto[sel.PNums()[i]][j];
                if(surfacefacs[sel.GetIndex()] > 0) Swap(points[0], points[2]);
                for(auto i : Range(points))
                  el[sel.PNums().Size() + i] = points[i];
                el.SetIndex(new_mat_nr);
                mesh.AddVolumeElement(el);
              }
            Element2d newel = sel;
            for(auto& p : newel.PNums())
              p = mapto[p].Last();
            newel.SetIndex(si_map[sel.GetIndex()]);
            mesh.AddSurfaceElement(newel);
          }
        else
          {
            bool has_moved = false;
            for(auto p : sel.PNums())
              if(mapto[p].Size())
                has_moved = true;
            if(has_moved)
              for(auto p : sel.PNums())
                {
                  if(!mapto[p].Size())
                    {
                      fixed_points.SetBit(p);
                      if(move_boundaries.Test(sel.GetIndex()))
                        moveboundarypoint.SetBit(p);
                    }
                }
          }
        if(move_boundaries.Test(sel.GetIndex()))
          {
            for(auto& p : mesh[si].PNums())
              if(mapto[p].Size())
                p = mapto[p].Last();
          }
      }

    for(SegmentIndex sei = 0; sei < nseg; sei++)
      {
        auto& seg = mesh[sei];
        if(move_boundaries.Test(seg.si))
          for(auto& p : seg.PNums())
            if(mapto[p].Size())
              p = mapto[p].Last();
      }

    for(ElementIndex ei = 0; ei < ne; ei++)
      {
        auto el = mesh[ei];
        ArrayMem<PointIndex,4> fixed;
        ArrayMem<PointIndex,4> moved;
        bool moved_bnd = false;
        for(const auto& p : el.PNums())
          {
            if(fixed_points.Test(p))
              fixed.Append(p);
            if(mapto[p].Size())
              moved.Append(p);
            if(moveboundarypoint.Test(p))
              moved_bnd = true;
          }

        bool do_move, do_insert;
        if(domains.Test(el.GetIndex()))
          {
            do_move = fixed.Size() && moved_bnd;
            do_insert = do_move;
          }
        else
          {
            do_move = !fixed.Size() || moved_bnd;
            do_insert = !do_move;
          }

        if(do_move)
          {
            for(auto& p : mesh[ei].PNums())
              if(mapto[p].Size())
                p = mapto[p].Last();
          }
        if(do_insert)
          {
            if(el.GetType() != TET)
              throw Exception("Boundarylayer only implemented for tets outside yet!");
            if(moved.Size() == 2)
              {
                if(fixed.Size() == 2)
                  throw Exception("This should not be possible!");
                PointIndex p1 = moved[0];
                PointIndex p2 = moved[1];
                for(auto i : Range(blp.heights))
                  {
                    PointIndex p3 = mapto[moved[1]][i];
                    PointIndex p4 = mapto[moved[0]][i];
                    Element nel(PYRAMID);
                    nel[0] = p1;
                    nel[1] = p2;
                    nel[2] = p3;
                    nel[3] = p4;
                    nel[4] = el[0] + el[1] + el[2] + el[3] - fixed[0] - moved[0] - moved[1];
                    if(Cross(mesh[p2]-mesh[p1], mesh[p4]-mesh[p1]) * (mesh[nel[4]]-mesh[nel[1]]) > 0)
                      Swap(nel[1], nel[3]);
                    nel.SetIndex(el.GetIndex());
                    mesh.AddVolumeElement(nel);
                    p1 = p4;
                    p2 = p3;
                  }
              }
            if(moved.Size() == 1 && fixed.Size() == 1)
              {
                PointIndex p1 = moved[0];
                for(auto i : Range(blp.heights))
                  {
                    Element nel = el;
                    PointIndex p2 = mapto[moved[0]][i];
                    for(auto& p : nel.PNums())
                      {
                        if(p == moved[0])
                          p = p1;
                        else if(p == fixed[0])
                          p = p2;
                      }
                    p1 = p2;
                    mesh.AddVolumeElement(nel);
                  }
              }
          }
      }

    for(auto i : Range(1, fd_old+1))
      if(si_map[i] != -1)
        {
          if(mesh.GetFaceDescriptor(mesh.GetNFD()).DomainIn() == new_mat_nr)
            mesh.GetFaceDescriptor(i).SetDomainOut(new_mat_nr);
          else
            mesh.GetFaceDescriptor(i).SetDomainIn(new_mat_nr);
        }
  }
}

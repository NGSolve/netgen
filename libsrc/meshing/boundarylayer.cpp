#include <mystdlib.h>
#include "meshing.hpp"
#include "debugging.hpp"
#include "global.hpp"

#include <set>
#include <regex>

namespace netgen
{
  // checks if a segment is intersecting a plane, spanned by three points, lam will be set s.t. p_intersect = seg[0] + lam * (seg[1]-seg[0])
  bool isIntersectingPlane ( const array<Point<3>, 2> & seg, const array<Point<3>, 3> & trig, double & lam)
  {
      auto n = Cross(trig[1]-trig[0], trig[2]-trig[0]);
      auto v0n = (seg[0]-trig[0])*n;
      auto v1n = (seg[1]-trig[0])*n;
      if(v0n * v1n >= 0)
          return false;

      lam = -v0n/(v1n-v0n);
      lam *= 0.9;
      if(lam < -1e-8 || lam>1+1e-8)
          return false;
      return true;
  }

  bool isIntersectingPlane ( const array<Point<3>, 2> & seg, const ArrayMem<Point<3>, 4> & face, double & lam)
  {
     lam = 1.0;
     bool intersect0 = isIntersectingPlane( seg, array<Point<3>, 3>{face[0], face[1], face[2]}, lam );
     if(face.Size()==3)
         return intersect0;

     double lam1 = 1.0;
     bool intersect1 = isIntersectingPlane( seg, array<Point<3>, 3>{face[2], face[3], face[0]}, lam1 );
     lam = min(lam, lam1);
     return intersect0 || intersect1;
  }

  bool isIntersectingTrig ( const array<Point<3>, 2> & seg, const array<Point<3>, 3> & trig, double & lam)
  {
      if(!isIntersectingPlane(seg, trig, lam))
          return false;


      //buffer enlargement of triangle
      auto pt0 = trig[0];
      auto pt1 = trig[1];
      auto pt2 = trig[2];
      Point<3> center = { (pt0[0] + pt1[0] + pt2[0]) / 3.0, (pt0[1] + pt1[1] + pt2[1]) / 3.0, (pt0[2] + pt1[2] + pt2[2]) / 3.0 };
      array<Point<3>, 3>  larger_trig = {
      center + (pt0 - center) * 1.1,
      center + (pt1 - center) * 1.1,
      center + (pt2 - center) * 1.1, };

      auto p = seg[0] + lam/0.9*(seg[1]-seg[0]);

      auto n_trig = Cross(trig[1]-trig[0], trig[2]-trig[0]).Normalize();
      for(auto i : Range(3))
      {
          // check if p0 and p are on same side of segment p1-p2
          auto p0 = larger_trig[i];
          auto p1 = larger_trig[(i+1)%3];
          auto p2 = larger_trig[(i+2)%3];
          auto n = Cross(p2-p1, n_trig);

          auto v0 = (p2-p1).Normalize();
          auto v1 = (p0-p1).Normalize();
          auto inside_dir = (v1 - (v1*v0) * v0).Normalize();
          auto v2 = (p-p1).Normalize();
          if(inside_dir * v1 < 0)
              inside_dir = -inside_dir;

          if( (inside_dir*v2) < 0 )
              return false;
      }
      return true;
  };

  bool isIntersectingFace( const array<Point<3>, 2> & seg,  const ArrayMem<Point<3>, 4> & face, double & lam )
  {
      lam = 1.0;
      double lam0 = 1.0;
      bool intersect0 = isIntersectingTrig( seg, {face[0], face[1], face[2]}, lam0 );
      if(intersect0)
          lam = min(lam, lam0);
      if(face.Size()==3)
          return intersect0;

      double lam1 = 1.0;
      bool intersect1 = isIntersectingTrig( seg, {face[2], face[3], face[0]}, lam1 );
      if(intersect1)
          lam = min(lam, lam1);
      return intersect0 || intersect1;
  }

  array<Point<3>, 2> BoundaryLayerTool :: GetMappedSeg( PointIndex pi )
  {
      return { mesh[pi], mesh[pi] + height*limits[pi]*growthvectors[pi] * 1.5 };
  }

  ArrayMem<Point<3>, 4> BoundaryLayerTool :: GetFace( SurfaceElementIndex sei )
  {
      const auto & sel = mesh[sei];
      ArrayMem<Point<3>, 4> points(sel.GetNP());
      for(auto i : Range(sel.GetNP()))
          points[i] = mesh[sel[i]];
      return points;
  }

  ArrayMem<Point<3>, 4> BoundaryLayerTool :: GetMappedFace( SurfaceElementIndex sei )
  {
      const auto & sel = mesh[sei];
      ArrayMem<Point<3>, 4> points(sel.GetNP());
      for(auto i : Range(sel.GetNP()))
          points[i] = mesh[sel[i]] + height * limits[sel[i]]*growthvectors[sel[i]];
      return points;
  }

  ArrayMem<Point<3>, 4> BoundaryLayerTool :: GetMappedFace( SurfaceElementIndex sei, int face )
  {
      if(face == -1) return GetFace(sei);
      if(face == -2) return GetMappedFace(sei);
      const auto & sel = mesh[sei];
      auto np = sel.GetNP();
      auto pi0 = sel[face % np];
      auto pi1 = sel[(face+1) % np];
      ArrayMem<Point<3>, 4> points(4);
      points[0] = points[3] = mesh[pi0];
      points[1] = points[2] = mesh[pi1];
      points[3] += height * limits[pi0]*growthvectors[pi0];
      points[2] += height * limits[pi1]*growthvectors[pi1];
      return points;
  }

  Vec<3> BoundaryLayerTool :: getEdgeTangent(PointIndex pi, int edgenr)
  {
      Vec<3> tangent = 0.0;
      ArrayMem<PointIndex,2> pts;
      for(auto segi : topo.GetVertexSegments(pi))
      {
          auto & seg = mesh[segi];
          if(seg.edgenr != edgenr+1)
              continue;
          PointIndex other = seg[0]+seg[1]-pi;
          if(!pts.Contains(other))
            pts.Append(other);
      }
      if(pts.Size() != 2)
        throw Exception("Something went wrong in getEdgeTangent!");
      tangent = mesh[pts[1]] - mesh[pts[0]];
      return tangent.Normalize();
  }

  void BoundaryLayerTool :: LimitGrowthVectorLengths()
  {
    static Timer tall("BoundaryLayerTool::LimitGrowthVectorLengths"); RegionTimer rtall(tall);

    limits.SetSize(np);
    limits = 1.0;
   
    // Function to calculate the dot product of two 3D vectors
    // Is there netgen native function for this?
    const auto Dot = [](Vec<3> a, Vec<3> b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    auto parallel_limiter = [&](PointIndex pi1, PointIndex pi2, SurfaceElementIndex si) {
        MeshPoint& a_base = mesh[pi1];
        MeshPoint& b_base = mesh[pi2];
        MeshPoint a_end = mesh[pi1] + height * limits[pi1] * growthvectors[pi1];
        MeshPoint b_end = mesh[pi2] + height * limits[pi2] * growthvectors[pi2];

        double ab_base = (b_base - a_base).Length();
        Vec<3> a_vec = (a_end - a_base);
        Vec<3> b_vec = (b_end - b_base);

        // Calculate parallel projections
        Vec<3> ab_base_norm = (b_base - a_base).Normalize();
        double a_vec_x = Dot(a_vec, ab_base_norm);
        double b_vec_x = Dot(b_vec, -ab_base_norm);
        double ratio_parallel = (a_vec_x + b_vec_x) / ab_base;

        double PARALLEL_RATIO_LIMIT = 0.85;
        if (ratio_parallel > PARALLEL_RATIO_LIMIT) {
            // Adjust limits, vectors, and projections if parallel ratio exceeds the limit
            double corrector = PARALLEL_RATIO_LIMIT / ratio_parallel;
            limits[pi1] *= corrector;
            limits[pi2] *= corrector;
        }
    };      
    
    auto perpendicular_limiter = [&](PointIndex pi1, PointIndex pi2, SurfaceElementIndex si) {
        // this part is same as in parallel limiter, but note that limits contents are already changed
        MeshPoint& a_base = mesh[pi1];
        MeshPoint& b_base = mesh[pi2];
        MeshPoint a_end = mesh[pi1] + height * limits[pi1] * growthvectors[pi1];
        MeshPoint b_end = mesh[pi2] + height * limits[pi2] * growthvectors[pi2];

        double ab_base = (b_base - a_base).Length();
        Vec<3> a_vec = (a_end - a_base);
        Vec<3> b_vec = (b_end - b_base);

        // Calculate parallel projections
        Vec<3> ab_base_norm = (b_base - a_base).Normalize();
        double a_vec_x = Dot(a_vec, ab_base_norm);
        double b_vec_x = Dot(b_vec, -ab_base_norm);
        double ratio_parallel = (a_vec_x + b_vec_x) / ab_base;

        // Calculate surface normal at point si
        Vec<3> surface_normal = getNormal(mesh[si]);

        double a_vec_y = abs(Dot(a_vec, surface_normal));
        double b_vec_y = abs(Dot(b_vec, surface_normal));
        double diff_perpendicular = abs(a_vec_y - b_vec_y);
        double tan_alpha = diff_perpendicular / (ab_base - a_vec_x - b_vec_x);

        double TAN_ALPHA_LIMIT = 0.36397; // Approximately 20 degrees in radians
        if (tan_alpha > TAN_ALPHA_LIMIT) {
            if (a_vec_y > b_vec_y) {
                double correction = (TAN_ALPHA_LIMIT / tan_alpha * diff_perpendicular + b_vec_y) / a_vec_y;
                limits[pi1] *= correction;
            }
            else {
                double correction = (TAN_ALPHA_LIMIT / tan_alpha * diff_perpendicular + a_vec_y) / b_vec_y;
                limits[pi2] *= correction;
            }
        }
    };

    auto neighbour_limiter = [&](PointIndex pi1, PointIndex pi2, SurfaceElementIndex si) {
        parallel_limiter(pi1, pi2, si);
        perpendicular_limiter(pi1, pi2, si);
    };
    
    auto modifiedsmooth = [&](size_t nsteps) {
        for (auto i : Range(nsteps))
            for (SurfaceElementIndex sei : mesh.SurfaceElements().Range())
            {
            // assuming triangle
            neighbour_limiter(mesh[sei].PNum(1), mesh[sei].PNum(2), sei);
            neighbour_limiter(mesh[sei].PNum(2), mesh[sei].PNum(3), sei);
            neighbour_limiter(mesh[sei].PNum(3), mesh[sei].PNum(1), sei);  
        }
    };

    auto smooth = [&] (size_t nsteps) {
        for(auto i : Range(nsteps))
            for(const auto & sel : mesh.SurfaceElements())
            {
                double min_limit = 999;
                for(auto pi : sel.PNums())
                    min_limit = min(min_limit, limits[pi]);
                for(auto pi : sel.PNums())
                    limits[pi] = min(limits[pi], 1.4*min_limit);
            }
    };

    // check for self-intersection within new elements (prisms/hexes)
    auto self_intersection = [&] () {
        for(SurfaceElementIndex sei : mesh.SurfaceElements().Range())
          {
            auto facei = mesh[sei].GetIndex();
            if(facei < nfd_old && !params.surfid.Contains(facei))
                continue;

            auto sel = mesh[sei];
            auto np = sel.GetNP();
            // check if a new edge intesects the plane of any opposing face
            double lam;
            for(auto i : Range(np))
                for(auto fi : Range(np-2))
                    if(isIntersectingPlane(GetMappedSeg(sel[i]), GetMappedFace(sei, i+fi+1), lam))
                        if(lam < 1.0)
                            limits[sel[i]] *= lam;
          }
    };

    // first step: intersect with other surface elements that are boundary of domain the layer is grown into
    // second (and subsequent) steps: intersect with other boundary layers, allow restriction by 20% in each step
    auto changed_domains = domains;
      if(!params.outside)
          changed_domains.Invert();

    bool limit_reached = true;
    double lam_lower_limit = 1.0;
    int step = 0;

    while(limit_reached || step<3)
    {
        Array<double, PointIndex> new_limits;
        new_limits.SetSize(np);
        new_limits = 1.0;

        if(step>1)
            lam_lower_limit *= 0.8;
        limit_reached = false;

        // build search tree with all surface elements (bounding box of a surface element also covers the generated boundary layer)
        Box<3> bbox(Box<3>::EMPTY_BOX);
        for(auto pi : mesh.Points().Range())
        {
            bbox.Add(mesh[pi]);
            bbox.Add(mesh[pi]+limits[pi]*height*growthvectors[pi]);
        }
        BoxTree<3> tree(bbox);

        for(auto sei : mesh.SurfaceElements().Range())
        {
            const auto & sel = mesh[sei];
            Box<3> box(Box<3>::EMPTY_BOX);
            const auto& fd = mesh.GetFaceDescriptor(sel.GetIndex());
            if(!changed_domains.Test(fd.DomainIn()) &&
               !changed_domains.Test(fd.DomainOut()))
              continue;
            for(auto pi : sel.PNums())
                box.Add(mesh[pi]);
            // also add moved points to bounding box
            if(params.surfid.Contains(sel.GetIndex()))
                for(auto pi : sel.PNums())
                    box.Add(mesh[pi]+limits[pi]*height*growthvectors[pi]);
            tree.Insert(box, sei);
        }

        for(auto pi : mesh.Points().Range())
        {
            if(mesh[pi].Type() == INNERPOINT)
                continue;
            if(growthvectors[pi].Length2() == 0.0)
                continue;
            Box<3> box(Box<3>::EMPTY_BOX);
            auto seg = GetMappedSeg(pi);
            box.Add(seg[0]);
            box.Add(seg[1]);
            double lam = 1.0;
            tree.GetFirstIntersecting(box.PMin(), box.PMax(), [&](SurfaceElementIndex sei)
              {
                const auto & sel = mesh[sei];
                if(sel.PNums().Contains(pi))
                    return false;
                auto face = GetMappedFace(sei, -2);
                double lam_ = 999;
                bool is_bl_sel = params.surfid.Contains(sel.GetIndex());

                if (step == 0)
                {
                    face = GetMappedFace(sei, -1);
                    if (isIntersectingFace(seg, face, lam_))
                    {
                        if (is_bl_sel)
                            lam_ *= params.limit_safety;
                        lam = min(lam, lam_);
                    }
                }

                if(step==1)
                {
                    if(isIntersectingFace(seg, face, lam_))
                    {
                        if(is_bl_sel) // allow only half the distance if the opposing surface element has a boundary layer too
                            lam_ *= params.limit_safety;
                        lam = min(lam, lam_);
                    }
                }
                // if the opposing surface element has a boundary layer, we need to additionally intersect with the new faces 
                if(step>1 && is_bl_sel)
                {
                    for(auto facei : Range(-1, sel.GetNP()))
                    {
                        auto face = GetMappedFace(sei, facei);
                        if(isIntersectingFace(seg, face, lam_)) // && lam_ > other_limit)
                        {
                            lam = min(lam, lam_);
                        }
                    }
                }
                return false;
              });
            if(lam<1)
            {
                if(lam<lam_lower_limit && step>1)
                {
                    limit_reached = true;
                    lam = lam_lower_limit;
                }
            }

            new_limits[pi] = min(limits[pi], lam* limits[pi]);
        }
        step++;
        limits = new_limits;
        if (step > 0)
           modifiedsmooth(1);
    }

    self_intersection();
    modifiedsmooth(1);

    for(auto pi : Range(growthvectors))
        growthvectors[pi] *= limits[pi];

  }


  // depending on the geometry type, the mesh contains segments multiple times (once for each face)
  bool HaveSingleSegments( const Mesh & mesh )
  {
      auto& topo = mesh.GetTopology();
      NgArray<SurfaceElementIndex> surf_els;

      for(auto segi : Range(mesh.LineSegments()))
      {
          mesh.GetTopology().GetSegmentSurfaceElements(segi+1, surf_els);
          if(surf_els.Size()<2)
              continue;

          auto seg = mesh[segi];
          auto pi0 = min(seg[0], seg[1]);
          auto pi1 = max(seg[0], seg[1]);
          auto p0_segs = topo.GetVertexSegments(seg[0]);

          for(auto segi_other : p0_segs)
          {
              if(segi_other == segi)
                  continue;

              auto seg_other = mesh[segi_other];
              auto pi0_other = min(seg_other[0], seg_other[1]);
              auto pi1_other = max(seg_other[0], seg_other[1]);
              if( pi0_other == pi0 && pi1_other == pi1 )
                  return false;
          }

          // found segment with multiple adjacent surface elements but no other segments with same points -> have single segments
          return true;
      }

      return true;
  }

  // duplicates segments (and sets seg.si accordingly) to have a unified data structure for all geometry types
  Array<Segment> BuildSegments( Mesh & mesh )
  {
      Array<Segment> segments;
      auto& topo = mesh.GetTopology();

      NgArray<SurfaceElementIndex> surf_els;

      for(auto segi : Range(mesh.LineSegments()))
      {
          auto seg = mesh[segi];
          mesh.GetTopology().GetSegmentSurfaceElements(segi+1, surf_els);
          for(auto seli : surf_els)
          {
              const auto & sel = mesh[seli];
              seg.si = sel.GetIndex();

              auto np = sel.GetNP();
              for(auto i : Range(np))
              {
                  if(sel[i] == seg[0])
                  {
                      if(sel[(i+1)%np] != seg[1])
                          swap(seg[0], seg[1]);
                      break;
                  }
              }

              segments.Append(seg);
          }
      }
      return segments;
  }

  void MergeAndAddSegments( Mesh & mesh, FlatArray<Segment> new_segments)
  {
      INDEX_2_HASHTABLE<bool> already_added( mesh.LineSegments().Size() + 2*new_segments.Size() );

      for(auto & seg : mesh.LineSegments())
      {
          INDEX_2 i2 (seg[0], seg[1]);
          i2.Sort();
          if(!already_added.Used(i2))
              already_added.Set(i2, true);
      }

      for(auto & seg : new_segments)
      {
          INDEX_2 i2 (seg[0], seg[1]);
          i2.Sort();

          if(!already_added.Used(i2))
          {
              mesh.AddSegment(seg);
              already_added.Set(i2, true);
          }
      }
  }

  void BoundaryLayerTool :: InterpolateSurfaceGrowthVectors()
  {
    static Timer tall("InterpolateSurfaceGrowthVectors"); RegionTimer rtall(tall);
    static Timer tsmooth("InterpolateSurfaceGrowthVectors-Smoothing");
    auto np = mesh.GetNP();
    BitArray is_point_on_bl_surface(np+1);
    is_point_on_bl_surface.Clear();
    BitArray is_point_on_other_surface(np+1);
    is_point_on_other_surface.Clear();
    Array<Vec<3>, PointIndex> normals(np);
    for(auto pi : Range(growthvectors))
        normals[pi] = growthvectors[pi];

    ParallelForRange( mesh.SurfaceElements().Range(), [&] ( auto myrange )
      {
        for(SurfaceElementIndex sei : myrange)
          {
            auto facei = mesh[sei].GetIndex();
            if(facei < nfd_old && !params.surfid.Contains(facei))
              {
                for(auto pi : mesh[sei].PNums())
                  if(mesh[pi].Type() == SURFACEPOINT)
                    is_point_on_other_surface.SetBitAtomic(pi);
              }
            else
              {
                for(auto pi : mesh[sei].PNums())
                  if(mesh[pi].Type() == SURFACEPOINT)
                    is_point_on_bl_surface.SetBitAtomic(pi);
              }
          }
      });

    Array<PointIndex> points;
    for(PointIndex pi : mesh.Points().Range())
      {
        if(is_point_on_bl_surface[pi])
        {
            points.Append(pi);
            growthvectors[pi] = 0.0;
        }
        if(is_point_on_other_surface[pi])
          {
            points.Append(pi);
          }
      }

    // smooth tangential part of growth vectors from edges to surface elements
    RegionTimer rtsmooth(tsmooth);
    for(auto i : Range(10))
    {
        for(auto pi : points)
        {
            auto sels = p2sel[pi];
            Vec<3> new_gw = growthvectors[pi];
            int cnt = 1;
            std::set<PointIndex> suround;
            suround.insert(pi);
            auto normal = normals[pi];
            for(auto sei: sels)
            {
                const auto & sel = mesh[sei];
                for(auto pi1 : sel.PNums())
                    if(suround.count(pi1)==0)
                    {
                        suround.insert(pi1);
                        auto gw_other = growthvectors[pi1];
                        auto normal_other = getNormal(mesh[sei]);
                        auto tangent_part = gw_other - (gw_other*normal_other)*normal_other;
                        if(is_point_on_bl_surface[pi])
                          new_gw += tangent_part;
                        else
                          new_gw += gw_other;
                    }
            }

            growthvectors[pi] = 1.0/suround.size() * new_gw;
        }
    }

    for(auto pi : points)
        growthvectors[pi] += normals[pi];
  }


  BoundaryLayerTool::BoundaryLayerTool(Mesh & mesh_, const BoundaryLayerParameters & params_)
      : mesh(mesh_), topo(mesh_.GetTopology()), params(params_)
  {
    static Timer timer("BoundaryLayerTool::ctor");
    RegionTimer regt(timer);

    //for(auto & seg : mesh.LineSegments())
        //seg.edgenr = seg.epgeominfo[1].edgenr;

    height = 0.0;
    for (auto h : params.heights)
      height += h;

    max_edge_nr = -1;
    for(const auto& seg : mesh.LineSegments())
      if(seg.edgenr > max_edge_nr)
        max_edge_nr = seg.edgenr;

    int ndom = mesh.GetNDomains();
    ndom_old = ndom;

    new_mat_nrs.SetSize(mesh.FaceDescriptors().Size() + 1);
    new_mat_nrs = -1;
    for(auto [bcname, matname] : params.new_mat)
      {
        mesh.SetMaterial(++ndom, matname);
        regex pattern(bcname);
        for(auto i : Range(1, mesh.GetNFD()+1))
          {
            auto& fd = mesh.GetFaceDescriptor(i);
            if(regex_match(fd.GetBCName(), pattern))
              new_mat_nrs[i] = ndom;
          }
      }

    domains = params.domains;
    if(!params.outside)
      domains.Invert();

    topo.SetBuildVertex2Element(true);
    mesh.UpdateTopology();

    have_single_segments = HaveSingleSegments(mesh);
    if(have_single_segments)
        segments = BuildSegments(mesh);
    else
        segments = mesh.LineSegments();

    np = mesh.GetNP();
    ne = mesh.GetNE();
    nse = mesh.GetNSE();
    nseg = segments.Size();

    p2sel = mesh.CreatePoint2SurfaceElementTable();

    nfd_old = mesh.GetNFD();
    moved_surfaces.SetSize(nfd_old+1);
    moved_surfaces.Clear();
    si_map.SetSize(nfd_old+1);
    for(auto i : Range(nfd_old+1))
      si_map[i] = i;
  }

  void BoundaryLayerTool :: CreateNewFaceDescriptors()
  {
    surfacefacs.SetSize(nfd_old+1);
    surfacefacs = 0.0;
    // create new FaceDescriptors
    for(auto i : Range(1, nfd_old+1))
      {
        const auto& fd = mesh.GetFaceDescriptor(i);
        string name = fd.GetBCName();
        if(params.surfid.Contains(i))
          {
            if(auto isIn = domains.Test(fd.DomainIn()); isIn != domains.Test(fd.DomainOut()))
              {
                int new_si = mesh.GetNFD()+1;
                surfacefacs[i] = isIn ? 1. : -1.;
                // -1 surf nr is so that curving does not do anything
                FaceDescriptor new_fd(-1, isIn ? new_mat_nrs[i] : fd.DomainIn(),
                                      isIn ? fd.DomainOut() : new_mat_nrs[i], -1);
                new_fd.SetBCProperty(new_si);
                mesh.AddFaceDescriptor(new_fd);
                si_map[i] = new_si;
                moved_surfaces.SetBit(i);
                mesh.SetBCName(new_si-1, "mapped_" + name);
              }
          }
      }

    for(auto si : params.surfid)
      if(surfacefacs[si] == 0.0)
        throw Exception("Surface " + to_string(si) + " is not a boundary of the domain to be grown into!");
  }

  void BoundaryLayerTool ::CreateFaceDescriptorsSides()
  {
    BitArray face_done(mesh.GetNFD()+1);
    face_done.Clear();
    for(const auto& sel : mesh.SurfaceElements())
      {
        auto facei = sel.GetIndex();
        if(face_done.Test(facei))
          continue;
        bool point_moved = false;
        bool point_fixed = false;
        for(auto pi : sel.PNums())
          {
            if(growthvectors[pi].Length() > 0)
              point_moved = true;
            else
              point_fixed = true;
          }
        if(point_moved && !moved_surfaces.Test(facei))
          {
            int new_si = mesh.GetNFD()+1;
            const auto& fd = mesh.GetFaceDescriptor(facei);
            auto isIn = domains.Test(fd.DomainIn());
            auto isOut = domains.Test(fd.DomainOut());
            int si = params.sides_keep_surfaceindex ? facei : -1;
            // domin and domout can only be set later
            FaceDescriptor new_fd(si, -1,
                                  -1, si);
            new_fd.SetBCProperty(new_si);
            mesh.AddFaceDescriptor(new_fd);
            si_map[facei] = new_si;
            mesh.SetBCName(new_si-1, fd.GetBCName());
            face_done.SetBit(facei);
          }
      }
  }

  void BoundaryLayerTool :: CalculateGrowthVectors()
  {
    growthvectors.SetSize(np);
    growthvectors = 0.;

    for(auto pi : mesh.Points().Range())
    {
        const auto & p = mesh[pi];
        if(p.Type() == INNERPOINT)
            continue;

        std::map<int, Vec<3>> normals;

        // calculate one normal vector per face (average with angles as weights for multiple surface elements within a face)
        for(auto sei : p2sel[pi])
        {
            const auto & sel = mesh[sei];
            auto facei = sel.GetIndex();
            if(!params.surfid.Contains(facei))
                continue;

            auto n = surfacefacs[sel.GetIndex()] * getNormal(sel);

            int itrig = sel.PNums().Pos(pi);
            itrig += sel.GetNP();
            auto v0 = (mesh[sel.PNumMod(itrig+1)] - mesh[pi]).Normalize();
            auto v1 = (mesh[sel.PNumMod(itrig-1)] - mesh[pi]).Normalize();
            if(normals.count(facei)==0)
                normals[facei] = {0.,0.,0.};
            normals[facei] += acos(v0*v1)*n;
        }

        for(auto & [facei, n] : normals)
            n *= 1.0/n.Length();

        // combine normal vectors for each face to keep uniform distances
        auto & np = growthvectors[pi];
        ArrayMem<Vec<3>, 3> ns;
        for (auto &[facei, n] : normals) {
          ns.Append(n);
        }

        ArrayMem<Vec<3>, 3> removed;
        // reduce to full rank of max 3
        while(true)
          {
            if(ns.Size() <= 1)
              break;
            if(ns.Size() == 2 && ns[0] * ns[1] < 1 - 1e-6)
              break;
            if (ns.Size() == 3)
              {
                DenseMatrix mat(3,3);
                for(auto i : Range(3))
                  for(auto j : Range(3))
                    mat(i,j) = ns[i][j];
                if(fabs(mat.Det()) > 1e-6)
                  break;
              }
            int maxpos1;
            int maxpos2;
            double val = 0;
            for (auto i : Range(ns))
              {
                for (auto j : Range(i + 1, ns.Size()))
                  {
                    double ip = ns[i] * ns[j];
                    if(ip > val)
                      {
                        val = ip;
                        maxpos1 = i;
                        maxpos2 = j;
                      }
                  }
              }
            removed.Append(ns[maxpos1]);
            removed.Append(ns[maxpos2]);
            ns[maxpos1] = 0.5 * (ns[maxpos1] + ns[maxpos2]);
            ns.DeleteElement(maxpos2);
          }

        if(ns.Size() == 0)
          continue;
        if(ns.Size() == 1)
          np = ns[0];
        else if(ns.Size() == 2)
          {
            np = ns[0];
            auto n = ns[1];
            auto npn = np * n;
            auto npnp = np * np;
            auto nn = n * n;
            if(nn-npn*npn/npnp == 0) { np = n; continue; }
            np += (nn - npn)/(nn - npn*npn/npnp) * (n - npn/npnp * np);
          }
        else // ns.Size() == 3
          {
            DenseMatrix mat(3,3);
            for(auto i : Range(3))
              for(auto j : Range(3))
                mat(i, j) = ns[i] * ns[j];
            Vector rhs(3);
            rhs = 1.;
            Vector res(3);
            DenseMatrix inv(3, ns.Size());
            CalcInverse(mat, inv);
            inv.Mult(rhs, res);
            for(auto i : Range(ns))
              np += res[i] * ns[i];
          }
        for(auto& n : removed)
          if(n * np < 0)
            cout << "WARNING: Growth vector at point " << pi << " in opposite direction to face normal!" << endl << "Growthvector = " << np << ", face normal = " << n << endl;
    }
  }

  Array<Array<pair<SegmentIndex, int>>, SegmentIndex> BoundaryLayerTool :: BuildSegMap()
  {
    // Bit array to keep track of segments already processed
    BitArray segs_done(nseg+1);
    segs_done.Clear();

    // map for all segments with same points
    // points to pair of SegmentIndex, int
    // int is type of other segment, either:
    // 0 == adjacent surface grows layer
    // 1 == adjacent surface doesn't grow layer, but layer ends on it
    // 2 == adjacent surface is interior surface that ends on layer
    // 3 == adjacent surface is exterior surface that ends on layer (not allowed yet)
    Array<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap(segments.Size());

    // moved segments
    is_edge_moved.SetSize(max_edge_nr+1);
    is_edge_moved = false;

    // boundaries to project endings to
    is_boundary_projected.SetSize(nfd_old+1);
    is_boundary_projected.Clear();
    is_boundary_moved.SetSize(nfd_old+1);
    is_boundary_moved.Clear();

    for(auto si : Range(segments))
      {
        if(segs_done[si]) continue;
        const auto& segi = segments[si];
        if(!moved_surfaces.Test(segi.si)) continue;
        segs_done.SetBit(si);
        segmap[si].Append(make_pair(si, 0));
        moved_segs.Append(si);
        is_edge_moved.SetBit(segi.edgenr);
        for(auto sj : Range(segments))
          {
            if(segs_done.Test(sj)) continue;
            const auto& segj = segments[sj];
            if((segi[0] == segj[0] && segi[1] == segj[1]) ||
               (segi[0] == segj[1] && segi[1] == segj[0]))
              {
                segs_done.SetBit(sj);
                int type;
                if(moved_surfaces.Test(segj.si))
                  type = 0;
                else if(const auto& fd = mesh.GetFaceDescriptor(segj.si); domains.Test(fd.DomainIn()) && domains.Test(fd.DomainOut()))
                  {
                    type = 2;
                    if(fd.DomainIn() == 0 || fd.DomainOut() == 0)
                        is_boundary_projected.SetBit(segj.si);
                  }
                else if(const auto& fd = mesh.GetFaceDescriptor(segj.si); !domains.Test(fd.DomainIn()) && !domains.Test(fd.DomainOut()))
                  {
                    type = 3;
                    is_boundary_moved.SetBit(segj.si);
                  }
                else
                  {
                    type = 1;
                    // in case 1 we project the growthvector onto the surface
                    is_boundary_projected.SetBit(segj.si);
                  }
                segmap[si].Append(make_pair(sj, type));
              }
          }
      }

    return segmap;
  }

  BitArray BoundaryLayerTool :: ProjectGrowthVectorsOnSurface()
  {
    BitArray in_surface_direction(nfd_old+1);
    in_surface_direction.Clear();
    // project growthvector on surface for inner angles
    if(params.grow_edges)
      {
        for(const auto& sel : mesh.SurfaceElements())
          if(is_boundary_projected.Test(sel.GetIndex()))
            {
              auto n = getNormal(sel);
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
                  auto tol = v1.Length() * 1e-12;
                  if((v1 * v3 > -tol) && (v2 * v3 > -tol))
                    in_surface_direction.SetBit(sel.GetIndex());
                  else
                    continue;

                  if(!params.project_boundaries.Contains(sel.GetIndex()))
                    continue;
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
        for(const auto& seg : segments)
          {
            int count = 0;
            for(const auto& seg2 : segments)
              if(((seg[0] == seg2[0] && seg[1] == seg2[1]) || (seg[0] == seg2[1] && seg[1] == seg2[0])) && params.surfid.Contains(seg2.si))
                count++;
            if(count == 1)
              {
                growthvectors[seg[0]] = {0., 0., 0.};
                growthvectors[seg[1]] = {0., 0., 0.};
              }
          }
      }

        return in_surface_direction;
    }

  void BoundaryLayerTool :: InterpolateGrowthVectors()
  {
    // interpolate tangential component of growth vector along edge
    for(auto edgenr : Range(max_edge_nr))
      {
        // if(!is_edge_moved[edgenr+1]) continue;

        // build sorted list of edge
        Array<PointIndex> points;
        // find first vertex on edge
        double edge_len = 0.;
        auto is_end_point = [&] (PointIndex pi)
        {
          // if(mesh[pi].Type() == FIXEDPOINT)
          //   return true;
          // return false;
            auto segs = topo.GetVertexSegments(pi);
            auto first_edgenr = mesh[segs[0]].edgenr;
            for(auto segi : segs)
                if(mesh[segi].edgenr != first_edgenr)
                    return true;
            return false;
        };

        bool any_grows = false;
        
        for(const auto& seg : segments)
          {
            if(seg.edgenr-1 == edgenr)
              {
                if(growthvectors[seg[0]].Length2() != 0 ||
                   growthvectors[seg[1]].Length2() != 0)
                  any_grows = true;
                if(points.Size() == 0 && is_end_point(seg[0]))
                  {
                    points.Append(seg[0]);
                    points.Append(seg[1]);
                    edge_len += (mesh[seg[1]] - mesh[seg[0]]).Length();
                  }
              }
          }

        if(!any_grows)
          continue;

        if(!points.Size())
          throw Exception("Could not find startpoint for edge " + ToString(edgenr));

        while(true)
          {
            bool point_found = false;
            for(auto si : topo.GetVertexSegments(points.Last()))
              {
                const auto& seg = mesh[si];
                if(seg.edgenr-1 != edgenr)
                    continue;
                if(seg[0] == points.Last() && points[points.Size()-2] !=seg[1])
                  {
                    edge_len += (mesh[points.Last()] - mesh[seg[1]]).Length();
                    points.Append(seg[1]);
                    point_found = true;
                    break;
                  }
                else if(seg[1] == points.Last() &&
                        points[points.Size()-2] != seg[0])
                  {
                    edge_len += (mesh[points.Last()] - mesh[seg[0]]).Length();
                    points.Append(seg[0]);
                    point_found = true;
                    break;
                  }
              }
            if(is_end_point(points.Last()))
              break;
            if(!point_found)
              {
                throw Exception(string("Could not find connected list of line segments for edge ") + edgenr);
              }
          }

        if(growthvectors[points[0]].Length2() == 0 &&
           growthvectors[points.Last()].Length2() == 0)
          continue;

        // tangential part of growth vectors
        auto t1 = (mesh[points[1]]-mesh[points[0]]).Normalize();
        auto gt1 = growthvectors[points[0]] * t1 * t1;
        auto t2 = (mesh[points.Last()]-mesh[points[points.Size()-2]]).Normalize();
        auto gt2 = growthvectors[points.Last()] * t2 * t2;

        if(!is_edge_moved[edgenr+1])
          {
            if(growthvectors[points[0]] * (mesh[points[1]] - mesh[points[0]]) < 0)
              gt1 = 0.;
            if(growthvectors[points.Last()] * (mesh[points[points.Size()-2]] - mesh[points.Last()]) < 0)
              gt2 = 0.;
          }

        double len = 0.;
        for(size_t i = 1; i < points.Size()-1; i++)
          {
            auto pi = points[i];
            len += (mesh[pi] - mesh[points[i-1]]).Length();
            auto t = getEdgeTangent(pi, edgenr);
            auto lam = len/edge_len;
            auto interpol = (1-lam) * (gt1 * t) * t + lam * (gt2 * t) * t;
            growthvectors[pi] += interpol;
          }
      }

    InterpolateSurfaceGrowthVectors();
  }

  void BoundaryLayerTool :: InsertNewElements( FlatArray<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap, const BitArray & in_surface_direction )
  {
    static Timer timer("BoundaryLayerTool::InsertNewElements"); RegionTimer rt(timer);
    Array<Array<PointIndex>, PointIndex> mapto(np);
    // insert new points
    for (PointIndex pi = 1; pi <= np; pi++)
      if (growthvectors[pi].Length2() != 0)
        {
          Point<3> p = mesh[pi];
          for(auto i : Range(params.heights))
            {
              p += params.heights[i] * growthvectors[pi];
              mapto[pi].Append(mesh.AddPoint(p));
            }
        }

    // add 2d quads on required surfaces
    map<pair<PointIndex, PointIndex>, int> seg2edge;
    if(params.grow_edges)
      {
        for(auto sei : moved_segs)
          {
            // copy here since we will add segments and this would
            // invalidate a reference!
            auto segi = segments[sei];
            for(auto [sej, type] : segmap[sei])
              {
                auto segj = segments[sej];
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
                    new_segments.Append(s);
                  }
                // here we need to grow the quad elements
                else if(type == 1)
                  {
                    PointIndex pp1 = segj[1];
                    PointIndex pp2 = segj[0];
                    if(in_surface_direction.Test(segj.si))
                      {
                        Swap(pp1, pp2);
                        is_boundary_moved.SetBit(segj.si);
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
                    new_segments.Append(s0);

                    for(auto i : Range(params.heights))
                      {
                        Element2d sel(QUAD);
                        p3 = mapto[pp2][i];
                        p4 = mapto[pp1][i];
                        sel[0] = p1;
                        sel[1] = p2;
                        sel[2] = p3;
                        sel[3] = p4;
                        for(auto i : Range(4))
                        {
                            sel.GeomInfo()[i].u = 0.0;
                            sel.GeomInfo()[i].v = 0.0;
                        }
                        sel.SetIndex(si_map[segj.si]);
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
                        new_segments.Append(s1);
                        Segment s2;
                        s2[0] = p4;
                        s2[1] = p1;
                        s2[2] = PointIndex::INVALID;
                        pair = make_pair(p1, p4);
                        if(seg2edge.find(pair) == seg2edge.end())
                          seg2edge[pair] = ++max_edge_nr;
                        s2.edgenr = seg2edge[pair];
                        s2.si = segj.si;
                        new_segments.Append(s2);
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
                    new_segments.Append(s3);
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
        if(moved_surfaces.Test(sel.GetIndex()))
          {
            Array<PointIndex> points(sel.PNums());
            if(surfacefacs[sel.GetIndex()] > 0) Swap(points[0], points[2]);
            for(auto j : Range(params.heights))
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
                el.SetIndex(new_mat_nrs[sel.GetIndex()]);
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
                      if(is_boundary_moved.Test(sel.GetIndex()))
                        moveboundarypoint.SetBit(p);
                    }
                }
          }
        if(is_boundary_moved.Test(sel.GetIndex()))
          {
            for(auto& p : mesh[si].PNums())
              if(mapto[p].Size())
                p = mapto[p].Last();
          }
      }

    for(SegmentIndex sei = 0; sei < nseg; sei++)
      {
        auto& seg = segments[sei];
        if(is_boundary_moved.Test(seg.si))
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
            if(el.GetType() == TET)
              {
                if(moved.Size() == 3) // inner corner
                  {
                    PointIndex p1 = moved[0];
                    PointIndex p2 = moved[1];
                    PointIndex p3 = moved[2];
                    auto v1 = mesh[p1];
                    auto n = Cross(mesh[p2]-v1, mesh[p3]-v1);
                    auto d = mesh[mapto[p1][0]] - v1;
                    if(n*d > 0)
                      Swap(p2,p3);
                    PointIndex p4 = p1;
                    PointIndex p5 = p2;
                    PointIndex p6 = p3;
                    for(auto i : Range(params.heights))
                      {
                        Element nel(PRISM);
                        nel[0] = p4; nel[1] = p5; nel[2] = p6;
                        p4 = mapto[p1][i]; p5 = mapto[p2][i]; p6 = mapto[p3][i];
                        nel[3] = p4; nel[4] = p5; nel[5] = p6;
                        nel.SetIndex(el.GetIndex());
                        mesh.AddVolumeElement(nel);
                      }
                  }
                if(moved.Size() == 2)
                  {
                    if(fixed.Size() == 1)
                      {
                        PointIndex p1 = moved[0];
                        PointIndex p2 = moved[1];
                        for(auto i : Range(params.heights))
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
                  }
                if(moved.Size() == 1 && fixed.Size() == 1)
                  {
                    PointIndex p1 = moved[0];
                    for(auto i : Range(params.heights))
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
            else if(el.GetType() == PYRAMID)
              {
                if(moved.Size() == 2)
                  {
                    if(fixed.Size() != 2)
                      throw Exception("This case is not implemented yet! Fixed size = " + ToString(fixed.Size()));
                    PointIndex p1 = moved[0];
                    PointIndex p2 = moved[1];
                    for(auto i : Range(params.heights))
                      {
                        PointIndex p3 = mapto[moved[1]][i];
                        PointIndex p4 = mapto[moved[0]][i];
                        Element nel(PYRAMID);
                        nel[0] = p1;
                        nel[1] = p2;
                        nel[2] = p3;
                        nel[3] = p4;
                        nel[4] = el[0] + el[1] + el[2] + el[3] + el[4] - fixed[0] - fixed[1] - moved[0] - moved[1];
                        if(Cross(mesh[p2] - mesh[p1], mesh[p4]-mesh[p1]) * (mesh[nel[4]]-mesh[nel[1]]) > 0)
                          Swap(nel[1], nel[3]);
                        nel.SetIndex(el.GetIndex());
                        mesh.AddVolumeElement(nel);
                        p1 = p4;
                        p2 = p3;
                      }
                  }
                else if(moved.Size() == 1)
                  throw Exception("This case is not implemented yet!");
              }
            else
              throw Exception("Boundarylayer only implemented for tets and pyramids outside yet!");
          }
      }
  }

  void BoundaryLayerTool :: SetDomInOut()
  {
    for(auto i : Range(1, nfd_old+1))
      if(moved_surfaces.Test(i))
        {
          if(auto dom = mesh.GetFaceDescriptor(si_map[i]).DomainIn(); dom > ndom_old)
            mesh.GetFaceDescriptor(i).SetDomainOut(dom);
          else
            mesh.GetFaceDescriptor(i).SetDomainIn(mesh.GetFaceDescriptor(si_map[i]).DomainOut());
        }
  }

  void BoundaryLayerTool :: SetDomInOutSides()
  {
    BitArray done(mesh.GetNFD()+1);
    done.Clear();
    for(auto sei : Range(mesh.SurfaceElements()))
      {
        auto& sel = mesh[sei];
        auto index = sel.GetIndex();
        if(done.Test(index))
          continue;
        done.SetBit(index);
        auto& fd = mesh.GetFaceDescriptor(index);
        if(fd.DomainIn() != -1)
          continue;
        int e1, e2;
        mesh.GetTopology().GetSurface2VolumeElement(sei+1, e1, e2);
        if(e1 == 0)
          fd.SetDomainIn(0);
        else
          fd.SetDomainIn(mesh.VolumeElement(e1).GetIndex());
        if(e2 == 0)
          fd.SetDomainOut(0);
        else
          fd.SetDomainOut(mesh.VolumeElement(e2).GetIndex());
      }
  }

  void BoundaryLayerTool :: AddSegments()
  {
    if(have_single_segments)
        MergeAndAddSegments(mesh, new_segments);
    else
    {
      for(auto & seg : new_segments)
        mesh.AddSegment(seg);
    }
  }

  void BoundaryLayerTool :: FixVolumeElements()
  {
      static Timer timer("BoundaryLayerTool::FixVolumeElements"); RegionTimer rt(timer);
      BitArray is_inner_point(mesh.GetNP()+1);
      is_inner_point.Clear();

      auto changed_domains = domains;
      if(!params.outside)
          changed_domains.Invert();

      for(ElementIndex ei : Range(ne))
          if(changed_domains.Test(mesh[ei].GetIndex()))
              for(auto pi : mesh[ei].PNums())
                  if(mesh[pi].Type() == INNERPOINT)
                      is_inner_point.SetBit(pi);

      Array<PointIndex> points;
      for(auto pi : mesh.Points().Range())
          if(is_inner_point.Test(pi))
              points.Append(pi);

      auto p2el = mesh.CreatePoint2ElementTable(is_inner_point);

      // smooth growth vectors to shift additional element layers to the inside and fix flipped tets
      for(auto step : Range(10))
      {
          for(auto pi : points)
          {
              Vec<3> average_gw = 0.0;
              auto & els = p2el[pi];
              size_t cnt = 0;
              for(auto ei : els)
                  if(ei<ne)
                      for(auto pi1 : mesh[ei].PNums())
                          if(pi1<=np)
                          {
                              average_gw += growthvectors[pi1];
                              cnt++;
                          }
              growthvectors[pi] = 1.0/cnt * average_gw;
          }
      }
      
      for(auto pi : points)
      {
          mesh[pi] += height * growthvectors[pi];
          growthvectors[pi] = 0.0;
      }
  }

  void BoundaryLayerTool :: Perform()
  {
      CreateNewFaceDescriptors();
      CalculateGrowthVectors();
      CreateFaceDescriptorsSides();
      auto segmap = BuildSegMap();

      auto in_surface_direction = ProjectGrowthVectorsOnSurface();

      if(params.limit_growth_vectors)
        LimitGrowthVectorLengths();

      InterpolateGrowthVectors();
      FixVolumeElements();
      InsertNewElements(segmap, in_surface_direction);
      SetDomInOut();
      AddSegments();
      mesh.GetTopology().ClearEdges();
      mesh.SetNextMajorTimeStamp();
      mesh.UpdateTopology();
      SetDomInOutSides();
      MeshingParameters mp;
      mp.optimize3d ="m";
      mp.optsteps3d = 4;
      OptimizeVolume(mp, mesh);
  }

  void GenerateBoundaryLayer(Mesh& mesh, const BoundaryLayerParameters& blp)
  {
    static Timer timer("Create Boundarylayers");
    RegionTimer regt(timer);

    BoundaryLayerTool tool(mesh, blp);
    tool.Perform();
  }

} // namespace netgen

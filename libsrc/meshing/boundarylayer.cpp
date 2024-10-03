#include "boundarylayer.hpp"
#include "boundarylayer_limiter.hpp"

#include <regex>
#include <set>

#include "debugging.hpp"
#include "global.hpp"
#include "meshfunc.hpp"

namespace netgen {

struct SpecialPointException : public Exception {
  SpecialPointException() : Exception("") {}
};

std::tuple<int, int> FindCloseVectors(FlatArray<Vec<3>> ns,
                                      bool find_max = true) {
  int maxpos1;
  int maxpos2;

  double val = find_max ? -1e99 : 1e99;
  for (auto i : Range(ns))
    for (auto j : Range(i + 1, ns.Size())) {
      double ip = ns[i] * ns[j];
      if ((find_max && (ip > val)) || (!find_max && (ip < val))) {
        val = ip;
        maxpos1 = i;
        maxpos2 = j;
      }
    }
  return {maxpos1, maxpos2};
}

Vec<3> CalcGrowthVector(FlatArray<Vec<3>> ns) {
  if (ns.Size() == 0)
    return {0, 0, 0};
  if (ns.Size() == 1)
    return ns[0];
  if (ns.Size() == 2) {
    auto gw = ns[0];
    auto n = ns[1];
    auto npn = gw * n;
    auto npnp = gw * gw;
    auto nn = n * n;
    if (fabs(nn - npn * npn / npnp) < 1e-6)
      return n;
    gw += (nn - npn) / (nn - npn * npn / npnp) * (n - npn / npnp * gw);
    return gw;
  }
  if (ns.Size() == 3) {
    DenseMatrix mat(3, 3);
    for (auto i : Range(3))
      for (auto j : Range(3))
        mat(i, j) = ns[i][j];

    if (fabs(mat.Det()) > 1e-6) {
      DenseMatrix mat(3, 3);
      for (auto i : Range(3))
        for (auto j : Range(3))
          mat(i, j) = ns[i] * ns[j];
      Vector rhs(3);
      rhs = 1.;
      Vector res(3);
      DenseMatrix inv(3, ns.Size());
      CalcInverse(mat, inv);
      inv.Mult(rhs, res);
      Vec<3> growth = 0.;
      for (auto i : Range(ns))
        growth += res[i] * ns[i];
      return growth;
    }
  }
  auto [maxpos1, maxpos2] = FindCloseVectors(ns);
  Array<Vec<3>> new_normals;
  new_normals = ns;
  const auto dot = ns[maxpos1] * ns[maxpos2];
  auto average = 0.5 * (ns[maxpos1] + ns[maxpos2]);
  average.Normalize();
  new_normals[maxpos1] = average;
  new_normals.DeleteElement(maxpos2);
  auto gw = CalcGrowthVector(new_normals);

  for (auto n : ns)
    if (n * gw < 0)
      throw SpecialPointException();
  return gw;
}

SpecialBoundaryPoint ::GrowthGroup ::GrowthGroup(FlatArray<int> faces_,
                                                 FlatArray<Vec<3>> normals) {
  faces = faces_;
  growth_vector = CalcGrowthVector(normals);
}

SpecialBoundaryPoint ::SpecialBoundaryPoint(
    const std::map<int, Vec<3>> &normals) {
  // find opposing face normals
  Array<Vec<3>> ns;
  Array<int> faces;
  for (auto [face, normal] : normals) {
    ns.Append(normal);
    faces.Append(face);
  }

  auto [minface1, minface2] = FindCloseVectors(ns, false);
  minface1 = faces[minface1];
  minface2 = faces[minface2];
  Array<int> g1_faces;
  g1_faces.Append(minface1);
  Array<int> g2_faces;
  g2_faces.Append(minface2);
  auto n1 = normals.at(minface1);
  auto n2 = normals.at(minface2);
  separating_direction = 0.5 * (n2 - n1);

  Array<Vec<3>> normals1, normals2;
  for (auto [facei, normali] : normals)
    if (facei != minface1 && facei != minface2) {
      g1_faces.Append(facei);
      g2_faces.Append(facei);
    }
  for (auto fi : g1_faces)
    normals1.Append(normals.at(fi));
  for (auto fi : g2_faces)
    normals2.Append(normals.at(fi));
  growth_groups.Append(GrowthGroup(g1_faces, normals1));
  growth_groups.Append(GrowthGroup(g2_faces, normals2));
}

Vec<3> BoundaryLayerTool ::getEdgeTangent(PointIndex pi, int edgenr,
                                          FlatArray<Segment *> segs) {
  Vec<3> tangent = 0.0;
  ArrayMem<PointIndex, 2> pts;
  for (auto *p_seg : segs) {
    auto &seg = *p_seg;
    if (seg.edgenr != edgenr)
      continue;
    PointIndex other = seg[0] + seg[1] - pi;
    if (!pts.Contains(other))
      pts.Append(other);
  }
  if (pts.Size() != 2) {
    cout << "getEdgeTangent pi = " << pi << ", edgenr = " << edgenr << endl;
    cout << pts << endl;
    for (auto *p_seg : segs)
      cout << *p_seg << endl;
    throw NG_EXCEPTION("Something went wrong in getEdgeTangent!");
  }
  tangent = mesh[pts[1]] - mesh[pts[0]];
  return tangent.Normalize();
}

void BoundaryLayerTool ::LimitGrowthVectorLengths() {
  static Timer tall("BoundaryLayerTool::LimitGrowthVectorLengths");
  RegionTimer rtall(tall);

  GrowthVectorLimiter limiter(*this);
  limiter.Perform();
}

// depending on the geometry type, the mesh contains segments multiple times
// (once for each face)
bool HaveSingleSegments(const Mesh &mesh) {
  auto &topo = mesh.GetTopology();
  NgArray<SurfaceElementIndex> surf_els;

  for (auto segi : Range(mesh.LineSegments())) {
    mesh.GetTopology().GetSegmentSurfaceElements(segi + 1, surf_els);
    if (surf_els.Size() < 2)
      continue;

    auto seg = mesh[segi];
    auto pi0 = min(seg[0], seg[1]);
    auto pi1 = max(seg[0], seg[1]);
    auto p0_segs = topo.GetVertexSegments(seg[0]);

    for (auto segi_other : p0_segs) {
      if (segi_other == segi)
        continue;

      auto seg_other = mesh[segi_other];
      auto pi0_other = min(seg_other[0], seg_other[1]);
      auto pi1_other = max(seg_other[0], seg_other[1]);
      if (pi0_other == pi0 && pi1_other == pi1)
        return false;
    }

    // found segment with multiple adjacent surface elements but no other
    // segments with same points -> have single segments
    return true;
  }

  return true;
}

// duplicates segments (and sets seg.si accordingly) to have a unified data
// structure for all geometry types
Array<Segment> BuildSegments(Mesh &mesh) {
  Array<Segment> segments;
  // auto& topo = mesh.GetTopology();

  NgArray<SurfaceElementIndex> surf_els;

  for (auto segi : Range(mesh.LineSegments())) {
    auto seg = mesh[segi];
    mesh.GetTopology().GetSegmentSurfaceElements(segi + 1, surf_els);
    for (auto seli : surf_els) {
      const auto &sel = mesh[seli];
      seg.si = sel.GetIndex();

      auto np = sel.GetNP();
      for (auto i : Range(np)) {
        if (sel[i] == seg[0]) {
          if (sel[(i + 1) % np] != seg[1])
            swap(seg[0], seg[1]);
          break;
        }
      }

      segments.Append(seg);
    }
  }
  return segments;
}

void MergeAndAddSegments(Mesh &mesh, FlatArray<Segment> segments,
                         FlatArray<Segment> new_segments) {
  INDEX_2_HASHTABLE<bool> already_added(segments.Size() +
                                        2 * new_segments.Size());

  mesh.LineSegments().SetSize0();

  auto addSegment = [&](const auto &seg) {
    INDEX_2 i2(seg[0], seg[1]);
    i2.Sort();
    if (!already_added.Used(i2)) {
      mesh.AddSegment(seg);
      already_added.Set(i2, true);
    }
  };

  for (const auto &seg : segments)
    addSegment(seg);

  for (const auto &seg : new_segments)
    addSegment(seg);
}

BoundaryLayerTool::BoundaryLayerTool(Mesh &mesh_,
                                     const BoundaryLayerParameters &params_)
    : mesh(mesh_), topo(mesh_.GetTopology()), params(params_) {
  static Timer timer("BoundaryLayerTool::ctor");
  RegionTimer regt(timer);
  ProcessParameters();
  if (domains.NumSet() == 0)
    return;

  topo.SetBuildVertex2Element(true);
  mesh.UpdateTopology();

  have_single_segments = HaveSingleSegments(mesh);

  if (have_single_segments)
    segments = BuildSegments(mesh);
  else
    segments = mesh.LineSegments();

  np = mesh.GetNP();
  ne = mesh.GetNE();
  nse = mesh.GetNSE();
  nseg = segments.Size();

  p2sel = mesh.CreatePoint2SurfaceElementTable();

  nfd_old = mesh.GetNFD();
  moved_surfaces.SetSize(nfd_old + 1);
  moved_surfaces.Clear();
  si_map.SetSize(nfd_old + 1);
  for (auto i : Range(nfd_old + 1))
    si_map[i] = i;
}

void BoundaryLayerTool ::CreateNewFaceDescriptors() {
  surfacefacs.SetSize(nfd_old + 1);
  surfacefacs = 0.0;
  // create new FaceDescriptors
  for (auto i : Range(1, nfd_old + 1)) {
    const auto &fd = mesh.GetFaceDescriptor(i);
    string name = fd.GetBCName();
    if (par_surfid.Contains(i)) {
      if (auto isIn = domains.Test(fd.DomainIn());
          isIn != domains.Test(fd.DomainOut())) {
        int new_si = mesh.GetNFD() + 1;
        surfacefacs[i] = isIn ? 1. : -1.;
        moved_surfaces.SetBit(i);
        if (!insert_only_volume_elements) {
          // -1 surf nr is so that curving does not do anything
          FaceDescriptor new_fd(-1, isIn ? new_mat_nrs[i] : fd.DomainIn(),
                                isIn ? fd.DomainOut() : new_mat_nrs[i], -1);
          new_fd.SetBCProperty(new_si);
          new_fd.SetSurfColour(fd.SurfColour());
          mesh.AddFaceDescriptor(new_fd);
          si_map[i] = new_si;
          mesh.SetBCName(new_si - 1, "mapped_" + name);
        }
      }
      // curving of surfaces with boundary layers will often
      // result in pushed through elements, since we do not (yet)
      // curvature through layers.
      // Therefore we disable curving for these surfaces.
      if (!params.keep_surfaceindex)
        mesh.GetFaceDescriptor(i).SetSurfNr(-1);
    }
  }

  for (auto si : par_surfid)
    if (surfacefacs[si] == 0.0)
      throw Exception("Surface " + to_string(si) +
                      " is not a boundary of the domain to be grown into!");
}

void BoundaryLayerTool ::CreateFaceDescriptorsSides() {
  if (insert_only_volume_elements)
    return;
  BitArray face_done(mesh.GetNFD() + 1);
  face_done.Clear();
  for (const auto &sel : mesh.SurfaceElements()) {
    auto facei = sel.GetIndex();
    if (face_done.Test(facei))
      continue;
    bool point_moved = false;
    // bool point_fixed = false;
    for (auto pi : sel.PNums()) {
      if (growthvectors[pi].Length() > 0)
        point_moved = true;
      /*
      else
        point_fixed = true;
      */
    }
    if (point_moved && !moved_surfaces.Test(facei)) {
      int new_si = mesh.GetNFD() + 1;
      const auto &fd = mesh.GetFaceDescriptor(facei);
      // auto isIn = domains.Test(fd.DomainIn());
      // auto isOut = domains.Test(fd.DomainOut());
      int si = params.sides_keep_surfaceindex ? facei : -1;
      // domin and domout can only be set later
      FaceDescriptor new_fd(si, -1, -1, si);
      new_fd.SetBCProperty(new_si);
      mesh.AddFaceDescriptor(new_fd);
      si_map[facei] = new_si;
      mesh.SetBCName(new_si - 1, fd.GetBCName());
      face_done.SetBit(facei);
    }
  }
}

void BoundaryLayerTool ::CalculateGrowthVectors() {
  growthvectors.SetSize(np);
  growthvectors = 0.;

  for (auto pi : mesh.Points().Range()) {
    const auto &p = mesh[pi];
    if (p.Type() == INNERPOINT)
      continue;

    std::map<int, Vec<3>> normals;

    // calculate one normal vector per face (average with angles as weights for
    // multiple surface elements within a face)
    for (auto sei : p2sel[pi]) {
      const auto &sel = mesh[sei];
      auto facei = sel.GetIndex();
      if (!par_surfid.Contains(facei))
        continue;

      auto n = surfacefacs[sel.GetIndex()] * getNormal(sel);

      int itrig = sel.PNums().Pos(pi);
      itrig += sel.GetNP();
      auto v0 = (mesh[sel.PNumMod(itrig + 1)] - mesh[pi]).Normalize();
      auto v1 = (mesh[sel.PNumMod(itrig - 1)] - mesh[pi]).Normalize();
      if (normals.count(facei) == 0)
        normals[facei] = {0., 0., 0.};
      normals[facei] += acos(v0 * v1) * n;
    }

    for (auto &[facei, n] : normals)
      n *= 1.0 / n.Length();

    // combine normal vectors for each face to keep uniform distances
    ArrayMem<Vec<3>, 5> ns;
    for (auto &[facei, n] : normals) {
      ns.Append(n);
    }

    try {
      growthvectors[pi] = CalcGrowthVector(ns);
    } catch (const SpecialPointException &e) {
      special_boundary_points.emplace(pi, normals);
      growthvectors[pi] =
          special_boundary_points[pi].growth_groups[0].growth_vector;
    }
  }
}

Array<Array<pair<SegmentIndex, int>>, SegmentIndex>
BoundaryLayerTool ::BuildSegMap() {
  // Bit array to keep track of segments already processed
  BitArray segs_done(nseg + 1);
  segs_done.Clear();

  // map for all segments with same points
  // points to pair of SegmentIndex, int
  // int is type of other segment, either:
  // 0 == adjacent surface grows layer
  // 1 == adjacent surface doesn't grow layer, but layer ends on it
  // 2 == adjacent surface is interior surface that ends on layer
  // 3 == adjacent surface is exterior surface that ends on layer (not allowed
  // yet)
  Array<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap(segments.Size());

  // moved segments
  is_edge_moved.SetSize(max_edge_nr + 1);
  is_edge_moved = false;

  // boundaries to project endings to
  is_boundary_projected.SetSize(nfd_old + 1);
  is_boundary_projected.Clear();
  is_boundary_moved.SetSize(nfd_old + 1);
  is_boundary_moved.Clear();

  for (auto si : Range(segments)) {
    if (segs_done[si])
      continue;
    const auto &segi = segments[si];
    if (!moved_surfaces.Test(segi.si))
      continue;
    segs_done.SetBit(si);
    segmap[si].Append(make_pair(si, 0));
    moved_segs.Append(si);
    is_edge_moved.SetBit(segi.edgenr);
    for (auto sj : Range(segments)) {
      if (segs_done.Test(sj))
        continue;
      const auto &segj = segments[sj];
      if ((segi[0] == segj[0] && segi[1] == segj[1]) ||
          (segi[0] == segj[1] && segi[1] == segj[0])) {
        segs_done.SetBit(sj);
        int type;
        if (moved_surfaces.Test(segj.si)) {
          type = 0;
          moved_segs.Append(sj);
        } else if (const auto &fd = mesh.GetFaceDescriptor(segj.si);
                   domains.Test(fd.DomainIn()) &&
                   domains.Test(fd.DomainOut())) {
          type = 2;
          if (fd.DomainIn() == 0 || fd.DomainOut() == 0)
            is_boundary_projected.SetBit(segj.si);
        } else if (const auto &fd = mesh.GetFaceDescriptor(segj.si);
                   !domains.Test(fd.DomainIn()) &&
                   !domains.Test(fd.DomainOut())) {
          type = 3;
          is_boundary_moved.SetBit(segj.si);
        } else {
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

BitArray BoundaryLayerTool ::ProjectGrowthVectorsOnSurface() {
  BitArray in_surface_direction(nfd_old + 1);
  in_surface_direction.Clear();
  // project growthvector on surface for inner angles
  if (params.grow_edges) {
    for (const auto &sel : mesh.SurfaceElements())
      if (is_boundary_projected.Test(sel.GetIndex())) {
        auto n = getNormal(sel);
        for (auto i : Range(sel.PNums())) {
          auto pi = sel.PNums()[i];
          if (growthvectors[pi].Length2() == 0.)
            continue;
          auto next = sel.PNums()[(i + 1) % sel.GetNV()];
          auto prev = sel.PNums()[i == 0 ? sel.GetNV() - 1 : i - 1];
          auto v1 = (mesh[next] - mesh[pi]).Normalize();
          auto v2 = (mesh[prev] - mesh[pi]).Normalize();
          auto v3 = growthvectors[pi];
          v3.Normalize();
          auto tol = v1.Length() * 1e-12;
          if ((v1 * v3 > -tol) && (v2 * v3 > -tol))
            in_surface_direction.SetBit(sel.GetIndex());
          else
            continue;

          if (!par_project_boundaries.Contains(sel.GetIndex()))
            continue;
          auto &g = growthvectors[pi];
          auto ng = n * g;
          auto gg = g * g;
          auto nn = n * n;
          // if(fabs(ng*ng-nn*gg) < 1e-12 || fabs(ng) < 1e-12) continue;
          auto a = -ng * ng / (ng * ng - nn * gg);
          auto b = ng * gg / (ng * ng - nn * gg);
          g += a * g + b * n;
        }
      }
  } else {
    for (const auto &seg : segments) {
      int count = 0;
      for (const auto &seg2 : segments)
        if (((seg[0] == seg2[0] && seg[1] == seg2[1]) ||
             (seg[0] == seg2[1] && seg[1] == seg2[0])) &&
            par_surfid.Contains(seg2.si))
          count++;
      if (count == 1) {
        growthvectors[seg[0]] = {0., 0., 0.};
        growthvectors[seg[1]] = {0., 0., 0.};
      }
    }
  }

  return in_surface_direction;
}

void BoundaryLayerTool ::InsertNewElements(
    FlatArray<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap,
    const BitArray &in_surface_direction) {
  static Timer timer("BoundaryLayerTool::InsertNewElements");
  RegionTimer rt(timer);
  mapto.SetSize(0);
  mapto.SetSize(np);
  mapfrom.SetSize(mesh.GetNP());
  mapfrom = PointIndex::INVALID;

  auto changed_domains = domains;
  if (!params.outside)
    changed_domains.Invert();

  auto &identifications = mesh.GetIdentifications();
  const int identnr = identifications.GetNr("boundarylayer");

  auto add_points = [&](PointIndex pi, Vec<3> &growth_vector,
                        Array<PointIndex> &new_points) {
    Point<3> p = mesh[pi];
    PointIndex pi_last = pi;
    double height = 0.0;
    for (auto i : Range(par_heights)) {
      height += par_heights[i];
      auto pi_new = mesh.AddPoint(p);
      mapfrom.Append(pi);
      new_points.Append(pi_new);
      growth_vector_map[pi_new] = {&growth_vector, height};
      if (special_boundary_points.count(pi) > 0)
        mesh.AddLockedPoint(pi_new);
      pi_last = pi_new;
    }
  };

  // insert new points
  for (PointIndex pi = 1; pi <= np; pi++) {
    if (growthvectors[pi].Length2() != 0) {
      if (special_boundary_points.count(pi)) {
        for (auto &group : special_boundary_points[pi].growth_groups)
          add_points(pi, group.growth_vector, group.new_points);
      } else
        add_points(pi, growthvectors[pi], mapto[pi]);
    }
  }

  // get point from mapto (or the group if point is mapped to multiple new
  // points) layer = -1 means last point (top of boundary layer)
  auto newPoint = [&](PointIndex pi, int layer = -1, int group = 0) {
    if (layer == -1)
      layer = par_heights.Size() - 1;
    if (special_boundary_points.count(pi))
      return special_boundary_points[pi].growth_groups[group].new_points[layer];
    else
      return mapto[pi][layer];
  };

  auto hasMoved = [&](PointIndex pi) {
    return mapto[pi].Size() > 0 || special_boundary_points.count(pi);
  };

  auto numGroups = [&](PointIndex pi) -> size_t {
    if (special_boundary_points.count(pi))
      return special_boundary_points[pi].growth_groups.Size();
    else
      return 1;
  };

  auto getGroups = [&](PointIndex pi, int face_index) -> Array<int> {
    auto n = numGroups(pi);
    Array<int> groups;
    if (n == 1) {
      groups.Append(0);
      return groups;
    }
    const auto &all_groups = special_boundary_points[pi].growth_groups;
    for (auto i : Range(n))
      if (all_groups[i].faces.Contains(face_index))
        groups.Append(i);
    // cout << "groups " << pi << ", " << face_index << endl << groups;
    return groups;
  };

  // add 2d quads on required surfaces
  map<pair<PointIndex, PointIndex>, int> seg2edge;
  map<int, int> edge_map;
  int edge_nr = max_edge_nr;
  auto getEdgeNr = [&](int ei) {
    if (edge_map.count(ei) == 0)
      edge_map[ei] = ++edge_nr;
    return edge_map[ei];
  };
  if (params.grow_edges) {
    for (auto sei : moved_segs) {
      // copy here since we will add segments and this would
      // invalidate a reference!
      // auto segi = segments[sei];
      for (auto [sej, type] : segmap[sei]) {
        auto segj = segments[sej];
        if (type == 0) {
          auto addSegment = [&](PointIndex p0, PointIndex p1,
                                bool extra_edge_nr = false) {
            Segment s;
            s[0] = p0;
            s[1] = p1;
            s[2] = PointIndex::INVALID;
            auto pair =
                s[0] < s[1] ? make_pair(s[0], s[1]) : make_pair(s[1], s[0]);
            if (extra_edge_nr)
              s.edgenr = ++edge_nr;
            else
              s.edgenr = getEdgeNr(segj.edgenr);
            s.si = si_map[segj.si];
            new_segments.Append(s);
            // cout << __LINE__ <<"\t" << s << endl;
            return s;
          };

          auto p0 = segj[0], p1 = segj[1];
          auto g0 = getGroups(p0, segj.si);
          auto g1 = getGroups(p1, segj.si);

          if (g0.Size() == 1 && g1.Size() == 1)
            auto s =
                addSegment(newPoint(p0, -1, g0[0]), newPoint(p1, -1, g1[0]));
          else {
            if (g0.Size() == 2)
              addSegment(newPoint(p0, -1, g0[0]), newPoint(p0, -1, g0[1]));
            if (g1.Size() == 2)
              addSegment(newPoint(p1, -1, g1[0]), newPoint(p1, -1, g1[1]));
          }
        }
        // here we need to grow the quad elements
        else if (type == 1) {
          PointIndex pp1 = segj[1];
          PointIndex pp2 = segj[0];
          if (in_surface_direction.Test(segj.si)) {
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

          for (auto i : Range(par_heights)) {
            Element2d sel(QUAD);
            p3 = newPoint(pp2, i);
            p4 = newPoint(pp1, i);
            sel[0] = p1;
            sel[1] = p2;
            sel[2] = p3;
            sel[3] = p4;
            for (auto i : Range(4)) {
              sel.GeomInfo()[i].u = 0.0;
              sel.GeomInfo()[i].v = 0.0;
            }
            sel.SetIndex(si_map[segj.si]);
            new_sels.Append(sel);
            new_sels_on_moved_bnd.Append(sel);

            // TODO: Too many, would be enough to only add outermost ones
            Segment s1;
            s1[0] = p2;
            s1[1] = p3;
            s1[2] = PointIndex::INVALID;
            auto pair = make_pair(p2, p3);
            s1.edgenr = getEdgeNr(segj.edgenr);
            s1.si = segj.si;
            // new_segments.Append(s1);
            Segment s2;
            s2[0] = p4;
            s2[1] = p1;
            s2[2] = PointIndex::INVALID;
            pair = make_pair(p1, p4);
            s2.edgenr = getEdgeNr(segj.edgenr);
            s2.si = segj.si;
            // new_segments.Append(s2);
            p1 = p4;
            p2 = p3;
          }
          Segment s3;
          s3[0] = p3;
          s3[1] = p4;
          s3[2] = PointIndex::INVALID;
          auto pair = p3 < p4 ? make_pair(p3, p4) : make_pair(p4, p3);
          s3.edgenr = getEdgeNr(segj.edgenr);
          s3.si = segj.si;
          new_segments.Append(s3);
        }
      }
    }
  }

  auto getClosestGroup = [&](PointIndex pi, SurfaceElementIndex sei) {
    auto n = numGroups(pi);
    if (n == 1)
      return 0;
    const auto &sel = mesh[sei];
    auto groups = getGroups(pi, sel.GetIndex());
    if (groups.Size() == 1)
      return groups[0];

    auto &growth_groups = special_boundary_points[pi].growth_groups;

    auto vdir = Center(mesh[sel[0]], mesh[sel[1]], mesh[sel[2]]) - mesh[pi];
    auto dot = vdir * special_boundary_points[pi].separating_direction;

    return dot > 0 ? 1 : 0;
  };

  BitArray fixed_points(np + 1);
  fixed_points.Clear();
  auto p2el = mesh.CreatePoint2ElementTable();
  for (SurfaceElementIndex si = 0; si < nse; si++) {
    // copy because surfaceels array will be resized!
    const auto sel = mesh[si];
    if (moved_surfaces.Test(sel.GetIndex())) {
      Array<PointIndex> points(sel.PNums());
      if (surfacefacs[sel.GetIndex()] > 0)
        Swap(points[0], points[2]);
      ArrayMem<int, 4> groups(points.Size());
      for (auto i : Range(points))
        groups[i] = getClosestGroup(sel[i], si);
      bool add_volume_element = true;
      for (auto pi : sel.PNums())
        if (numGroups(pi) > 1)
          add_volume_element = false;
      for (auto j : Range(par_heights)) {
        auto eltype = points.Size() == 3 ? PRISM : HEX;
        Element el(eltype);
        for (auto i : Range(points))
          el[i] = points[i];
        for (auto i : Range(points))
          points[i] = newPoint(sel.PNums()[i], j, groups[i]);
        if (surfacefacs[sel.GetIndex()] > 0)
          Swap(points[0], points[2]);
        for (auto i : Range(points))
          el[sel.PNums().Size() + i] = points[i];
        auto new_index = new_mat_nrs[sel.GetIndex()];
        if (new_index == -1)
          throw Exception("Boundary " + ToString(sel.GetIndex()) +
                          " with name " + mesh.GetBCName(sel.GetIndex() - 1) +
                          " extruded, but no new material specified for it!");
        el.SetIndex(new_mat_nrs[sel.GetIndex()]);
        if (add_volume_element)
          mesh.AddVolumeElement(el);
        else {
          // Let the volume mesher fill the hole with pyramids/tets
          // To insert pyramids, we need close surface identifications on open
          // quads
          for (auto i : Range(points))
            if (numGroups(sel[i]) == 1)
              identifications.Add(el[i], el[i + points.Size()], identnr);
        }
      }
      Element2d newel = sel;
      for (auto i : Range(points))
        newel[i] = newPoint(sel[i], -1, groups[i]);
      newel.SetIndex(si_map[sel.GetIndex()]);
      new_sels.Append(newel);
    }
    if (is_boundary_moved.Test(sel.GetIndex())) {
      for (auto &p : mesh[si].PNums())
        if (hasMoved(p))
          p = newPoint(p);
    }
  }

  // for (SegmentIndex sei = 0; sei < nseg; sei++) {
  //   auto &seg = segments[sei];
  //   if (is_boundary_moved.Test(seg.si)) {
  //     if(insert_only_volume_elements) throw
  //     Exception("insert_only_volume_elements and is_boundary_moved are
  //     incompatible"); for (auto &p : seg.PNums())
  //       if (hasMoved(p))
  //         p = newPoint(p);
  //   }
  // }

  // fill holes in surface mesh at special boundary points (i.e. points with >=4
  // adjacent boundary faces)
  auto p2sel = ngcore::CreateSortedTable<SurfaceElementIndex, PointIndex>(
      new_sels.Range(),
      [&](auto &table, SurfaceElementIndex ei) {
        for (PointIndex pi : new_sels[ei].PNums())
          table.Add(pi, ei);
      },
      mesh.GetNP());

  for (auto &[special_pi, special_point] : special_boundary_points) {
    if (special_point.growth_groups.Size() != 2)
      throw Exception("special_point.growth_groups.Size() != 2");

    // Special points are split into two new points, when mapping a surface
    // element, we choose the closer one to the center. Now, find points which
    // are mapped to both new points (for different surface elements they belong
    // to). At exactly these points we need to insert new surface elements to
    // fill the hole.
    std::map<int, std::array<std::set<PointIndex>, 2>> close_group;
    for (auto sei : p2sel[special_pi]) {
      const auto &sel = mesh[sei];
      for (auto p : sel.PNums())
        if (p != special_pi)
          close_group[sel.GetIndex()][getClosestGroup(special_pi, sei)].insert(
              p);
    }

    for (auto [fi, groups] : close_group) {
      const auto mapped_fi = si_map[fi];
      std::set<PointIndex> common_points;
      for (auto pi : groups[0])
        if (groups[1].count(pi) == 1)
          common_points.insert(pi);
      if (common_points.size() > 0) {
        auto pi_common = mapto[*common_points.begin()].Last();
        auto new_special_pi0 = special_point.growth_groups[0].new_points.Last();
        auto new_special_pi1 = special_point.growth_groups[1].new_points.Last();
        for (auto sei : p2sel[pi_common]) {
          if (mesh[sei].GetIndex() == mapped_fi &&
              mesh[sei].PNums().Contains(new_special_pi0)) {
            auto sel = mesh[sei];
            sel.Invert();
            for (auto &pi : sel.PNums())
              if (pi != pi_common && pi != new_special_pi0)
                pi = new_special_pi1;
            new_sels.Append(sel);
          }
        }
      }
    }
  }

  for (auto &[pi, special_point] : special_boundary_points) {
    if (special_point.growth_groups.Size() != 2)
      throw Exception("special_point.growth_groups.Size() != 2");
    for (auto igroup : Range(2)) {
      auto &group = special_point.growth_groups[igroup];
      std::set<int> faces;
      for (auto face : group.faces)
        faces.insert(si_map[face]);
      auto pi_new = group.new_points.Last();
      auto pi_new_other =
          special_point.growth_groups[1 - igroup].new_points.Last();
      for (auto sei : p2sel[pi_new])
        faces.erase(mesh[sei].GetIndex());
      for (auto face : faces)
        for (auto seg : new_segments) {
          if ( // seg.si == face
              (seg[0] == pi_new || seg[1] == pi_new) &&
              (seg[0] != pi_new_other && seg[1] != pi_new_other)) {
            bool is_correct_face = false;
            auto pi_other = seg[0] == pi_new ? seg[1] : seg[0];
            for (auto sei : p2sel[pi_other]) {
              if (mesh[sei].GetIndex() == face) {
                is_correct_face = true;
                break;
              }
            }
            if (is_correct_face) {
              Element2d sel;
              sel[0] = seg[1];
              sel[1] = seg[0];
              sel[2] = pi_new_other;
              sel.SetIndex(face);
              new_sels.Append(sel);
            }
          }
        }
    }
  }
}

void BoundaryLayerTool ::SetDomInOut() {
  if (insert_only_volume_elements)
    return;
  for (auto i : Range(1, nfd_old + 1))
    if (moved_surfaces.Test(i)) {
      if (auto dom = mesh.GetFaceDescriptor(si_map[i]).DomainIn();
          dom > ndom_old)
        mesh.GetFaceDescriptor(i).SetDomainOut(dom);
      else
        mesh.GetFaceDescriptor(i).SetDomainIn(
            mesh.GetFaceDescriptor(si_map[i]).DomainOut());
    }
}

void BoundaryLayerTool ::SetDomInOutSides() {
  if (insert_only_volume_elements)
    return;
  BitArray done(mesh.GetNFD() + 1);
  done.Clear();
  for (auto sei : Range(mesh.SurfaceElements())) {
    auto &sel = mesh[sei];
    auto index = sel.GetIndex();
    if (done.Test(index))
      continue;
    done.SetBit(index);
    auto &fd = mesh.GetFaceDescriptor(index);
    if (fd.DomainIn() != -1)
      continue;
    int e1, e2;
    mesh.GetTopology().GetSurface2VolumeElement(sei + 1, e1, e2);
    if (e1 == 0)
      fd.SetDomainIn(0);
    else
      fd.SetDomainIn(mesh.VolumeElement(e1).GetIndex());
    if (e2 == 0)
      fd.SetDomainOut(0);
    else
      fd.SetDomainOut(mesh.VolumeElement(e2).GetIndex());
  }
}

void BoundaryLayerTool ::AddSegments() {
  if (insert_only_volume_elements)
    return;
  if (have_single_segments)
    MergeAndAddSegments(mesh, segments, new_segments);
  else {
    mesh.LineSegments() = segments;
    for (auto &seg : new_segments)
      mesh.AddSegment(seg);
  }
}

void BoundaryLayerTool ::AddSurfaceElements() {
  for (auto &sel :
       insert_only_volume_elements ? new_sels_on_moved_bnd : new_sels) {
    cout << "add surface element " << sel << endl;
    mesh.AddSurfaceElement(sel);
  }
}

void BoundaryLayerTool ::ProcessParameters() {
  if (int *bc = get_if<int>(&params.boundary); bc) {
    for (int i = 1; i <= mesh.GetNFD(); i++)
      if (mesh.GetFaceDescriptor(i).BCProperty() == *bc)
        par_surfid.Append(i);
  } else if (string *s = get_if<string>(&params.boundary); s) {
    regex pattern(*s);
    BitArray boundaries(mesh.GetNFD() + 1);
    boundaries.Clear();
    for (int i = 1; i <= mesh.GetNFD(); i++) {
      auto &fd = mesh.GetFaceDescriptor(i);
      if (regex_match(fd.GetBCName(), pattern)) {
        boundaries.SetBit(i);
        auto dom_pattern = get_if<string>(&params.domain);
        // only add if adjacent to domain
        if (dom_pattern) {
          regex pattern(*dom_pattern);
          bool mat1_match =
              fd.DomainIn() > 0 &&
              regex_match(mesh.GetMaterial(fd.DomainIn()), pattern);
          bool mat2_match =
              fd.DomainOut() > 0 &&
              regex_match(mesh.GetMaterial(fd.DomainOut()), pattern);
          // if boundary is inner or outer remove from list
          if (mat1_match == mat2_match)
            boundaries.Clear(i);
          // if((fd.DomainIn() > 0 &&
          // regex_match(mesh.GetMaterial(fd.DomainIn()), pattern)) ||
          // (fd.DomainOut() > 0 &&
          // regex_match(self.GetMaterial(fd.DomainOut()), pattern)))
          // boundaries.Clear(i);
          // par_surfid.Append(i);
        }
        // else
        //   par_surfid.Append(i);
      }
    }
    for (int i = 1; i <= mesh.GetNFD(); i++)
      if (boundaries.Test(i))
        par_surfid.Append(i);
  } else {
    auto &surfids = *get_if<std::vector<int>>(&params.boundary);
    for (auto id : surfids)
      par_surfid.Append(id);
  }

  insert_only_volume_elements = !params.new_material.has_value();
  if (params.new_material) {
    if (string *mat = get_if<string>(&*params.new_material); mat)
      par_new_mat = {{".*", *mat}};
    else
      par_new_mat = *get_if<map<string, string>>(&*params.new_material);
  }

  if (params.project_boundaries.has_value()) {
    auto proj_bnd = *params.project_boundaries;
    if (string *s = get_if<string>(&proj_bnd); s) {
      regex pattern(*s);
      for (int i = 1; i <= mesh.GetNFD(); i++)
        if (regex_match(mesh.GetFaceDescriptor(i).GetBCName(), pattern))
          par_project_boundaries.Append(i);
    } else {
      for (auto id : *get_if<std::vector<int>>(&proj_bnd))
        par_project_boundaries.Append(id);
    }
  }

  if (double *height = get_if<double>(&params.thickness); height) {
    par_heights.Append(*height);
  } else {
    auto &heights = *get_if<std::vector<double>>(&params.thickness);
    for (auto val : heights)
      par_heights.Append(val);
  }

  int nr_domains = mesh.GetNDomains();
  domains.SetSize(nr_domains + 1); // one based
  domains.Clear();
  if (string *pdomain = get_if<string>(&params.domain); pdomain) {
    regex pattern(*pdomain);
    for (auto i : Range(1, nr_domains + 1))
      if (regex_match(mesh.GetMaterial(i), pattern))
        domains.SetBit(i);
  } else if (int *idomain = get_if<int>(&params.domain); idomain) {
    domains.SetBit(*idomain);
  } else {
    for (auto i : *get_if<std::vector<int>>(&params.domain))
      domains.SetBit(i);
  }
  if (domains.NumSet() == 0)
    return;
  total_height = 0.0;
  for (auto h : par_heights)
    total_height += h;

  max_edge_nr = -1;
  for (const auto &seg : mesh.LineSegments())
    if (seg.edgenr > max_edge_nr)
      max_edge_nr = seg.edgenr;

  int ndom = mesh.GetNDomains();
  ndom_old = ndom;

  new_mat_nrs.SetSize(mesh.FaceDescriptors().Size() + 1);
  new_mat_nrs = -1;
  if (insert_only_volume_elements) {
    for (auto i : Range(1, mesh.GetNFD() + 1)) {
      auto &fd = mesh.GetFaceDescriptor(i);
      auto domin = fd.DomainIn();
      auto domout = fd.DomainOut();
      for (int dom : {domin, domout})
        if (domains.Test(dom)) {
          if (params.outside) {
            dom = domin + domout - dom;
            if (dom == 0)
              throw NG_EXCEPTION("No new material specified for boundarylayer "
                                 "on the outside of domain");
          }
          new_mat_nrs[i] = dom;
        }
    }
  } else {
    for (auto [bcname, matname] : par_new_mat) {
      mesh.SetMaterial(++ndom, matname);
      regex pattern(bcname);
      for (auto i : Range(1, mesh.GetNFD() + 1)) {
        auto &fd = mesh.GetFaceDescriptor(i);
        if (regex_match(fd.GetBCName(), pattern))
          new_mat_nrs[i] = ndom;
      }
    }
  }
  cout << "new mat numbers " << endl << new_mat_nrs << endl;

  if (!params.outside)
    domains.Invert();
}

void BoundaryLayerTool ::Perform() {
  CreateNewFaceDescriptors();
  CalculateGrowthVectors();
  CreateFaceDescriptorsSides();
  auto segmap = BuildSegMap();

  auto in_surface_direction = ProjectGrowthVectorsOnSurface();

  InsertNewElements(segmap, in_surface_direction);

  SetDomInOut();
  AddSegments();
  AddSurfaceElements();

  mesh.CalcSurfacesOfNode();
  topo.SetBuildVertex2Element(true);
  mesh.UpdateTopology();
  InterpolateGrowthVectors();
  InterpolateSurfaceGrowthVectors();

  if (params.limit_growth_vectors)
    LimitGrowthVectorLengths();

  for (auto [pi, data] : growth_vector_map) {
    auto [gw, height] = data;
    mesh[pi] += height * (*gw);
  }

  mesh.CalcSurfacesOfNode();
  mesh.GetTopology().ClearEdges();
  mesh.SetNextMajorTimeStamp();
  mesh.UpdateTopology();
  SetDomInOutSides();
}

void GenerateBoundaryLayer(Mesh &mesh, const BoundaryLayerParameters &blp) {
  static Timer timer("Create Boundarylayers");
  RegionTimer regt(timer);

  BoundaryLayerTool tool(mesh, blp);
  tool.Perform();
}

} // namespace netgen

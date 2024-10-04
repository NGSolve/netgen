#include "boundarylayer.hpp"

namespace netgen {

void BoundaryLayerTool ::InterpolateGrowthVectors() {
  int new_max_edge_nr = max_edge_nr;
  for (const auto &seg : segments)
    if (seg.edgenr > new_max_edge_nr)
      new_max_edge_nr = seg.edgenr;
  for (const auto &seg : new_segments)
    if (seg.edgenr > new_max_edge_nr)
      new_max_edge_nr = seg.edgenr;

  auto getGW = [&](PointIndex pi) -> Vec<3> {
    if (growth_vector_map.count(pi) == 0)
      growth_vector_map[pi] = {&growthvectors[pi], total_height};
    auto [gw, height] = growth_vector_map[pi];
    return height * (*gw);
  };
  auto addGW = [&](PointIndex pi, Vec<3> vec) {
    if (growth_vector_map.count(pi) == 0)
      growth_vector_map[pi] = {&growthvectors[pi], total_height};
    auto [gw, height] = growth_vector_map[pi];
    *gw += 1.0 / height * vec;
  };

  // interpolate tangential component of growth vector along edge
  if (max_edge_nr >= new_max_edge_nr)
    return;

  auto edgenr2seg = ngcore::CreateSortedTable<Segment *, int>(
      Range(segments.Size() + new_segments.Size()),
      [&](auto &table, size_t segi) {
        auto &seg = segi < segments.Size()
                        ? segments[segi]
                        : new_segments[segi - segments.Size()];
        table.Add(seg.edgenr, &seg);
      },
      new_max_edge_nr + 1);
  auto point2seg = ngcore::CreateSortedTable<Segment *, PointIndex>(
      Range(segments.Size() + new_segments.Size()),
      [&](auto &table, size_t segi) {
        auto &seg = segi < segments.Size()
                        ? segments[segi]
                        : new_segments[segi - segments.Size()];
        table.Add(seg[0], &seg);
        table.Add(seg[1], &seg);
      },
      mesh.GetNP());

  for (auto edgenr : Range(max_edge_nr + 1, new_max_edge_nr + 1)) {
    double edge_len = 0.;

    auto is_end_point = [&](PointIndex pi) {
      auto segs = point2seg[pi];
      if (segs.Size() == 1)
        return true;
      auto first_edgenr = (*segs[0]).edgenr;
      for (auto *p_seg : segs)
        if (p_seg->edgenr != first_edgenr)
          return true;
      return false;
    };

    bool any_grows = false;

    Array<PointIndex> points;
    for (auto *p_seg : edgenr2seg[edgenr]) {
      auto &seg = *p_seg;

      if (getGW(seg[0]).Length2() != 0 || getGW(seg[1]).Length2() != 0)
        any_grows = true;

      if (points.Size() == 0)
        for (auto i : Range(2))
          if (is_end_point(seg[i])) {
            points.Append(seg[i]);
            points.Append(seg[1 - i]);
            edge_len += (mesh[seg[1]] - mesh[seg[0]]).Length();
            break;
          }
    }

    if (!any_grows) {
      PrintMessage(1, "BLayer: skip interpolating growth vectors at edge ",
                   edgenr + 1);
      continue;
    }

    if (!points.Size()) {
      cerr << "Could not find startpoint for edge " << edgenr << endl;
      continue;
    }

    std::set<PointIndex> points_set;
    points_set.insert(points[0]);
    points_set.insert(points[1]);

    bool point_found = true;
    while (point_found) {
      if (is_end_point(points.Last()))
        break;
      point_found = false;
      for (auto *p_seg : point2seg[points.Last()]) {
        const auto &seg = *p_seg;
        if (seg.edgenr != edgenr)
          continue;
        auto plast = points.Last();
        if (plast != seg[0] && plast != seg[1])
          continue;
        auto pnew = plast == seg[0] ? seg[1] : seg[0];
        if (pnew == points[0] && points.Size() > 1) {
        }
        if (points_set.count(pnew) > 0 &&
            (pnew != points[0] || points.Size() == 2))
          continue;
        edge_len += (mesh[points.Last()] - mesh[pnew]).Length();
        points.Append(pnew);
        points_set.insert(pnew);
        point_found = true;
        break;
      }
    }
    if (!point_found) {
      cerr << "Could not find connected list of line segments for edge "
           << edgenr << endl;
      cerr << "current points: " << endl << points << endl;
      continue;
    }

    if (getGW(points[0]).Length2() == 0 && getGW(points.Last()).Length2() == 0)
      continue;

    // tangential part of growth vectors
    auto t1 = (mesh[points[1]] - mesh[points[0]]).Normalize();
    auto gt1 = getGW(points[0]) * t1 * t1;
    auto t2 =
        (mesh[points.Last()] - mesh[points[points.Size() - 2]]).Normalize();
    auto gt2 = getGW(points.Last()) * t2 * t2;

    double len = 0.;
    for (auto i : IntRange(1, points.Size() - 1)) {
      auto pi = points[i];
      len += (mesh[pi] - mesh[points[i - 1]]).Length();
      auto t = getEdgeTangent(pi, edgenr, point2seg[pi]);
      auto lam = len / edge_len;
      auto interpol = (1 - lam) * (gt1 * t) * t + lam * (gt2 * t) * t;
      addGW(pi, interpol);
    }
  }
}

void BoundaryLayerTool ::InterpolateSurfaceGrowthVectors() {
  static Timer tall("InterpolateSurfaceGrowthVectors");
  RegionTimer rtall(tall);
  static Timer tsmooth("InterpolateSurfaceGrowthVectors-Smoothing");
  auto np_old = this->np;
  auto np = mesh.GetNP();

  non_bl_growth_vectors.clear();
  // for(const auto & sel : new_sels) {
  //   for(auto pi : sel.PNums()) {
  //     if(mesh[pi].Type() == INNERPOINT)
  //       mesh[pi].SetType(SURFACEPOINT);
  //   }
  // }
  // cout << __FILE__ << ":" << __LINE__ << endl;

  auto getGW = [&](PointIndex pi) -> Vec<3> {
    return growthvectors[pi];
    // if (growth_vector_map.count(pi) == 0) {
    //   non_bl_growth_vectors[pi] = .0;
    //   growth_vector_map[pi] = {&non_bl_growth_vectors[pi], 1.0};
    // }
    // auto [gw, height] = growth_vector_map[pi];
    // return height * (*gw);
  };
  auto addGW = [&](PointIndex pi, Vec<3> vec) {
    growthvectors[pi] += vec;
    // // cout << "add gw " << pi << "\t" << vec << endl;
    // if (growth_vector_map.count(pi) == 0) {
    //   // cout << "\t make new entry" << endl;
    //   non_bl_growth_vectors[pi] = .0;
    //   growth_vector_map[pi] = {&non_bl_growth_vectors[pi], 1.0};
    // }
    // auto [gw, height] = growth_vector_map[pi];
    // // cout << "\tcurrent gw " << *gw << "\t" << height << endl;
    // *gw += 1.0 / height * vec;
  };

  auto hasMoved = [&](PointIndex pi) {
    return (pi - PointIndex::BASE >= np_old) || mapto[pi].Size() > 0 ||
           special_boundary_points.count(pi);
  };

  // cout << __FILE__ << ":" << __LINE__ << endl;
  std::set<PointIndex> points_set;
  for (const auto &sel: mesh.SurfaceElements())
    for (auto pi : sel.PNums())
      if(mesh[pi].Type() == SURFACEPOINT && hasMoved(pi))
        points_set.insert(pi);

  // Array<bool> has_moved_points(max_edge_nr + 1);
  // has_moved_points = false;
  // std::set<PointIndex> moved_edge_points;

  // for (auto seg : segments) {
  //   if (hasMoved(seg[0]) != hasMoved(seg[1]))
  //     has_moved_points[seg.edgenr] = true;
  // }
  // cout << __FILE__ << ":" << __LINE__ << endl;

  // for (auto seg : segments)
  //   if (has_moved_points[seg.edgenr])
  //     for (auto pi : seg.PNums())
  //       if (mesh[pi].Type() == EDGEPOINT)
  //         points_set.insert(pi);

  Array<PointIndex> points;
  for (auto pi : points_set)
    points.Append(pi);
  QuickSort(points);

  // cout << __FILE__ << ":" << __LINE__ << endl;
  // cout << "points to interpolate " << endl << points << endl;

  // cout << __FILE__ << ":" << __LINE__ << endl;
  auto p2sel = mesh.CreatePoint2SurfaceElementTable();
  // cout << __FILE__ << ":" << __LINE__ << endl;
  // auto p2sel = ngcore::CreateSortedTable<SurfaceElementIndex, PointIndex>(
  //       new_sels.Range(),
  //       [&](auto &table, SurfaceElementIndex ei) {
  //         for (PointIndex pi : new_sels[ei].PNums())
  //           table.Add(pi, ei);
  //       },
  //       mesh.GetNP());
  // smooth tangential part of growth vectors from edges to surface elements
  Array<Vec<3>, PointIndex> corrections(mesh.GetNP());
  corrections = 0.0;
  RegionTimer rtsmooth(tsmooth);
  struct Neighbor {
    PointIndex pi;
    SurfaceElementIndex sei;
    double weight;
  };
  Array<ArrayMem<Neighbor, 20>> neighbors(points.Size());
  // cout << __FILE__ << ":" << __LINE__ << endl;

  ArrayMem<double, 20> angles;
  ArrayMem<double, 20> inv_dists;
  for (auto i : points.Range()) {
    auto &p_neighbors = neighbors[i];
    auto pi = points[i];
    angles.SetSize(0);
    inv_dists.SetSize(0);
    for (auto sei : p2sel[pi]) {
      const auto &sel = mesh[sei];
      for (auto pi1 : sel.PNums()) {
        if (pi1 == pi)
          continue;
        auto pi2 = pi1;
        for (auto pi_ : sel.PNums()) {
          if (pi_ != pi && pi_ != pi1) {
            pi2 = pi_;
            break;
          }
        }
        p_neighbors.Append({pi1, sei, 0.0});
        // if((mesh[pi1]-mesh[pi]).Length() < 1e-10) {
        //   cout << "close points " << pi << "\t" << pi1 << endl;
        // }
        inv_dists.Append(1.0 / (mesh[pi1] - mesh[pi]).Length());
        auto dot = (mesh[pi1] - mesh[pi]).Normalize() *
                   (mesh[pi2] - mesh[pi]).Normalize();
        angles.Append(acos(dot));
      }
    }
    double sum_inv_dist = 0.0;
    for (auto inv_dist : inv_dists)
      sum_inv_dist += inv_dist;
    double sum_angle = 0.0;
    for (auto angle : angles)
      sum_angle += angle;

    // cout << "angles " << angles << endl;
    // cout << "inv_dists " << inv_dists << endl;

    double sum_weight = 0.0;
    for (auto i : Range(inv_dists)) {
      p_neighbors[i].weight =
          inv_dists[i] * angles[i] / sum_inv_dist / sum_angle;
      sum_weight += p_neighbors[i].weight;
    }
    for (auto i : Range(inv_dists))
      p_neighbors[i].weight /= sum_weight;

  //   if(pi == 19911) {
  //     cout << "pi " << pi << endl;
  //     for(auto & nb : p_neighbors) {
  //       cout << "neighbor " << nb.pi << "\t" << nb.weight << endl;
  //     }
  //     cout << "inv_dist " << inv_dists << endl;
  //     cout << "angles " << angles << endl;
  //   }
  }
  // cout << __FILE__ << ":" << __LINE__ << endl;

  Array<Vec<3>, SurfaceElementIndex> surf_normals(mesh.GetNSE());
  for (auto sei : mesh.SurfaceElements().Range())
    surf_normals[sei] = getNormal(mesh[sei]);
  // cout << __FILE__ << ":" << __LINE__ << endl;

  BitArray interpolate_tangent(mesh.GetNP()+1);
  // cout << __FILE__ << ":" << __LINE__ << endl;
  interpolate_tangent = false;
  for(auto pi : points) {
    for (auto sei : p2sel[pi])
      if(is_boundary_moved[mesh[sei].GetIndex()])
        interpolate_tangent.SetBit(pi);
  }
  // cout << __FILE__ << ":" << __LINE__ << endl;

  constexpr int N_STEPS = 64;
  for ([[maybe_unused]] auto i : Range(N_STEPS)) {
    for (auto i : points.Range()) {
      auto pi = points[i];
      // cout << "AVERAGE " << pi << endl;
      auto &p_neighbors = neighbors[i];

      ArrayMem<Vec<3>, 20> g_vectors;
      double max_len = 0.0;
      double sum_len = 0.0;

      // average only tangent component on new bl points, average whole growth
      // vector otherwise
      bool do_average_tangent = true;
      for (const auto &s : p_neighbors) {
        auto gw_other = getGW(s.pi) + corrections[s.pi];
        // if(pi == 19911) cout << "neighbor " << s.pi << "\t" << s.weight << '\t' << gw_other << "\tdo avg: " << do_average_tangent << endl;
        // if(pi == 19911) cout << "\tneighbor gw" << gw_other << endl;
        if (do_average_tangent) {
          auto n = surf_normals[s.sei];
          gw_other = gw_other - (gw_other * n) * n;
        }
        // if(pi == 19911) cout << "\tneighbor gw" << gw_other << endl;
        auto v = gw_other;
        auto len = v.Length2();
        sum_len += len;
        max_len = max(max_len, len);
        // if(pi == 19911) cout << "\tneighbor v" << v << endl;
        g_vectors.Append(v);
      }

      if (max_len == 0.0)
        continue;

      double lambda = 0;
      if (i > N_STEPS / 4.)
        lambda = 2.0 * (i - N_STEPS / 4.) / (N_STEPS / 2.);
      lambda = min(1.0, lambda);

      auto &correction = corrections[pi];
      correction = 0.0;
      for (const auto i : p_neighbors.Range()) {
        auto v = g_vectors[i];
        double weight = lambda * p_neighbors[i].weight +
                        (1.0 - lambda) * v.Length2() / sum_len;
        // if(pi == 19911) cout << "pi " << pi << "\tneighbor " << p_neighbors[i].pi << "\tweight " << weight << endl;
        correction += weight * v;
      }

      if (!do_average_tangent)
        correction -= getGW(pi);
      // if(pi == 19911) cout << "pi " << pi << "\tcorrection " << correction << endl;
    }
  }

  for (auto pi : points)
    addGW(pi, corrections[pi]);
}

} // namespace netgen

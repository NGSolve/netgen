#include "boundarylayer.hpp"

namespace netgen
{

namespace detail
{
struct Neighbor
{
  PointIndex pi;
  SurfaceElementIndex sei;
  double weight;
};
} // namespace detail

Array<ArrayMem<detail::Neighbor, 20>>
BuildNeighbors (FlatArray<PointIndex> points, const Mesh& mesh)
{
  auto p2sel = mesh.CreatePoint2SurfaceElementTable();

  Array<ArrayMem<detail::Neighbor, 20>> neighbors(points.Size());

  ArrayMem<double, 20> angles;
  ArrayMem<double, 20> inv_dists;
  for (auto i : points.Range())
    {
      auto& p_neighbors = neighbors[i];
      auto pi = points[i];
      angles.SetSize(0);
      inv_dists.SetSize(0);
      for (auto sei : p2sel[pi])
        {
          const auto& sel = mesh[sei];
          for (auto pi1 : sel.PNums())
            {
              if (pi1 == pi)
                continue;
              auto pi2 = pi1;
              for (auto pi_ : sel.PNums())
                {
                  if (pi_ != pi && pi_ != pi1)
                    {
                      pi2 = pi_;
                      break;
                    }
                }
              p_neighbors.Append({pi1, sei, 0.0});
              inv_dists.Append(1.0 / (mesh[pi1] - mesh[pi]).Length());
              auto dot = (mesh[pi1] - mesh[pi]).Normalize() * (mesh[pi2] - mesh[pi]).Normalize();
              angles.Append(acos(dot));
            }
        }
      double sum_inv_dist = 0.0;
      for (auto inv_dist : inv_dists)
        sum_inv_dist += inv_dist;
      double sum_angle = 0.0;
      for (auto angle : angles)
        sum_angle += angle;

      double sum_weight = 0.0;
      for (auto i : Range(inv_dists))
        {
          p_neighbors[i].weight =
            inv_dists[i] * angles[i] / sum_inv_dist / sum_angle;
          sum_weight += p_neighbors[i].weight;
        }
      for (auto i : Range(inv_dists))
        p_neighbors[i].weight /= sum_weight;
    }
  return neighbors;
}

void BoundaryLayerTool ::InterpolateGrowthVectors()
{
  point_types.SetSize(mesh.GetNP());
  for (auto p : mesh.Points().Range())
    point_types[p] = mesh[p].Type();

  int new_max_edge_nr = max_edge_nr;
  for (const auto& seg : segments)
    if (seg.edgenr > new_max_edge_nr)
      new_max_edge_nr = seg.edgenr;
  for (const auto& seg : new_segments)
    if (seg.edgenr > new_max_edge_nr)
      new_max_edge_nr = seg.edgenr;

  auto getGW = [&] (PointIndex pi) -> Vec<3> {
    if (growth_vector_map.count(pi) == 0)
      growth_vector_map[pi] = {&growthvectors[pi], total_height};
    auto [gw, height] = growth_vector_map[pi];
    return height * (*gw);
  };
  auto addGW = [&] (PointIndex pi, Vec<3> vec) {
    if (growth_vector_map.count(pi) == 0)
      growth_vector_map[pi] = {&growthvectors[pi], total_height};
    auto [gw, height] = growth_vector_map[pi];
    *gw += 1.0 / height * vec;
  };

  // interpolate tangential component of growth vector along edge
  if (max_edge_nr >= new_max_edge_nr)
    return;

  auto edgenr2seg = ngcore::CreateSortedTable<Segment*, int>(
    Range(segments.Size() + new_segments.Size()),
    [&] (auto& table, size_t segi) {
      auto& seg = segi < segments.Size()
                    ? segments[segi]
                    : new_segments[segi - segments.Size()];
      table.Add(seg.edgenr, &seg);
    },
    new_max_edge_nr + 1);
  auto point2seg = ngcore::CreateSortedTable<Segment*, PointIndex>(
    Range(segments.Size() + new_segments.Size()),
    [&] (auto& table, size_t segi) {
      auto& seg = segi < segments.Size()
                    ? segments[segi]
                    : new_segments[segi - segments.Size()];
      table.Add(seg[0], &seg);
      table.Add(seg[1], &seg);
    },
    mesh.GetNP());

  for (auto edgenr : Range(1, new_max_edge_nr + 1))
    {
      // "inner" edges between two flat faces are not treated as edges for interpolation
      bool no_angles = true;
      ArrayMem<SurfaceElementIndex, 4> faces;

      for (auto* p_seg : edgenr2seg[edgenr])
        {
          auto& seg = *p_seg;
          faces.SetSize(0);
          // if (seg[0] <= p2sel.Size())
          if (seg[0] < IndexBASE<PointIndex>() + p2sel.Size())
            {
              for (auto sei : p2sel[seg[0]])
                if (moved_surfaces.Test(mesh[sei].GetIndex()) && p2sel[seg[1]].Contains(sei))
                  faces.Append(sei);
            }

          if (faces.Size() == 2)
            {
              auto n0 = getNormal(mesh[faces[0]]);
              auto n1 = getNormal(mesh[faces[1]]);
              if (n0 * n1 < 0.9)
                no_angles = false;
            }
          else
            {
              no_angles = false;
            }
        }
      if (no_angles)
        {
          for (auto* p_seg : edgenr2seg[edgenr])
            for (auto pi : p_seg->PNums())
              {
                if (pi >= first_new_pi)
                  continue;
                if (point_types[pi] == EDGEPOINT)
                  point_types[pi] = SURFACEPOINT;
                else if (point_types[pi] == FIXEDPOINT)
                  {
                    // Check at edge corners if all adjacent surface elements have roughly the same normal.
                    // If so, also treat this point as surface point for growth vector interpolation
                    Vec<3> n = 0.0;
                    for (auto si : p2sel[pi])
                      n += getNormal(mesh[si]);
                    n.Normalize();
                    bool is_corner = false;
                    for (auto si : p2sel[pi])
                      if (getNormal(mesh[si]) * n < 0.9)
                        is_corner = true;
                    if (!is_corner)
                      point_types[pi] = SURFACEPOINT;
                  }
              }
          continue;
        }
    }

  for (auto edgenr : Range(max_edge_nr + 1, new_max_edge_nr + 1))
    {
      double edge_len = 0.;
      bool any_grows = false;

      auto is_end_point = [&] (PointIndex pi) {
        auto segs = point2seg[pi];
        if (segs.Size() == 1)
          return true;
        auto first_edgenr = (*segs[0]).edgenr;
        for (auto* p_seg : segs)
          if (p_seg->edgenr != first_edgenr)
            return true;
        return false;
      };

      Array<PointIndex> points;
      for (auto* p_seg : edgenr2seg[edgenr])
        {
          auto& seg = *p_seg;

          if (getGW(seg[0]).Length2() != 0 || getGW(seg[1]).Length2() != 0)
            any_grows = true;

          if (points.Size() == 0)
            for (auto i : Range(2))
              if (is_end_point(seg[i]))
                {
                  points.Append(seg[i]);
                  points.Append(seg[1 - i]);
                  edge_len += (mesh[seg[1]] - mesh[seg[0]]).Length();
                  break;
                }
        }

      if (!any_grows)
        {
          PrintMessage(1, "BLayer: skip interpolating growth vectors at edge ", edgenr + 1);
          continue;
        }

      if (!points.Size())
        {
          if (debugparam.debugoutput)
            cerr << "Could not find startpoint for edge " << edgenr << endl;
          continue;
        }

      std::set<PointIndex> points_set;
      points_set.insert(points[0]);
      points_set.insert(points[1]);

      bool point_found = true;
      while (point_found)
        {
          if (is_end_point(points.Last()))
            break;
          point_found = false;
          for (auto* p_seg : point2seg[points.Last()])
            {
              const auto& seg = *p_seg;
              if (seg.edgenr != edgenr)
                continue;
              auto plast = points.Last();
              if (plast != seg[0] && plast != seg[1])
                continue;
              auto pnew = plast == seg[0] ? seg[1] : seg[0];
              if (pnew == points[0] && points.Size() > 1)
                {
                }
              if (points_set.count(pnew) > 0 && (pnew != points[0] || points.Size() == 2))
                continue;
              edge_len += (mesh[points.Last()] - mesh[pnew]).Length();
              points.Append(pnew);
              points_set.insert(pnew);
              point_found = true;
              break;
            }
        }
      if (!point_found)
        {
          if (debugparam.debugoutput)
            {
              cerr << "Could not find connected list of line segments for edge "
                   << edgenr << endl;
              cerr << "current points: " << endl
                   << points << endl;
            }
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
      for (auto i : IntRange(1, points.Size() - 1))
        {
          auto pi = points[i];
          len += (mesh[pi] - mesh[points[i - 1]]).Length();
          auto t = getEdgeTangent(pi, edgenr, point2seg[pi]);
          auto lam = len / edge_len;
          auto interpol = (1 - lam) * (gt1 * t) * t + lam * (gt2 * t) * t;
          addGW(pi, interpol);
        }
    }
}

void BoundaryLayerTool ::InterpolateSurfaceGrowthVectors()
{
  static Timer tall("InterpolateSurfaceGrowthVectors");
  RegionTimer rtall(tall);
  static Timer tsmooth("InterpolateSurfaceGrowthVectors-Smoothing");
  auto np_old = this->np;
  [[maybe_unused]] auto np = mesh.GetNP();

  auto hasMoved = [&] (PointIndex pi) {
    return (pi - IndexBASE<PointIndex>() >= np_old) || mapto[pi].Size() > 0 || special_boundary_points.count(pi);
  };

  std::set<PointIndex> points_set;
  for (const auto& sel : mesh.SurfaceElements())
    {
      for (auto pi : sel.PNums())
        if (point_types[pi] == SURFACEPOINT && hasMoved(pi))
          points_set.insert(pi);
    }

  Array<PointIndex> points;
  for (auto pi : points_set)
    points.Append(pi);
  QuickSort(points);

  // smooth tangential part of growth vectors from edges to surface elements
  Array<Vec<3>, PointIndex> corrections(mesh.GetNP());
  corrections = 0.0;
  RegionTimer rtsmooth(tsmooth);
  auto neighbors = BuildNeighbors(points, mesh);

  Array<Vec<3>, SurfaceElementIndex> surf_normals(mesh.GetNSE());
  for (auto sei : mesh.SurfaceElements().Range())
    surf_normals[sei] = getNormal(mesh[sei]);

  BitArray interpolate_tangent(mesh.GetNP() + 1);
  interpolate_tangent = false;
  for (auto pi : points)
    {
      for (auto sei : p2sel[pi])
        if (is_boundary_moved[mesh[sei].GetIndex()])
          interpolate_tangent.SetBit(pi);
    }

  constexpr int N_STEPS = 64;
  for ([[maybe_unused]] auto i : Range(N_STEPS))
    {
      for (auto i : points.Range())
        {
          auto pi = points[i];
          auto& p_neighbors = neighbors[i];

          ArrayMem<Vec<3>, 20> g_vectors;
          double max_len = 0.0;
          double sum_len = 0.0;

          // average only tangent component on new bl points, average whole growth
          // vector otherwise
          bool do_average_tangent = true;
          for (const auto& s : p_neighbors)
            {
              auto gw_other = growthvectors[s.pi] + corrections[s.pi];
              if (do_average_tangent)
                {
                  auto n = surf_normals[s.sei];
                  gw_other = gw_other - (gw_other * n) * n;
                }
              auto v = gw_other;
              auto len = v.Length2();
              sum_len += len;
              max_len = max(max_len, len);
              g_vectors.Append(v);
            }

          if (max_len == 0.0)
            continue;

          double lambda = 0;
          if (i > N_STEPS / 4.)
            lambda = 2.0 * (i - N_STEPS / 4.) / (N_STEPS / 2.);
          lambda = min(1.0, lambda);

          auto& correction = corrections[pi];
          correction = 0.0;
          for (const auto i : p_neighbors.Range())
            {
              auto v = g_vectors[i];
              double weight = lambda * p_neighbors[i].weight + (1.0 - lambda) * v.Length2() / sum_len;
              correction += weight * v;
            }

          if (!do_average_tangent)
            correction -= growthvectors[pi];
        }
    }

  for (auto pi : points)
    growthvectors[pi] += corrections[pi];
}

void BoundaryLayerTool ::FixSurfaceElements()
{
  static Timer tall("FixSurfaceElements");
  RegionTimer rtall(tall);
  [[maybe_unused]] auto np_old = this->np;
  [[maybe_unused]] auto np = mesh.GetNP();

  non_bl_growth_vectors.clear();

  auto getGW = [&] (PointIndex pi) -> Vec<3> {
    // return growthvectors[pi];
    if (growth_vector_map.count(pi) == 0)
      {
        non_bl_growth_vectors[pi] = .0;
        growth_vector_map[pi] = {&non_bl_growth_vectors[pi], 1.0};
      }
    auto [gw, height] = growth_vector_map[pi];
    return height * (*gw);
  };

  auto addGW = [&] (PointIndex pi, Vec<3> vec) {
    if (growth_vector_map.count(pi) == 0)
      {
        non_bl_growth_vectors[pi] = .0;
        growth_vector_map[pi] = {&non_bl_growth_vectors[pi], 1.0};
      }
    auto [gw, height] = growth_vector_map[pi];
    *gw += 1.0 / height * vec;
  };

  std::set<PointIndex> points_set;
  // only smooth over old surface elements
  for (SurfaceElementIndex sei : Range(nse))
    {
      const auto& sel = mesh[sei];
      if (sel.GetNP() == 3 && is_boundary_moved[sel.GetIndex()])
        for (auto pi : sel.PNums())
          if (point_types[pi] == SURFACEPOINT)
            points_set.insert(pi);
    }

  Array<PointIndex> points;
  for (auto pi : points_set)
    points.Append(pi);
  QuickSort(points);

  Array<Vec<3>, PointIndex> corrections(mesh.GetNP());
  corrections = 0.0;

  auto neighbors = BuildNeighbors(points, mesh);

  constexpr int N_STEPS = 32;
  for ([[maybe_unused]] auto i : Range(N_STEPS))
    {
      for (auto i : points.Range())
        {
          auto pi = points[i];
          auto& p_neighbors = neighbors[i];

          ArrayMem<Vec<3>, 20> g_vectors;
          double max_len = 0.0;
          double sum_len = 0.0;

          for (const auto& s : p_neighbors)
            {
              auto v = getGW(s.pi) + corrections[s.pi];
              auto len = v.Length2();
              sum_len += len;
              max_len = max(max_len, len);
              g_vectors.Append(v);
            }

          if (max_len == 0.0)
            continue;

          double lambda = 0;
          if (i > N_STEPS / 4.)
            lambda = 2.0 * (i - N_STEPS / 4.) / (N_STEPS / 2.);
          lambda = min(1.0, lambda);

          auto& correction = corrections[pi];
          correction = 0.0;
          for (const auto i : p_neighbors.Range())
            {
              auto v = g_vectors[i];
              double weight = lambda * p_neighbors[i].weight + (1.0 - lambda) * v.Length2() / sum_len;
              correction += weight * v;
            }
        }
    }

  for (auto pi : points)
    addGW(pi, corrections[pi]);
}

} // namespace netgen

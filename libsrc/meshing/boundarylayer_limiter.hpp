#include "boundarylayer.hpp"

namespace netgen {

struct Intersection_ {
  bool is_intersecting = false;
  double lam0 = -1, lam1 = -1;
  Point<3> p;
  double bary[3];
  operator bool() const { return is_intersecting; }
};

struct GrowthVectorLimiter {
  BoundaryLayerTool &tool;
  const BoundaryLayerParameters &params;
  Mesh &mesh;
  double height;
  FlatArray<double, PointIndex> limits;
  FlatArray<Vec<3>, PointIndex> growthvectors;
  BitArray changed_domains;
  unique_ptr<BoxTree<3>> tree;
  Array<PointIndex, PointIndex> map_from;
  ofstream debug;

  GrowthVectorLimiter(BoundaryLayerTool &tool_)
      : tool(tool_), params(tool_.params), mesh(tool_.mesh),
        height(tool_.total_height), limits(tool_.limits),
        growthvectors(tool_.growthvectors), map_from(mesh.Points().Size()),
        debug("debug.txt") {
    changed_domains = tool.domains;
    if (!params.outside)
      changed_domains.Invert();

    map_from = tool.mapfrom;
  }

  double GetLimit(PointIndex pi) {
    if (pi <= tool.np)
      return limits[pi];
    return limits[map_from[pi]];
  }

  bool SetLimit(PointIndex pi, double new_limit) {
    double &limit = (pi <= tool.np) ? limits[pi] : limits[map_from[pi]];
    if (limit <= new_limit)
      return false;
    limit = new_limit;
    return true;
  }

  bool ScaleLimit(PointIndex pi, double factor) {
    double &limit = (pi <= tool.np) ? limits[pi] : limits[map_from[pi]];
    return SetLimit(pi, limit * factor);
  }

  Point<3> GetPoint(PointIndex pi_to, double shift = 1.,
                    bool apply_limit = false) {
    if (tool.growth_vector_map.count(pi_to) == 0)
      return mesh[pi_to];

    auto [gw, height] = tool.growth_vector_map[pi_to];
    if (apply_limit)
      shift *= GetLimit(pi_to);
    return mesh[pi_to] + shift * height * (*gw);
  }

  Point<3> GetMappedPoint(PointIndex pi_from, double shift = 1.) {
    auto pi_to = tool.mapto[pi_from].Last();
    return GetPoint(pi_to, shift);
  }

  std::array<Point<3>, 2> GetMappedSeg(PointIndex pi_from, double shift = 1.) {
    return {mesh[pi_from], GetMappedPoint(pi_from, shift)};
  }

  std::array<Point<3>, 2> GetSeg(PointIndex pi_to, double shift = 1.,
                                 bool apply_limit = false) {
    return {GetPoint(pi_to, 0), GetPoint(pi_to, shift, apply_limit)};
  }

  auto GetTrig(SurfaceElementIndex sei, double shift = 0.0,
               bool apply_limit = false) {
    auto sel = mesh[sei];
    std::array<Point<3>, 3> trig;
    for (auto i : Range(3))
      trig[i] = GetPoint(sel[i], shift, apply_limit);
    return trig;
  }

  auto GetMappedTrig(SurfaceElementIndex sei, double shift = 0.0) {
    auto sel = mesh[sei];
    std::array<Point<3>, 3> trig;
    for (auto i : Range(3))
      trig[i] = GetMappedPoint(sel[i], shift);
    return trig;
  }

  auto GetSideTrig(SurfaceElementIndex sei, int index, double shift = 0.0,
                   bool grow_first_vertex = true) {
    auto trig = GetMappedTrig(sei, 0.0);
    auto sel = mesh[sei];
    auto index1 = (index + 1) % 3;
    if (!grow_first_vertex)
      index1 = (index + 2) % 3;
    trig[index] = GetMappedPoint(sel[index1], shift);
    return trig;
  }

  static constexpr double INTERSECTION_SAFETY = .9;
  bool LimitGrowthVector(PointIndex pi_to, SurfaceElementIndex sei,
                         double trig_shift, double seg_shift,
                         bool check_prism_sides = false) {
    auto pi_from = map_from[pi_to];
    if (!pi_from.IsValid())
      return false;

    auto seg = GetSeg(pi_to, seg_shift, true);

    for (auto pi : mesh[sei].PNums()) {
      if (pi == pi_from)
        return false;
      if (map_from[pi] == pi_from)
        return false;
    }

    if (check_prism_sides) {
      for (auto i : Range(3)) {
        auto side = GetSideTrig(sei, i, trig_shift, true);
        auto intersection = isIntersectingTrig(seg, side);
        if (intersection)
          return ScaleLimit(pi_to, intersection.lam0 * INTERSECTION_SAFETY);
      }
      return false;
    } else if (trig_shift > 0) {
      auto intersection =
          isIntersectingTrig(seg, GetTrig(sei, trig_shift, true));
      if (!intersection)
        return false;

      double scaling_factor = 0.9;
      double s = 1.0;

      while (true) {
        s *= scaling_factor;
        auto reduced_intersection =
            isIntersectingTrig(GetSeg(pi_to, s * seg_shift, true),
                               GetTrig(sei, s * trig_shift, true));
        if (!reduced_intersection)
          break;
      }

      // cout << "Scale limits " << s << endl;

      bool result = false;
      result |= ScaleLimit(pi_to, s);
      for (auto pi : mesh[sei].PNums())
        result |= ScaleLimit(pi, s);
      return result;

      double dshift = trig_shift;
      double lam0 = intersection.lam0 * seg_shift * GetLimit(pi_from);
      while (dshift / trig_shift > lam0) {
        dshift *= 0.9;
        auto reduced_intersection =
            isIntersectingTrig(seg, GetTrig(sei, dshift, true));
        if (!reduced_intersection)
          break;
        // cout << "still intersecting " << dshift*trig_shift << " > " << lam0
        // << endl;
        intersection = reduced_intersection;
      }
      lam0 = intersection.lam0 * seg_shift;
      double max_trig_limit = 1e99;
      auto sel = mesh[sei];
      for (auto i : Range(3))
        max_trig_limit = min(max_trig_limit, GetLimit(sel[i]));

      double new_seg_limit = lam0 * INTERSECTION_SAFETY;
      double new_trig_limit = dshift * trig_shift * INTERSECTION_SAFETY;

      if (new_trig_limit >= max_trig_limit &&
          new_seg_limit >= GetLimit(pi_from))
        return false; // nothing to do

      result = false;
      result |= SetLimit(pi_from, new_seg_limit);
      for (auto pi : sel.PNums())
        result |= SetLimit(pi, new_trig_limit);

      return result;
    } else {
      auto trig = GetTrig(sei, 0.0);
      auto intersection = isIntersectingTrig(seg, trig);
      // checking with original surface elements -> allow only half the distance
      auto new_seg_limit = 0.40 * intersection.lam0 * seg_shift;
      if (intersection && new_seg_limit < GetLimit(pi_from)) {
        auto p0 = seg[0];
        auto p1 = seg[1];
        auto d = Dist(p0, p1);
        auto [gw, height] = tool.growth_vector_map[pi_to];
        return SetLimit(pi_from, new_seg_limit);
      }
      return false;
    }
  }

  void LimitSelfIntersection() {
    // check for self-intersection within new elements (prisms/hexes)
    bool found_debug_element = false;
    auto isIntersecting = [&](SurfaceElementIndex sei, double shift) {
      // checks if surface element is self intersecting when growing with factor
      // shift

      // ignore new surface elements, side trigs are only built
      // from original surface elements
      if (sei >= tool.nse)
        return false;
      const auto sel = mesh[sei];
      auto np = sel.GetNP();
      for (auto i : Range(np)) {
        if (sel[i] > tool.np)
          return false;
        if (tool.mapto[sel[i]].Size() == 0)
          return false;
      }
      for (auto i : Range(np)) {
        auto seg = GetMappedSeg(sel[i], shift * limits[sel[i]]);
        for (auto fi : Range(np - 2)) {
          for (auto side : {true, false}) {
            auto trig = GetSideTrig(sei, i + fi, 1.0, side);
            if (isIntersectingPlane(seg, trig))
              return true;
          }
        }
      }
      return false;
    };

    auto equalizeLimits = [&](SurfaceElementIndex sei) {
      const auto sel = mesh[sei];
      auto np = sel.GetNP();
      double max_limit = 0;
      double min_limit = 1e99;
      for (auto i : Range(np)) {
        max_limit = max(max_limit, limits[sel[i]]);
        min_limit = min(min_limit, limits[sel[i]]);
      }
      // equalize
      if (max_limit / min_limit > 1.2) {
        max_limit = min_limit * 1.2;
        for (auto i : Range(np))
          SetLimit(sel[i], min(limits[sel[i]], max_limit));
      }
    };

    for (SurfaceElementIndex sei : mesh.SurfaceElements().Range()) {
      auto sel = mesh[sei];
      const auto &fd = mesh.GetFaceDescriptor(sel.GetIndex());
      if (sei >= tool.nse)
        continue;
      if (sel.GetNP() == 4)
        continue;
      // if(sei >= tool.nse || (!changed_domains.Test(fd.DomainIn()) &&
      //   !changed_domains.Test(fd.DomainOut())))
      //   continue;

      auto np = sel.GetNP();
      // ArrayMem<double, 4> ori_limits;
      // ori_limits.SetSize(np);
      // for(auto i : Range(np))
      //     ori_limits[i] = limits[sel[i]];

      equalizeLimits(sei);

      double shift = 1.0;
      double safety = 1.4;
      const double step_factor = 0.9;
      while (isIntersecting(sei, shift * safety)) {
        shift *= step_factor;
        double max_limit = 0;
        for (auto i : Range(np))
          max_limit = max(max_limit, limits[sel[i]]);
        for (auto i : Range(np))
          if (max_limit == limits[sel[i]])
            ScaleLimit(sel[i], step_factor);
        // if (max_limit < 0.01) break;
      }
    }
  }

  // checks if a segment is intersecting a plane, spanned by three points, lam
  // will be set s.t. p_intersect = seg[0] + lam * (seg[1]-seg[0])
  Intersection_ isIntersectingPlane(std::array<Point<3>, 2> seg,
                                    std::array<Point<3>, 3> trig) {
    auto t1 = trig[1] - trig[0];
    auto t2 = trig[2] - trig[0];
    auto n = Cross(t1, t2);
    auto v0n = (seg[0] - trig[0]) * n;
    auto v1n = (seg[1] - trig[0]) * n;

    Intersection_ intersection;
    intersection.lam0 = -v0n / (v1n - v0n);
    intersection.p = seg[0] + intersection.lam0 * (seg[1] - seg[0]);
    intersection.is_intersecting = (v0n * v1n < 0) &&
                                   (intersection.lam0 > -1e-8) &&
                                   (intersection.lam0 < 1 + 1e-8);

    return intersection;
  }

  // Intersection_ isIntersectingPlane(PointIndex pi, PointIndex pi_to,
  //                                   SurfaceElementIndex sei,
  //                                   double shift = 0.0) {
  //   return isIntersectingPlane(GetSeg(pi, pi_to), GetTrig(sei, shift));
  // }

  Intersection_ isIntersectingTrig(std::array<Point<3>, 2> seg,
                                   std::array<Point<3>, 3> trig) {
    auto intersection = isIntersectingPlane(seg, trig);
    if (!intersection)
      return intersection;

    auto p = seg[0] + intersection.lam0 * (seg[1] - seg[0]) - trig[0];

    Vec3d col1 = trig[1] - trig[0];
    Vec3d col2 = trig[2] - trig[0];
    Vec3d col3 = Cross(col1, col2);
    Vec3d rhs = p;
    Vec3d bary;
    SolveLinearSystem(col1, col2, col3, rhs, bary);

    intersection.lam1 = 0;
    double eps = 0.1;
    if (bary.X() >= -eps && bary.Y() >= -eps &&
        bary.X() + bary.Y() <= 1 + eps) {
      intersection.bary[0] = bary.X();
      intersection.bary[1] = bary.Y();
      intersection.bary[2] = 1.0 - bary.X() - bary.Y();
    } else
      intersection.is_intersecting = false;
    return intersection;
  }

  Intersection_ isIntersectingTrig(PointIndex pi_from, PointIndex pi_to,
                                   SurfaceElementIndex sei,
                                   double shift = 0.0) {
    return isIntersectingTrig(GetSeg(pi_from, pi_to), GetTrig(sei, shift));
  }

  void BuildSearchTree(double trig_shift) {
    Box<3> bbox(Box<3>::EMPTY_BOX);
    for (PointIndex pi : mesh.Points().Range()) {
      bbox.Add(mesh[pi]);
      bbox.Add(GetPoint(pi, 1.1));
      // if(tool.mapto[pi].Size() >0)
      // bbox.Add(mesh[tool.mapto[pi].Last()]);
    }

    tree = make_unique<BoxTree<3>>(bbox);

    for (auto sei : mesh.SurfaceElements().Range()) {
      const auto &sel = mesh[sei];
      auto sel_index = mesh[sei].GetIndex();

      Box<3> box(Box<3>::EMPTY_BOX);
      for (auto pi : sel.PNums()) {
        box.Add(GetPoint(pi, 0.));
        box.Add(GetPoint(pi, trig_shift * GetLimit(pi)));
      }
      tree->Insert(box, sei);
    }
  }

  template <typename TFunc>
  void FindTreeIntersections(double trig_shift, double seg_shift, TFunc f) {
    BuildSearchTree(trig_shift);
    auto np_new = mesh.Points().Size();
    int counter = 0;
    for (auto i : IntRange(tool.np, np_new)) {
      PointIndex pi_to = i + PointIndex::BASE;
      PointIndex pi_from = map_from[pi_to];
      if (!pi_from.IsValid())
        throw Exception("Point not mapped");

      // if(mesh[pi_to].Type() == INNERPOINT)
      // continue;
      // if(growthvectors[pi_to].Length2() == 0.0)
      // continue;
      Box<3> box(Box<3>::EMPTY_BOX);
      auto seg = GetSeg(pi_to, seg_shift);

      box.Add(GetPoint(pi_to, 0));
      box.Add(GetPoint(pi_to, GetLimit(pi_from)));
      tree->GetFirstIntersecting(box.PMin(), box.PMax(),
                                 [&](SurfaceElementIndex sei) {
                                   const auto &sel = mesh[sei];
                                   if (sel.PNums().Contains(pi_from))
                                     return false;
                                   if (sel.PNums().Contains(pi_to))
                                     return false;
                                   counter++;
                                   f(pi_to, sei);
                                   return false;
                                 });
    }
  }

  void Perform() {
    tool.limits.SetSize(mesh.Points().Size());
    tool.limits = 1.0;

    // limit to not intersect with other (original) surface elements
    double trig_shift = 0;
    double seg_shift = 2.1;
    FindTreeIntersections(
        trig_shift, seg_shift, [&](PointIndex pi_to, SurfaceElementIndex sei) {
          if (sei >= tool.nse)
            return; // ignore new surface elements in first pass
          LimitGrowthVector(pi_to, sei, trig_shift, seg_shift);
        });

    LimitSelfIntersection();

    // for(auto i : Range(growthvectors))
    //   growthvectors[i] *= limits[i];
    // limits = 1.0;

    // now limit again with shifted surface elements
    trig_shift = 1.1;
    seg_shift = 1.1;
    size_t limit_counter = 1;

    while (limit_counter) {
      limit_counter = 0;
      FindTreeIntersections(
          trig_shift, seg_shift,
          [&](PointIndex pi_to, SurfaceElementIndex sei) {
            if (LimitGrowthVector(pi_to, sei, trig_shift, seg_shift))
              limit_counter++;
            auto sel = mesh[sei];
            bool is_mapped = true;
            for (auto pi : sel.PNums()) {
              if (pi >= tool.np)
                return;
              if (tool.mapto[pi].Size() == 0)
                return;
            }
            if (LimitGrowthVector(pi_to, sei, trig_shift, seg_shift, true))
              limit_counter++;
          });
    }

    // check if surface trigs are intersecting each other
    {
      Point3d pmin, pmax;
      mesh.GetBox(pmin, pmax);
      BoxTree<3, SurfaceElementIndex> setree(pmin, pmax);

      for (auto sei : mesh.SurfaceElements().Range()) {
        const Element2d &tri = mesh[sei];

        Box<3> box(Box<3>::EMPTY_BOX);
        for (PointIndex pi : tri.PNums())
          box.Add(GetPoint(pi, 1.0, true));

        box.Increase(1e-3 * box.Diam());
        setree.Insert(box, sei);
      }
      for (auto sei : mesh.SurfaceElements().Range()) {
        const Element2d &tri = mesh[sei];

        Box<3> box(Box<3>::EMPTY_BOX);
        for (PointIndex pi : tri.PNums())
          box.Add(GetPoint(pi, 1.0, true));

        setree.GetFirstIntersecting(
            box.PMin(), box.PMax(), [&](SurfaceElementIndex sej) {
              const Element2d &tri2 = mesh[sej];

              if (mesh[tri[0]].GetLayer() != mesh[tri2[0]].GetLayer())
                return false;

              netgen::Point<3> tri1_points[3], tri2_points[3];
              const netgen::Point<3> *trip1[3], *trip2[3];
              for (int k = 0; k < 3; k++) {
                trip1[k] = &tri1_points[k];
                trip2[k] = &tri2_points[k];
              }
              auto set_points = [&]() {
                for (int k = 0; k < 3; k++) {
                  tri1_points[k] = GetPoint(tri[k], 1.0, true);
                  tri2_points[k] = GetPoint(tri2[k], 1.0, true);
                }
              };

              set_points();

              int counter = 0;
              while (IntersectTriangleTriangle(&trip1[0], &trip2[0])) {
                PointIndex pi_max_limit = PointIndex::INVALID;
                for (PointIndex pi :
                     {tri[0], tri[1], tri[2], tri2[0], tri2[1], tri2[2]})
                  if (pi > tool.np &&
                      (!pi_max_limit.IsValid() ||
                       limits[tool.mapfrom[pi]] > limits[pi_max_limit]))
                    pi_max_limit = tool.mapfrom[pi];

                if (!pi_max_limit.IsValid())
                  break;

                limits[pi_max_limit] *= 0.9;
                set_points();
                counter++;
                if (counter > 20) {
                  cerr << "Limit intersecting sourface elements: too many "
                          "limitation steps"
                       << endl;
                  break;
                }
              }
              return false;
            });
      }
    }

    for (auto i : Range(growthvectors))
      growthvectors[i] *= limits[i];
    for (auto &[special_pi, special_point] : tool.special_boundary_points) {
      for (auto &group : special_point.growth_groups) {
        group.growth_vector *= limits[special_pi];
      }
    }
  }
};

} // namespace netgen

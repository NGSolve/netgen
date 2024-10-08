#include "boundarylayer.hpp"
#include <core/array.hpp>

namespace netgen {

struct Intersection_ {
  bool is_intersecting = false;
  double lam0 = -1, lam1 = -1;
  Point<3> p;
  double bary[3];
  operator bool() const { return is_intersecting; }
};

struct GrowthVectorLimiter {
  typedef std::array<Point<3>, 2> Seg;
  typedef std::array<Point<3>, 3> Trig;

  BoundaryLayerTool &tool;
  const BoundaryLayerParameters &params;
  Mesh &mesh;
  double height;
  Array<double, PointIndex> limits;
  FlatArray<Vec<3>, PointIndex> growthvectors;
  BitArray changed_domains;
  unique_ptr<BoxTree<3>> tree;
  Array<PointIndex, PointIndex> map_from;
  Table<SurfaceElementIndex, PointIndex> p2sel;

  GrowthVectorLimiter(BoundaryLayerTool &tool_)
      : tool(tool_), params(tool_.params), mesh(tool_.mesh),
        height(tool_.total_height), growthvectors(tool_.growthvectors),
        map_from(mesh.Points().Size()) {
    changed_domains = tool.domains;
    if (!params.outside)
      changed_domains.Invert();

    map_from = tool.mapfrom;
    p2sel = ngcore::CreateSortedTable<SurfaceElementIndex, PointIndex>(
        tool.new_sels.Range(),
        [&](auto &table, SurfaceElementIndex ei) {
          for (PointIndex pi : tool.new_sels[ei].PNums())
            table.Add(pi, ei);
        },
        mesh.GetNP());
  }

  auto SurfaceElementsRange() { return Range(tool.nse + tool.new_sels.Size()); }

  const auto &Get(SurfaceElementIndex sei) {
    if (sei < tool.nse)
      return mesh[sei];
    return tool.new_sels[sei - tool.nse];
  }

  std::pair<double, double> GetMinMaxLimit(SurfaceElementIndex sei) {
    const auto &sel = Get(sei);
    double min_limit = GetLimit(sel[0]);
    double max_limit = min_limit;
    for (auto i : IntRange(1, sel.GetNP())) {
      auto limit = GetLimit(sel[i]);
      min_limit = min(min_limit, limit);
      max_limit = max(max_limit, limit);
    }
    return {min_limit, max_limit};
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

  Vec<3> GetVector(PointIndex pi_to, double shift = 1.,
                   bool apply_limit = false) {
    auto [gw, height] = tool.growth_vector_map[pi_to];
    if (apply_limit)
      shift *= GetLimit(pi_to);
    return shift * height * (*gw);
  }

  Point<3> GetPoint(PointIndex pi_to, double shift = 1.,
                    bool apply_limit = false) {
    if (pi_to <= tool.np || tool.growth_vector_map.count(pi_to) == 0)
      return mesh[pi_to];

    return mesh[pi_to] + GetVector(pi_to, shift, apply_limit);
  }

  Point<3> GetMappedPoint(PointIndex pi_from, double shift = 1.) {
    auto pi_to = tool.mapto[pi_from].Last();
    return GetPoint(pi_to, shift);
  }

  Seg GetMappedSeg(PointIndex pi_from, double shift = 1.) {
    return {mesh[pi_from], GetMappedPoint(pi_from, shift)};
  }

  Seg GetSeg(PointIndex pi_to, double shift = 1., bool apply_limit = false) {
    return {GetPoint(pi_to, 0), GetPoint(pi_to, shift, apply_limit)};
  }

  Trig GetTrig(SurfaceElementIndex sei, double shift = 0.0,
               bool apply_limit = false) {
    auto sel = Get(sei);
    Trig trig;
    for (auto i : Range(3))
      trig[i] = GetPoint(sel[i], shift, apply_limit);
    return trig;
  }

  Trig GetMappedTrig(SurfaceElementIndex sei, double shift = 0.0) {
    auto sel = Get(sei);
    Trig trig;
    for (auto i : Range(3))
      trig[i] = GetMappedPoint(sel[i], shift);
    return trig;
  }

  Trig GetSideTrig(SurfaceElementIndex sei, int index, double shift = 0.0,
                   bool grow_first_vertex = true) {
    auto trig = GetMappedTrig(sei, 0.0);
    auto sel = Get(sei);
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

    for (auto pi : Get(sei).PNums()) {
      if (pi == pi_from)
        return false;
      if (map_from[pi] == pi_from)
        return false;
    }

    if (check_prism_sides || trig_shift > .0) {
      auto [trig_min_limit, trig_max_limit] = GetMinMaxLimit(sei);
      if (GetLimit(pi_to) < trig_min_limit)
        return false;

      auto getTrigs = [&](double scaling = 1.0) -> ArrayMem<Trig, 3> {
        ArrayMem<Trig, 3> trigs;
        if (check_prism_sides)
          for (auto i : Range(3))
            trigs.Append(GetSideTrig(sei, i, scaling * trig_shift, true));
        else
          trigs.Append(GetTrig(sei, scaling * trig_shift, true));
        return trigs;
      };

      double scaling = 1.0;
      while (true) {
        bool have_intersection = false;
        auto seg = GetSeg(pi_to, scaling * seg_shift, true);
        for (auto trig : getTrigs(scaling))
          have_intersection |= isIntersectingTrig(seg, trig);
        if (!have_intersection)
          break;
        scaling *= 0.9;
      }
      if (scaling == 1.0)
        return false;

      double new_limit = scaling * max(GetLimit(pi_to), trig_max_limit);
      SetLimit(pi_to, new_limit);
      for (auto pi : Get(sei).PNums())
        SetLimit(pi, new_limit);
      return true;
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

  void EqualizeLimits(double factor = .5) {
    static Timer t("GrowthVectorLimiter::EqualizeLimits");
    PrintMessage(5, "GrowthVectorLimiter - equalize limits");
    RegionTimer reg(t);
    if (factor == 0.0)
      return;
    for (PointIndex pi : IntRange(tool.np, mesh.GetNP())) {
      auto pi_from = map_from[pi];
      std::set<PointIndex> pis;
      for (auto sei : p2sel[pi])
        for (auto pi_ : tool.new_sels[sei].PNums())
          pis.insert(pi_);
      ArrayMem<double, 20> limits;
      for (auto pi1 : pis) {
        auto limit = GetLimit(pi1);
        if (limit > 0.0)
          limits.Append(GetLimit(pi1));
      }
      if (limits.Size() == 0)
        continue;
      QuickSort(limits);

      double mean_limit = limits[limits.Size() / 2];
      // if mean limit is the maximum limit, take the average of second-highest
      // and highest value
      if (mean_limit > limits[0] && mean_limit == limits.Last()) {
        auto i = limits.Size() - 1;
        while (limits[i] == limits.Last())
          i--;
        mean_limit = 0.5 * (limits[i] + limits.Last());
      }

      if (limits.Size() % 2 == 0)
        mean_limit = 0.5 * (mean_limit + limits[(limits.Size() - 1) / 2]);

      SetLimit(pi, factor * mean_limit + (1.0 - factor) * GetLimit(pi));
    }
  }

  void LimitSelfIntersection(double safety = 1.4) {
    static Timer t("GrowthVectorLimiter::LimitSelfIntersection");
    PrintMessage(5, "GrowthVectorLimiter - self intersection");
    RegionTimer reg(t);
    // check for self-intersection within new elements (prisms/hexes)
    auto isIntersecting = [&](SurfaceElementIndex sei, double shift) {
      // checks if surface element is self intersecting when growing with factor
      // shift

      // ignore new surface elements, side trigs are only built
      // from original surface elements
      if (sei >= tool.nse)
        return false;
      const auto sel = Get(sei);
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

    for (SurfaceElementIndex sei : mesh.SurfaceElements().Range()) {
      auto sel = mesh[sei];
      const auto &fd = mesh.GetFaceDescriptor(sel.GetIndex());
      if (sel.GetNP() == 4)
        continue;

      auto np = sel.GetNP();

      double shift = 1.0;
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
  Intersection_ isIntersectingPlane(const Seg & seg,
                                    const Trig & trig) {
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

  Intersection_ isIntersectingTrig(const Seg &seg, const Trig &trig) {
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
    static Timer t("BuildSearchTree");
    RegionTimer rt(t);
    Box<3> bbox(Box<3>::EMPTY_BOX);
    for (PointIndex pi : mesh.Points().Range()) {
      bbox.Add(mesh[pi]);
      bbox.Add(GetPoint(pi, 1.1));
    }

    tree = make_unique<BoxTree<3>>(bbox);

    for (auto sei : SurfaceElementsRange()) {
      const auto &sel = Get(sei);
      auto sel_index = sel.GetIndex();

      Box<3> box(Box<3>::EMPTY_BOX);
      for (auto pi : sel.PNums()) {
        box.Add(GetPoint(pi, 0.));
        box.Add(GetPoint(pi, trig_shift * GetLimit(pi)));
      }
      tree->Insert(box, sei);
    }
  }

  template <typename TFunc>
  void FindTreeIntersections(double trig_shift, double seg_shift, TFunc f,
                             BitArray *relevant_points = nullptr) {
    static Timer t("GrowthVectorLimiter::FindTreeIntersections");
    RegionTimer rt(t);
    BuildSearchTree(trig_shift);
    auto np_new = mesh.Points().Size();
    int counter = 0;
    for (auto i : IntRange(tool.np, np_new)) {
      PointIndex pi_to = i + PointIndex::BASE;
      PointIndex pi_from = map_from[pi_to];
      if (!pi_from.IsValid())
        throw Exception("Point not mapped");

      if (relevant_points && !relevant_points->Test(pi_to) &&
          !relevant_points->Test(pi_from))
        continue;

      Box<3> box(Box<3>::EMPTY_BOX);
      auto seg = GetSeg(pi_to, seg_shift);

      box.Add(GetPoint(pi_to, 0));
      box.Add(GetPoint(pi_to, GetLimit(pi_from)));
      tree->GetFirstIntersecting(box.PMin(), box.PMax(),
                                 [&](SurfaceElementIndex sei) {
                                   const auto &sel = Get(sei);
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

  void FixIntersectingSurfaceTrigs() {
    static Timer t("GrowthVectorLimiter::FixIntersectingSurfaceTrigs");
    RegionTimer reg(t);
    // check if surface trigs are intersecting each other
    bool changed = true;
    while (changed) {
      changed = false;
      Point3d pmin, pmax;
      mesh.GetBox(pmin, pmax);
      BoxTree<3, SurfaceElementIndex> setree(pmin, pmax);

      for (auto sei : SurfaceElementsRange()) {
        const Element2d &tri = Get(sei);

        Box<3> box(Box<3>::EMPTY_BOX);
        for (PointIndex pi : tri.PNums())
          box.Add(GetPoint(pi, 1.0, true));

        box.Increase(1e-3 * box.Diam());
        setree.Insert(box, sei);
      }

      for (auto sei : SurfaceElementsRange()) {
        const Element2d &tri = Get(sei);

        Box<3> box(Box<3>::EMPTY_BOX);
        for (PointIndex pi : tri.PNums())
          box.Add(GetPoint(pi, 1.0, true));

        setree.GetFirstIntersecting(box.PMin(), box.PMax(), [&](size_t sej) {
          const Element2d &tri2 = Get(sej);

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
            changed = true;
            PointIndex pi_max_limit = PointIndex::INVALID;
            for (PointIndex pi :
                 {tri[0], tri[1], tri[2], tri2[0], tri2[1], tri2[2]})
              if (pi > tool.np && (!pi_max_limit.IsValid() ||
                                   GetLimit(pi) > GetLimit(pi_max_limit)))
                pi_max_limit = map_from[pi];

            if (!pi_max_limit.IsValid())
              break;

            ScaleLimit(pi_max_limit, 0.9);
            set_points();
            counter++;
            if (counter > 20) {
              cerr << "Limit intersecting surface elements: too many "
                      "limitation steps, sels: "
                   << Get(sei) << '\t' << Get(sej) << endl;
              for (auto si : {sei, sej}) {
                auto sel = Get(si);
                cerr << "Limits: ";
                for (auto pi : sel.PNums())
                  cerr << GetLimit(pi) << ",\t";
                cerr << endl;
                for (auto pi : sel.PNums())
                  cerr << GetPoint(pi, 1.0, true) << "\t";
                cerr << endl;
              }
              cerr << "pi_max_limit " << pi_max_limit << endl;
              break;
            }
          }
          return false;
        });
      }
    }
  }

  void LimitOriginalSurface(double safety) {
    static Timer t("GrowthVectorLimiter::LimitOriginalSurface");
    RegionTimer reg(t);
    PrintMessage(5, "GrowthVectorLimiter - original surface");
    // limit to not intersect with other (original) surface elements
    double trig_shift = 0;
    double seg_shift = safety;
    FindTreeIntersections(
        trig_shift, seg_shift, [&](PointIndex pi_to, SurfaceElementIndex sei) {
          if (sei >= tool.nse)
            return; // ignore new surface elements in first pass
          LimitGrowthVector(pi_to, sei, trig_shift, seg_shift);
        });
  }

  void LimitBoundaryLayer(double safety = 1.1) {
    static Timer t("GrowthVectorLimiter::LimitBoundaryLayer");
    PrintMessage(5, "GrowthVectorLimiter - boundary layer");
    // now limit again with shifted surface elements
    double trig_shift = safety;
    double seg_shift = safety;
    size_t limit_counter = 1;

    BitArray relevant_points, relevant_points_next;
    relevant_points.SetSize(mesh.Points().Size() + 1);
    relevant_points_next.SetSize(mesh.Points().Size() + 1);
    relevant_points.Set();

    while (limit_counter) {
      RegionTimer reg(t);
      size_t find_counter = 0;
      limit_counter = 0;
      relevant_points_next.Clear();
      FindTreeIntersections(
          trig_shift, seg_shift,
          [&](PointIndex pi_to, SurfaceElementIndex sei) {
            find_counter++;
            auto sel = Get(sei);

            if (LimitGrowthVector(pi_to, sei, trig_shift, seg_shift)) {
              limit_counter++;
              relevant_points_next.SetBit(pi_to);
              relevant_points_next.SetBit(map_from[pi_to]);
              for (auto pi : sel.PNums()) {
                relevant_points_next.SetBit(pi);
                if (pi >= tool.np)
                  relevant_points_next.SetBit(map_from[pi]);
                else
                  relevant_points_next.SetBit(map_from[pi]);
              }
            }

            for (auto pi : sel.PNums()) {
              if (pi >= tool.np)
                return;
              if (tool.mapto[pi].Size() == 0)
                return;
            }
            if (LimitGrowthVector(pi_to, sei, trig_shift, seg_shift, true))
              limit_counter++;
          },
          &relevant_points);
      relevant_points = relevant_points_next;
    }
  }

  void Perform() {
    limits.SetSize(mesh.Points().Size());
    limits = 1.0;

    std::array safeties = {0.5, 1.1, 1.5, 1.5};

    // No smoothing in the last pass, to avoid generating new intersections
    std::array smoothing_factors = {0.8, 0.7, 0.5, 0.0};

    for (auto i_pass : Range(safeties.size())) {
      PrintMessage(4, "GrowthVectorLimiter pass ", i_pass);
      double safety = safeties[i_pass];
      // intersect segment with original surface elements
      LimitOriginalSurface(2.1);
      // intersect prisms with themself
      LimitSelfIntersection(1.3 * safety);
      // intesect segment with prism
      LimitBoundaryLayer(safety);

      for (auto i : Range(3))
        EqualizeLimits(smoothing_factors[i_pass]);

      if (i_pass == safeties.size() - 1)
        FixIntersectingSurfaceTrigs();
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

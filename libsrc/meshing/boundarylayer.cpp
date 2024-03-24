
#include <set>
#include <regex>
#include <mystdlib.h>

#include "global.hpp"
#include "debugging.hpp"

#include "boundarylayer.hpp"
#include "meshfunc.hpp"


namespace netgen
{
  struct Face {
    ArrayMem<Point<3>, 4> p;
    ArrayMem<double, 4> lam;
  };

  struct Intersection_ {
    bool is_intersecting=false;
    double lam0=-1, lam1=-1;
    Point<3> p;
    double bary[3];
    operator bool() const { return is_intersecting; }
  };

std::tuple<int, int> FindCloseVectors(FlatArray<Vec<3>> ns, bool find_max = true) {
  int maxpos1;
  int maxpos2;
    
  double val = find_max ? -1e99 : 1e99;
  for (auto i : Range(ns))
      for (auto j : Range(i + 1, ns.Size()))
        {
          double ip = ns[i] * ns[j];
          if((find_max && (ip > val)) || (!find_max && (ip < val)))
            {
              val = ip;
              maxpos1 = i;
              maxpos2 = j;
            }
        }
  return {maxpos1, maxpos2};
}

Vec<3> CalcGrowthVector(FlatArray<Vec<3>> ns) {
  if(ns.Size() == 0) return {0,0,0};
  if(ns.Size() == 1) return ns[0];
  if(ns.Size() == 2)
    {
      auto gw = ns[0];
      auto n = ns[1];
      auto npn = gw * n;
      auto npnp = gw * gw;
      auto nn = n * n;
      if(fabs(nn-npn*npn/npnp) < 1e-6) 
        return n;
      gw += (nn - npn)/(nn - npn*npn/npnp) * (n - npn/npnp * gw);
      return gw;
    }
  if (ns.Size() == 3)
    {
      DenseMatrix mat(3,3);
      for(auto i : Range(3))
        for(auto j : Range(3))
          mat(i,j) = ns[i][j];

      if(fabs(mat.Det()) > 1e-6) {
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
        Vec<3> growth = 0.;
        for(auto i : Range(ns))
          growth += res[i] * ns[i];
        return growth;
      }
    }
  auto [maxpos1, maxpos2] = FindCloseVectors(ns);
  Array<Vec<3>> new_normals;
  new_normals = ns;
  const auto dot = ns[maxpos1] * ns[maxpos2];
  auto average = 0.5*(ns[maxpos1] + ns[maxpos2]);
  average.Normalize();
  new_normals[maxpos1] = average;
  new_normals.DeleteElement(maxpos2);
  auto gw = CalcGrowthVector(new_normals);

  for(auto n : ns)
    if(n*gw < 0)
      throw Exception("Normals not pointing in same direction as growth vector");
  return gw;
}

SpecialBoundaryPoint :: GrowthGroup :: GrowthGroup(FlatArray<int> faces_, FlatArray<Vec<3>> normals)
{
  faces = faces_;
  growth_vector = CalcGrowthVector(normals);
}

SpecialBoundaryPoint :: SpecialBoundaryPoint( const std::map<int, Vec<3>> & normals )
{
  // find opposing face normals
  Array<Vec<3>> ns;
  Array<int> faces;
  for(auto [face, normal] : normals){
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

  Array<Vec<3>> normals1, normals2;
  for(auto [facei, normali] : normals)
    if(facei != minface1 && facei != minface2)
    {
      g1_faces.Append(facei);
      g2_faces.Append(facei);
    }
  for(auto fi : g1_faces)
    normals1.Append(normals.at(fi));
  for(auto fi : g2_faces)
    normals2.Append(normals.at(fi));
  growth_groups.Append(GrowthGroup(g1_faces, normals1));
  growth_groups.Append(GrowthGroup(g2_faces, normals2));
}


struct GrowthVectorLimiter {
  BoundaryLayerTool & tool;
  const BoundaryLayerParameters & params;
  Mesh & mesh;
  double height;
  FlatArray<double, PointIndex> limits;
  FlatArray<Vec<3>, PointIndex> growthvectors;
  BitArray changed_domains;
  unique_ptr<BoxTree<3>> tree;
  Array<PointIndex, PointIndex> map_from;
  ofstream debug;

  GrowthVectorLimiter( BoundaryLayerTool & tool_) :
    tool(tool_),
    params(tool_.params), mesh(tool_.mesh), height(tool_.total_height), limits(tool_.limits), growthvectors(tool_.growthvectors), map_from(mesh.Points().Size()), debug("debug.txt") {
    changed_domains = params.domains;
    if(!params.outside)
      changed_domains.Invert();

    map_from = PointIndex::INVALID;
    for(auto pi : tool.mapto.Range())
      for (auto pi_to : tool.mapto[pi])
        map_from[pi_to] = pi;
  }

  // Face GetFace( SurfaceElementIndex sei );
  // Face GetMappedFace( SurfaceElementIndex sei );
  // Face GetMappedFace( SurfaceElementIndex sei, int face );

  Point<3> GetPoint(PointIndex pi_to, double shift=1.) {
    if(tool.growth_vector_map.count(pi_to) == 0)
      return mesh[pi_to];

    auto [gw, height] = tool.growth_vector_map[pi_to];
    return mesh[pi_to] + shift * height * (*gw);

    // if(shift == 1.) return mesh[pi_to];
    // auto pi_from = map_from[pi_to];
    // if(!pi_from.IsValid()) return mesh[pi_to];
    // if(shift == 0.) return mesh[pi_from];
    // Point<3> p_from = mesh[pi_from];
    // auto [gw_ptr, h] = tool.growth_vector_map.at(pi_to);
    // Point<3> p_to = mesh[pi_to] + h * limits[pi_to] * (*gw_ptr);
    // return p_from + shift * (p_to - p_from);
  }

  Point<3> GetMappedPoint(PointIndex pi_from, double shift=1.) {
    auto pi_to = tool.mapto[pi_from].Last();
    return GetPoint(pi_to, shift);
  }


  std::array<Point<3>, 2> GetMappedSeg(PointIndex pi_from, double shift=1.) {
    return {mesh[pi_from], GetMappedPoint(pi_from, shift)};
  }

  std::array<Point<3>, 2> GetSeg(PointIndex pi_to, double shift=1.) {
    return {GetPoint(pi_to, 0), GetPoint(pi_to, shift)};
  }

  auto GetTrig(SurfaceElementIndex sei, double shift = 0.0) {
    auto sel = mesh[sei];
    std::array<Point<3>, 3> trig;
    for (auto i : Range(3))
      trig[i] = GetPoint(sel[i], shift);
    return trig;
  }

  auto GetMappedTrig(SurfaceElementIndex sei, double shift = 0.0) {
    auto sel = mesh[sei];
    std::array<Point<3>, 3> trig;
    for (auto i : Range(3))
      trig[i] = GetMappedPoint(sel[i], shift);
    return trig;
  }

  auto GetSideTrig(SurfaceElementIndex sei, int index, double shift = 0.0, bool grow_first = true) {
    // cout << "get side trig " << sei << "/" << mesh.SurfaceElements().Size() << endl;
    auto trig = GetMappedTrig(sei, 0.0);
    // cout << "get side trig " << sei << "/" << mesh.SurfaceElements().Size() << endl;
    auto sel = mesh[sei];
    // cout << "get side trig " << sei << "/" << mesh.SurfaceElements().Size() << endl;
    auto index1 = (index+1)%3;
    if(!grow_first) index1 = (index+2)%3;
    trig[index] = GetMappedPoint(sel[index1], shift);
    // auto [gw_ptr, h] = tool.growth_vector_map.at(sel[index1]);
    // trig[index] += h * shift * (*gw_ptr);
    return trig;
  }

  string DebugPoint(PointIndex pi, string col, bool shift=false ) {
    stringstream ss;
    stringstream pos;
    auto p = mesh[pi];
    auto p1 = mesh[pi]+height*growthvectors[pi];
    pos << "\"position\": [" << p[0] << "," << p[1] << "," << p[2];
    if(shift)
      pos << "," << p1[0] << "," << p1[1] << "," << p1[2];
    pos << "]}";

    ss <<  R"({"name": "point", "color": ")" + col + R"(", "type": "points", )";
    ss << pos.str();
    if(shift) {
      ss <<  R"(, {"name": "line", "color": ")" + col + R"(", "type": "lines", )";
      ss << pos.str();
    }
    return ss.str();
  }

  void DebugOut(PointIndex pi, SurfaceElementIndex sei) {
    auto sel = mesh[sei];
    auto trig = GetTrig(sei, 1.0);
    debug << "[";
    for(auto i : Range(3))
      if(pi!=sel[i])
        debug << DebugPoint(sel[i], "blue", true) << ',';
    debug << DebugPoint(pi, "red", true);
    debug << "]" << endl;
  }

  static constexpr double INTERSECTION_SAFETY = .99;
  void LimitGrowthVector(PointIndex pi_to, SurfaceElementIndex sei, double trig_shift, double seg_shift) {
    auto pi_from = map_from[pi_to];
    if(!pi_from.IsValid()) return;

    if(trig_shift > 0) {
      auto intersection = isIntersectingTrig(pi_from, pi_to, sei, trig_shift);
      if(!intersection) return;
      double dshift = trig_shift;
      while(intersection && dshift > 0.01 && dshift > intersection.lam0) {
        dshift *= 0.9;
        intersection = isIntersectingTrig(pi_from, pi_to, sei, dshift);
      }
      dshift /= 0.9;
      intersection = isIntersectingTrig(pi_from, pi_to, sei, dshift);
      // if(dshift < 1) cout << "final dshift " << dshift << "\t" << intersection.lam0 << endl;
      limits[pi_from] *= intersection.lam0;
      auto sel = mesh[sei];
      for(auto i : Range(3))
        limits[sel[i]] *= dshift;
    }
    else {
      auto seg = GetSeg(pi_to, seg_shift);
      auto trig = GetTrig(sei, 0.0);
      // cout << mesh[sei] << " " << pi_from << endl;
      auto intersection = isIntersectingTrig(seg, trig);
      auto lam = intersection.lam0;
      // cout << "lam " << intersection.is_intersecting <<"\t" << lam << endl;
      if(intersection) {
        // check with original surface elements
        limits[pi_from] = min(limits[pi_from], seg_shift*0.45*INTERSECTION_SAFETY*lam);
        // cout << "set limit " << pi_from << '\t' << limits[pi_from] << endl;
        // if(limits[pi_from] < 0.1) DebugOut(pi_from, sei);
        return;
      }
      // cout << "lam " << lam << endl;
      // if (lam > 0 && lam < 2./INTERSECTION_SAFETY) {
      //   limits[pi_from] = min(limits[pi_from], 0.45*INTERSECTION_SAFETY*lam);
      //   cout << "limit by 2 safety " << limits[pi_from] << endl;
      // }
    }

    // check with shifted surface elements using growthvectors
    // return;
  }

  void LimitSelfIntersection () {
    // check for self-intersection within new elements (prisms/hexes)
    bool found_debug_element = false;
    auto isIntersecting = [&](SurfaceElementIndex sei, double shift) {
      if(sei>=tool.nse) return false; // ignore new surface elements, side trigs are only built from original surface elements
      const auto sel = mesh[sei];
      // cout << "check " << sei << " " << mesh[sei] << " shift = " << shift << endl;
      auto np = sel.GetNP();
      for(auto i : Range(np))
      {
        if(sel[i] > tool.np) return false;
        if(tool.mapto[sel[i]].Size() == 0) return false;
      }
      for(auto i : Range(np)) {
        auto seg = GetMappedSeg(sel[i], shift*limits[sel[i]]);
        for(auto fi : Range(np-2)) {
          for(auto side : {true, false}) {
            auto trig = GetSideTrig(sei, i+fi, 1.0, side);
            // if(sel.PNums().Contains(1)) {
            //   cout << "check " << seg[0] << seg[1] << endl;
            //   cout << "\t" << trig[0] << trig[1] << trig[2] << endl;
            // }
            if(false && sei==104) {
              auto isect = isIntersectingPlane(seg, trig);

              Array<Point<3>> points;
              if(i==0 && side) {
                points.Append(mesh[sel[0]]);
                points.Append(mesh[sel[1]]);
                points.Append(mesh[sel[2]]);
                // debug_gui.DrawTrigs(points, "sel_"+ToString(sei), "green");

                points.Append(GetMappedPoint(sel[0], shift*limits[sel[0]]));
                points.Append(GetMappedPoint(sel[1], shift*limits[sel[1]]));
                points.Append(GetMappedPoint(sel[2], shift*limits[sel[2]]));
                // debug_gui.DrawTrigs(points, "sel_"+ToString(sei), "green");
              }

              string name = "face_"+ToString(i);

              points.SetSize(0);
              points.Append(trig[0]);
              points.Append(trig[1]);
              points.Append(trig[2]);
              // debug_gui.DrawTrigs(points, name, "blue");

              points.SetSize(0);
              points.Append(seg[0]);
              points.Append(seg[1]);
              // debug_gui.DrawLines(points, name, "red");

              if(isect) {
              points.SetSize(0);
              points.Append(isect.p);
              // debug_gui.DrawPoints(points, name, "red");

                              


              }
              // debug_gui.DrawLines(points, "points_sel_"+ToString(sei)+"_"+ToString(i)+"_"+ToString(side), "blue");


              found_debug_element = true;

            }
            if(auto isect = isIntersectingPlane(seg, trig); isect) {
              // if(sel.PNums().Contains(1))
              //   cout << "found intersection " << shift << "\t"<< sei << " " << mesh[sei] << endl;
              // cout << isect.lam0 << ',' << isect.p << endl;
              //   cout << "\t" << seg[0] << seg[1] << endl;
              //   cout << "\t" << trig[0] << trig[1] << trig[2] << endl;
              return true;
            }
          }
        }
      }

      if(false && sei==104) {
        // int i_;
        // cin >> i_;
        // debug_gui.ResetObjects();
      }
      return false;
    };

    auto equalizeLimits = [&](SurfaceElementIndex sei) {
        const auto sel = mesh[sei];
        auto np = sel.GetNP();
        double max_limit = 0;
        double min_limit = 1e99;
        for(auto i : Range(np)) {
          max_limit = max(max_limit, limits[sel[i]]);
          min_limit = min(min_limit, limits[sel[i]]);
        }
        // equalize
        if(max_limit/min_limit > 1.2) {
          max_limit = min_limit * 1.2;
          for(auto i : Range(np))
            limits[sel[i]] = min(limits[sel[i]], max_limit);
        }
    };

    for(SurfaceElementIndex sei : mesh.SurfaceElements().Range()) {
      auto sel = mesh[sei];
      const auto& fd = mesh.GetFaceDescriptor(sel.GetIndex());
      if(sei >= tool.nse) continue;
      if(sel.GetNP()==4) continue;
      // if(sei >= tool.nse || (!changed_domains.Test(fd.DomainIn()) &&
      //   !changed_domains.Test(fd.DomainOut())))
      //   continue;

      auto np = sel.GetNP();
      // ArrayMem<double, 4> ori_limits;
      // ori_limits.SetSize(np);
      // for(auto i : Range(np))
      //     ori_limits[i] = limits[sel[i]];

      // equalizeLimits(sei);
      
      double shift = 1.0;
      double safety = 1.1;
      const double step_factor = 0.9;
      while(shift>0.01 && isIntersecting(sei, shift*safety)) {
        shift *= step_factor;
        // cout << "reduce " << sel << "\tshift = " << shift << endl;
        double max_limit = 0;
        // cout << "\tlimits before " << limits[sel[0]] << '\t' <<  limits[sel[1]] << '\t' <<  limits[sel[2]] << endl;
        for(auto i : Range(np))
          max_limit = max(max_limit, limits[sel[i]]);
        for(auto i : Range(np))
          if(max_limit == limits[sel[i]])
            limits[sel[i]] *= step_factor;
        // cout << "\tlimits after  " << limits[sel[0]] << '\t' <<  limits[sel[1]] << '\t' <<  limits[sel[2]] << endl;
        if(max_limit < 0.01) {
          // if(sel.PNums().Contains(1))
            // cout << "self intersection" << endl;
          break;
        }
      }
      // if(shift < 1) {
      //   if(shift < 0.3) {
      //     cout << "self intersection " << sel << "\t" << shift << endl;
      //     cout << "\t" << limits[sel[0]] << '\t' <<  limits[sel[1]] << '\t' <<  limits[sel[2]] << endl;
      //   }
      //   for(auto pi : sel.PNums())
      //     limits[pi] *= INTERSECTION_SAFETY*shift;
      // }

      // cout << 
      // auto np = sel.GetNP();
      // // check if a new edge intesects the plane of any opposing face
      // for(auto i : Range(np)) {
      //   auto seg = GetSeg(sel[i]);
      //   for(auto fi : Range(np-2)) {
      //     auto trig = GetSideTrig(sei, i+fi, 1.0);
      //     auto intersection = isIntersectingPlane(seg, trig);
      //     if(intersection) {
      //       if(intersection.lam0 < 0.2) {
      //          DebugOut(sel[i], sei);
      //       }
      //       limits[sel[i]] *= INTERSECTION_SAFETY*intersection.lam0;
      //     }
      //   }
      // }
    }
  }

  // checks if a segment is intersecting a plane, spanned by three points, lam will be set s.t. p_intersect = seg[0] + lam * (seg[1]-seg[0])
  Intersection_ isIntersectingPlane ( std::array<Point<3>, 2> seg, std::array<Point<3>, 3> trig )
  {
    auto t1 = trig[1]-trig[0];
    auto t2 = trig[2]-trig[0];
    auto n = Cross(trig[1]-trig[0], trig[2]-trig[0]);
    auto v0n = (seg[0]-trig[0])*n;
    auto v1n = (seg[1]-trig[0])*n;

    Intersection_ intersection;
    intersection.lam0 = -v0n/(v1n-v0n);
    // cout << "lam0 " << intersection.lam0 << ", " << v0n << ", " << v1n << endl;
    intersection.p = seg[0] + intersection.lam0*(seg[1]-seg[0]);
    intersection.is_intersecting = (v0n * v1n < 0) && (intersection.lam0 > -1e-8) && (intersection.lam0<1+1e-8);

    return intersection;
  }

  Intersection_ isIntersectingPlane ( PointIndex pi, PointIndex pi_to, SurfaceElementIndex sei, double shift = 0.0 )
  {
    return isIntersectingPlane(GetSeg(pi, pi_to), GetTrig(sei, shift));
  }

  Intersection_ isIntersectingTrig ( std::array<Point<3>, 2> seg, std::array<Point<3>, 3> trig )
  {
    auto intersection = isIntersectingPlane(seg, trig);
    if(!intersection)
      return intersection;

    auto p = seg[0] + intersection.lam0*(seg[1]-seg[0]) - trig[0];

    Vec3d col1 = trig[1]-trig[0];
    Vec3d col2 = trig[2]-trig[0];
    Vec3d col3 = Cross(col1,col2);
    Vec3d rhs = p;
    Vec3d bary;
    SolveLinearSystem (col1, col2, col3, rhs, bary);

    intersection.lam1 = 0;
    double eps = 0;
    if (bary.X() >= -eps && bary.Y() >= -eps && 
      bary.X() + bary.Y() <= 1+eps) 
    {
      intersection.bary[0] = bary.X();
      intersection.bary[1] = bary.Y();
      intersection.bary[2] = 1.0-bary.X()-bary.Y();
      // cout << "\tFOUND INTERSECTION " << intersection.bary[0] << ", " << intersection.bary[1] << ", " << intersection.bary[2] << endl;
    }
    else
      intersection.is_intersecting = false;
    // cout << "return intersection " << intersection.is_intersecting << endl;
    return intersection;
  }

  Intersection_ isIntersectingTrig ( PointIndex pi_from, PointIndex pi_to, SurfaceElementIndex sei, double shift=0.0)
  {
    return isIntersectingTrig(GetSeg(pi_from, pi_to), GetTrig(sei, shift));
  }

  void BuildSearchTree(double trig_shift) {
    Box<3> bbox(Box<3>::EMPTY_BOX);
    for(PointIndex pi : mesh.Points().Range())
    {
      bbox.Add(mesh[pi]);
      bbox.Add(GetPoint(pi, 1.1));
      // if(tool.mapto[pi].Size() >0)
        // bbox.Add(mesh[tool.mapto[pi].Last()]);
    }

    tree = make_unique<BoxTree<3>>(bbox);

    for(auto sei : mesh.SurfaceElements().Range())
    {
      // cout << "handle sei " << sei << ", old nse " << tool.nse << endl;
      const auto & sel = mesh[sei];
      auto sel_index = mesh[sei].GetIndex();
      // if(sel_index)
      // {
      //   cout << "index " << sel_index << endl;
      //   const auto& fd = mesh.GetFaceDescriptor(sel_index);
      //   if( !changed_domains.Test(fd.DomainIn()) &&
      //     !changed_domains.Test(fd.DomainOut()))
      //     continue;
      // }

      Box<3> box(Box<3>::EMPTY_BOX);
      for(auto p : GetTrig(sei, trig_shift))
        box.Add(p);
      // for(auto pi : sel.PNums())
      // box.Add(mesh[pi]);
      // also add moved points to bounding box
      // for(auto p : GetTrig(sei, trig_shift))
          // box.Add(mesh[pi]+trig_shift*limits[pi]*height*growthvectors[pi]);
      tree->Insert(box, sei);
    }
  }

  template<typename TFunc>
  void FindTreeIntersections( double trig_shift, double seg_shift, TFunc f) {
    BuildSearchTree(trig_shift);
    auto np_new = mesh.Points().Size();
    // cout << "np_new " << np_new << endl;
    // cout << "too.np " << tool.np << endl;
    int counter = 0;
    for(auto i : IntRange(tool.np, np_new)) {
      // cout << "handle point " << i << endl;
      PointIndex pi_to = i+PointIndex::BASE;
      PointIndex pi_from = map_from[pi_to];
      if(!pi_from.IsValid())
        throw Exception("Point not mapped");

      // if(mesh[pi_to].Type() == INNERPOINT)
        // continue;
      // if(growthvectors[pi_to].Length2() == 0.0)
        // continue;
      Box<3> box(Box<3>::EMPTY_BOX);
      auto seg = GetSeg(pi_to, seg_shift);

      box.Add(GetPoint(pi_to, 0));
      box.Add(GetPoint(pi_to, limits[pi_from]));
      tree->GetFirstIntersecting(box.PMin(), box.PMax(), [&](SurfaceElementIndex sei) {
        const auto & sel = mesh[sei];
        // cout << "got intersecting " << sei << endl;
        if(sel.PNums().Contains(pi_from))
          return false;
        counter++;
        f(pi_to, sei);
        return false;
      });
    }
    cout << "intersections counter " << counter << endl;
  }
};


  // Face GrowthVectorLimiter :: GetMappedFace( SurfaceElementIndex sei, int face_nr )
  // {
  //     if(face_nr == -1) return GetFace(sei);
  //     if(face_nr == -2) return GetMappedFace(sei);
  //     const auto & sel = mesh[sei];
  //     auto np = sel.GetNP();
  //     Face face;
  //     face.p.SetSize(4);
  //     face.lam.SetSize(4);
  //     face.lam = 0;
  //     auto pi0 = sel[face_nr % np];
  //     auto pi1 = sel[(face_nr+1) % np];
  //     face.p[0] = face.p[3] = mesh[pi0];
  //     face.p[1] = face.p[2] = mesh[pi1];
  //     face.p[3] += height * limits[pi0]*growthvectors[pi0];
  //     face.p[2] += height * limits[pi1]*growthvectors[pi1];
  //     face.lam[2] = 1.0;
  //     face.lam[3] = 1.0;
  //     return face;
  // }

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
      if(pts.Size() != 2) {
        cout << "getEdgeTangent pi = " << pi << ", edgenr = " << edgenr << endl;
        for(auto segi : topo.GetVertexSegments(pi))
          cout  << mesh[segi] << endl;
        throw Exception("Something went wrong in getEdgeTangent!");
      }
      tangent = mesh[pts[1]] - mesh[pts[0]];
      return tangent.Normalize();
  }

  void BoundaryLayerTool :: LimitGrowthVectorLengths()
  {
    static Timer tall("BoundaryLayerTool::LimitGrowthVectorLengths"); RegionTimer rtall(tall);
    mesh.Save("mesh_before_limit.vol");
  
    limits.SetSize(mesh.Points().Size());
    limits = 1.0;

    GrowthVectorLimiter limiter(*this); //, mesh, params, limits, growthvectors, total_height);

    // limit to not intersect with other (original) surface elements
    double trig_shift = 0;
    double seg_shift = 2.1;
    limiter.FindTreeIntersections(trig_shift, seg_shift, [&] (PointIndex pi_to, SurfaceElementIndex sei) {
    // cout << "found intersection 1 " << pi_to << ", " << sei << endl;
      if(sei >= nse) return; // ignore new surface elements in first pass
      limiter.LimitGrowthVector(pi_to, sei, trig_shift, seg_shift);
    });

    limiter.LimitSelfIntersection();

    // for(auto i : Range(growthvectors))
    //   growthvectors[i] *= limits[i];
    // limits = 1.0;

    // // now limit again with shifted surface elements
    trig_shift = 1.1;
    seg_shift = 1.1;
    limiter.FindTreeIntersections(trig_shift, seg_shift, [&] (PointIndex pi_to, SurfaceElementIndex sei) {
      // cout << "Should not have intersection with original surface element anymore" << endl;
      // cout << "found intersection 2 " << pi_to << ", " << sei << endl;
      limiter.LimitGrowthVector(pi_to, sei, trig_shift, seg_shift);
    });

  // for (auto i : Range(limits))
  // if(limits[i] < 1.0)
  //   cout << i << ": " << limits[i] << endl;

    // for (auto [pi_to, data] : growth_vector_map) {
    //   auto pi_from = limiter.map_from[pi_to];
    //   if(pi_from.IsValid())
    //     limits[pi_from] = min(limits[pi_from], limits[pi_to]);
    // }

  // cout << "limits " << endl << limits << endl;
    for(auto i : Range(growthvectors))
      growthvectors[i] *= limits[i];

}

  // void BoundaryLayerTool :: LimitGrowthVectorLengths()
  // {
  //   static Timer tall("BoundaryLayerTool::LimitGrowthVectorLengths"); RegionTimer rtall(tall);

  //   limits.SetSize(np);
  //   limits = 1.0;
  //  
  //   // Function to calculate the dot product of two 3D vectors
  //   // Is there netgen native function for this?
  //   const auto Dot = [](Vec<3> a, Vec<3> b) {
  //       return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  //   };

  //   auto parallel_limiter = [&](PointIndex pi1, PointIndex pi2, SurfaceElementIndex si) {
  //       MeshPoint& a_base = mesh[pi1];
  //       MeshPoint& b_base = mesh[pi2];
  //       MeshPoint a_end = mesh[pi1] + height * limits[pi1] * growthvectors[pi1];
  //       MeshPoint b_end = mesh[pi2] + height * limits[pi2] * growthvectors[pi2];

  //       double ab_base = (b_base - a_base).Length();
  //       Vec<3> a_vec = (a_end - a_base);
  //       Vec<3> b_vec = (b_end - b_base);

  //       // Calculate parallel projections
  //       Vec<3> ab_base_norm = (b_base - a_base).Normalize();
  //       double a_vec_x = Dot(a_vec, ab_base_norm);
  //       double b_vec_x = Dot(b_vec, -ab_base_norm);
  //       // double ratio_parallel = (a_vec_x + b_vec_x) / ab_base;

  //       double PARALLEL_RATIO_LIMIT = 0.85;
  //       if (ratio_parallel > PARALLEL_RATIO_LIMIT) {
  //           // Adjust limits, vectors, and projections if parallel ratio exceeds the limit
  //           double corrector = PARALLEL_RATIO_LIMIT / ratio_parallel;
  //           limits[pi1] *= corrector;
  //           limits[pi2] *= corrector;
  //       }
  //   };      
  //   
  //   auto perpendicular_limiter = [&](PointIndex pi1, PointIndex pi2, SurfaceElementIndex si) {
  //       // this part is same as in parallel limiter, but note that limits contents are already changed
  //       MeshPoint& a_base = mesh[pi1];
  //       MeshPoint& b_base = mesh[pi2];
  //       MeshPoint a_end = mesh[pi1] + height * limits[pi1] * growthvectors[pi1];
  //       MeshPoint b_end = mesh[pi2] + height * limits[pi2] * growthvectors[pi2];

  //       double ab_base = (b_base - a_base).Length();
  //       Vec<3> a_vec = (a_end - a_base);
  //       Vec<3> b_vec = (b_end - b_base);

  //       // Calculate parallel projections
  //       Vec<3> ab_base_norm = (b_base - a_base).Normalize();
  //       double a_vec_x = Dot(a_vec, ab_base_norm);
  //       double b_vec_x = Dot(b_vec, -ab_base_norm);
  //       double ratio_parallel = (a_vec_x + b_vec_x) / ab_base;

  //       // Calculate surface normal at point si
  //       Vec<3> surface_normal = getNormal(mesh[si]);

  //       double a_vec_y = abs(Dot(a_vec, surface_normal));
  //       double b_vec_y = abs(Dot(b_vec, surface_normal));
  //       double diff_perpendicular = abs(a_vec_y - b_vec_y);
  //       double tan_alpha = diff_perpendicular / (ab_base - a_vec_x - b_vec_x);

  //       double TAN_ALPHA_LIMIT = 0.36397; // Approximately 20 degrees in radians
  //       if (tan_alpha > TAN_ALPHA_LIMIT) {
  //           if (a_vec_y > b_vec_y) {
  //               double correction = (TAN_ALPHA_LIMIT / tan_alpha * diff_perpendicular + b_vec_y) / a_vec_y;
  //               limits[pi1] *= correction;
  //           }
  //           else {
  //               double correction = (TAN_ALPHA_LIMIT / tan_alpha * diff_perpendicular + a_vec_y) / b_vec_y;
  //               limits[pi2] *= correction;
  //           }
  //       }
  //   };

  //   auto neighbour_limiter = [&](PointIndex pi1, PointIndex pi2, SurfaceElementIndex si) {
  //       parallel_limiter(pi1, pi2, si);
  //       perpendicular_limiter(pi1, pi2, si);
  //   };
  //   
  //   auto modifiedsmooth = [&](size_t nsteps) {
  //       for ([[maybe_unused]] auto i : Range(nsteps))
  //           for (SurfaceElementIndex sei : mesh.SurfaceElements().Range())
  //           {
  //           // assuming triangle
  //           neighbour_limiter(mesh[sei].PNum(1), mesh[sei].PNum(2), sei);
  //           neighbour_limiter(mesh[sei].PNum(2), mesh[sei].PNum(3), sei);
  //           neighbour_limiter(mesh[sei].PNum(3), mesh[sei].PNum(1), sei);  
  //       }
  //   };

  //   auto smooth = [&] (size_t nsteps) {
  //       for([[maybe_unused]] auto i : Range(nsteps))
  //           for(const auto & sel : mesh.SurfaceElements())
  //           {
  //               double min_limit = 999;
  //               for(auto pi : sel.PNums())
  //                   min_limit = min(min_limit, limits[pi]);
  //               for(auto pi : sel.PNums())
  //                   limits[pi] = min(limits[pi], 1.4*min_limit);
  //           }
  //   };

  //   // check for self-intersection within new elements (prisms/hexes)
  //   auto self_intersection = [&] () {
  //       for(SurfaceElementIndex sei : mesh.SurfaceElements().Range())
  //         {
  //           auto facei = mesh[sei].GetIndex();
  //           if(facei < nfd_old && !params.surfid.Contains(facei))
  //               continue;

  //           auto sel = mesh[sei];
  //           auto np = sel.GetNP();
  //           // check if a new edge intesects the plane of any opposing face
  //           double lam0, lam1;
  //           for(auto i : Range(np))
  //               for(auto fi : Range(np-2))
  //                   if(isIntersectingPlane(GetMappedSeg(sel[i]), GetMappedFace(sei, i+fi+1), lam0, lam1))
  //                       if(lam0 < 1.0)
  //                           limits[sel[i]] *= lam0;
  //         }
  //   };

  //   // first step: intersect with other surface elements that are boundary of domain the layer is grown into
  //   // second (and subsequent) steps: intersect with other boundary layers, allow restriction by 20% in each step
  //   auto changed_domains = domains;
  //     if(!params.outside)
  //         changed_domains.Invert();

  //   bool limit_reached = true;
  //   double lam_lower_limit = 1.0;
  //   int step = 0;

  //   while(step<2)
  //   {
  //       Array<double, PointIndex> new_limits;
  //       new_limits.SetSize(np);
  //       new_limits = 1.0;

  //   // if(step==1) break;

  //       if(step>1)
  //           lam_lower_limit *= 0.8;
  //       limit_reached = false;

  //       // build search tree with all surface elements (bounding box of a surface element also covers the generated boundary layer)

  //       for(auto pi : mesh.Points().Range())
  //       {
  //           if(mesh[pi].Type() == INNERPOINT)
  //               continue;
  //           if(growthvectors[pi].Length2() == 0.0)
  //               continue;
  //           const auto debug = false;
  //           Box<3> box(Box<3>::EMPTY_BOX);
  //           auto seg = GetMappedSeg(pi);
  //           box.Add(seg[0]);
  //           box.Add(seg[1]);
  //           double lam = 1.0;
  //           tree.GetFirstIntersecting(box.PMin(), box.PMax(), [&](SurfaceElementIndex sei)
  //             {
  //               const auto & sel = mesh[sei];
  //               if(sel.PNums().Contains(pi))
  //                   return false;
  //                                     cout << "LIMIT STEP " << step << endl;
  //               if(step == 0)
  //                 LimitGrowthVector(pi, sei, new_limits, {-1, 0});
  //               else
  //                 LimitGrowthVector(pi, sei, new_limits, {-2, -1});
  //               // auto face = GetMappedFace(sei, -2);
  //               // double lam0_ = 999;
  //               // double lam1_ = 999;
  //               // bool is_bl_sel = params.surfid.Contains(sel.GetIndex());

  //               // if (step == 0)
  //               // {
  //               //     face = GetFace(sei);
  //               //     if (isIntersectingFace(seg, face, lam0_, lam1_))
  //               //     {
  //               //         if(lam0_ < lam) {
  //               //         if(debug) cout << "intersecting face " << sei << endl;
  //               //         if(debug) cout << "\t" << lam0_ << endl;
  //               //                       }
  //               //         // if (is_bl_sel)
  //               //         //     lam_ *= params.limit_safety;
  //               //         lam = min(lam, lam0_);
  //               //     }
  //               // }

  //               // if(step==1)
  //               // {
  //               //     if(isIntersectingFace(seg, face, lam0_, lam1_))
  //               //     {
  //               //         // if(is_bl_sel) // allow only half the distance if the opposing surface element has a boundary layer too
  //               //         //     lam_ *= params.limit_safety;
  //               //         lam = min(lam, lam0_);
  //               //     }
  //               // }
  //               // // if the opposing surface element has a boundary layer, we need to additionally intersect with the new faces 
  //               // if(step>1 && is_bl_sel)
  //               // {
  //               //     for(auto facei : Range(-1, sel.GetNP()))
  //               //     {
  //               //         auto face = GetMappedFace(sei, facei);
  //               //         if(isIntersectingFace(seg, face, lam0_, lam1_)) // && lam_ > other_limit)
  //               //         {
  //               //             lam = min(lam, lam0_);
  //               //         }
  //               //     }
  //               // }
  //               return false;
  //             });
  //           // if(lam<1)
  //           // {
  //           //     if(lam<lam_lower_limit && step>1)
  //           //     {
  //           //         limit_reached = true;
  //           //         lam = lam_lower_limit;
  //           //     }
  //           // }

  //           // new_limits[pi] = min(limits[pi], lam* limits[pi]);
  //       }
  //       // cout << "new limits " << endl;
  //       // cout << new_limits << endl;
  //       // for(auto pi : mesh.Points().Range())
  //       //   if(growthvectors[pi].Length2())
  //       //   {
  //       //     if(new_limits[pi] < 0.001) {
  //       //       cout << pi << " " << new_limits[pi] << endl;
  //       //       new_limits[pi] = 0.001;
  //       //     }
  //       //   }

  //       if(step == 0) {
  //     cout << "limit with 0.9" << endl;
  //         for(auto & v : new_limits)
  //           v *= 0.9;
  //       }
  //       for(auto i : Range(limits))
  //         limits[i] *= new_limits[i];
  //   cout << "new limits " << endl << limits << endl;
  //       // for(auto pi : mesh.Points().Range())
  //       //   if(growthvectors[pi].Length2())
  //       //     cout << "apply limit " << pi << " \t" << limits[pi] << endl;


  //       // break;
  //       // if (step > 0)
  //       //    modifiedsmooth(1);
  //       step++;
  //   }

  //   // self_intersection();
  //   // modifiedsmooth(1);

  // cout << "final limits " << limits << endl;
  //   for(auto pi : Range(growthvectors))
  //       growthvectors[pi] *= limits[pi];

  // }


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
      // auto& topo = mesh.GetTopology();

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

  void MergeAndAddSegments( Mesh & mesh, FlatArray<Segment> segments, FlatArray<Segment> new_segments)
  {
      INDEX_2_HASHTABLE<bool> already_added( segments.Size() + 2*new_segments.Size() );

      mesh.LineSegments().SetSize0();

      auto addSegment = [&] (const auto & seg) {
          INDEX_2 i2 (seg[0], seg[1]);
          i2.Sort();
          if(!already_added.Used(i2)) {
              mesh.AddSegment(seg);
              already_added.Set(i2, true);
          }
      };

      for(const auto & seg : segments)
        addSegment(seg);

      for(const auto & seg : new_segments)
        addSegment(seg);
  }

  void BoundaryLayerTool :: InterpolateSurfaceGrowthVectors()
  {
    static Timer tall("InterpolateSurfaceGrowthVectors"); RegionTimer rtall(tall);
    static Timer tsmooth("InterpolateSurfaceGrowthVectors-Smoothing");
    auto np_old = this->np;
    auto np = mesh.GetNP();
    BitArray is_point_on_bl_surface(np+1);
    is_point_on_bl_surface.Clear();
    BitArray is_point_on_other_surface(np+1);
    is_point_on_other_surface.Clear();

    auto getGW = [&] (PointIndex pi) -> Vec<3> {
        if(pi-PointIndex::BASE<np_old && mapto[pi].Size())
          return {.0};
        if(growth_vector_map.count(pi)==0)
          growth_vector_map[pi] = { &growthvectors[pi], total_height };
        auto [gw, height] = growth_vector_map[pi];
        return height * (*gw);
    };
    auto setGW = [&] (PointIndex pi, Vec<3> vec) {
        if(growth_vector_map.count(pi)==0)
          growth_vector_map[pi] = { &growthvectors[pi], total_height };
        auto [gw, height] = growth_vector_map[pi];
        *gw = 1.0/height*vec;
    };
    // auto getGW = [&] (PointIndex pi) -> Vec<3>& {
    //     if(pi<=np_old)
    //      return growthvectors[pi];
    //    return *get<0>(growth_vector_map[pi]);
    // };

    Array<Vec<3>, PointIndex> normals(np);
    for(auto pi = np_old; pi<np; pi++){
        normals[pi+PointIndex::BASE] = getGW(pi+PointIndex::BASE);
    }

    auto hasMoved =  [&](PointIndex pi) {
      return (pi-PointIndex::BASE>=np_old) || mapto[pi].Size() > 0 || special_boundary_points.count(pi);
    };


    ParallelForRange( Range(nse), [&] ( auto myrange )
      {
        for(SurfaceElementIndex sei : myrange)
          {
            auto facei = mesh[sei].GetIndex();
            if(facei < nfd_old && !params.surfid.Contains(facei))
              {
                for(auto pi : mesh[sei].PNums()) 
                  if(auto & p = mesh[pi]; p.Type() == SURFACEPOINT)
                    is_point_on_other_surface.SetBitAtomic(pi);
              }
            else
              {
                for(auto pi : mesh[sei].PNums())
                  if(pi >= np_old + PointIndex::BASE && mesh[pi].Type() == SURFACEPOINT)
                    is_point_on_bl_surface.SetBitAtomic(pi);
              }
          }
      });

    Array<PointIndex> points;
    Array<bool> has_moved_points(max_edge_nr+1);
    has_moved_points = false;
    std::set<PointIndex> moved_edge_points;

    for(auto seg : segments) {
      if(hasMoved(seg[0])!= hasMoved(seg[1]))
        has_moved_points[seg.edgenr] = true;
  }


    for(auto seg : segments)
      if(has_moved_points[seg.edgenr])
        for(auto pi : seg.PNums())
          if(mesh[pi].Type() == EDGEPOINT)
            moved_edge_points.insert(pi);

    for(auto pi: moved_edge_points)
      points.Append(pi);

    for(PointIndex pi : mesh.Points().Range())
      if(is_point_on_bl_surface[pi] || is_point_on_other_surface[pi])
        points.Append(pi);


    auto p2sel = mesh.CreatePoint2SurfaceElementTable();
    // smooth tangential part of growth vectors from edges to surface elements
    RegionTimer rtsmooth(tsmooth);
    for([[maybe_unused]] auto i : Range(10))
    {
        for(auto pi : points)
        {
            auto sels = p2sel[pi];
            Vec<3> new_gw = getGW(pi);
            new_gw = 0.;
            // int cnt = 1;
            std::set<PointIndex> suround;
            suround.insert(pi);
            double total_weight = 0;
            // auto normal = normals[pi];
            for(auto sei: sels)
            {
                const auto & sel = mesh[sei];
                for(auto pi1 : sel.PNums())
                    if(suround.count(pi1)==0)
                    {
                        suround.insert(pi1);
                        auto gw_other = getGW(pi1);
                        auto normal_other = getNormal(mesh[sei]);
                        auto tangent_part = gw_other - (gw_other*normal_other)*normal_other;
                        double weight = 1.0;
                        // cout << "tangent part " << pi1 << tangent_part << endl;
          
                        if(is_point_on_bl_surface[pi]) {
                          if(mesh[pi1].Type() == FIXEDPOINT)
                            weight *= 1.0; //13-i;
                          else
                            weight = 1.0;
                          new_gw += weight * tangent_part;
                        }
                        else {
                          new_gw += weight * gw_other;
                        }
                        total_weight += weight;
                    }
            }
            // total_weight += suround.size();

            setGW(pi, 1.0/total_weight * new_gw);
            // cout << "average " << pi << " " << getGW(pi) << endl;
        }
    }

    // for(auto pi : points) 
    //     getGW(pi) += normals[pi];
    // for(auto pi : mesh.Points().Range()) 
    //   cout << "point " << pi << " has type " << (int)(mesh[pi].Type()) << endl;
  }


  BoundaryLayerTool::BoundaryLayerTool(Mesh & mesh_, const BoundaryLayerParameters & params_)
      : mesh(mesh_), topo(mesh_.GetTopology()), params(params_)
  {
    static Timer timer("BoundaryLayerTool::ctor");
    RegionTimer regt(timer);

    //for(auto & seg : mesh.LineSegments())
        //seg.edgenr = seg.epgeominfo[1].edgenr;

    total_height = 0.0;
    for (auto h : params.heights)
      total_height += h;

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
            // curving of surfaces with boundary layers will often
            // result in pushed through elements, since we do not (yet)
            // curvature through layers.
            // Therefore we disable curving for these surfaces.
            if(!params.keep_surfaceindex)
              mesh.GetFaceDescriptor(i).SetSurfNr(-1);
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
        // bool point_fixed = false;
        for(auto pi : sel.PNums())
          {
            if(growthvectors[pi].Length() > 0)
              point_moved = true;
            /*
            else
              point_fixed = true;
            */
          }
        if(point_moved && !moved_surfaces.Test(facei))
          {
            int new_si = mesh.GetNFD()+1;
            const auto& fd = mesh.GetFaceDescriptor(facei);
            // auto isIn = domains.Test(fd.DomainIn());
            // auto isOut = domains.Test(fd.DomainOut());
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
        ArrayMem<Vec<3>, 5> ns;
        for (auto &[facei, n] : normals) {
          ns.Append(n);
        }

        try {
          growthvectors[pi] = CalcGrowthVector(ns);
        }
        catch(const Exception & e) {
          cout << "caught exception for point " << pi << ":\t" << e.what() << endl;
          special_boundary_points.emplace(pi, normals);
          growthvectors[pi] = special_boundary_points[pi].growth_groups[0].growth_vector;
      }

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
                if(moved_surfaces.Test(segj.si)) {
                  type = 0;
                  moved_segs.Append(sj);
                }
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
    // mesh.Save("interpolate.vol");
  // cout << "new number of line segments " << mesh.LineSegments().Size() << endl;
    int new_max_edge_nr = max_edge_nr;
    for(const auto& seg : segments)
      if(seg.edgenr > new_max_edge_nr)
        new_max_edge_nr = seg.edgenr;
    for(const auto& seg : new_segments)
      if(seg.edgenr > new_max_edge_nr)
        new_max_edge_nr = seg.edgenr;

    auto getGW = [&] (PointIndex pi) -> Vec<3> {
        if(growth_vector_map.count(pi)==0)
          growth_vector_map[pi] = { &growthvectors[pi], total_height };
        auto [gw, height] = growth_vector_map[pi];
        return height * (*gw);
    };
    auto addGW = [&] (PointIndex pi, Vec<3> vec) {
        if(growth_vector_map.count(pi)==0)
          growth_vector_map[pi] = { &growthvectors[pi], total_height };
        auto [gw, height] = growth_vector_map[pi];
        *gw += 1.0/height*vec;
    };
  // cout << "edge range " << max_edge_nr << ", " << new_max_edge_nr << endl;

    // interpolate tangential component of growth vector along edge
    for(auto edgenr : Range(max_edge_nr+1, new_max_edge_nr))
      {
    // cout << "SEARCH EDGE " << edgenr +1 << endl;
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
            if(segs.Size() == 1)
              return true;
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
                if(getGW(seg[0]).Length2() != 0 ||
                   getGW(seg[1]).Length2() != 0)
                  any_grows = true;
                if(points.Size() == 0 && (is_end_point(seg[0]) || is_end_point(seg[1])))
                  {
                    PointIndex seg0 = seg[0], seg1 = seg[1];
                    if(is_end_point(seg[1]))
                      Swap(seg0, seg1);
                    points.Append(seg0);
                    points.Append(seg1);
                    edge_len += (mesh[seg[1]] - mesh[seg[0]]).Length();
                  }
              }
          }

        if(!any_grows) {
          // cout << "skip edge " << edgenr+1 << endl;
          continue;
        }

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

        if(getGW(points[0]).Length2() == 0 &&
           getGW(points.Last()).Length2() == 0)
          continue;
        // cout << "Points to average " << endl << points << endl;

        // tangential part of growth vectors
        auto t1 = (mesh[points[1]]-mesh[points[0]]).Normalize();
        auto gt1 = getGW(points[0]) * t1 * t1;
        auto t2 = (mesh[points.Last()]-mesh[points[points.Size()-2]]).Normalize();
        auto gt2 = getGW(points.Last()) * t2 * t2;

        // if(!is_edge_moved[edgenr+1])
          // {
          //   if(getGW(points[0]) * (mesh[points[1]] - mesh[points[0]]) < 0)
          //     gt1 = 0.;
          //   if(getGW(points.Last()) * (mesh[points[points.Size()-2]] - mesh[points.Last()]) < 0)
          //     gt2 = 0.;
          // }

        // cout << "edgenr " << edgenr << endl;
        // cout << "points " << endl << points << endl;

        double len = 0.;
        for(auto i : IntRange(1, points.Size()-1))
          {
            auto pi = points[i];
            len += (mesh[pi] - mesh[points[i-1]]).Length();
            auto t = getEdgeTangent(pi, edgenr);
            auto lam = len/edge_len;
            auto interpol = (1-lam) * (gt1 * t) * t + lam * (gt2 * t) * t;
      // if(pi==89) {
      //   cout << "points " << points << endl;
      //   cout << "INTERPOL" << len << ',' << t << ',' << lam << ',' << interpol << endl;
      //   cout << gt1 << endl;
      //   cout << gt2 << endl;
      //   cout << getGW(pi) << endl;

      // }
            addGW(pi, interpol);
          }
      }

      InterpolateSurfaceGrowthVectors();
  }

  void BoundaryLayerTool :: InsertNewElements( FlatArray<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap, const BitArray & in_surface_direction )
  {
    static Timer timer("BoundaryLayerTool::InsertNewElements"); RegionTimer rt(timer);
    mapto.SetSize(0);
    mapto.SetSize(np);

    auto changed_domains = domains;
    if(!params.outside)
        changed_domains.Invert();


    auto & identifications = mesh.GetIdentifications();
    const int identnr = identifications.GetNr("boundarylayer");

    auto add_points = [&](PointIndex pi, Vec<3> & growth_vector, Array<PointIndex> & new_points)
    {
      Point<3> p = mesh[pi];
      PointIndex pi_last = pi;
      for(auto i : Range(params.heights))
        {
          // p += params.heights[i] * growth_vector;
          auto pi_new = mesh.AddPoint(p);
          new_points.Append(pi_new);
          growth_vector_map[pi_new] = { &growth_vector, params.heights[i] };
          if(special_boundary_points.count(pi) > 0)
            mesh.AddLockedPoint(pi_new);
          pi_last = pi_new;
        }
    };

    // insert new points
    for (PointIndex pi = 1; pi <= np; pi++) {
      if (growthvectors[pi].Length2() != 0) {
        if(special_boundary_points.count(pi))
        {
          for(auto & group : special_boundary_points[pi].growth_groups)
            add_points(pi, group.growth_vector, group.new_points);
        }
        else
          add_points(pi, growthvectors[pi], mapto[pi]);
       }
    }

    // get point from mapto (or the group if point is mapped to multiple new points)
    // layer = -1 means last point (top of boundary layer)
    auto newPoint = [&](PointIndex pi, int layer = -1, int group = 0) {
      if(layer == -1) layer = params.heights.Size()-1;
      if(special_boundary_points.count(pi))
          return special_boundary_points[pi].growth_groups[group].new_points[layer];
      else
          return mapto[pi][layer];
    };

    auto hasMoved =  [&](PointIndex pi) {
      return mapto[pi].Size() > 0 || special_boundary_points.count(pi);
    };

    auto numGroups = [&](PointIndex pi) -> size_t {
        if(special_boundary_points.count(pi))
            return special_boundary_points[pi].growth_groups.Size();
        else
            return 1;
    };


    auto getGroups = [&] (PointIndex pi, int face_index) -> Array<int> {
      auto n = numGroups(pi);
      Array<int> groups;
      if(n == 1) {
        groups.Append(0);
        return groups;
      }
      const auto & all_groups = special_boundary_points[pi].growth_groups;
      for(auto i : Range(n))
          if(all_groups[i].faces.Contains(face_index))
            groups.Append(i);
      // cout << "groups " << pi << ", " << face_index << endl << groups;
      return groups;
    };


    // add 2d quads on required surfaces
    map<pair<PointIndex, PointIndex>, int> seg2edge;
    map<int, int> edge_map;
    int edge_nr = max_edge_nr; 
    auto getEdgeNr = [&] (int ei) {
      if(edge_map.count(ei) == 0)
        edge_map[ei] = ++edge_nr;
      return edge_map[ei];
    };
    if(params.grow_edges)
      {
        for(auto sei : moved_segs)
          {
            // copy here since we will add segments and this would
            // invalidate a reference!
            // auto segi = segments[sei];
            for(auto [sej, type] : segmap[sei])
              {
                auto segj = segments[sej];
                if(type == 0)
                  {
                    auto addSegment = [&] (PointIndex p0, PointIndex p1, bool extra_edge_nr = false) {
                      Segment s;
                      s[0] = p0;
                      s[1] = p1;
                      s[2] = PointIndex::INVALID;
                      auto pair = s[0] < s[1] ? make_pair(s[0], s[1]) : make_pair(s[1], s[0]);
                      if(extra_edge_nr)
                        s.edgenr = ++edge_nr;
                      else
                        s.edgenr = getEdgeNr(segj.edgenr);
                      s.si = si_map[segj.si];
                      new_segments.Append(s);
                      // cout << __LINE__ <<"\t" << s << endl;
                      return s;
                    };

                    auto p0=segj[0], p1=segj[1];
                    auto g0 = getGroups(p0, segj.si);
                    auto g1 = getGroups(p1, segj.si);

                    if(g0.Size() == 1 && g1.Size() == 1)
                      auto s = addSegment( newPoint(p0, -1, g0[0]), newPoint(p1, -1, g1[0]) );
                    else {
                      if(g0.Size() == 2)
                        addSegment( newPoint(p0, -1, g0[0]), newPoint(p0, -1, g0[1]) );
                      if(g1.Size() == 2)
                        addSegment( newPoint(p1, -1, g1[0]), newPoint(p1, -1, g1[1]) );
                    }
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
                    // cout << __LINE__ <<"\t" << s0 << endl;

                    for(auto i : Range(params.heights))
                      {
                        Element2d sel(QUAD);
                        p3 = newPoint(pp2, i);
                        p4 = newPoint(pp1, i);
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
                        s1.edgenr = getEdgeNr(segj.edgenr);
                        s1.si = segj.si;
                        // new_segments.Append(s1);
                        // cout << __LINE__ <<"\t" << s1 << endl;
                        Segment s2;
                        s2[0] = p4;
                        s2[1] = p1;
                        s2[2] = PointIndex::INVALID;
                        pair = make_pair(p1, p4);
                        s2.edgenr = getEdgeNr(segj.edgenr);
                        s2.si = segj.si;
                        // new_segments.Append(s2);
                        // cout << __LINE__ <<"\t" << s2 << endl;
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
                    // cout << __LINE__ << "\t" << s3 << endl;
                  }
              }
          }
      }

    auto getClosestGroup = [&] (PointIndex pi, SurfaceElementIndex sei)
    {
      auto n = numGroups(pi);
      if(n == 1) return 0;
      const auto & sel = mesh[sei];
      auto igroup = 0;
      double distance = 1e99;
      for(auto j : Range(n)) {
        // auto g = getGroups(pi, sel.GetIndex());
        auto vcenter = Center(mesh[sel[0]], mesh[sel[1]], mesh[sel[2]]);
        auto dist = (vcenter-(mesh[pi]+special_boundary_points[pi].growth_groups[j].growth_vector)).Length2();
        if(dist < distance) {
          distance = dist;
          igroup = j;
        }
      }
      return getGroups(pi, sel.GetIndex())[igroup];
    };

    BitArray fixed_points(np+1);
    fixed_points.Clear();
    BitArray moveboundarypoint(np+1);
    moveboundarypoint.Clear();
    auto p2el = mesh.CreatePoint2ElementTable();
    for(SurfaceElementIndex si = 0; si < nse; si++)
      {
        // copy because surfaceels array will be resized!
        const auto sel = mesh[si];
        if(moved_surfaces.Test(sel.GetIndex()))
          {
            Array<PointIndex> points(sel.PNums());
            if(surfacefacs[sel.GetIndex()] > 0) Swap(points[0], points[2]);
            ArrayMem<int, 4> groups(points.Size());
            for(auto i : Range(points))
              groups[i] = getClosestGroup(sel[i], si);
            bool add_volume_element = true;
            for(auto pi : sel.PNums())
              if(numGroups(pi)>1)
                add_volume_element = false;
            for(auto j : Range(params.heights))
              {
                auto eltype = points.Size() == 3 ? PRISM : HEX;
                Element el(eltype);
                for(auto i : Range(points))
                    el[i] = points[i];
                for(auto i : Range(points))
                  points[i] = newPoint(sel.PNums()[i], j, groups[i]);
                if(surfacefacs[sel.GetIndex()] > 0) Swap(points[0], points[2]);
                for(auto i : Range(points))
                  el[sel.PNums().Size() + i] = points[i];
                auto new_index = new_mat_nrs[sel.GetIndex()];
                if(new_index == -1)
                  throw Exception("Boundary " + ToString(sel.GetIndex()) + " with name " + mesh.GetBCName(sel.GetIndex()-1) + " extruded, but no new material specified for it!");
                el.SetIndex(new_mat_nrs[sel.GetIndex()]);
                if(add_volume_element)
                  mesh.AddVolumeElement(el);
                else
                {
                  // Let the volume mesher fill the hole with pyramids/tets
                  // To insert pyramids, we need close surface identifications on open quads
                  for(auto i : Range(points))
                    if(numGroups(sel[i]) == 1)
                      identifications.Add(el[i], el[i+points.Size()], identnr);
                }
              }
            Element2d newel = sel;
            for(auto i: Range(points))
              newel[i] = newPoint(sel[i], -1, groups[i]);
            newel.SetIndex(si_map[sel.GetIndex()]);
            mesh.AddSurfaceElement(newel);

            // also move volume element adjacent to this surface element accordingly
            ElementIndex ei = -1;
            // if(groups[0] || groups[1] || groups[2])
            // for(auto ei_ : p2el[sel.PNums()[0]])
            //   {
            //     const auto & el = mesh[ei_];
            //     // if(!domains.Test(el.GetIndex())) continue;
            //     cout << "check " << ei_ << "\t" << el << "\t" << sel << endl;
            //     auto pnums = el.PNums();
            //     if(pnums.Contains(sel[1]) && pnums.Contains(sel[2])) {
            //       ei = ei_;
            //       break;
            //     }
            //   }
            if(ei != -1) {
              auto & el = mesh[ei];
              for (auto i : Range(el.GetNP()))
                for (auto j : Range(3)) 
                {
                  if(groups[j] && el[i] == sel[j]) {
                    el[i] = newel[j];
                    break;
                  }
                }
            }
          }
        else
          {
            bool has_moved = false;
            for(auto p : sel.PNums())
              has_moved |= hasMoved(p);
            if(has_moved)
              for(auto p : sel.PNums())
                {
                  if(hasMoved(p))
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
              if(hasMoved(p)) 
                p = newPoint(p);
          }
      }

    for(SegmentIndex sei = 0; sei < nseg; sei++)
      {
        auto& seg = segments[sei];
        if(is_boundary_moved.Test(seg.si))
          for(auto& p : seg.PNums())
            if(hasMoved(p))
              p = newPoint(p);
        // else if(hasMoved(seg[0]) || hasMoved(seg[1]))
        // {
        //   auto tangent = mesh[seg[1]] - mesh[seg[0]];
        //   if(hasMoved(seg[0]) && growthvectors[seg[0]] * tangent > 0)
        //     seg[0] = newPoint(seg[0]);
        //   if(hasMoved(seg[1]) && growthvectors[seg[1]] * tangent < 0)
        //     seg[1] = newPoint(seg[1]);
        // }
      }

    // fill holes in surface mesh at special boundary points (with >=4 adjacent boundary faces)
    auto p2sel = mesh.CreatePoint2SurfaceElementTable();
    for(auto & [pi, special_point] : special_boundary_points) {
    if(special_point.growth_groups.Size() != 2)
      throw Exception("special_point.growth_groups.Size() != 2");
    for(auto igroup : Range(2)) {
      auto & group = special_point.growth_groups[igroup];
      std::set<int> faces;
      for(auto face : group.faces)
        faces.insert(si_map[face]);
      auto pi_new = group.new_points.Last();
      auto pi_new_other = special_point.growth_groups[1-igroup].new_points.Last();
      for(auto sei : p2sel[pi_new])
        faces.erase(mesh[sei].GetIndex()); 
      for(auto face : faces)
        for(auto seg : new_segments) {
          if(//seg.si == face
            (seg[0] == pi_new || seg[1] == pi_new)
            && (seg[0] != pi_new_other && seg[1] != pi_new_other)
          ) {
            bool is_correct_face = false;
            auto pi_other = seg[0] == pi_new ? seg[1] : seg[0];
            for(auto sei : p2sel[pi_other]) {
              if(mesh[sei].GetIndex() == face) {
                is_correct_face = true;
                break;
              }
            }
            if(is_correct_face) {
              Element2d sel;
              sel[0] = seg[1];
              sel[1] = seg[0];
              sel[2] = pi_new_other;
              sel.SetIndex(face);
              mesh.AddSurfaceElement(sel);
            }
          }
        }
      }
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
            if(hasMoved(p))
              moved.Append(p);
            if(moveboundarypoint.Test(p))
              moved_bnd = true;
          }

        bool do_move, do_insert;
        if(changed_domains.Test(el.GetIndex()))
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
              if(hasMoved(p)) {
                if(special_boundary_points.count(p)) {
                  auto & special_point = special_boundary_points[p];
                  auto & group = special_point.growth_groups[0];
                  p = group.new_points.Last();
                }
                else
                  p = newPoint(p);
              }
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
                    auto d = mesh[newPoint(p1,0)] - v1;
                    if(n*d > 0)
                      Swap(p2,p3);
                    PointIndex p4 = p1;
                    PointIndex p5 = p2;
                    PointIndex p6 = p3;
                    for(auto i : Range(params.heights))
                      {
                        Element nel(PRISM);
                        nel[0] = p4; nel[1] = p5; nel[2] = p6;
                        p4 = newPoint(p1, i); p5 = newPoint(p2, i); p6 = newPoint(p3, i);
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
                            PointIndex p3 = newPoint(moved[1], i);
                            PointIndex p4 = newPoint(moved[0], i);
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
                        PointIndex p2 = newPoint(moved[0], i);
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
                        PointIndex p3 = newPoint(moved[1], i);
                        PointIndex p4 = newPoint(moved[0], i);
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
              else if(do_move) {
                throw Exception("Boundarylayer only implemented for tets and pyramids outside yet!");
              }
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
        MergeAndAddSegments(mesh, segments, new_segments);
    else
    {
      mesh.LineSegments() = segments;
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
      for([[maybe_unused]] auto step : Range(0))
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
      
      // for(auto pi : points)
      // {
      //     mesh[pi] += height * growthvectors[pi];
      //     growthvectors[pi] = 0.0;
      // }
  }

  void BoundaryLayerTool :: Perform()
  {
      CreateNewFaceDescriptors();
      CalculateGrowthVectors();
      // cout << "growthvectors " << __LINE__ << endl << growthvectors << endl;
      CreateFaceDescriptorsSides();
      auto segmap = BuildSegMap();

      auto in_surface_direction = ProjectGrowthVectorsOnSurface();
      // cout << "growthvectors " << __LINE__ << endl << growthvectors << endl;

      // auto fout = ofstream("growthvectors.txt");
      // for (auto pi : Range(mesh.Points()))
      // {
      //   for(auto i : Range(3))
      //     fout << mesh[pi][i] << " ";
      //   for(auto i : Range(3))
      //     fout << mesh[pi][i]+height*growthvectors[pi][i] << " ";
      // }
      // fout << endl;


      // FixVolumeElements();
      // mesh.Save("before_insert.vol");

      InsertNewElements(segmap, in_surface_direction);

      SetDomInOut();
      AddSegments();

      mesh.CalcSurfacesOfNode();
      topo.SetBuildVertex2Element(true);
      mesh.UpdateTopology();
      InterpolateGrowthVectors();


      // cout << "growthvectors before "  << endl<< growthvectors << endl;
      // cout << "growthvectors after " << endl << growthvectors << endl;

      if(params.limit_growth_vectors)
        LimitGrowthVectorLengths();
      // cout << "growthvectors " << __LINE__ << endl << growthvectors << endl;

      // for(PointIndex pi : Range(PointIndex::BASE, this->np + PointIndex::BASE))
      // {
      //   cout << "move " << pi << "\tby " << total_height << " * " << growthvectors[pi] << endl;
      //   mesh[pi] += total_height * growthvectors[pi];
      // }
      for (auto [pi, data] : growth_vector_map) {
        auto [gw, height] = data;
        mesh[pi] += height * (*gw);
      }


      
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

  cout << "generate boundary layers" << endl;
  cout << "points before " << mesh.GetNP() << endl;
  cout << "sels before " << mesh.GetNSE() << endl;

    // debug_gui.DrawMesh(mesh);

    BoundaryLayerTool tool(mesh, blp);
  cout << tool.nse << ',' << tool.np << endl;
    tool.Perform();
    mesh.Save("layer.vol");
  }

} // namespace netgen

#ifndef NETGEN_BOUNDARYLAYER_HPP
#define NETGEN_BOUNDARYLAYER_HPP

#include <core/array.hpp>
#include <mystdlib.h>
#include <meshing.hpp>

namespace netgen
{

///
DLL_HEADER extern void InsertVirtualBoundaryLayer (Mesh& mesh);

/// Create a typical prismatic boundary layer on the given
/// surfaces

struct SpecialBoundaryPoint
{
  struct GrowthGroup
  {
    Array<int> faces;
    Vec<3> growth_vector;
    Array<PointIndex> new_points;

    GrowthGroup (FlatArray<int> faces_, FlatArray<Vec<3>> normals);
    GrowthGroup (const GrowthGroup&) = default;
    GrowthGroup () = default;
  };
  Array<GrowthGroup> growth_groups;
  Vec<3> separating_direction;

  SpecialBoundaryPoint (const std::map<int, Vec<3>>& normals);
  SpecialBoundaryPoint () = default;
};

DLL_HEADER void GenerateBoundaryLayer (Mesh& mesh,
                                       const BoundaryLayerParameters& blp);

DLL_HEADER int /* new_domain_number */ GenerateBoundaryLayer2 (Mesh& mesh, int domain, const Array<double>& thicknesses, bool should_make_new_domain = true, const Array<int>& boundaries = Array<int>{});

class BoundaryLayerTool
{
public:
  BoundaryLayerTool (Mesh& mesh_, const BoundaryLayerParameters& params_);
  void ProcessParameters ();
  void Perform ();

  Mesh& mesh;
  MeshTopology& topo;
  BoundaryLayerParameters params;
  Array<Vec<3>, PointIndex> growthvectors;
  std::map<PointIndex, Vec<3>> non_bl_growth_vectors;
  Table<SurfaceElementIndex, PointIndex> p2sel;

  BitArray domains, is_edge_moved, is_boundary_projected, is_boundary_moved;
  Array<SegmentIndex> moved_segs;
  int max_edge_nr, nfd_old, ndom_old;
  Array<int> new_mat_nrs;
  BitArray moved_surfaces;
  int np, nseg, nse, ne;
  PointIndex first_new_pi;
  double total_height;
  Array<POINTTYPE, PointIndex> point_types;

  // These parameters are derived from given BoundaryLayerParameters and the Mesh
  Array<double> par_heights;
  Array<int> par_surfid;
  bool insert_only_volume_elements;
  map<string, string> par_new_mat;
  bool have_material_map = false;
  Array<size_t> par_project_boundaries;

  bool have_single_segments;
  Array<Segment> old_segments, free_segments, segments, new_segments, new_segments_on_moved_bnd;
  Array<Element2d, SurfaceElementIndex> new_sels, new_sels_on_moved_bnd;
  Array<Array<PointIndex>, PointIndex> mapto;
  Array<PointIndex, PointIndex> mapfrom;

  Array<double> surfacefacs;
  Array<int> si_map;

  std::map<PointIndex, SpecialBoundaryPoint> special_boundary_points;
  std::map<PointIndex, std::tuple<Vec<3>*, double>> growth_vector_map;

  // major steps called in Perform()
  void CreateNewFaceDescriptors ();
  void CreateFaceDescriptorsSides ();
  void CalculateGrowthVectors ();
  Array<Array<pair<SegmentIndex, int>>, SegmentIndex> BuildSegMap ();

  BitArray ProjectGrowthVectorsOnSurface ();
  void InterpolateSurfaceGrowthVectors ();
  void InterpolateGrowthVectors ();
  void LimitGrowthVectorLengths ();
  void FixSurfaceElements ();

  void InsertNewElements (FlatArray<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap, const BitArray& in_surface_direction);
  void SetDomInOut ();
  void SetDomInOutSides ();
  void AddSegments ();
  void AddSurfaceElements ();

  Vec<3> getNormal (const Element2d& el)
  {
    auto v0 = mesh[el[0]];
    return Cross(mesh[el[1]] - v0, mesh[el[2]] - v0).Normalize();
  }

  Vec<3> getEdgeTangent (PointIndex pi, int edgenr, FlatArray<Segment*> segs);
};

} // namespace netgen
#endif // NETGEN_BOUNDARYLAYER_HPP

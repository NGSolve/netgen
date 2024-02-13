#ifndef NETGEN_BOUNDARYLAYER_HPP
#define NETGEN_BOUNDARYLAYER_HPP

namespace netgen
{

///
DLL_HEADER extern void InsertVirtualBoundaryLayer (Mesh & mesh);

/// Create a typical prismatic boundary layer on the given 
/// surfaces

class BoundaryLayerParameters
{
public:
  // parameters by Philippose ..
  Array<int> surfid;
  Array<double> heights;
  map<string, string> new_mat;
  BitArray domains;
  bool outside = false; // set the boundary layer on the outside
  bool grow_edges = false;
  bool limit_growth_vectors = true;
  double limit_safety = 0.3; // alloow only 30% of the growth vector length
  bool sides_keep_surfaceindex = false;
  bool keep_surfaceindex = false;
  Array<size_t> project_boundaries;
};

DLL_HEADER void GenerateBoundaryLayer (Mesh & mesh,
                                       const BoundaryLayerParameters & blp);

DLL_HEADER int /* new_domain_number */ GenerateBoundaryLayer2 (Mesh & mesh, int domain, const Array<double> & thicknesses, bool should_make_new_domain=true, const Array<int> & boundaries=Array<int>{});

class BoundaryLayerTool
{
  public:
    BoundaryLayerTool(Mesh & mesh_, const BoundaryLayerParameters & params_);
    void Perform();

  protected:
    Mesh & mesh;
    MeshTopology & topo;
    BoundaryLayerParameters params;
    Array<Vec<3>, PointIndex> growthvectors;
    Table<SurfaceElementIndex, PointIndex> p2sel;

    BitArray domains, is_edge_moved, is_boundary_projected, is_boundary_moved;
    Array<SegmentIndex> moved_segs;
    int max_edge_nr, nfd_old, ndom_old;
    Array<int> new_mat_nrs;
    BitArray moved_surfaces;
    int np, nseg, nse, ne;
    double height;

    bool have_single_segments;
    Array<Segment> segments, new_segments;

    Array<double> surfacefacs;
    Array<int> si_map;
    Array<double, PointIndex> limits;

    // major steps called in Perform()
    void CreateNewFaceDescriptors();
    void CreateFaceDescriptorsSides();
    void CalculateGrowthVectors();
    Array<Array<pair<SegmentIndex, int>>, SegmentIndex> BuildSegMap();

    BitArray ProjectGrowthVectorsOnSurface();
    void InterpolateSurfaceGrowthVectors();
    void InterpolateGrowthVectors();
    void LimitGrowthVectorLengths();

    void InsertNewElements(FlatArray<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap, const BitArray & in_surface_direction);
    void SetDomInOut();
    void SetDomInOutSides();
    void AddSegments();
    void FixVolumeElements();

    // utility functions
    array<Point<3>, 2> GetMappedSeg( PointIndex pi );
    ArrayMem<Point<3>, 4> GetFace( SurfaceElementIndex sei );
    ArrayMem<Point<3>, 4> GetMappedFace( SurfaceElementIndex sei );
    ArrayMem<Point<3>, 4> GetMappedFace( SurfaceElementIndex sei, int face );

    Vec<3> getNormal(const Element2d & el)
    {
        auto v0 = mesh[el[0]];
        return Cross(mesh[el[1]]-v0, mesh[el[2]]-v0).Normalize();
    }

    Vec<3> getEdgeTangent(PointIndex pi, int edgenr);
};

} // namespace netgen
#endif // NETGEN_BOUNDARYLAYER_HPP

#ifndef FILE_BOUNDARYLAYER
#define FILE_BOUNDARYLAYER


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
  string new_mat;
  BitArray domains;
  bool outside = false; // set the boundary layer on the outside
  bool grow_edges = false;
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
    int max_edge_nr, new_mat_nr, nfd_old;
    int np, nseg, nse, ne;

    bool have_single_segments;
    Array<Segment> segments, new_segments;

    Array<double> surfacefacs;
    Array<int> si_map;

    // major steps called in Perform()
    void CreateNewFaceDescriptors();
    void CalculateGrowthVectors();
    Array<Array<pair<SegmentIndex, int>>, SegmentIndex> BuildSegMap();

    BitArray ProjectGrowthVectorsOnSurface();
    void InterpolateGrowthVectors();

    void InsertNewElements(FlatArray<Array<pair<SegmentIndex, int>>, SegmentIndex> segmap, const BitArray & in_surface_direction);
    void SetDomInOut();
    void AddSegments();

    // utility functions
    array<Point<3>, 2> GetGrowSeg( PointIndex pi0 );

    ArrayMem<Point<3>, 4> GetGrowFace( SurfaceElementIndex sei );
    ArrayMem<Point<3>, 4> GetGrowFace( SurfaceElementIndex sei, int face );

    bool IsSegIntersectingPlane ( array<Point<3>, 2> seg, array<Point<3>, 3> trig, double & lam);
    bool IsIntersectingTrig ( array<Point<3>, 2> seg, array<Point<3>, 3> trig, double & lam);

    Vec<3> getNormal(const Element2d & el)
    {
        auto v0 = mesh[el[0]];
        return Cross(mesh[el[1]]-v0, mesh[el[2]]-v0).Normalize();
    }

    Vec<3> getEdgeTangent(PointIndex pi)
    {
        Vec<3> tangent = 0.0;
        for(auto segi : topo.GetVertexSegments(pi))
        {
            auto & seg = mesh[segi];
            tangent += (mesh[seg[1]] - mesh[seg[0]]);
        }
        return tangent.Normalize();
    }

};

#endif

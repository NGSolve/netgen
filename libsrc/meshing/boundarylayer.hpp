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
  Array<size_t> new_matnrs;
  BitArray domains;
  bool outside = false; // set the boundary layer on the outside
  bool grow_edges = false;
};

DLL_HEADER void GenerateBoundaryLayer (Mesh & mesh,
                                       const BoundaryLayerParameters & blp);


#endif

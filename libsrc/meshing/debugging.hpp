#include "meshfunc.hpp"


namespace netgen
{
    unique_ptr<Mesh> GetOpenElements( const Mesh & m, int dom = 0, bool only_quads = false );

    unique_ptr<Mesh> FilterMesh( const Mesh & m, FlatArray<PointIndex> points, FlatArray<SurfaceElementIndex> sels = Array<SurfaceElementIndex>{}, FlatArray<ElementIndex> els = Array<ElementIndex>{} );

    // Checks if the mesh is valid. This is called automatically on various places if debugparam.slowchecks is set
    void CheckMesh( const Mesh & m, MESHING_STEP meshing_step, const char *filename = "", int line = -1 );

    // Sometimes during SwapImprove we discover topological errors in the mesh. For instance, an edge is adjacent to 8 tets around it, but
    // the 8 "other" points of the tets don't form a closed path around the edge. Instead there are 2 sets of 4 points/tets each, which are not connected.
    // This function checks for such errors and returns true if any are found.
    void CheckElementsAroundEdges( const Mesh & m );
}

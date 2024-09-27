#include "meshclass.hpp"


namespace netgen
{
    unique_ptr<Mesh> GetOpenElements( const Mesh & m, int dom = 0 );

    unique_ptr<Mesh> FilterMesh( const Mesh & m, FlatArray<PointIndex> points, FlatArray<SurfaceElementIndex> sels = Array<SurfaceElementIndex>{}, FlatArray<ElementIndex> els = Array<ElementIndex>{} );



}

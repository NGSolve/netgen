#include "meshing.hpp"

namespace netgen
{

  static inline Point<2> P2( Point<3> p )
  {
    return {p[0], p[1]};
  }

  static inline Point<3> P3( Point<2> p )
  {
    return {p[0], p[1], 0};
  }

  class DelaunayTrig
  {
    PointIndex pnums[3];
    Point<2> c;

  public:
    double r;
    double rad2;
    DelaunayTrig () = default;
    DelaunayTrig (int p1, int p2, int p3)
    {
      pnums[0] = p1;
      pnums[1] = p2;
      pnums[2] = p3;
    }

    PointIndex & operator[] (int j) { return pnums[j]; }
    const PointIndex & operator[] (int j) const { return pnums[j]; }

    void CalcCenter (FlatArray<Point<2>, PointIndex> points);

    Point<2> Center() const { return c; }
    double Radius2() const { return rad2; }
    Box<2> BoundingBox() const { return Box<2> (c-Vec<2>(r,r), c+Vec<2>(r,r)); }

    mutable PointIndex visited_pi = -1;
  };

  class DelaunayMesh
  {
    ngcore::ClosedHashTable<IVec<2>, IVec<2>> edge_to_trig;
    Array<DelaunayTrig> trigs;
    unique_ptr<DelaunayTree<2>> tree;
    Array<Point<2>, PointIndex> & points;

    Array<int> closeels;
    Array<int> intersecting;
    Array<IVec<2>> edges;

    int GetNeighbour( int eli, int edge );

    void SetNeighbour( int eli, int edge );

    void UnsetNeighbours( int eli );

    void AppendTrig( int pi0, int pi1, int pi2 );

    public:
    DelaunayMesh( Array<Point<2>, PointIndex> & points_, Box<2> box  );

    void CalcIntersecting( PointIndex pi_new );
    void CalcWeights( PointIndex pi_new, std::map<PointIndex, double> & weights );
    void AddPoint( PointIndex pi_new );
    Array<DelaunayTrig> & GetElements() { return trigs; }
    unique_ptr<Mesh> GetMesh(PointIndex pi_new); // for debugging purposes
  };

} // namespace netgen

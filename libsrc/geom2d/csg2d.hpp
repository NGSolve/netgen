#ifndef NETGEN_CSG2D_HPP_INCLUDED
#define NETGEN_CSG2D_HPP_INCLUDED

#include <variant>

#include "geometry2d.hpp"

namespace netgen
{
using namespace std;
using namespace ngcore;
using netgen::Point;
using netgen::Vec;
using Spline = SplineSeg3<2>;
using netgen::Box;

inline double Area(const Point<2>& P, const Point<2>& Q, const Point<2>& R)
{
  return (Q[0]-P[0]) * (R[1]-P[1]) - (Q[1]-P[1]) * (R[0]-P[0]);
}

// compute weight of spline such that p lies on it
void ComputeWeight( Spline & s, Point<2> p );

enum IntersectionType
{                 // types of intersection (detected in the first phase)
  NO_INTERSECTION = 0,
  X_INTERSECTION,
  T_INTERSECTION_Q,
  T_INTERSECTION_P,
  V_INTERSECTION,
  X_OVERLAP,      // Q0 -- P1 -- Q1 -- P0   (different direction)
  T_OVERLAP_Q,    // same direction or P inside Q
  T_OVERLAP_P,    // same direction or Q inside P
  V_OVERLAP       // one common point
};

enum IntersectionLabel
{      // for the classification of intersection vertices in the second phase
  NONE,
  CROSSING,
  BOUNCING,
  LEFT_ON,
  RIGHT_ON,
  ON_ON,
  ON_LEFT,
  ON_RIGHT,
  DELAYED_CROSSING,
  DELAYED_BOUNCING
};

enum EntryExitLabel
{         // for marking intersection vertices as "entry" or "exit"
  EXIT,
  ENTRY,
  NEITHER
};

enum IteratorType
{
  SOURCE,
  INTERSECTION,
  CROSSING_INTERSECTION,
  ALL
};

inline constexpr const double MAXH_DEFAULT{1e99};
inline const string POINT_NAME_DEFAULT{""};
inline const string BC_DEFAULT{""};
inline const string MAT_DEFAULT{""};

struct EdgeInfo
{
  optional<Point<2>> control_point = nullopt; // for spline segments
  double maxh = MAXH_DEFAULT;
  string bc = BC_DEFAULT;

  EdgeInfo() = default;
  EdgeInfo(Point<2> p) : control_point(p) {}
  EdgeInfo(double h) : maxh(h) {}
  EdgeInfo(string s) : bc(s) {}
  EdgeInfo(optional<Point<2>> p, double h, string s)
      : control_point(p), maxh(h), bc(s)
    {}

  void Assign( EdgeInfo other )
    {
      if(other.control_point != nullopt)
          control_point = other.control_point;
      if(other.bc != BC_DEFAULT)
          bc = other.bc;
      if(other.maxh != MAXH_DEFAULT)
          maxh = min(maxh, other.maxh);
    }
};

struct PointInfo
{
  double maxh = MAXH_DEFAULT;
  string name = POINT_NAME_DEFAULT;
  PointInfo() = default;
  PointInfo(const PointInfo& other) = default;
  PointInfo(double amaxh) : maxh(amaxh) {}
  PointInfo(string aname) : name(aname) {}
  PointInfo(double amaxh, string aname) : maxh(amaxh), name(aname) {}

  void Assign(const PointInfo& other)
  {
    maxh = min(maxh, other.maxh);
    if(other.name != POINT_NAME_DEFAULT)
      name = other.name;
  }
};

struct Vertex : Point<2>
{
  Vertex (Point<2> p) : Point<2>(p) {}
  Vertex (const Vertex & v) : Point<2>(v)
  {
    spline = v.spline;
    info = v.info;
    pinfo = v.pinfo;
    is_source = true;
  }

  Vertex * prev = nullptr;
  Vertex * next = nullptr;
  unique_ptr<Vertex> pnext = nullptr;
  Vertex * neighbour = nullptr; // same vertex in other polygon (at intersections)
  double lam = -1.0;
  bool is_intersection = false;
  bool is_source = false;

  IntersectionLabel label = NONE;    // type of intersection vertex
  EntryExitLabel enex = NEITHER;        // entry/exit "flag"

  // In case the edge this - next is curved, store the spline information here
  optional<Spline> spline = nullopt;
  EdgeInfo info;
  PointInfo pinfo;

  DLL_HEADER Vertex * Insert(Point<2> p, double lam = -1.0);

  void Link( Vertex * v )
  {
    neighbour = v;
    v->neighbour = this;
    is_intersection = true;
    v->is_intersection = true;
  }
};

struct VertexIterator
{
  struct iterator
  {
    iterator(Vertex* root, IteratorType IterType) :
      root(root), V(NULL), iterType(IterType)
    {
      if (root == NULL)
        return;

      if (nextVertex() == NULL)         // no (source/intersection) vertex found
        root = V = NULL;                // -> mark iterator as "end"
    }

    const iterator& operator++()
    {
      nextVertex();
      return *this;
    }

    Vertex* operator*()
    {
      return V;
    }

    bool operator!=(const iterator& other) const
    {
      return (root != other.root) || (V != other.V);
    }

    private:
    Vertex* root;
    Vertex* V;
    IteratorType iterType;

    //
    // find the next vertex
    // if iterType is ALL, then it is just the next vertex
    // if iterType is SOURCE, then it is the next source vertex
    // if iterType is INTERSECTION, then it is the next intersection vertex
    // if iterType is CROSSING_INTERSECTION, then it is the next intersection vertex with CROSSING label
    //
    Vertex* nextVertex()
    {
      bool nextFound = false;

      if (V == NULL)
      {                  // find first (source/intersection) vertex
        V = root;
        switch(iterType)
        {
          case ALL:
            nextFound = true;
            break;
          case SOURCE:
            if (V->is_source)
              nextFound = true;
            break;
          case INTERSECTION:
            if (V->is_intersection)
              nextFound = true;
            break;
          case CROSSING_INTERSECTION:
            if (V->is_intersection && (V->label == CROSSING))
              nextFound = true;
            break;
        }
      }

      while (!nextFound)
      {              // find next (source/intersection) vertex
        switch(iterType)
        {
          case ALL:
            V = V->next;
            break;
          case SOURCE:
            do {
              V = V->next;
            } while (!V->is_source && V != root);
            break;
          case INTERSECTION:
            do {
              V = V->next;
            } while (!V->is_intersection && V != root);
            break;
          case CROSSING_INTERSECTION:
            do {
              V = V->next;
            } while ( ( !V->is_intersection || (V->label != CROSSING) ) && V != root);
            break;
        }

        if (V == root)
        {                // back at the root vertex?
          root = V = NULL;              // -> mark iterator as "end"
          return(V);
        }

        switch(iterType)
        {
          case ALL:
            nextFound = true;
            break;
          case SOURCE:
            if (V->is_source)
              nextFound = true;
            break;
          case INTERSECTION:
            if (V->is_intersection)
              nextFound = true;
            break;
          case CROSSING_INTERSECTION:
            if (V->is_intersection && (V->label == CROSSING))
              nextFound = true;
            break;
        }
      }
      return(V);
    }
  };

  public:
  VertexIterator() : root(NULL) {};

  iterator begin() { return iterator(root, iterType); }
  iterator end()   { return iterator(NULL, iterType); }

  Vertex* root;
  IteratorType iterType;
};


struct Edge
{
  Vertex * v0 = nullptr;
  Vertex * v1 = nullptr;

  Edge (Vertex* v, Vertex* w) : v0(v), v1(w) { };
};

struct EdgeIterator
{
  struct iterator
  {
    iterator(Vertex* root, IteratorType IterType) :
      root(root), one(NULL), two(NULL), iterType(IterType)
    {
      if (root == NULL)
        return;

      if (nextEdge() == NULL)           // no source edge found
        root = one = two = NULL;        // -> mark iterator as "end"
    }

    const iterator& operator++() { nextEdge(); return *this; }

    Edge operator*()
    {
      return Edge(one,two);
    }

    bool operator!=(const iterator& other) const
    {
      return (root != other.root) || (one != other.one) || (two != other.two);
    }

    private:
    Vertex* root;
    Vertex* one;
    Vertex* two;
    IteratorType iterType;

    //
    // find the next vertex, starting at curr
    // if iterType is ALL, then it is just the next vertex
    // if iterType is SOURCE, then it is the next source vertex
    //
    Vertex* nextVertex(Vertex* curr)
    {
      if (curr == NULL)
        return(NULL);

      switch(iterType)
      {
        case ALL:
          curr = curr->next;
          break;

        case SOURCE:
          do {
            curr = curr->next;
          } while (!curr->is_source);
          break;
        default:
          ;
      }

      return(curr);
    }

    //
    // find the next edge
    //
    Vertex* nextEdge()
    {
      if (root == NULL)                 // empty polygon?
        return (NULL);

      if (one == NULL)
      {                // find one (source) vertex
        one = root;                     // note: root is always a (source) vertex
        two = nextVertex(one);
        if (two == one)                 // just one (source) vertex
          return(NULL);                 // -> no (source) edges
        return(one);
      }

      if (two == root)
      {                // back at the root vertex?
        root = one = two = NULL;        // -> mark iterator as "end"
        return(NULL);
      }

      one = two;
      two = nextVertex(one);

      return (one);
    }
  };

  public:
  EdgeIterator() : root(NULL) {};

  iterator begin() { return iterator(root, iterType); }
  iterator end()   { return iterator(NULL, iterType); }

  Vertex* root;
  IteratorType iterType;
};


inline int CalcSide( const Point<2> & p0, const Point<2> & p1, const Point<2> & r )
{
  if ( (p0[1] < r[1]) != (p1[1] < r[1]) )
  {
    if (p0[0] >= r[0])
    {
      if (p1[0] > r[0])
        return 2 * (p1[1] > p0[1]) - 1;
      else
        if ( (Area(p0,p1,r) > 0) == (p1[1] > p0[1]) )
          return 2 * (p1[1] > p0[1]) - 1;
    }
    else
    {
      if (p1[0] > r[0])
        if ( (Area(p0,p1,r) > 0) == (p1[1] > p0[1]) )
          return 2 * (p1[1] > p0[1]) - 1;
    }
  }
  return 0;
}

struct Loop
{
  unique_ptr<Vertex> first = nullptr;
  unique_ptr<Box<2>> bbox = nullptr;

  Loop() = default;

  Loop(const Loop & p)
    : first(nullptr)
  {
    for(auto v : p.Vertices(ALL))
      AppendVertex(*v);
  }

  Loop(Loop && p) = default;

  Loop & operator=(Loop && p) = default;

  Loop & operator=(const Loop & p)
  {
    // static Timer t("Loop::operator="); RegionTimer rt(t);
    first = nullptr;
    if(p.first)
    {
      size_t n = p.Size();
      Array<unique_ptr<Vertex>> new_verts(n);
      {
        size_t i = 0;
        for(const auto v : p.Vertices(ALL))
          new_verts[i++] = make_unique<Vertex>(*v);
      }

      for(auto i : IntRange(n-1))
      {
        Vertex * v  = new_verts[i].get();
        Vertex * vn = new_verts[i+1].get();
        v->next = vn;
        vn->prev = v;
      }
      Vertex * vfirst = new_verts[0].get();
      Vertex * vlast = new_verts[n-1].get();
      vfirst->prev = vlast;
      vlast->next = vfirst;

      for(auto i : IntRange(1,n))
        new_verts[n-1-i]->pnext = std::move(new_verts[n-i]);

      first = std::move(new_verts[0]);
    }
    bbox = nullptr;
    return *this;
  }

  void Clear()
  {
    first = nullptr;
  }

  Vertex & AppendVertex(const Vertex & v)
  {
    auto & vnew = Append( static_cast<Point<2>>(v), true );
    vnew.info = v.info;
    vnew.pinfo = v.pinfo;
    if(v.spline)
      vnew.spline = *v.spline;
    if(bbox)
      bbox->Add(v);
    return vnew;
  }

  Vertex & Append(Point<2> p, bool source = false)
  {
    Vertex * vnew;
    if(first==nullptr)
    {
      first = make_unique<Vertex>(p);
      first->next = first.get();
      first->prev = first.get();
      vnew = first.get();
    }
    else
    {
      vnew = first->prev->Insert(p);
    }

    vnew->is_source = source;
    //     cout << "size after " << Size() << endl;
    if(bbox)
      bbox->Add(p);
    return *vnew;
  }

  void Remove (Vertex* v)
  {
    v->prev->next = v->next;
    v->next->prev = v->prev;
    if(first.get() == v)
      first = std::move(v->pnext);
    else
      v->prev->pnext = std::move(v->pnext);
    bbox.reset();
  }

  bool IsInside( Point<2> r ) const;
  bool IsLeftInside( const Vertex & p0 );
  bool IsRightInside( const Vertex & p0 );

  EdgeIterator Edges(IteratorType iterType) const
  {
    EdgeIterator it;
    it.iterType = iterType;
    it.root = first.get();
    return it;
  }

  VertexIterator Vertices(IteratorType iterType, Vertex* first_ = nullptr) const
  {
    VertexIterator it;
    it.iterType = iterType;
    it.root = (first_ == nullptr) ? first.get() : first_;
    return it;
  }

  //
  // check, if all vertices have the ON_ON label
  //
  bool allOnOn()
  {
    for (Vertex* v : Vertices(ALL))
      if (v->label != ON_ON)
        return(false);
    return(true);
  }

  //
  // check, if the polygon does not contain any crossing intersection vertex
  // or crossing intersection chain or (if we want to compute the union instead
  // of the intersection) a bouncing vertex or a bouncing intersection chain
  //
  bool noCrossingVertex(bool union_case = false)
  {
    for (Vertex* v : Vertices(ALL))
      if (v->is_intersection)
      {
        if ( (v->label == CROSSING) || (v->label == DELAYED_CROSSING) )
          return(false);

        if (union_case && ( (v->label == BOUNCING) || (v->label == DELAYED_BOUNCING) ) )
          return(false);
      }
    return(true);
  }

  //
  // return a non-intersection point
  //
  Point<2> getNonIntersectionPoint()
  {
    for (Vertex* v : Vertices(ALL))
      if (!v->is_intersection)
        return *v;

    // no non-intersection vertex found -> find suitable edge midpoint
    for (Vertex* v : Vertices(ALL))
      // make sure that edge from V to V->next is not collinear with other polygon
      if ( (v->next->neighbour != v->neighbour->prev) && (v->next->neighbour != v->neighbour->next) )
        // return edge midpoint
        return Center(*v, *v->next);
    throw Exception("no point found");
  }

  //
  // return and insert a non-intersection vertex
  //
  Vertex* getNonIntersectionVertex();

  void SetBC(string bc)
  {
    for(auto v : Vertices(ALL))
      v->info.bc = bc;
  }

  size_t Size() const
  {
    if(first==nullptr) return 0;

    size_t cnt = 0;

    for([[maybe_unused]] auto v : Vertices(ALL))
      cnt++;

    return cnt;
  }

  const Box<2> & GetBoundingBox()
  {
    if(bbox==nullptr)
    {
      static Timer tall("Loop::GetBoundingBox"); RegionTimer rt(tall);
      bbox = make_unique<Box<2>>(Box<2>::EMPTY_BOX);
      for(auto v : Vertices(ALL))
      {
        bbox->Add(*v);
        if(v->spline)
          bbox->Add(v->spline->TangentPoint());
      }
    }
    return *bbox;
  }
};


struct Solid2d
{
  Array<Loop> polys;

  int layer = 1;
  string name = MAT_DEFAULT;
  double maxh = MAXH_DEFAULT;

  Solid2d() = default;
  Solid2d(string name_) : name(name_) {}
  DLL_HEADER Solid2d(const Array<std::variant<Point<2>, EdgeInfo, PointInfo>> & points, string name_=MAT_DEFAULT, string bc_=BC_DEFAULT);
  Solid2d(Solid2d && other) = default;
  Solid2d(const Solid2d & other) = default;

  DLL_HEADER Solid2d operator+(const Solid2d & other) const;
  DLL_HEADER Solid2d operator*(const Solid2d & other) const;
  DLL_HEADER Solid2d operator-(const Solid2d & other) const;

  Solid2d& operator=(Solid2d && other) = default;
  Solid2d& operator=(const Solid2d & other) = default;
  DLL_HEADER Solid2d& operator+=(const Solid2d & other);
  DLL_HEADER Solid2d& operator*=(const Solid2d & other);
  DLL_HEADER Solid2d& operator-=(const Solid2d & other);

  void Append( const Loop & poly )
  {
    polys.Append(poly);
  }

  bool IsInside( Point<2> r ) const;
  bool IsLeftInside( const Vertex & p0 );
  bool IsRightInside( const Vertex & p0 );

  template<typename TFunc>
  Solid2d & Transform( const TFunc & func )
    {
      for(auto & poly : polys)
          for(auto v : poly.Vertices(ALL))
            {
              auto p = func(*v);
              (*v)[0] = p[0];
              (*v)[1] = p[1];
              if(v->spline)
                {
                  auto &s = *v->spline;
                  auto pmid = func(s.GetPoint(0.5));
                  s = Spline(func(s.StartPI()), func(s.TangentPoint()), func(s.EndPI()));
                  ComputeWeight(s, pmid);
                }
            }
      return *this;
    }

  Solid2d & Move( Vec<2> v );
  Solid2d & Scale( double s );
  Solid2d & Scale( Vec<2> s );
  Solid2d & RotateRad( double ang, Point<2> center = {0,0} );
  Solid2d & RotateDeg( double ang, Point<2> center = {0,0} )
    {
      return RotateRad( ang/180.*M_PI, center );
    }

  Solid2d & BC(string bc)
  {
    for(auto & p : polys)
      for(auto v : p.Vertices(ALL))
        v->info.bc = bc;
    return *this;
  }

  Solid2d & Maxh(double maxh)
  {
    this->maxh = maxh;
    for(auto & p : polys)
      for(auto v : p.Vertices(ALL))
        v->info.maxh = maxh;
    return *this;
  }

  Solid2d & Mat(string mat)
  {
    name = mat;
    return *this;
  }

  Solid2d & Layer(int layer_)
  {
    layer = layer_;
    return *this;
  }

  Box<2> GetBoundingBox() const;
};


class CSG2d
{
  public:
    Array<Solid2d> solids;

    void Add ( Solid2d s )
    {
      solids.Append(s);
    }

    DLL_HEADER shared_ptr<netgen::SplineGeometry2d> GenerateSplineGeometry();
    DLL_HEADER shared_ptr<netgen::Mesh> GenerateMesh(MeshingParameters & mp);
};

DLL_HEADER Solid2d Circle( Point<2> center, double r, string name="", string bc="");
DLL_HEADER Solid2d Rectangle( Point<2> p0, Point<2> p1, string mat=MAT_DEFAULT, string bc=BC_DEFAULT );

DLL_HEADER void AddIntersectionPoints ( Solid2d & s1, Solid2d & s2 );
DLL_HEADER Solid2d ClipSolids ( const Solid2d & s1, const Solid2d & s2, char op);
DLL_HEADER Solid2d ClipSolids ( const Solid2d & s1, Solid2d && s2, char op);
DLL_HEADER Solid2d ClipSolids ( Solid2d && s1, const Solid2d & s2, char op);
DLL_HEADER Solid2d ClipSolids ( Solid2d && s1, Solid2d && s2, char op);

DLL_HEADER IntersectionType intersect(const Point<2> P1, const Point<2> P2, const Point<2> Q1, const Point<2> Q2, double& alpha, double& beta);

}
#endif // NETGEN_CSG2D_HPP_INCLUDED

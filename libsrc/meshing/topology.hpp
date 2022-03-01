#ifndef TOPOLOGY
#define TOPOLOGY

/**************************************************************************/
/* File:   topology.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Apr. 01                                                    */
/**************************************************************************/

/*
    Mesh topology
    (Elements, Faces, Edges, Vertices
*/

namespace netgen
{
  
  /*
    struct T_EDGE
    {
    // int orient:1;
    int nr;    // 0-based
    };
    
    struct T_FACE
    {
    // int forient:3;
    int fnr;    // 0-based
    };
  */
  
  typedef int T_EDGE;
  typedef int T_FACE;
  

  /*
  template <typename T, int S>
  struct FixArray
  {
    T vals[S];
    T & operator[] (size_t i) { return vals[i]; }
    T operator[] (size_t i) const { return vals[i]; }
  };
  */

  template <typename T, int S>
  using FixArray = std::array<T,S>;

class MeshTopology
{
  const Mesh * mesh;
  bool buildvertex2element;
  bool buildedges;
  bool buildfaces;
  bool build_parent_edges = false; // may be changed to default = false
  bool build_parent_faces = false; // may be changed to default = false
  static bool static_buildedges, static_buildfaces, static_buildvertex2element;

  NgArray<INDEX_2> edge2vert;
  NgArray<INDEX_4> face2vert;
  /*
  NgArray<T_EDGE[12]> edges;
  NgArray<T_FACE[6]> faces;
  NgArray<T_EDGE[4]> surfedges;
  */
  NgArray<FixArray<T_EDGE,12>> edges;
  NgArray<FixArray<T_FACE,6>> faces;
  NgArray<FixArray<T_EDGE,4>> surfedges;
  
  NgArray<T_EDGE> segedges;
  NgArray<T_FACE> surffaces;
  NgArray<INDEX_2> surf2volelement;
  NgArray<int> face2surfel;
  Array<SegmentIndex> edge2segment;
  Table<ElementIndex, PointIndex> vert2element;
  Table<SurfaceElementIndex, PointIndex> vert2surfelement;
  Table<SegmentIndex,PointIndex> vert2segment;
  Table<int,PointIndex> vert2pointelement;
  int timestamp;
public:
  int GetNSurfedges() const {return surfedges.Size();}

  MeshTopology () = default;
  MeshTopology (MeshTopology && top) = default;
  MeshTopology (const Mesh & amesh);
  ~MeshTopology ();
  MeshTopology & operator= (MeshTopology && top) = default;

  void SetBuildVertex2Element (bool bv2e) { buildvertex2element = bv2e; }  
  void SetBuildEdges (bool be) { buildedges = be; }
  void SetBuildFaces (bool bf) { buildfaces = bf; }
  void SetBuildParentEdges (bool bh) { build_parent_edges = bh; }
  void SetBuildParentFaces (bool bh) { build_parent_faces = bh; }

  void EnableTable (string name, bool set);
  static void EnableTableStatic (string name, bool set);

  bool HasEdges () const  { return buildedges; }
  bool HasFaces () const  { return buildfaces; }
  bool HasParentEdges () const { return build_parent_edges; }

  void Update(NgTaskManager tm = &DummyTaskManager, NgTracer tracer = &DummyTracer);
  bool NeedsUpdate() const;


  int GetNEdges () const { return edge2vert.Size(); }
  int GetNFaces () const { return face2vert.Size(); }

  static inline short int GetNVertices (ELEMENT_TYPE et);
  static inline short int GetNPoints (ELEMENT_TYPE et);
  static inline short int GetNEdges (ELEMENT_TYPE et);
  static inline short int GetNFaces (ELEMENT_TYPE et);

  static const Point3d * GetVertices (ELEMENT_TYPE et);
  inline static const ELEMENT_EDGE * GetEdges1 (ELEMENT_TYPE et);
  inline static const ELEMENT_EDGE * GetEdges0 (ELEMENT_TYPE et);
  inline static const ELEMENT_FACE * GetFaces1 (ELEMENT_TYPE et);
  inline static const ELEMENT_FACE * GetFaces0 (ELEMENT_TYPE et);
  
  int GetSegmentEdge (int segnr) const { return segedges[segnr-1]+1; }
  int GetEdge (SegmentIndex segnr) const { return segedges[segnr]; }

  void GetSegmentEdge (int segnr, int & enr, int & orient) const
  {
    enr = segedges.Get(segnr)+1;
    // orient = segedges.Get(segnr).orient;
    orient = GetSegmentEdgeOrientation(segnr);
  }

  void GetElementEdges (int elnr, NgArray<int> & edges) const;
  void GetElementFaces (int elnr, NgArray<int> & faces, bool withorientation = false) const;
  void GetElementEdgeOrientations (int elnr, NgArray<int> & eorient) const;
  void GetElementFaceOrientations (int elnr, NgArray<int> & forient) const;

  int GetElementEdges (int elnr, int * edges, int * orient) const;
  int GetElementFaces (int elnr, int * faces, int * orient) const;

  int GetElementEdgeOrientation (int elnr, int locedgenr) const; // old style
  int GetElementFaceOrientation (int elnr, int locfacenr) const; // old style
  int GetSurfaceElementEdgeOrientation (int elnr, int locedgenr) const; // old style
  int GetSurfaceElementFaceOrientation2 (int elnr) const; // old style
  int GetSegmentEdgeOrientation (int elnr) const; // old style
  
  
  void GetFaceVertices (int fnr, NgArray<int> & vertices) const;
  void GetFaceVertices (int fnr, int * vertices) const;
  void GetEdgeVertices (int enr, int & v1, int & v2) const;
  void GetEdgeVertices (int enr, PointIndex & v1, PointIndex & v2) const;
  const int * GetEdgeVerticesPtr (int enr) const { return &edge2vert[enr][0]; }
  const int * GetFaceVerticesPtr (int fnr) const { return &face2vert[fnr][0]; }
  void GetFaceEdges (int fnr, NgArray<int> & edges, bool withorientation = false) const;

  ELEMENT_TYPE GetFaceType (int fnr) const
  { return (face2vert.Get(fnr)[3] == 0) ? TRIG : QUAD; }

  void GetSurfaceElementEdges (int elnr, NgArray<int> & edges) const;
  int GetSurfaceElementFace (int elnr) const;
  void GetSurfaceElementEdgeOrientations (int elnr, NgArray<int> & eorient) const;
  int GetSurfaceElementFaceOrientation (int elnr) const;
  void GetEdges (SurfaceElementIndex elnr, NgArray<int> & edges) const;
  int GetFace (SurfaceElementIndex elnr) const
  { return surffaces[elnr]; }

  int GetSurfaceElementEdges (int elnr, int * edges, int * orient) const;

  const T_EDGE * GetElementEdgesPtr (int elnr) const { return &edges[elnr][0]; }
  const T_EDGE * GetSurfaceElementEdgesPtr (int selnr) const { return &surfedges[selnr][0]; }
  const T_EDGE * GetSegmentElementEdgesPtr (int selnr) const { return &segedges[selnr]; }

  const T_FACE * GetElementFacesPtr (int elnr) const { return &faces[elnr][0]; }
  const T_FACE * GetSurfaceElementFacesPtr (int selnr) const { return &surffaces[selnr]; }


  void GetSurface2VolumeElement (int selnr, int & elnr1, int & elnr2) const
  { 
    elnr1 = surf2volelement.Get(selnr)[0];
    elnr2 = surf2volelement.Get(selnr)[1];
  }

  int GetFace2SurfaceElement (int fnr) const { return face2surfel[fnr-1]; }

  SegmentIndex GetSegmentOfEdge(int edgenr) const { return edge2segment[edgenr-1]; }
  
  void GetVertexElements (int vnr, Array<ElementIndex> & elements) const;
  FlatArray<ElementIndex> GetVertexElements (int vnr) const
  { return vert2element[vnr]; }

  void GetVertexSurfaceElements( int vnr, Array<SurfaceElementIndex>& elements ) const;
  const auto & GetVertexSurfaceElements( ) const { return vert2surfelement; }
  FlatArray<SurfaceElementIndex> GetVertexSurfaceElements(PointIndex vnr) const
  { return vert2surfelement[vnr]; }

  FlatArray<SegmentIndex> GetVertexSegments (int vnr) const
  { return vert2segment[vnr]; }

  FlatArray<int> GetVertexPointElements (int vnr) const
  { return vert2pointelement[vnr]; }
  
  int GetVerticesEdge ( int v1, int v2) const;
  void GetSegmentVolumeElements ( int segnr, NgArray<ElementIndex> & els ) const;
  void GetSegmentSurfaceElements ( int segnr, NgArray<SurfaceElementIndex> & els ) const;

  // Call this before Update() to discard old edges
  void ClearEdges() { edge2vert.SetSize(0); }

private:
  Array<std::tuple<int, std::array<int,3>>> parent_edges;
  void BuildParentEdges ();

  Array<std::tuple<int, std::array<int,4>>> parent_faces;
  void BuildParentFaces ();
public:
  auto GetParentEdges (int enr) const { return parent_edges[enr]; }
  auto GetParentFaces (int fnr) const { return parent_faces[fnr]; }
};










inline short int MeshTopology :: GetNVertices (ELEMENT_TYPE et)
{
  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return 2;

    case TRIG:
    case TRIG6:
      return 3;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return 4;

    case TET:
    case TET10:
      return 4;

    case PYRAMID:
    case PYRAMID13:
      return 5;

    case PRISM:
    case PRISM12:
    case PRISM15:
      return 6;

    case HEX:
    case HEX20:
      return 8;

      // default:
      // cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
  return 0;
}


inline short int MeshTopology :: GetNPoints (ELEMENT_TYPE et)
{
  switch (et)
    {
    case SEGMENT:
      return 2;
    case SEGMENT3:
      return 3;

    case TRIG:
      return 3;
    case TRIG6:
      return 6;

    case QUAD:
    case QUAD6:
      return 4;

    case QUAD8:
      return 8;

    case TET:
      return 4;
    case TET10:
      return 10;

    case PYRAMID:
      return 5;
    case PYRAMID13:
      return 13;

    case PRISM:
      return 6;
    case PRISM12:
      return 12;
    case PRISM15:
      return 15;

    case HEX:
      return 8;

    case HEX20:
      return 20;
      // default:
      // cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
  return 0;
}



inline short int MeshTopology :: GetNEdges (ELEMENT_TYPE et)
{
  __assume(et >= SEGMENT && et <= PYRAMID13);
  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return 1;

    case TRIG:
    case TRIG6:
      return 3;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return 4;

    case TET:
    case TET10:
      return 6;

    case PYRAMID:
    case PYRAMID13:
      return 8;

    case PRISM:
    case PRISM12:
    case PRISM15:
      return 9;

    case HEX:
    case HEX20:
      return 12;
    default:
      return 0;
      // default:
      // cerr << "Ng_ME_GetNEdges, illegal element type " << et << endl;
    }
  // return 0;
}


inline short int MeshTopology :: GetNFaces (ELEMENT_TYPE et)
{
  __assume(et >= SEGMENT && et <= PYRAMID13);
  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return 0;

    case TRIG:
    case TRIG6:
      return 1;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return 1;

    case TET:
    case TET10:
      return 4;

    case PYRAMID:
    case PYRAMID13:
      return 5;

    case PRISM:
    case PRISM12:
    case PRISM15:
      return 5;

    case HEX:
    case HEX20:
      return 6;

    default:
      return -99;
      // default:
      // cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
}






const ELEMENT_EDGE * MeshTopology :: GetEdges1 (ELEMENT_TYPE et)
{
  static ELEMENT_EDGE segm_edges[1] =
    { { 1, 2 }};

  static ELEMENT_EDGE trig_edges[3] =
    { { 3, 1 },
      { 2, 3 },        
      { 1, 2 }};

  static ELEMENT_EDGE quad_edges[4] =
    { { 1, 2 },
      { 3, 4 },
      { 4, 1 },
      { 2, 3 }};


  static ELEMENT_EDGE tet_edges[6] =
    { { 4, 1 },
      { 4, 2 },
      { 4, 3 }, 
      { 1, 2 },
      { 1, 3 },
      { 2, 3 }};

  static ELEMENT_EDGE prism_edges[9] =
    { { 3, 1 },
      { 1, 2 },
      { 3, 2 },
      { 6, 4 },
      { 4, 5 },
      { 6, 5 },
      { 3, 6 },
      { 1, 4 },
      { 2, 5 }};

  static ELEMENT_EDGE pyramid_edges[8] =
    { { 1, 2 },
      { 2, 3 },
      { 1, 4 },
      { 4, 3 },
      { 1, 5 },
      { 2, 5 },
      { 3, 5 },
      { 4, 5 }};

  static ELEMENT_EDGE hex_edges[12] =
    {
      { 1, 2 },
      { 3, 4 },
      { 4, 1 },
      { 2, 3 },
      { 5, 6 },
      { 7, 8 },
      { 8, 5 },
      { 6, 7 },
      { 1, 5 },
      { 2, 6 },
      { 3, 7 },
      { 4, 8 },
    };

  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return segm_edges;

    case TRIG:
    case TRIG6:
      return trig_edges;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return quad_edges;

    case TET:
    case TET10:
      return tet_edges;

    case PYRAMID:
    case PYRAMID13:
      return pyramid_edges;

    case PRISM:
    case PRISM12:
    case PRISM15:
      return prism_edges;

    case HEX:
    case HEX20:
      return hex_edges;
      // default:
      // cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
    }
   return 0;  
}



const ELEMENT_EDGE * MeshTopology :: GetEdges0 (ELEMENT_TYPE et)
{
  static ELEMENT_EDGE segm_edges[1] =
    { { 0, 1 }};

  static ELEMENT_EDGE trig_edges[3] =
    { { 2, 0 },
      { 1, 2 },        
      { 0, 1 }};

  static ELEMENT_EDGE quad_edges[4] =
    { { 0, 1 },
      { 2, 3 },
      { 3, 0 },
      { 1, 2 }};


  static ELEMENT_EDGE tet_edges[6] =
    { { 3, 0 },
      { 3, 1 },
      { 3, 2 }, 
      { 0, 1 },
      { 0, 2 },
      { 1, 2 }};

  static ELEMENT_EDGE prism_edges[9] =
    { { 2, 0 },
      { 0, 1 },
      { 2, 1 },
      { 5, 3 },
      { 3, 4 },
      { 5, 4 },
      { 2, 5 },
      { 0, 3 },
      { 1, 4 }};

  static ELEMENT_EDGE pyramid_edges[8] =
    { { 0, 1 },
      { 1, 2 },
      { 0, 3 },
      { 3, 2 },
      { 0, 4 },
      { 1, 4 },
      { 2, 4 },
      { 3, 4 }};

  static ELEMENT_EDGE hex_edges[12] =
    {
      { 0, 1 },
      { 2, 3 },
      { 3, 0 },
      { 1, 2 },
      { 4, 5 },
      { 6, 7 },
      { 7, 4 },
      { 5, 6 },
      { 0, 4 },
      { 1, 5 },
      { 2, 6 },
      { 3, 7 },
    };
  
  switch (et)
    {
    case SEGMENT:
    case SEGMENT3:
      return segm_edges;

    case TRIG:
    case TRIG6:
      return trig_edges;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return quad_edges;

    case TET:
    case TET10:
      return tet_edges;

    case PYRAMID:
    case PYRAMID13:
      return pyramid_edges;

    case PRISM:
    case PRISM12:
    case PRISM15:
      return prism_edges;

    case HEX:
    case HEX20:
      return hex_edges;
      // default:
      // cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
    }
   return 0;  
}










inline const ELEMENT_FACE * MeshTopology :: GetFaces1 (ELEMENT_TYPE et)
{
  static const ELEMENT_FACE trig_faces[1] = 
    { { 1, 2, 3, 0 } };
  static const ELEMENT_FACE quad_faces[1] = 
    { { 1, 2, 3, 4 } };

  static const ELEMENT_FACE tet_faces[4] =
    { { 4, 2, 3, 0 },
      { 4, 3, 1, 0 },
      { 4, 1, 2, 0 },
      { 1, 3, 2, 0 } };
  
  static const ELEMENT_FACE prism_faces[5] =
    {
      { 1, 3, 2, 0 },
      { 4, 5, 6, 0 },
      { 3, 1, 4, 6 },
      { 1, 2, 5, 4 },
      { 2, 3, 6, 5 } 
    };

  static const ELEMENT_FACE pyramid_faces[5] =
    {
      { 1, 2, 5, 0 },
      { 2, 3, 5, 0 },
      { 3, 4, 5, 0 },
      { 4, 1, 5, 0 },
      { 1, 4, 3, 2 } 
    };

  static const ELEMENT_FACE hex_faces[6] =
    {
      { 1, 4, 3, 2 },
      { 5, 6, 7, 8 },
      { 1, 2, 6, 5 },
      { 2, 3, 7, 6 },
      { 3, 4, 8, 7 },
      { 4, 1, 5, 8 }
    };


  
  switch (et)
    {
    case TRIG:
    case TRIG6:
      return trig_faces;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return quad_faces;


    case TET:
    case TET10:
      return tet_faces;

    case PRISM:
    case PRISM12:
    case PRISM15:
      return prism_faces;

    case PYRAMID:
    case PYRAMID13:
      return pyramid_faces;

    case SEGMENT:
    case SEGMENT3:

    case HEX:
    case HEX20:
      return hex_faces;

      // default:
      // cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
    }
  return 0;
}





inline const ELEMENT_FACE * MeshTopology :: GetFaces0 (ELEMENT_TYPE et)
{
  static const ELEMENT_FACE trig_faces[1] = 
    { { 0, 1, 2, -1 } };
  static const ELEMENT_FACE quad_faces[1] = 
    { { 0, 1, 2, 3 } };

  static const ELEMENT_FACE tet_faces[4] =
    { { 3, 1, 2, -1 },
      { 3, 2, 0, -1 },
      { 3, 0, 1, -1 },
      { 0, 2, 1, -1 } };
  
  static const ELEMENT_FACE prism_faces[5] =
    {
      { 0, 2, 1, -1 },
      { 3, 4, 5, -1 },
      { 2, 0, 3, 5 },
      { 0, 1, 4, 3 },
      { 1, 2, 5, 4 } 
    };

  static const ELEMENT_FACE pyramid_faces[5] =
    {
      { 0, 1, 4, -1 },
      { 1, 2, 4, -1 },
      { 2, 3, 4, -1 },
      { 3, 0, 4, -1 },
      { 0, 3, 2, 1 } 
    };

  static const ELEMENT_FACE hex_faces[6] =
    {
      { 0, 3, 2, 1 },
      { 4, 5, 6, 7 },
      { 0, 1, 5, 4 },
      { 1, 2, 6, 5 },
      { 2, 3, 7, 6 },
      { 3, 0, 4, 7 }
    };


  
  switch (et)
    {
    case TRIG:
    case TRIG6:
      return trig_faces;

    case QUAD:
    case QUAD6:
    case QUAD8:
      return quad_faces;


    case TET:
    case TET10:
      return tet_faces;

    case PRISM:
    case PRISM12:
    case PRISM15:
      return prism_faces;

    case PYRAMID:
    case PYRAMID13:
      return pyramid_faces;

    case SEGMENT:
    case SEGMENT3:

    case HEX:
    case HEX20:
      return hex_faces;

      // default:
      // cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
    }
  return 0;
}

}

#endif

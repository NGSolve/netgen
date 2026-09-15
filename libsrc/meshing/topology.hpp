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

#include "meshtype.hpp"

namespace netgen
{
  // typedef int T_EDGE;
  // typedef int T_FACE;

  class EdgeIndex : public Index<int,EdgeIndex,0>
  {
  public:
    using Index::Index;
  };

  class FaceIndex : public Index<int,FaceIndex,0>
  {
  public:
    using Index::Index;
  };

  typedef EdgeIndex T_EDGE;  
  typedef FaceIndex T_FACE;  
  
class MeshTopology
{
  const Mesh * mesh;
  bool buildvertex2element;
  bool buildedges;
  bool buildfaces;
  bool build_parent_edges = false; // may be changed to default = false
  bool build_parent_faces = false; // may be changed to default = false
  static bool static_buildedges, static_buildfaces, static_buildvertex2element;

  Array<std::array<PointIndex,2>> edge2vert;
  Array<std::array<PointIndex,4>> face2vert;

  Array<std::array<EdgeIndex,12>, ElementIndex> edges;
  Array<std::array<FaceIndex,6>, ElementIndex> faces;
  Array<std::array<EdgeIndex,4>, SurfaceElementIndex> surfedges;
  
  Array<EdgeIndex,SegmentIndex> segedges;
  Array<FaceIndex,SurfaceElementIndex> surffaces;
  // Array<INDEX_2, SurfaceElementIndex> surf2volelement;
  Array<std::array<ElementIndex,2>, SurfaceElementIndex> surf2volelement;
  Array<SurfaceElementIndex> face2surfel;
  
  Array<SegmentIndex> edge2segment;
  Table<ElementIndex, PointIndex> vert2element;
  Table<SurfaceElementIndex, PointIndex> vert2surfelement;
  Table<SegmentIndex,PointIndex> vert2segment;
  Table<int,PointIndex> vert2pointelement;
  int timestamp;
public:
  MeshTopology () = default;
  MeshTopology (MeshTopology && top) = default;
  DLL_HEADER MeshTopology (const Mesh & amesh);
  DLL_HEADER ~MeshTopology ();
  MeshTopology & operator= (MeshTopology && top) = default;

  void SetBuildVertex2Element (bool bv2e) { buildvertex2element = bv2e; }  
  void SetBuildEdges (bool be) { buildedges = be; }
  void SetBuildFaces (bool bf) { buildfaces = bf; }
  void SetBuildParentEdges (bool bh) { build_parent_edges = bh; }
  void SetBuildParentFaces (bool bh) { build_parent_faces = bh; }

  DLL_HEADER void EnableTable (string name, bool set);
  static void EnableTableStatic (string name, bool set);

  bool HasEdges () const  { return buildedges; }
  bool HasFaces () const  { return buildfaces; }
  bool HasParentEdges () const { return build_parent_edges; }
  bool HasParentFaces () const { return build_parent_faces; }

  void Update(NgTaskManager tm = &DummyTaskManager, NgTracer tracer = &DummyTracer);
  bool NeedsUpdate() const;


  size_t GetNEdges () const { return edge2vert.Size(); }
  size_t GetNFaces () const { return face2vert.Size(); }

  static inline short int GetNVertices (ELEMENT_TYPE et);
  static inline short int GetNPoints (ELEMENT_TYPE et);
  static inline short int GetNEdges (ELEMENT_TYPE et);
  static inline short int GetNFaces (ELEMENT_TYPE et);

  DLL_HEADER static const Point<3> * GetVertices (ELEMENT_TYPE et);
  inline static const ELEMENT_EDGE * GetEdges1 (ELEMENT_TYPE et);
  inline static const ELEMENT_EDGE * GetEdges0 (ELEMENT_TYPE et);
  inline static FlatArray<ELEMENT_EDGE> GetEdges (ELEMENT_TYPE et);
  inline static const ELEMENT_FACE * GetFaces1 (ELEMENT_TYPE et);
  inline static const ELEMENT_FACE * GetFaces0 (ELEMENT_TYPE et);

  [[deprecated("use GetEdge(SegmentIndex) instead")]]                    
  EdgeIndex GetSegmentEdge (int segnr) const { return segedges[IndexBASE<SegmentIndex>()+(segnr-1)]+1; }
  
  EdgeIndex GetEdge (SegmentIndex segnr) const { return segedges[segnr]; }
  inline FlatArray<EdgeIndex> GetEdges (SegmentIndex segnr) const;

  [[deprecated("use GetEdge(SegmentIndex) instead")]]                      
  void GetSegmentEdge (int segnr, int & enr, int & orient) const;

  void GetElementFaces (int elnr, Array<int> & faces, bool withorientation) const;  

  // definition in meshclass.hpp 
  inline FlatArray<EdgeIndex> GetEdges (ElementIndex elnr) const;
  inline FlatArray<FaceIndex> GetFaces (ElementIndex elnr) const;    

  
  // [[deprecated("use GetElementEdge instead")]]                        
  void GetElementEdgeOrientations (int elnr, Array<int> & eorient) const;
  // [[deprecated("use GetElementEdge instead")]]                        
  void GetElementFaceOrientations (int elnr, Array<int> & forient) const;

  [[deprecated("use GetEdges (ElementIndex) -> FlatArray")]]                            
  int GetElementEdges (int elnr, int * edges, int * orient) const;

  // [[deprecated("use GetFaces (ElementIndex) -> FlatArray")]]                              
  int GetElementFaces (int elnr, int * faces, int * orient) const;

  // [[deprecated("use GetElementEdge instead")]]                      
  int GetElementEdgeOrientation (int elnr, int locedgenr) const; // old style
  // [[deprecated("use GetElementEdge instead")]]                        
  int GetElementFaceOrientation (int elnr, int locfacenr) const; // old style
  // [[deprecated("use GetElementEdge instead")]]                        
  int GetSurfaceElementEdgeOrientation (int elnr, int locedgenr) const; // old style
  // [[deprecated("use GetElementEdge instead")]]                        
  int GetSurfaceElementFaceOrientation2 (int elnr) const; // old style
  // [[deprecated("use GetElementEdge instead")]]                        
  int GetSegmentEdgeOrientation (int elnr) const; // old style
  
  DLL_HEADER void GetFaceVertices (int fnr, Array<int> & vertices) const;
  DLL_HEADER void GetFaceVertices (int fnr, int * vertices) const;
  auto GetFaceVertices (int fnr) const
  { return FlatArray (face2vert[fnr][3].IsValid() ? 4 : 3, &face2vert[fnr][0]); }
  auto GetEdgeVertices (int enr) const { return std::array{edge2vert[enr][0], edge2vert[enr][1]}; }
  auto GetEdgeVerticesPtr (int enr) const { return &edge2vert[enr][0]; }
  auto GetFaceVerticesPtr (int fnr) const { return &face2vert[fnr][0]; }
  DLL_HEADER void GetFaceEdges (int fnr, Array<int> & edges, bool withorientation = false) const;

  // ELEMENT_TYPE GetFaceType (int fnr) const
  // { return (!face2vert[fnr-1][3].IsValid()) ? TRIG : QUAD; }    
  ELEMENT_TYPE GetFaceType0 (int fnr) const   // a face number, not a surface element
  { return (!face2vert[fnr][3].IsValid()) ? TRIG : QUAD; }    

  // [[deprecated("orientation is outdated")]]                            
  int GetSurfaceElementFaceOrientation (int elnr) const;


  inline FlatArray<EdgeIndex> GetEdges (SurfaceElementIndex elnr) const;
  inline FlatArray<FaceIndex> GetFaces (SurfaceElementIndex elnr) const;
  // { return FlatArray<EdgeIndex>(GetNEdges ( (*mesh)[elnr].GetType()), &surfedges[elnr][0]); }
  
  int GetFace (SurfaceElementIndex elnr) const
  { return surffaces[elnr]; }

  int GetSurfaceElementEdges (int elnr, int * edges, int * orient) const;

  int GetNSurfedges() const {return surfedges.Size();}



  void GetSurface2VolumeElement (SurfaceElementIndex selnr, ElementIndex & elnr1, ElementIndex & elnr2) const
  { 
    elnr1 = surf2volelement[selnr][0];
    elnr2 = surf2volelement[selnr][1];
  }

  std::array<ElementIndex,2> GetSurface2VolumeElement (SurfaceElementIndex sei) 
  {
    return surf2volelement[sei];
  }

  SurfaceElementIndex GetFace2SurfaceElement (int fnr) const { return face2surfel[fnr]; }

  SegmentIndex GetSegmentOfEdge(int edgenr) const { return edge2segment[edgenr-1]; }

  
  FlatArray<ElementIndex> GetVertexElements (PointIndex vnr) const
  { return vert2element[vnr]; }

  const auto & GetVertexSurfaceElements( ) const { return vert2surfelement; }
  
  FlatArray<SurfaceElementIndex> GetVertexSurfaceElements(PointIndex vnr) const
  { return vert2surfelement[vnr]; }

  FlatArray<SegmentIndex> GetVertexSegments (PointIndex vnr) const
  { return vert2segment[vnr]; }

  FlatArray<int> GetVertexPointElements (PointIndex vnr) const
  { return vert2pointelement[vnr]; }
  
  DLL_HEADER int GetVerticesEdge ( PointIndex v1, PointIndex v2) const;
  void GetSegmentVolumeElements ( int segnr, Array<ElementIndex> & els ) const;
  void GetSegmentSurfaceElements ( int segnr, Array<SurfaceElementIndex> & els ) const;

  // Call this before Update() to discard old edges/faces (e.g. after Compress)
  void ClearEdges() { edge2vert.SetSize0(); }
  void ClearFaces() { face2vert.SetSize0(); }

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

    case HEX7:
      return 7;
      
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

    case HEX7:
      return 7;

    case HEX:
      return 8;

    case HEX20:
      return 20;
      // default:
      // cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
    }
  return -99;
}



inline short int MeshTopology :: GetNEdges (ELEMENT_TYPE et)
{
  // __assume(et >= SEGMENT && et <= PYRAMID13);
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

    case HEX7:
      return 11;
      
    case HEX:
    case HEX20:
      return 12;
      // default:
      // cerr << "Ng_ME_GetNEdges, illegal element type " << et << endl;
    }
  return 0;
}


inline short int MeshTopology :: GetNFaces (ELEMENT_TYPE et)
{
  // __assume(et >= SEGMENT && et <= PYRAMID13);
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
    case HEX7:      
      return 6;

    default:
      return 0;
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

  static ELEMENT_EDGE hex7_edges[11] =
    {
      { 1, 2 },
      { 3, 4 },
      { 4, 1 },
      { 2, 3 },
      { 5, 6 },
      { 7, 5 },
      { 6, 7 },
      { 1, 5 },
      { 2, 6 },
      { 3, 7 },
      { 4, 7 },
    };

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

    case HEX7:
      return hex7_edges;
      
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

  static ELEMENT_EDGE hex7_edges[11] =
    {
      { 0, 1 },
      { 2, 3 },
      { 3, 0 },
      { 1, 2 },
      { 4, 5 },
      { 6, 4 },
      { 5, 6 },
      { 0, 4 },
      { 1, 5 },
      { 2, 6 },
      { 3, 6 },
    };

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

    case HEX7:
      return hex7_edges;
      
    case HEX:
    case HEX20:
      return hex_edges;
      // default:
      // cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
    }
   return 0;  
}


FlatArray<ELEMENT_EDGE> MeshTopology :: GetEdges (ELEMENT_TYPE et)
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

  static ELEMENT_EDGE hex7_edges[11] =
    {
      { 0, 1 },
      { 2, 3 },
      { 3, 0 },
      { 1, 2 },
      { 4, 5 },
      { 6, 4 },
      { 5, 6 },
      { 0, 4 },
      { 1, 5 },
      { 2, 6 },
      { 3, 6 },
    };

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
      return { 1, segm_edges };

    case TRIG:
    case TRIG6:
      return { 3, trig_edges };

    case QUAD:
    case QUAD6:
    case QUAD8:
      return { 4, quad_edges };

    case TET:
    case TET10:
      return { 6, tet_edges };

    case PYRAMID:
    case PYRAMID13:
      return { 8, pyramid_edges };

    case PRISM:
    case PRISM12:
    case PRISM15:
      return { 9, prism_edges };

    case HEX7:
      return { 11, hex7_edges };

    case HEX:
    case HEX20:
      return { 12, hex_edges };
      // default:
      // cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
    }
  return { 0, nullptr };  
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

  static const ELEMENT_FACE hex7_faces[6] =
    {
      { 1, 4, 3, 2 },
      { 5, 6, 7, 0  },
      { 1, 2, 6, 5 },
      { 2, 3, 7, 6 },
      { 3, 4, 7, 0 },
      { 4, 1, 5, 7 }
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

    case HEX7:
      return hex7_faces;
    
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

  static const ELEMENT_FACE hex7_faces[6] =
    {
      { 0, 3, 2, 1 },
      { 4, 5, 6, -1},
      { 0, 1, 5, 4 },
      { 1, 2, 6, 5 },
      { 2, 3, 6, -1},
      { 3, 0, 4, 6 }
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

    case HEX7:
      return hex7_faces;

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

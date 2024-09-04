#ifndef NGINTERFACE_V2
#define NGINTERFACE_V2


/**************************************************************************/
/* File:   nginterface_v2.hpp                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   May  09                                                        */
/**************************************************************************/

#include "mydefs.hpp"
#include <core/mpi_wrapper.hpp>

/*
  C++ interface to Netgen
*/

#ifndef NGINTERFACE
  // implemented element types:
enum NG_ELEMENT_TYPE { 
  NG_PNT = 0,
  NG_SEGM = 1, NG_SEGM3 = 2,
  NG_TRIG = 10, NG_QUAD=11, NG_TRIG6 = 12, NG_QUAD6 = 13, NG_QUAD8 = 14,
  NG_TET = 20, NG_TET10 = 21, 
  NG_PYRAMID = 22, NG_PRISM = 23, NG_PRISM12 = 24, NG_PRISM15 = 27, NG_PYRAMID13 = 28,
  NG_HEX = 25, NG_HEX20 = 26, NG_HEX7 = 29
};

enum NG_REFINEMENT_TYPE { NG_REFINE_H = 0, NG_REFINE_P = 1, NG_REFINE_HP = 2 };
#endif

// #ifndef PARALLEL
// typedef int MPI_Comm;
// #endif


namespace netgen
{
  using namespace std;
  using namespace ngcore;
  
  // extern DLL_HEADER NgMPI_Comm ng_comm;
  
  static constexpr int POINTINDEX_BASE = 1;

  /*
  struct T_EDGE2
  {
    // int orient:1;
    // int nr:31;    // 0-based
    int nr;    // 0-based
  };
  struct T_FACE2
  {
    // int orient:3;
    // int nr:29;    // 0-based
    int nr;    // 0-based
  };
  */
  typedef int T_EDGE2; 
  typedef int T_FACE2; 

  template <typename T>
  class Ng_Buffer
  {
    size_t s;
    T * data;
  public:
    Ng_Buffer (size_t as, T * adata)
      : s(as), data(adata) { ; }
    Ng_Buffer (Ng_Buffer && buffer)
      : s(buffer.Size()), data(buffer.Release()) { ; }
    ~Ng_Buffer () { delete [] data; }
    size_t Size() const { return s; }
    T * Release() { T * hd = data; data = nullptr; return hd; }
  };

  template <typename T, int S>
  class Ng_BufferMS
  {
    size_t s;
    T data[S];
  public:
    Ng_BufferMS (size_t as) : s(as) { ; } 
    size_t Size() const { return s; }
    T & operator[] (size_t i) { return data[i]; }
    T operator[] (size_t i) const { return data[i]; }
  };

  
  class Ng_Element
  {

    class Ng_Points
    {
    public:
      size_t num;
      const int * ptr;
  
      size_t Size() const { return num; }
      int operator[] (size_t i) const { return ptr[i]-POINTINDEX_BASE; }
    };


    class Ng_Vertices
    {
    public:
      size_t num;
      const int * ptr;
  
      size_t Size() const { return num; }
      int operator[] (size_t i) const { return ptr[i]-POINTINDEX_BASE; }
    };

    /*
    class Ng_Edges
    {
    public:
      size_t num;
      const T_EDGE2 * ptr;
  
      size_t Size() const { return num; }
      int operator[] (size_t i) const { return ptr[i]; }
    };

    class Ng_Faces
    {
    public:
      size_t num;
      const T_FACE2 * ptr;
  
      size_t Size() const { return num; }
      int operator[] (size_t i) const { return ptr[i]; }
    };
    */

    class Ng_Facets
    {
    public:
      size_t num;
      int base;
      const int * ptr;
      
      size_t Size() const { return num; }
      int operator[] (size_t i) const { return ptr[i]-base; }
    };

    
  public:
    NG_ELEMENT_TYPE type;
    int index;           // material / boundary condition 
    string_view mat;   // material / boundary label
    NG_ELEMENT_TYPE GetType() const { return type; }
    int GetIndex() const { return index-1; }
    Ng_Points points;      // all points
    Ng_Vertices vertices;
    // Ng_Edges edges;
    FlatArray<T_EDGE2> edges;
    // Ng_Faces faces;
    FlatArray<T_FACE2> faces;    
    Ng_Facets facets;
    bool is_curved;
    int8_t newest_vertex;
  };

  
  class Ng_Point
  {
    double * pt;
  public:
    Ng_Point (double * apt) : pt(apt) { ; }
    double operator[] (size_t i)
    { return pt[i]; }
    operator const double * () { return pt; }
  };




  template <int DIM> class Ng_Node;

  template <>
  class Ng_Node<0>
  {
    class Ng_Elements
    {
    public:
      size_t ne;
      const int * ptr;
  
      size_t Size() const { return ne; }
      int operator[] (size_t i) const { return ptr[i]; }
    };


  public:
    Ng_Elements elements;
    Ng_Elements bnd_elements;
  };



  
  template <>
  class Ng_Node<1>
  {
    class Ng_Vertices
    {
    public:
      const int * ptr;
  
      size_t Size() const { return 2; }
      int operator[] (size_t i) const { return ptr[i]-POINTINDEX_BASE; }
    };


  public:
    Ng_Vertices vertices;
  };



  template <>
  class Ng_Node<2>
  {
    class Ng_Vertices
    {
    public:
      size_t nv;
      const int * ptr;
  
      size_t Size() const { return nv; }
      int operator[] (size_t i) const { return ptr[i]-POINTINDEX_BASE; }
    };

    /*
    class Ng_Edges
    {
    public:
      size_t ned;
      const int * ptr;
  
      size_t Size() const { return ned; }
      int operator[] (size_t i) const { return ptr[i]-1; }
    };
    */

  public:
    Ng_Vertices vertices;
    // Ng_Edges edges;
    int surface_el;  // -1 if face not on surface
  };



    




  class Mesh;


  inline void DummyTaskManager2 (function<void(int,int)> func)
  { func(0,1); }
  inline void DummyTracer2 (string, bool) { ; } 
  
  class DLL_HEADER Ngx_Mesh
  {
  private:
    shared_ptr<Mesh> mesh;
    
  public:
    // Ngx_Mesh () { ; }
    // Ngx_Mesh(class Mesh * amesh) : mesh(amesh) { ; }

    /** reuse a netgen-mesh **/
    Ngx_Mesh (shared_ptr<Mesh> amesh); 
    /** load a new mesh **/
    Ngx_Mesh (string filename, NgMPI_Comm acomm = NgMPI_Comm{});
    
    void LoadMesh (const string & filename, NgMPI_Comm comm = NgMPI_Comm{});

    void LoadMesh (istream & str, NgMPI_Comm comm = NgMPI_Comm{});
    void SaveMesh (ostream & str) const;
    void UpdateTopology ();
    void DoArchive (Archive & archive);

    const NgMPI_Comm & GetCommunicator() const;
    
    virtual ~Ngx_Mesh();

    bool Valid () const { return mesh != NULL; }
    
    int GetDimension() const;
    int GetNLevels() const;
    size_t GetNVLevel (int level) const;
    
    int GetNElements (int dim) const;
    int GetNNodes (int nt) const;

    Ng_Point GetPoint (int nr) const;

    template <int DIM> 
    Ng_Element GetElement (size_t nr) const;

    template <int DIM> 
    int GetElementIndex (size_t nr) const;

    /// material/boundary label of region, template argument is co-dimension
    template <int DIM> 
    string_view GetMaterialCD (int region_nr) const;

    /// Curved Elements:
    /// elnr .. element nr
    /// xi..... DIM_EL local coordinates
    /// x ..... DIM_SPACE global coordinates 
    /// dxdxi...DIM_SPACE x DIM_EL Jacobian matrix (row major storage)
    template <int DIM_EL, int DIM_SPACE> 
    void ElementTransformation (int elnr,
                                const double * xi, 
                                double * x, 
                                double * dxdxi) const;
    
    
    /// Curved Elements:
    /// elnr .. element nr
    /// npts .. number of points
    /// xi..... DIM_EL local coordinates
    /// sxi ... step xi
    /// x ..... DIM_SPACE global coordinates
    /// dxdxi...DIM_SPACE x DIM_EL Jacobian matrix (row major storage)
    template <int DIM_EL, int DIM_SPACE, typename T> 
    void MultiElementTransformation (int elnr, int npts,
                                     const T * xi, size_t sxi,
                                     T * x, size_t sx,
                                     T * dxdxi, size_t sdxdxi) const;
    

    template <int DIM>
    const Ng_Node<DIM> GetNode (int nr) const;

    Ng_BufferMS<int,4> GetFaceEdges (int fnr) const;
    
    template <int DIM>
    int GetNNodes ();

    // returns domain numbers of domains next to boundary bnr -> (domin, domout)
    // 3D only
    // std::pair<int,int> GetBoundaryNeighbouringDomains (int bnr);

    template <int DIM>
      void SetRefinementFlag (size_t elnr, bool flag);
    
    void Curve (int order);
    int GetCurveOrder ();

    void EnableTable (string name, bool set);

    void Refine (NG_REFINEMENT_TYPE reftype, bool onlyonce,
                 void (*taskmanager)(function<void(int,int)>) = &DummyTaskManager2,
                 void (*tracer)(string, bool) = &DummyTracer2);

    int GetHPElementLevel (int ei, int dir) const;
  
    void GetParentNodes (int ni, int * parents) const;
    int GetParentElement (int ei) const;
    int GetParentSElement (int ei) const;

    bool HasParentEdges() const;
    std::tuple<int, std::array<int,3>> GetParentEdges (int enr) const;
    std::tuple<int, std::array<int,4>> GetParentFaces (int fnr) const;
    
    int GetNIdentifications() const;
    int GetIdentificationType(int idnr) const;
    Ng_Buffer<int[2]> GetPeriodicVertices(int idnr) const;

    // Find element of point, returns local coordinates
    template <int DIM>
    int FindElementOfPoint 
    (double * p, double * lami,
     bool build_searchtrees = false, 
     int * const indices = NULL, int numind = 0) const;
    

    // for MPI-parallel
    FlatArray<int> GetDistantProcs (int nodetype, int locnum) const;
    size_t GetGlobalVertexNum (int locnum) const;
                               
    shared_ptr<Mesh> GetMesh () const { return mesh; } 
    shared_ptr<Mesh> SelectMesh () const;
    inline auto GetTimeStamp() const;


    // also added from nginterface.h, still 1-based, need redesign
    void HPRefinement (int levels, double parameter = 0.125,
                       bool setorders = true,bool ref_level = false);
    void SplitAlfeld ();
    
    size_t GetNP() const;
    int GetSurfaceElementSurfaceNumber (size_t ei) const;
    int GetSurfaceElementFDNumber (size_t ei) const;

    int GetElementOrder (int enr) const;
    void GetElementOrders (int enr, int * ox, int * oy, int * oz) const;
    void SetElementOrder (int enr, int order);
    void SetElementOrders (int enr, int ox, int oy, int oz);
    int GetSurfaceElementOrder (int enr) const;
    void GetSurfaceElementOrders (int enr, int * ox, int * oy) const;
    void SetSurfaceElementOrder (int enr, int order);
    void SetSurfaceElementOrders (int enr, int ox, int oy);
    int GetClusterRepVertex (int vi) const;
    int GetClusterRepEdge (int edi) const;
    int GetClusterRepFace (int fai) const;
    int GetClusterRepElement (int eli) const;
  };



  DLL_HEADER Ngx_Mesh * LoadMesh (const string & filename);
}


#ifdef HAVE_NETGEN_SOURCES
#include <meshing.hpp>

namespace netgen
{
#ifdef __GNUC__
#define NGX_INLINE  __attribute__ ((__always_inline__)) inline
#else
#define NGX_INLINE inline
#endif  
#include <nginterface_v2_impl.hpp>
}

#endif


#endif


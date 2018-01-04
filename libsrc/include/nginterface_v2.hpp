#ifndef NGINTERFACE_V2
#define NGINTERFACE_V2


/**************************************************************************/
/* File:   nginterface_v2.hpp                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   May  09                                                        */
/**************************************************************************/

/*
  C++ interface to Netgen
*/

namespace netgen
{

  static constexpr int POINTINDEX_BASE = 1;
  
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

    class Ng_Edges
    {
    public:
      size_t num;
      const T_EDGE2 * ptr;
  
      size_t Size() const { return num; }
      int operator[] (size_t i) const { return ptr[i].nr; }
    };

    class Ng_Faces
    {
    public:
      size_t num;
      const T_FACE2 * ptr;
  
      size_t Size() const { return num; }
      int operator[] (size_t i) const { return ptr[i].nr; }
    };

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
    const string * mat;   // material / boundary label
    NG_ELEMENT_TYPE GetType() const { return type; }
    int GetIndex() const { return index-1; }
    Ng_Points points;      // all points
    Ng_Vertices vertices;
    Ng_Edges edges;
    Ng_Faces faces;
    Ng_Facets facets;
    bool is_curved;
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

    class Ng_Edges
    {
    public:
      size_t ned;
      const int * ptr;
  
      size_t Size() const { return ned; }
      int operator[] (size_t i) const { return ptr[i]-1; }
    };


  public:
    Ng_Vertices vertices;
    Ng_Edges edges;
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
    Ngx_Mesh(shared_ptr<Mesh> amesh = NULL);
    void LoadMesh (const string & filename);

    void LoadMesh (istream & str);
    void SaveMesh (ostream & str) const;
    void UpdateTopology ();
    void DoArchive (ngstd::Archive & archive);

    virtual ~Ngx_Mesh();

    bool Valid () { return mesh != NULL; }
    
    int GetDimension() const;
    int GetNLevels() const;

    int GetNElements (int dim) const;
    int GetNNodes (int nt) const;

    Ng_Point GetPoint (int nr) const;

    template <int DIM> 
    Ng_Element GetElement (size_t nr) const;

    template <int DIM> 
    int GetElementIndex (size_t nr) const;

    /// material/boundary label of region, template argument is co-dimension
    template <int DIM> 
    const string & GetMaterialCD (int region_nr) const;

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
    
    
    template <int DIM>
    int GetNNodes ();

    // returns domain numbers of domains next to boundary bnr -> (domin, domout)
    // 3D only
    // std::pair<int,int> GetBoundaryNeighbouringDomains (int bnr);

    void Refine (NG_REFINEMENT_TYPE reftype,
                 void (*taskmanager)(function<void(int,int)>) = &DummyTaskManager2,
                 void (*tracer)(string, bool) = &DummyTracer2);


    // Find element of point, returns local coordinates
    template <int DIM>
    int FindElementOfPoint 
    (double * p, double * lami,
     bool build_searchtrees = false, 
     int * const indices = NULL, int numind = 0) const;
    

#ifdef PARALLEL
    std::tuple<int,int*> GetDistantProcs (int nodetype, int locnum) const;
#endif

    shared_ptr<Mesh> GetMesh () const { return mesh; } 
    shared_ptr<Mesh> SelectMesh () const;
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


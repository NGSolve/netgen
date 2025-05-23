#ifndef NGINTERFACE
#define NGINTERFACE





/**************************************************************************/
/* File:   nginterface.h                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

#include "mydefs.hpp"
#include <meshing/visual_interface.hpp>

/*
  Application program interface to Netgen

*/

// max number of nodes per element
#define NG_ELEMENT_MAXPOINTS 20

// max number of nodes per surface element
#define NG_SURFACE_ELEMENT_MAXPOINTS 8

// #ifndef PARALLEL
// typedef int MPI_Comm;
// #endif
// namespace netgen { extern DLL_HEADER ngcore::NgMPI_Comm ng_comm; }


// implemented element types:
enum NG_ELEMENT_TYPE { 
  NG_PNT = 0,
  NG_SEGM = 1, NG_SEGM3 = 2,
  NG_TRIG = 10, NG_QUAD=11, NG_TRIG6 = 12, NG_QUAD6 = 13, NG_QUAD8 = 14,
  NG_TET = 20, NG_TET10 = 21, 
  NG_PYRAMID = 22, NG_PRISM = 23, NG_PRISM12 = 24, NG_PRISM15 = 27, NG_PYRAMID13 = 28,
  NG_HEX = 25, NG_HEX20 = 26, NG_HEX7 = 29
};

typedef double NG_POINT[3];  // coordinates
typedef int NG_EDGE[2];      // initial point, end point
typedef int NG_FACE[4];      // points, last one is 0 for trig


#ifdef __cplusplus
extern "C" {
#endif
  
  // load geometry from file 
  DLL_HEADER void Ng_LoadGeometry (const char * filename);

  // load netgen mesh
  DLL_HEADER void Ng_LoadMesh (const char * filename, ngcore::NgMPI_Comm comm = ngcore::NgMPI_Comm{});

  // load netgen mesh
  DLL_HEADER void Ng_LoadMeshFromString (const char * mesh_as_string);

  // space dimension (2 or 3)
  DLL_HEADER int Ng_GetDimension ();

  // number of mesh points
  DLL_HEADER int Ng_GetNP ();
  
  // number of mesh vertices (differs from GetNP for 2nd order elements)
  DLL_HEADER int Ng_GetNV ();
  
  // number of mesh elements
  DLL_HEADER int Ng_GetNE ();
  
  // number of surface triangles
  DLL_HEADER int Ng_GetNSE ();
  
  // Get Point coordinates, index from 1 .. np
  DLL_HEADER void Ng_GetPoint (int pi, double * p);
  
  // Get Element Points
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetElement (int ei, int * epi, int * np = 0);

  // Get Element Type
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetElementType (int ei);

  // Get sub-domain of element ei
  DLL_HEADER int Ng_GetElementIndex (int ei);

  DLL_HEADER void Ng_SetElementIndex(const int ei, const int index);

  // Get Material of element ei
  DLL_HEADER const char * Ng_GetElementMaterial (int ei);

  // Get Material of domain dom
  DLL_HEADER const char * Ng_GetDomainMaterial (int dom);
  
  // Get User Data
  DLL_HEADER int Ng_GetUserDataSize (char * id);
  DLL_HEADER void Ng_GetUserData (char * id, double * data);

  // Get Surface Element Points
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetSurfaceElement (int ei, int * epi, int * np = 0);

  // Get Surface Element Type
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetSurfaceElementType (int ei);

  // Get Surface Element Index
  DLL_HEADER int Ng_GetSurfaceElementIndex (int ei);

  // Get Surface Element Surface Number
  DLL_HEADER int Ng_GetSurfaceElementSurfaceNumber (int ei);
  
  // Get Surface Element Number
  DLL_HEADER int Ng_GetSurfaceElementFDNumber (int ei);

  // Get BCName for Surface Element  
  DLL_HEADER char * Ng_GetSurfaceElementBCName (int ei);
  //void Ng_GetSurfaceElementBCName (int ei, char * name);

  // Get BCName for bc-number
  DLL_HEADER char * Ng_GetBCNumBCName (int bcnr);
  //void Ng_GetBCNumBCName (int bcnr, char * name);

  // Get BCName for bc-number of co dim 2
  DLL_HEADER char * Ng_GetCD2NumCD2Name (int cd2nr);

  // Get normal vector of surface element node
  // DLL_HEADER void Ng_GetNormalVector (int sei, int locpi, double * nv);     
  

  DLL_HEADER void Ng_SetPointSearchStartElement(int el);
  
  // Find element of point, returns local coordinates
  DLL_HEADER int Ng_FindElementOfPoint (double * p, double * lami,
                                        int build_searchtrees = 0, 
                                        const int * const indices = NULL, const int numind = 0);
  
  // Find surface element of point, returns local coordinates
  DLL_HEADER int Ng_FindSurfaceElementOfPoint (double * p, double * lami,
                                               int build_searchtrees = 0, 
                                               const int * const indices = NULL, const int numind = 0);
  

  // is element ei curved ?
  DLL_HEADER int Ng_IsElementCurved (int ei);
  // is element sei curved ?
  DLL_HEADER int Ng_IsSurfaceElementCurved (int sei);

  /// Curved Elements:
  /// xi..local coordinates
  /// x ..global coordinates
  /// dxdxi...D x D Jacobian matrix (row major storage)
  DLL_HEADER void Ng_GetElementTransformation (int ei, const double * xi, 
                                               double * x, double * dxdxi);

  
  /// buffer must be at least 100 doubles, alignment of double
  DLL_HEADER void Ng_GetBufferedElementTransformation (int ei, const double * xi, 
                                                       double * x, double * dxdxi,
                                                       void * buffer, int buffervalid);
  


  /// Curved Elements:
  /// xi..local coordinates
  /// x ..global coordinates
  /// dxdxi...D x D-1 Jacobian matrix (row major storage)
  /// curved ...is element curved ?
  DLL_HEADER void Ng_GetSurfaceElementTransformation (int sei, const double * xi, 
                                                      double * x, double * dxdxi);
  
  /// Curved Elements:
  /// xi..local coordinates
  /// sxi..step xi
  /// x ..global coordinates
  /// dxdxi...D x D Jacobian matrix (row major storage)
  DLL_HEADER void Ng_GetMultiElementTransformation (int ei, int n,
                                                    const double * xi, size_t sxi,
                                                    double * x, size_t sx,
                                                    double * dxdxi, size_t sdxdxi);

  
  
  DLL_HEADER int Ng_GetSegmentIndex (int elnr);
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetSegment (int elnr, int * epi, int * np = 0);




  // Mark element for refinement
  DLL_HEADER void Ng_SetRefinementFlag (int ei, int flag);
  DLL_HEADER void Ng_SetSurfaceRefinementFlag (int sei, int flag);

  // Do local refinement
  enum NG_REFINEMENT_TYPE { NG_REFINE_H = 0, NG_REFINE_P = 1, NG_REFINE_HP = 2 };
  DLL_HEADER void Ng_Refine (NG_REFINEMENT_TYPE reftype);

  // Use second order elements
  DLL_HEADER void Ng_SecondOrder ();
  DLL_HEADER void Ng_HighOrder (int order, bool rational = false);
  //void Ng_HPRefinement (int levels, double parameter = 0.125);
  DLL_HEADER void Ng_HPRefinement (int levels, double parameter = 0.125,
                                   bool setorders = true,bool ref_level = false);
  // void Ng_HPRefinement (int levels);
  // void Ng_HPRefinement (int levels, double parameter);


  // Topology and coordinate information of master element:

  DLL_HEADER int Ng_ME_GetNVertices (NG_ELEMENT_TYPE et);
  DLL_HEADER int Ng_ME_GetNEdges (NG_ELEMENT_TYPE et);
  DLL_HEADER int Ng_ME_GetNFaces (NG_ELEMENT_TYPE et);

  DLL_HEADER const NG_POINT * Ng_ME_GetVertices (NG_ELEMENT_TYPE et);
  DLL_HEADER const NG_EDGE * Ng_ME_GetEdges (NG_ELEMENT_TYPE et);
  DLL_HEADER const NG_FACE * Ng_ME_GetFaces (NG_ELEMENT_TYPE et);

  DLL_HEADER void Ng_UpdateTopology();

  DLL_HEADER int Ng_GetNEdges();
  DLL_HEADER int Ng_GetNFaces();

  [[deprecated("orientation is not supported anymore")]]                          
  DLL_HEADER int Ng_GetElement_Edges (int elnr, int * edges, int * orient = 0);
  // [[deprecated("orientation is not supported anymore")]]                            
  DLL_HEADER int Ng_GetElement_Faces (int elnr, int * faces, int * orient = 0);
  [[deprecated("orientation is not supported anymore")]]                          
  DLL_HEADER int Ng_GetSurfaceElement_Edges (int selnr, int * edges, int * orient = 0);
  // [[deprecated("orientation is not supported anymore")]]                            
  DLL_HEADER int Ng_GetSurfaceElement_Face (int selnr, int * orient = 0);

  DLL_HEADER void Ng_GetSurfaceElementNeighbouringDomains(const int selnr, int & in, int & out);
       
  DLL_HEADER int Ng_GetFace_Vertices (int fnr, int * vert);
  DLL_HEADER void Ng_GetEdge_Vertices (int ednr, int * vert);
  DLL_HEADER int Ng_GetFace_Edges (int fnr, int * edge);

  DLL_HEADER int Ng_GetNVertexElements (int vnr);
  DLL_HEADER void Ng_GetVertexElements (int vnr, int * els);

  DLL_HEADER int Ng_GetElementOrder (int enr);
  DLL_HEADER void Ng_GetElementOrders (int enr, int * ox, int * oy, int * oz);

  DLL_HEADER void Ng_SetElementOrder (int enr, int order);
  DLL_HEADER void Ng_SetElementOrders (int enr, int ox, int oy, int oz);

  DLL_HEADER int Ng_GetSurfaceElementOrder (int enr);
  DLL_HEADER void Ng_GetSurfaceElementOrders (int enr, int * ox, int * oy);

  DLL_HEADER void Ng_SetSurfaceElementOrder (int enr, int order);
  DLL_HEADER void Ng_SetSurfaceElementOrders (int enr, int ox, int oy);

  // Multilevel functions:

  // number of levels:
  DLL_HEADER int Ng_GetNLevels ();
  // get two parent nodes (indeed vertices !) of node ni
  DLL_HEADER void Ng_GetParentNodes (int ni, int * parents);

  // get parent element (first child has always same number)
  DLL_HEADER int Ng_GetParentElement (int ei);

  // get parent surface element (first child has always same number)
  DLL_HEADER int Ng_GetParentSElement (int ei);

  // representant of anisotropic cluster
  DLL_HEADER int Ng_GetClusterRepVertex (int vi);
  DLL_HEADER int Ng_GetClusterRepEdge (int edi);
  DLL_HEADER int Ng_GetClusterRepFace (int fai);
  DLL_HEADER int Ng_GetClusterRepElement (int eli);


  void Ng_SurfaceElementTransformation (int eli, double x, double y, 
					double * p3d, double * jacobian);

#ifdef PARALLEL

  // the following functions are 0-base  !!

  // number on distant processor 
  // returns pairs  (dist_proc, num_on_dist_proc)
  int NgPar_GetDistantNodeNums ( int nodetype, int locnum, int * pnums );
  int NgPar_GetNDistantNodeNums ( int nodetype, int locnum );
  
  DLL_HEADER int NgPar_GetGlobalNodeNum (int nodetype, int locnum);

#endif
  
  namespace netgen {
  // #include "../visualization/soldata.hpp"
    class SolutionData;
    class MouseEventHandler;
    class UserVisualizationObject;
  }

  enum Ng_SolutionType
  { NG_SOLUTION_NODAL = 1, 
    NG_SOLUTION_ELEMENT = 2, 
    NG_SOLUTION_SURFACE_ELEMENT = 3, 
    NG_SOLUTION_NONCONTINUOUS = 4,
    NG_SOLUTION_SURFACE_NONCONTINUOUS = 5,
    NG_SOLUTION_VIRTUAL_FUNCTION = 6,
    NG_SOLUTION_MARKED_ELEMENTS = 10,
    NG_SOLUTION_ELEMENT_ORDER = 11
  };
  
  struct Ng_SolutionData
  {
    std::string name;      // name of gridfunction
    std::string title = ""; // name of gridfunction ( printed on top of window )
    std::string number_format = "%.3e"; // printf-style string to format colormap values
    std::string unit = ""; // string to append to last number in colormap (ASCII only)
    double * data;    // solution values
    int components;   // relevant (double) components in solution vector
    int dist;         // # doubles per entry alignment! 
    int iscomplex;    // complex vector ? 
    bool draw_surface;
    bool draw_volume;
    std::shared_ptr<netgen::BitArray> draw_surfaces, draw_volumes;
    int order;        // order of elements, only partially supported 
    Ng_SolutionType soltype;  // type of solution function
    netgen::SolutionData * solclass;
  };
  
  // initialize solution data with default arguments
  DLL_HEADER void Ng_InitSolutionData (Ng_SolutionData * soldata);
  // set solution data
  DLL_HEADER void Ng_SetSolutionData (Ng_SolutionData * soldata);
  /// delete gridfunctions
  DLL_HEADER void Ng_ClearSolutionData();
  // redraw 
  DLL_HEADER void Ng_Redraw(bool blocking = false);
  ///
  DLL_HEADER void Ng_TclCmd(std::string cmd);
  ///
  DLL_HEADER void Ng_SetMouseEventHandler (netgen::MouseEventHandler * handler);
  ///
  DLL_HEADER void Ng_SetUserVisualizationObject (netgen::UserVisualizationObject * vis);
  //
  DLL_HEADER void Ng_SetVisualizationParameter (const char * name, 
                                                const char * value);
  

  // number of periodic vertices  
  DLL_HEADER int Ng_GetNPeriodicVertices (int idnr);
  // pairs should be an integer array of 2*npairs
  DLL_HEADER void Ng_GetPeriodicVertices (int idnr, int * pairs); 

  // number of periodic edges  
  DLL_HEADER int Ng_GetNPeriodicEdges (int idnr);
  // pairs should be an integer array of 2*npairs
  DLL_HEADER void Ng_GetPeriodicEdges (int idnr, int * pairs); 

  DLL_HEADER void RunParallel ( void * (*fun)(void *), void * in);

  DLL_HEADER void Ng_PushStatus (const char * str);
  DLL_HEADER void Ng_PopStatus ();
  DLL_HEADER void Ng_SetThreadPercentage (double percent);
  DLL_HEADER void Ng_GetStatus (char ** str, double & percent);

  DLL_HEADER void Ng_SetTerminate(void);
  DLL_HEADER void Ng_UnSetTerminate(void);
  DLL_HEADER int Ng_ShouldTerminate(void);
  DLL_HEADER void Ng_SetRunning(int flag);
  DLL_HEADER int Ng_IsRunning();
  
  //// added by Roman Stainko ....
  DLL_HEADER int Ng_GetVertex_Elements( int vnr, int* elems);
  DLL_HEADER int Ng_GetVertex_SurfaceElements( int vnr, int* elems );
  DLL_HEADER int Ng_GetVertex_NElements( int vnr );
  DLL_HEADER int Ng_GetVertex_NSurfaceElements( int vnr );


#ifdef SOCKETS
  int Ng_SocketClientOpen( const int port, const char * host );
  void Ng_SocketClientWrite( const char * write, char ** reply);
  void Ng_SocketClientClose ( void );
  void Ng_SocketClientGetServerHost ( const int number, char ** host );
  void Ng_SocketClientGetServerPort ( const int number, int * port );
  void Ng_SocketClientGetServerClientID ( const int number, int * id );
#endif

  DLL_HEADER void Ng_InitPointCurve(double red, double green, double blue);
  DLL_HEADER void Ng_AddPointCurvePoint(const double * point);

  
  // #ifdef PARALLEL
  // void Ng_SetElementPartition ( int elnr, int part );
  // int  Ng_GetElementPartition ( int elnr );
  // #endif

  DLL_HEADER void Ng_SaveMesh ( const char * meshfile );
  DLL_HEADER void Ng_Bisect ( const char * refinementfile );

  // if qualityloss is not equal to NULL at input, a (1-based) list of qualitylosses (due to projection)
  // is saved in *qualityloss, its size is the return value
  DLL_HEADER int Ng_Bisect_WithInfo ( const char * refinementfile, double ** qualityloss);

  typedef void * Ng_Mesh;
  DLL_HEADER Ng_Mesh Ng_SelectMesh (Ng_Mesh mesh);

  DLL_HEADER void Ng_GetArgs (int & argc, char ** &argv);


#ifdef __cplusplus
}
#endif

#endif






/*
  The new node interface ...
  it is 0-based !
 */

extern "C" {
  
  /*
    number of nodes of type nt
    nt = 0 is Vertex
    nt = 1 is Edge
    nt = 2 is Face
    nt = 3 is Cell
   */
  DLL_HEADER int Ng_GetNNodes (int nt);

  /*
    closure nodes of node (nt, nodenr):
    nodeset is bit-coded, bit 0 includes Vertices, bit 1 edges, etc
    E.g., nodeset = 6 includes edge and face nodes
    nodes consists of pairs of integers (nodetype, nodenr) 
    return value is number of nodes
   */
  DLL_HEADER int Ng_GetClosureNodes (int nt, int nodenr, int nodeset, int * nodes);

  
  /*
    number of dim-dimensional elements 
    dim = 3 ... volume elements
    dim = 2 ... surface elements
    dim = 1 ... segments
    dim = 0 ... not available
  */
  DLL_HEADER int Ng_GetNElements (int dim);

  /*
    closure nodes of dim-dimensional element elmentnr:
    nodeset is bit-coded, bit 0 includes Vertices, bit 1 edges, etc
    E.g., nodeset = 6 includes edge and face nodes
    nodes consists of pairs of integers (nodetype, nodenr) 
    return value is number of nodes
   */
  DLL_HEADER int Ng_GetElementClosureNodes (int dim, int elementnr, int nodeset, int * nodes);
}





#ifdef __cplusplus
#include <iostream>
namespace ngcore
{
  NGCORE_API extern int printmessage_importance;
}

#endif


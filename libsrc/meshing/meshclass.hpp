#ifndef NETGEN_MESHCLASS_HPP
#define NETGEN_MESHCLASS_HPP

/**************************************************************************/
/* File:   meshclass.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

/*
  The mesh class
*/

#include<filesystem>

#include <gprim/adtree.hpp>
#include <gprim/transform3d.hpp>

#include "meshtype.hpp"
#include "localh.hpp"
#include "topology.hpp"

namespace netgen
{
  class NetgenGeometry;
  using namespace std;

  static constexpr int  NG_MPI_TAG_MESH = 210;
  

  enum resthtype { RESTRICTH_FACE, RESTRICTH_EDGE, 
		   RESTRICTH_SURFACEELEMENT, RESTRICTH_POINT, RESTRICTH_SEGMENT };

  class HPRefElement;
  class CurvedElements;
  class AnisotropicClusters;
  class ParallelMeshTopology;

  class MarkedTet;
  class MarkedPrism;
  class MarkedIdentification;
  class MarkedTri;
  class MarkedQuad;

  typedef Array<MarkedTet> T_MTETS;
  typedef NgArray<MarkedPrism> T_MPRISMS;
  typedef NgArray<MarkedIdentification> T_MIDS;
  typedef NgArray<MarkedTri> T_MTRIS;
  typedef NgArray<MarkedQuad> T_MQUADS;

  struct BisectionInfo
  {
    unique_ptr<T_MTETS> mtets;
    unique_ptr<T_MPRISMS> mprisms;
    unique_ptr<T_MIDS> mids;
    unique_ptr<T_MTRIS> mtris;
    unique_ptr<T_MQUADS> mquads;

    BisectionInfo();
    ~BisectionInfo();
  };
  
  /// 2d/3d mesh
  class Mesh
  {
  public:
    // typedef Array<MeshPoint, PointIndex> T_POINTS;
    typedef netgen::T_POINTS T_POINTS;

  private:
    /// point coordinates
    T_POINTS points;

    // The communicator for this mesh. Just a dummy if compiled without MPI.  
    NgMPI_Comm comm;
    
    /// line-segments at edges
    Array<Segment, SegmentIndex> segments;
    /// surface elements, 2d-inner elements
    Array<Element2d, SurfaceElementIndex> surfelements;
    /// volume elements
    Array<Element, ElementIndex> volelements;
    /// points will be fixed forever
    Array<PointIndex> lockedpoints;


    /// surface indices at boundary nodes
    // TABLE<int,PointIndex::BASE> surfacesonnode;
    /// boundary edges  (1..normal bedge, 2..segment)
    unique_ptr<INDEX_2_CLOSED_HASHTABLE<int>> boundaryedges;
    ///
    unique_ptr<INDEX_2_CLOSED_HASHTABLE<int>> segmentht;
    ///
    unique_ptr<INDEX_3_CLOSED_HASHTABLE<int>> surfelementht;
    unique_ptr<INDEX_3_CLOSED_HASHTABLE<int>> illegal_trigs;

    /// faces of rest-solid
    NgArray<Element2d> openelements;
    /// open segments for surface meshing
    NgArray<Segment> opensegments;

    Array<int> tets_in_qualclass;



    /**
       Representation of local mesh-size h (one function per mesh layer)
    */
    Array<shared_ptr<LocalH>> lochfunc;
    ///
    double hglob;
    ///
    double hmin;
    ///
    NgArray<double> maxhdomain;
  
    /**
       the face-index of the surface element maps into
       this table.
    */
    Array<FaceDescriptor> facedecoding;

  
    /**
       the edge-index of the line element maps into
       this table.
    */
    NgArray<EdgeDescriptor> edgedecoding;

    /// sub-domain materials 
    NgArray<string*> materials;

    /// labels for boundary conditions
    NgArray<string*> bcnames;

    /// labels for co dim 2 bboundary conditions
    NgArray<string*> cd2names;

    /// labels for co dim 3 bbboundary conditions
    NgArray<string*> cd3names;

    /// Periodic surface, close surface, etc. identifications
    unique_ptr<Identifications> ident;


    /// number of vertices (if < 0, use np)
    int numvertices;

    /// geometric search tree for interval intersection search
    unique_ptr<BoxTree<3>> elementsearchtree;
    /// time stamp for tree
    mutable int elementsearchtreets;

    /// element -> face, element -> edge etc ...
    MeshTopology topology;
    /// methods for high order elements
    unique_ptr<CurvedElements> curvedelems;

    /// nodes identified by close points 
    unique_ptr<AnisotropicClusters> clusters;

    /// space dimension (2 or 3)
    int dimension;
  
    /// changed by every minor modification (addpoint, ...)
    int timestamp;
    /// changed after finishing global algorithm (improve, ...)
    int majortimestamp;

    /// mesh access semaphores.
    NgMutex mutex;
    /// mesh access semaphores.
    NgMutex majormutex;

    SymbolTable< NgArray<int>* > userdata_int;
    SymbolTable< NgArray<double>* > userdata_double;


    mutable NgArray< Point3d > pointcurves;
    mutable NgArray<int> pointcurves_startpoint;
    mutable NgArray<double> pointcurves_red,pointcurves_green,pointcurves_blue;


    /// start element for point search (GetElementOfPoint)
    mutable int ps_startelement;


#ifdef PARALLEL
    /// connection to parallel meshes
    unique_ptr<ParallelMeshTopology> paralleltop;
#endif

    
    shared_ptr<NetgenGeometry> geometry;


  public:
    DLL_HEADER void BuildBoundaryEdges(bool rebuild=true);

    DLL_HEADER bool PointContainedIn2DElement(const Point3d & p,
				   double lami[3],
				   const int element,
				   bool consider3D = false) const;
    DLL_HEADER bool PointContainedIn3DElement(const Point3d & p,
				   double lami[3],
				   const int element) const;
    DLL_HEADER bool PointContainedIn3DElementOld(const Point3d & p,
				      double lami[3],
				      const int element) const;

  public:
    Signal<> updateSignal;
    BisectionInfo bisectioninfo;

    // store coarse mesh before hp-refinement
    unique_ptr<NgArray<HPRefElement>> hpelements;
    unique_ptr<Mesh> coarsemesh;
  
  
    /// number of refinement levels
    // int mglevels;
    // number of vertices on each refinement level:
    NgArray<size_t> level_nv;
    /// refinement hierarchy
    NgArray<PointIndices<2>,PointIndex::BASE> mlbetweennodes;
    /// parent element of volume element
    NgArray<int> mlparentelement;
    /// parent element of surface element
    NgArray<int> mlparentsurfaceelement;



    ///
    DLL_HEADER Mesh();
    ///
    DLL_HEADER ~Mesh();

    DLL_HEADER Mesh & operator= (const Mesh & mesh2);
  
    ///
    DLL_HEADER void DeleteMesh();
  
    ///
    void ClearSurfaceElements();

    ///
    DLL_HEADER void ClearVolumeElements()
    {
      volelements.SetSize(0); 
      timestamp = NextTimeStamp();
    }

    ///
    DLL_HEADER void ClearSegments()
    { 
      segments.SetSize(0); 
      timestamp = NextTimeStamp();
    }
    
    ///
    bool TestOk () const;

    void SetAllocSize(int nnodes, int nsegs, int nsel, int nel);
    

    DLL_HEADER PointIndex AddPoint (const Point3d & p, int layer = 1);
    DLL_HEADER PointIndex AddPoint (const Point3d & p, int layer, POINTTYPE type);

    auto GetNP () const { return points.Size(); }

    // [[deprecated("Use Point(PointIndex) instead of int !")]]        
    MeshPoint & Point(int i)
    {
      // return points.Elem(i);
      return Point (PointIndex(i+PointIndex::BASE-1));
    }
    MeshPoint & Point(PointIndex pi) { return points[pi]; }
    // [[deprecated("Use Point(PointIndex) instead of int !")]]            
    const MeshPoint & Point(int i) const
    {
      // return points.Get(i);
      return Point (PointIndex(i+PointIndex::BASE-1));      
    }
    const MeshPoint & Point(PointIndex pi) const { return points[pi]; }

    const MeshPoint & operator[] (PointIndex pi) const { return points[pi]; }
    MeshPoint & operator[] (PointIndex pi) { return points[pi]; }

    const T_POINTS & Points() const { return points; }
    T_POINTS & Points() { return points; }


    DLL_HEADER SegmentIndex AddSegment (const Segment & s);
    void DeleteSegment (int segnr)
    {
      segments[segnr-1][0].Invalidate();
      segments[segnr-1][1].Invalidate();
    }
    /*
    void FullDeleteSegment (int segnr)  // von wem ist das ???
    {
      segments.Delete(segnr-PointIndex::BASE);
    }
    */

    int GetNSeg () const { return segments.Size(); }
    // [[deprecated("Use LineSegment(SegmentIndex) instead of int !")]]                
    Segment & LineSegment(int i) { return segments[i-1]; }
    // [[deprecated("Use LineSegment(SegmentIndex) instead of int !")]]                    
    const Segment & LineSegment(int i) const { return segments[i-1]; }

    Segment & LineSegment(SegmentIndex si) { return segments[si]; }
    const Segment & LineSegment(SegmentIndex si) const { return segments[si]; }
    const Segment & operator[] (SegmentIndex si) const { return segments[si]; }
    Segment & operator[] (SegmentIndex si) { return segments[si]; }

    const auto & LineSegments() const { return segments; }
    auto & LineSegments() { return segments; }
    
    Array<Element0d> pointelements;  // only via python interface

    DLL_HEADER SurfaceElementIndex AddSurfaceElement (const Element2d & el);
    // write to pre-allocated container, thread-safe
    DLL_HEADER void SetSurfaceElement (SurfaceElementIndex sei, const Element2d & el);
    
    [[deprecated("Use Delete(SurfaceElementIndex) instead of int !")]]
    void DeleteSurfaceElement (int eli)
    {
      /*
      surfelements.Elem(eli).Delete();
      surfelements.Elem(eli).PNum(1).Invalidate();
      surfelements.Elem(eli).PNum(2).Invalidate();
      surfelements.Elem(eli).PNum(3).Invalidate();
      */
      surfelements[eli-1].Delete();
      /*
      surfelements[eli-1].PNum(1).Invalidate();
      surfelements[eli-1].PNum(2).Invalidate();
      surfelements[eli-1].PNum(3).Invalidate();
      */
      timestamp = NextTimeStamp();
    }

    [[deprecated("Use Delete(SurfaceElementIndex) instead !")]]        
    void DeleteSurfaceElement (SurfaceElementIndex eli)
    {
      // for (auto & p : surfelements[eli].PNums()) p.Invalidate();
      surfelements[eli].Delete();
      timestamp = NextTimeStamp();
    }
    
    void Delete (SurfaceElementIndex eli)
    {
      // for (auto & p : surfelements[eli].PNums()) p.Invalidate();
      surfelements[eli].Delete();
      timestamp = NextTimeStamp();
    }

    auto GetNSE () const { return surfelements.Size(); }

    // [[deprecated("Use SurfaceElement(SurfaceElementIndex) instead of int !")]]    
    Element2d & SurfaceElement(int i) { return surfelements[i-1]; }
    // [[deprecated("Use SurfaceElement(SurfaceElementIndex) instead of int !")]]        
    const Element2d & SurfaceElement(int i) const { return surfelements[i-1]; }
    // [[deprecated("Use mesh[](SurfaceElementIndex) instead !")]]
    Element2d & SurfaceElement(SurfaceElementIndex i) { return surfelements[i]; }
    // [[deprecated("Use mesh[](SurfaceElementIndex) instead !")]]
    const Element2d & SurfaceElement(SurfaceElementIndex i) const { return surfelements[i]; }

    const Element2d & operator[] (SurfaceElementIndex ei) const
    { return surfelements[ei]; }
    Element2d & operator[] (SurfaceElementIndex ei)
    { return surfelements[ei]; }

    const auto & SurfaceElements() const { return surfelements; }
    auto & SurfaceElements() { return surfelements; }

  
    DLL_HEADER void RebuildSurfaceElementLists ();
    DLL_HEADER void GetSurfaceElementsOfFace (int facenr, Array<SurfaceElementIndex> & sei) const;

    DLL_HEADER ElementIndex AddVolumeElement (const Element & el);
    // write to pre-allocated container, thread-safe
    DLL_HEADER void SetVolumeElement (ElementIndex sei, const Element & el);

    auto GetNE () const { return volelements.Size(); }

    // [[deprecated("Use VolumeElement(ElementIndex) instead of int !")]]    
    Element & VolumeElement(int i) { return volelements[i-1]; }
    // [[deprecated("Use VolumeElement(ElementIndex) instead of int !")]]        
    const Element & VolumeElement(int i) const { return volelements[i-1]; }
    // [[deprecated("Use mesh[](VolumeElementIndex) instead !")]]
    Element & VolumeElement(ElementIndex i) { return volelements[i]; }
    // [[deprecated("Use mesh[](VolumeElementIndex) instead !")]]
    const Element & VolumeElement(ElementIndex i) const { return volelements[i]; }

    const Element & operator[] (ElementIndex ei) const { return volelements[ei]; }
    Element & operator[] (ElementIndex ei) { return volelements[ei]; }

    ELEMENTTYPE ElementType (ElementIndex i) const 
    { return (volelements[i].Flags().fixed) ? FIXEDELEMENT : FREEELEMENT; }

    const auto & VolumeElements() const { return volelements; }
    auto & VolumeElements() { return volelements; }

    ///
    DLL_HEADER double ElementError (int eli, const MeshingParameters & mp) const;

    /// 
    DLL_HEADER void AddLockedPoint (PointIndex pi);
    ///
    void ClearLockedPoints ();

    const auto & LockedPoints() const { return lockedpoints; }

    /// Returns number of domains
    DLL_HEADER int GetNDomains() const;
    ///
    int GetDimension() const { return dimension; }
    DLL_HEADER void SetDimension (int dim); //  { dimension = dim; }

    /// sets internal tables
    DLL_HEADER void CalcSurfacesOfNode ();

    /// additional (temporarily) fix points 
    void FixPoints (const NgBitArray & fixpoints);

    /**
       finds elements without neighbour and
       boundary elements without inner element.
       Results are stored in openelements.
       if dom == 0, all sub-domains, else subdomain dom */
    DLL_HEADER void FindOpenElements (int dom = 0);

  
    /**
       finds segments without surface element,
       and surface elements without neighbours.
       store in opensegmentsy
    */
    DLL_HEADER void FindOpenSegments (int surfnr = 0);
    /**
       remove one layer of surface elements
    */
    DLL_HEADER void RemoveOneLayerSurfaceElements ();


    int GetNOpenSegments () { return opensegments.Size(); }
    const Segment & GetOpenSegment (int nr) { return opensegments.Get(nr); }
  
    /**
       Checks overlap of boundary
       return == 1, iff overlap
    */
    DLL_HEADER int CheckOverlappingBoundary ();
    /**
       Checks consistent boundary
       return == 0, everything ok
    */
    DLL_HEADER int CheckConsistentBoundary () const;

    /*
      checks element orientation
    */
    DLL_HEADER int CheckVolumeMesh () const;


    /**
       finds average h of surface surfnr if surfnr > 0,
       else of all surfaces.
    */
    DLL_HEADER double AverageH (int surfnr = 0) const;
    /// Calculates localh 
    DLL_HEADER void CalcLocalH (double grading, int layer=1);
    ///
    DLL_HEADER void SetLocalH (netgen::Point<3> pmin, netgen::Point<3> pmax, double grading, int layer=1);
    ///
    DLL_HEADER void RestrictLocalH (const Point3d & p, double hloc, int layer=1);
    ///
    DLL_HEADER void RestrictLocalHLine (const Point3d & p1, const Point3d & p2, 
			     double hloc, int layer=1);
    /// number of elements per radius
    DLL_HEADER void CalcLocalHFromSurfaceCurvature(double grading, double elperr, int layer=1);
    ///
    DLL_HEADER void CalcLocalHFromPointDistances(double grading, int layer=1);
    ///
    DLL_HEADER void RestrictLocalH (resthtype rht, int nr, double loch);
    ///
    DLL_HEADER void LoadLocalMeshSize (const filesystem::path & meshsizefilename);
    ///
    DLL_HEADER void SetGlobalH (double h);
    ///
	DLL_HEADER void SetMinimalH (double h);
    ///
	DLL_HEADER double MaxHDomain (int dom) const;
    ///
	DLL_HEADER void SetMaxHDomain (const NgArray<double> & mhd);
    ///
    DLL_HEADER double GetH (const Point3d & p, int layer=1) const;
    DLL_HEADER double GetH (PointIndex pi) const { return GetH(points[pi], points[pi].GetLayer()); }
    ///
    double GetMinH (const Point3d & pmin, const Point3d & pmax, int layer=1);
    ///
    bool HasLocalHFunction (int layer=1) { return lochfunc[layer-1] != nullptr; }
    ///
    LocalH & LocalHFunction (int layer=1) { return * lochfunc[layer-1]; }

    shared_ptr<LocalH> GetLocalH(int layer=1) const
    {
      if(lochfunc.Size() == 1)
        return lochfunc[0];
      return lochfunc[layer-1];
    }
    DLL_HEADER void SetLocalH(shared_ptr<LocalH> loch, int layer=1);

    ///
    bool LocalHFunctionGenerated(int layer=1) const { return (lochfunc[layer-1] != NULL); }

    /// Find bounding box
    DLL_HEADER void GetBox (Point3d & pmin, Point3d & pmax, int dom = -1) const;

    /// Find bounding box of points of typ ptyp or less
    DLL_HEADER void GetBox (Point3d & pmin, Point3d & pmax, POINTTYPE ptyp ) const;

    ///
    int GetNOpenElements() const
    { return openelements.Size(); }
    ///
    const Element2d & OpenElement(int i) const
    { return openelements.Get(i); }

    auto & OpenElements() const { return openelements; }

    auto & OpenElements() { return openelements; }
    
    /// are also quads open elements
    bool HasOpenQuads () const;

    /// split into connected pieces
	DLL_HEADER void SplitIntoParts ();

    /// 
	DLL_HEADER void SplitSeparatedFaces ();

    /// Refines mesh and projects points to true surface
    // void Refine (int levels, const CSGeometry * geom);

    void ZRefine(const string& name, const Array<double>& slices);
    
    bool BoundaryEdge (PointIndex pi1, PointIndex pi2) const
    {
      if(!boundaryedges)
	const_cast<Mesh *>(this)->BuildBoundaryEdges();

      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return boundaryedges->Used (i2);
    }

    void DeleteBoundaryEdges ()
    {
        boundaryedges = nullptr;
    }

    bool IsSegment (PointIndex pi1, PointIndex pi2) const
    {
      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return segmentht->Used (i2);
    }

    SegmentIndex SegmentNr (PointIndex pi1, PointIndex pi2) const
    {
      INDEX_2 i2 (pi1, pi2);
      i2.Sort();
      return segmentht->Get (i2);
    }


    /**
       Remove unused points. etc.
    */
    DLL_HEADER void Compress ();

    /// first vertex has lowest index
    void OrderElements(); 

    ///
	DLL_HEADER void Save (ostream & outfile) const;
    ///
	DLL_HEADER void Load (istream & infile);
    ///
	DLL_HEADER void Merge (istream & infile, const int surfindex_offset = 0);
    ///
	DLL_HEADER void Save (const filesystem::path & filename) const;
    ///
	DLL_HEADER void Load (const filesystem::path & filename);
    ///
	DLL_HEADER void Merge (const filesystem::path & filename, const int surfindex_offset = 0);


    DLL_HEADER void DoArchive (Archive & archive);
    ///
	DLL_HEADER void ImproveMesh (const MeshingParameters & mp, OPTIMIZEGOAL goal = OPT_QUALITY);

    ///
    void ImproveMeshJacobian (const MeshingParameters & mp, OPTIMIZEGOAL goal = OPT_QUALITY, const NgBitArray * usepoint = NULL);
    ///
    void ImproveMeshJacobianOnSurface (const MeshingParameters & mp,
				       const NgBitArray & usepoint, 
				       const NgArray< Vec<3>* > & nv,
				       OPTIMIZEGOAL goal = OPT_QUALITY,
				       const NgArray< NgArray<int,PointIndex::BASE>* > * idmaps = NULL);
    /**
       free nodes in environment of openelements 
       for optimiztion
    */
    void FreeOpenElementsEnvironment (int layers);


    DLL_HEADER double CalcTotalBad (const MeshingParameters & mp);
    FlatArray<int> GetQualityHistogram() { return tets_in_qualclass; }

    ///
    bool LegalTet (Element & el) const
    {
      if (el.IllegalValid())
	return !el.Illegal();
      return LegalTet2 (el);
    }
    ///
    bool LegalTet2 (Element & el) const;


    ///
    // Find trigs with same vertices
    // return: number of illegal trigs
    int FindIllegalTrigs ();

    bool LegalTrig (const Element2d & el) const;
    /**
       if values non-null, return values in 4-double array:
       triangle angles min/max, tetangles min/max
       if null, output results on cout
    */
	DLL_HEADER void CalcMinMaxAngle (double badellimit, double * retvalues = NULL);

    /*
      Marks elements which are dangerous to refine
      return: number of illegal elements
    */
	DLL_HEADER int MarkIllegalElements ();

    /// orient surface mesh, for one sub-domain only
	DLL_HEADER void SurfaceMeshOrientation ();

    /// convert mixed element mesh to tet-mesh
	DLL_HEADER void Split2Tets();


    /// build box-search tree
    DLL_HEADER void BuildElementSearchTree ();

    void SetPointSearchStartElement(const int el) const {ps_startelement = el;}

    /// gives element of point, barycentric coordinates
    DLL_HEADER int GetElementOfPoint (const netgen::Point<3> & p,
			   double * lami,
			   bool build_searchtree = 0,
			   const int index = -1,
			   const bool allowindex = true) const;
    DLL_HEADER int GetElementOfPoint (const netgen::Point<3> & p,
			   double * lami,
			   const NgArray<int> * const indices,
			   bool build_searchtree = 0,
			   const bool allowindex = true) const;
    DLL_HEADER int GetSurfaceElementOfPoint (const netgen::Point<3> & p,
				  double * lami,
				  bool build_searchtree = 0,
				  const int index = -1,
				  const bool allowindex = true) const;
    DLL_HEADER int GetSurfaceElementOfPoint (const netgen::Point<3> & p,
				  double * lami,
				  const NgArray<int> * const indices,
				  bool build_searchtree = 0,
				  const bool allowindex = true) const;

    /// give list of vol elements which are int the box(p1,p2)
    void GetIntersectingVolEls(const Point3d& p1, const Point3d& p2, 
			       NgArray<int> & locels) const;

    ///
    int AddFaceDescriptor(const FaceDescriptor& fd)
    { facedecoding.Append(fd); return facedecoding.Size(); }

    int AddEdgeDescriptor(const EdgeDescriptor & fd)
    { edgedecoding.Append(fd); return edgedecoding.Size() - 1; }

    auto & GetCommunicator() const { return this->comm; }
    void SetCommunicator(NgMPI_Comm acomm);
    
    DLL_HEADER void SplitFacesByAdjacentDomains();
    DLL_HEADER shared_ptr<Mesh> GetSubMesh(string domains="", string faces="") const;

    ///
    DLL_HEADER void SetMaterial (int domnr, const string & mat);
    ///
    DLL_HEADER const string & GetMaterial (int domnr) const;
    DLL_HEADER static string defaultmat;
    const string * GetMaterialPtr (int domnr) const // 1-based
    {
      return domnr <= materials.Size() ? materials.Get(domnr) : &defaultmat;
    }
    
    DLL_HEADER void SetNBCNames ( int nbcn );

    DLL_HEADER void SetBCName ( int bcnr, const string & abcname );

    DLL_HEADER const string & GetBCName ( int bcnr ) const;

    DLL_HEADER void SetNCD2Names (int ncd2n);
    DLL_HEADER void SetCD2Name (int cd2nr, const string & abcname);

    DLL_HEADER const string & GetCD2Name (int cd2nr ) const;
    DLL_HEADER static string cd2_default_name;
    string * GetCD2NamePtr (int cd2nr ) const
    {
      if (cd2nr < cd2names.Size() && cd2names[cd2nr]) return cd2names[cd2nr];
      return &cd2_default_name;
    }
    size_t GetNCD2Names() const { return cd2names.Size(); }

    DLL_HEADER void SetNCD3Names (int ncd3n);
    DLL_HEADER void SetCD3Name (int cd3nr, const string & abcname);
    DLL_HEADER int AddCD3Name (const string & aname);

    DLL_HEADER const string & GetCD3Name (int cd3nr ) const;
    DLL_HEADER static string cd3_default_name;
    string * GetCD3NamePtr (int cd3nr ) const
    {
      if (cd3nr < cd3names.Size() && cd3names[cd3nr]) return cd3names[cd3nr];
      return &cd3_default_name;
    }
    size_t GetNCD3Names() const { return cd3names.Size(); }

    DLL_HEADER static string default_bc;
    string * GetBCNamePtr (int bcnr) const
    { return (bcnr < bcnames.Size() && bcnames[bcnr]) ? bcnames[bcnr] : &default_bc; }


    DLL_HEADER NgArray<string*> & GetRegionNamesCD (int codim);
    
    ///
    void ClearFaceDescriptors()
    { facedecoding.SetSize(0); }

    ///
    int GetNFD () const
    { return facedecoding.Size(); }

    const FaceDescriptor & GetFaceDescriptor (const Element2d & el) const
    { return facedecoding[el.GetIndex()-1]; }
    
    const FaceDescriptor & GetFaceDescriptor (int i) const
    { return facedecoding[i-1]; }      
    // { return facedecoding.Get(i); }

    auto & FaceDescriptors () const { return facedecoding; }

    const EdgeDescriptor & GetEdgeDescriptor (int i) const
    { return edgedecoding[i]; }


    ///
    FaceDescriptor & GetFaceDescriptor (int i)
    { return facedecoding[i-1]; }      
    // { return facedecoding.Elem(i); }

    int IdentifyPeriodicBoundaries(const string& id_name,
                                   const string& s1,
                                   const Transformation<3>& mapping,
                                   double pointTolerance);

    // #ifdef NONE
    //   /*
    //     Identify points pi1 and pi2, due to
    //     identification nr identnr
    //   */
    //   void AddIdentification (int pi1, int pi2, int identnr);

    //   int GetIdentification (int pi1, int pi2) const;
    //   int GetIdentificationSym (int pi1, int pi2) const;
    //   ///
    //   INDEX_2_HASHTABLE<int> & GetIdentifiedPoints () 
    //   { 
    //     return *identifiedpoints; 
    //   }

    //   ///
    //   void GetIdentificationMap (int identnr, NgArray<int> & identmap) const;
    //   ///
    //   void GetIdentificationPairs (int identnr, NgArray<INDEX_2> & identpairs) const;
    //   ///
    //   int GetMaxIdentificationNr () const
    //   { 
    //     return maxidentnr; 
    //   }
    // #endif

    /// return periodic, close surface etc. identifications
    Identifications & GetIdentifications () { return *ident; }
    /// return periodic, close surface etc. identifications
    const Identifications & GetIdentifications () const { return *ident; }
    ///
    bool HasIdentifications() const { return ident != nullptr; }

    DLL_HEADER void InitPointCurve(double red = 1, double green = 0, double blue = 0) const;
    DLL_HEADER void AddPointCurvePoint(const Point3d & pt) const;
    DLL_HEADER int GetNumPointCurves(void) const;
    DLL_HEADER int GetNumPointsOfPointCurve(int curve) const;
    DLL_HEADER Point3d & GetPointCurvePoint(int curve, int n) const;
    DLL_HEADER void GetPointCurveColor(int curve, double & red, double & green, double & blue) const;




    /// find number of vertices
    DLL_HEADER void ComputeNVertices ();
    /// number of vertices (no edge-midpoints)
    DLL_HEADER int GetNV () const;
    /// remove edge points
    DLL_HEADER void SetNP (int np);

  

    DLL_HEADER Table<ElementIndex, PointIndex> CreatePoint2ElementTable(std::optional<BitArray> points = std::nullopt, int domain = 0) const;
    DLL_HEADER Table<SurfaceElementIndex, PointIndex> CreatePoint2SurfaceElementTable( int faceindex=0 ) const;

    DLL_HEADER bool PureTrigMesh (int faceindex = 0) const;
    DLL_HEADER bool PureTetMesh () const;


    const MeshTopology & GetTopology () const { return topology; }
    MeshTopology & GetTopology () { return topology; }

    DLL_HEADER void UpdateTopology (NgTaskManager tm = &DummyTaskManager,
                                    NgTracer tracer = &DummyTracer);
  
    class CurvedElements & GetCurvedElements () const
    { return *curvedelems; }
    
    DLL_HEADER void BuildCurvedElements  (const class Refinement * ref, int aorder, bool arational = false);
    DLL_HEADER void BuildCurvedElements  (int aorder);

    const class AnisotropicClusters & GetClusters () const
    { return *clusters; }


    class CSurfaceArea
    {
      const Mesh & mesh;
      bool valid;
      double area;
    public:
      CSurfaceArea (const Mesh & amesh) 
	: mesh(amesh), valid(false), area(0.) { ; }

      void Add (const Element2d & sel)
      {
	if (sel.GetNP() == 3)
	  area += Cross ( mesh[sel[1]]-mesh[sel[0]],
			  mesh[sel[2]]-mesh[sel[0]] ).Length() / 2;
	else
	  area += Cross (Vec3d (mesh[sel.PNum(1)], mesh[sel.PNum(3)]),
			 Vec3d (mesh[sel.PNum(1)], mesh[sel.PNum(4)])).Length() / 2;;
      }
      void ReCalc ()
      {
	area = 0;
        /*
	for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	  Add (mesh[sei]);
        */
        for (const Element2d & el : mesh.SurfaceElements())
          Add (el);
	valid = true;
      }

      operator double () const { return area; }
      bool Valid() const { return valid; }
    };

    CSurfaceArea surfarea;
    CSurfaceArea & SurfaceArea() { return surfarea; }
    const CSurfaceArea & SurfaceArea() const { return surfarea; }



    int GetTimeStamp() const { return timestamp; }
    void SetNextTimeStamp() 
    { timestamp = NextTimeStamp(); }

    int GetMajorTimeStamp() const { return majortimestamp; }
    void SetNextMajorTimeStamp() 
    { majortimestamp = timestamp = NextTimeStamp(); }


    /// return mutex
    NgMutex & Mutex ()   { return mutex; }
    NgMutex & MajorMutex ()   { return majormutex; }


    DLL_HEADER shared_ptr<NetgenGeometry> GetGeometry() const;
    void SetGeometry (shared_ptr<NetgenGeometry> geom) 
    {
      geometry = geom;
    }

    ///
    void SetUserData(const char * id, NgArray<int> & data);
    ///
    bool GetUserData(const char * id, NgArray<int> & data, int shift = 0) const;
    ///
    void SetUserData(const char * id, NgArray<double> & data);
    ///
    bool GetUserData(const char * id, NgArray<double> & data, int shift = 0) const;

    ///
    friend void OptimizeRestart (Mesh & mesh3d);
    ///
    void PrintMemInfo (ostream & ost) const;
    /// 
    friend class Meshing3;

    // only for saving the geometry
    enum GEOM_TYPE { NO_GEOM = 0, GEOM_2D = 1, GEOM_CSG = 10, GEOM_STL = 11, GEOM_OCC = 12, GEOM_ACIS = 13 };
    GEOM_TYPE geomtype;
  

#ifdef PARALLEL
    /// returns parallel topology
    class ParallelMeshTopology & GetParallelTopology () const
    { return *paralleltop; }

    /// distributes the master-mesh to local meshes
    DLL_HEADER void Distribute ();
    DLL_HEADER void Distribute (NgArray<int> & volume_weights, NgArray<int> & surface_weights,
		     NgArray<int> & segment_weights);


    /// find connection to parallel meshes
    //   void FindExchangePoints () ;

    //   void FindExchangeEdges ();
    //   void FindExchangeFaces ();

    /// use metis to decompose master mesh 
    DLL_HEADER void ParallelMetis (int nproc); //  NgArray<int> & neloc );
    DLL_HEADER void ParallelMetis (NgArray<int> & volume_weights, NgArray<int> & surface_weights,
			NgArray<int> & segment_weights); 

    void PartHybridMesh (); //  NgArray<int> & neloc );
    void PartDualHybridMesh (); //  NgArray<int> & neloc );
    void PartDualHybridMesh2D ();  // ( NgArray<int> & neloc );


    /// send mesh from master to local procs
    void SendRecvMesh ();

    /// send mesh to parallel machine, keep global mesh at master 
    void SendMesh ( ) const;   // Mesh * mastermesh, NgArray<int> & neloc) const;
    /// loads a mesh sent from master processor
    void ReceiveParallelMesh ();

    
#else
    void ParallelMetis (int /* nproc */) {}
    void Distribute () {}
    void SendRecvMesh () {}
    void Distribute (NgArray<int> & volume_weights, NgArray<int> & surface_weights, 
      NgArray<int> & segment_weights){ }
#endif

    NgArray<int> vol_partition;
    NgArray<int> surf_partition;
    NgArray<int> seg_partition;

    shared_ptr<Mesh> Mirror( netgen::Point<3> p, Vec<3> n );

    private:
    MemoryTracer mem_tracer = {"Mesh",
      points, "points",
      segments, "segments",
      surfelements, "surfelements",
      volelements, "volelements"
    };
    public:
    const MemoryTracer & GetMemoryTracer() { return mem_tracer; }
  };

  inline ostream& operator<<(ostream& ost, const Mesh& mesh)
  {
    ost << "mesh: " << endl;
    mesh.Save(ost);
    return ost;
  }



  FlatArray<T_EDGE> MeshTopology :: GetEdges (SurfaceElementIndex elnr) const
  {
    return FlatArray<T_EDGE>(GetNEdges ( (*mesh)[elnr].GetType()), &surfedges[elnr][0]);
  }

  FlatArray<T_EDGE> MeshTopology :: GetEdges (ElementIndex elnr) const
  {
    return FlatArray<T_EDGE>(GetNEdges ( (*mesh)[elnr].GetType()), &edges[elnr][0]);
  }
  
  FlatArray<T_FACE> MeshTopology :: GetFaces (ElementIndex elnr) const
  {
    return FlatArray<T_FACE>(GetNFaces ( (*mesh)[elnr].GetType()), &faces[elnr][0]);
  }

  
}

#endif // NETGEN_MESHCLASS_HPP

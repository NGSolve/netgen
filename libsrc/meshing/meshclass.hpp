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

#include "paralleltop.hpp"

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

  typedef Array<MarkedTet,ElementIndex> T_MTETS;
  typedef Array<MarkedPrism> T_MPRISMS;
  typedef Array<MarkedIdentification> T_MIDS;
  typedef Array<MarkedTri> T_MTRIS;
  typedef Array<MarkedQuad> T_MQUADS;

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
    T_SURFELEMENTS surfelements;
    /// volume elements
    T_VOLELEMENTS volelements;
    /// points will be fixed forever
    Array<PointIndex> lockedpoints;


    /// surface indices at boundary nodes
    // TABLE<int,PointIndex::BASE> surfacesonnode;
    /// boundary edges  (1..normal bedge, 2..segment)
    unique_ptr<ClosedHashTable<SortedPointIndices<2>, int>> boundaryedges;
    ///
    unique_ptr<ClosedHashTable<SortedPointIndices<2>, SegmentIndex>> segmentht;
    ///
    unique_ptr<ClosedHashTable<SortedPointIndices<3>, SurfaceElementIndex>> surfelementht;
    unique_ptr<ClosedHashTable<SortedPointIndices<3>, int>> illegal_trigs;

    /// faces of rest-solid
    Array<Element2d> openelements;
    /// open segments for surface meshing
    Array<Segment> opensegments;
    /// face descriptor index for each open segment (parallel to opensegments)
    Array<int> opensegment_faces;

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
    Array<double> maxhdomain;
  
    /**
       regions of dimension 0..3 (vertices, edges, faces, volumes).
       The index of an element of dimension D maps into std::get<D>(regions).
    */
    std::tuple<RegionArray<0>, RegionArray<1>, RegionArray<2>, RegionArray<3>> regions;

    /// Periodic surface, close surface, etc. identifications
    unique_ptr<Identifications> ident;


    /// number of vertices (if < 0, use np)
    int numvertices;

    /// geometric search tree for interval intersection search
    unique_ptr<BoxTree<3, ElementIndex>> elementsearchtree_vol;
    unique_ptr<BoxTree<3, SurfaceElementIndex>> elementsearchtree_surf;
    /// time stamp for tree
    mutable size_t elementsearchtreets[4];

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
    std::mutex mutex;
    /// mesh access semaphores.
    std::mutex majormutex;

    SymbolTable< Array<int>* > userdata_int;
    SymbolTable< Array<double>* > userdata_double;


    mutable Array< netgen::Point<3> > pointcurves;
    mutable Array<int> pointcurves_startpoint;
    mutable Array<double> pointcurves_red,pointcurves_green,pointcurves_blue;


    /// start element for point search (GetElementOfPoint)
    mutable int ps_startelement;


#ifdef PARALLEL
    /// connection to parallel meshes
    unique_ptr<ParallelMeshTopology> paralleltop;
#endif

    
    shared_ptr<NetgenGeometry> geometry;


  public:
    DLL_HEADER void BuildBoundaryEdges(bool rebuild=true);

    DLL_HEADER bool PointContainedIn2DElement(const netgen::Point<3> & p,
                                   double lami[3],
                                   SurfaceElementIndex element,
                                   bool consider3D = false) const;
    DLL_HEADER bool PointContainedIn3DElement(const netgen::Point<3> & p,
                                   double lami[3],
                                   ElementIndex element,
                                   double tol=1e-4) const;
    DLL_HEADER bool PointContainedIn3DElementOld(const netgen::Point<3> & p,
                                      double lami[3],
                                      ElementIndex element,
                                      double tol=1e-4) const;

  public:
    Signal<> updateSignal;
    BisectionInfo bisectioninfo;

    // store coarse mesh before hp-refinement
    unique_ptr<Array<HPRefElement>> hpelements;
    unique_ptr<Mesh> coarsemesh;

    /// per-element data of hp- and p-refinement; the arrays stay empty until something is set
    struct HPElementInfo
    {
      int hp_elnr = -1;
      uint8_t orderx = 1, ordery = 1, orderz = 1;
    };
  private:
    Array<HPElementInfo, ElementIndex> hp_volinfo;
    Array<HPElementInfo, SurfaceElementIndex> hp_surfinfo;
    Array<HPElementInfo, SegmentIndex> hp_seginfo;

    template <typename TIndex>
    const Array<HPElementInfo,TIndex> & HPInfo () const
    {
      if constexpr (std::is_same_v<TIndex, ElementIndex>) return hp_volinfo;
      else if constexpr (std::is_same_v<TIndex, SurfaceElementIndex>) return hp_surfinfo;
      else return hp_seginfo;
    }
    template <typename TIndex>
    Array<HPElementInfo,TIndex> & HPInfo ()
    { return const_cast<Array<HPElementInfo,TIndex>&> (std::as_const(*this).HPInfo<TIndex>()); }
    template <typename TIndex>
    HPElementInfo GetHPInfo (TIndex i) const
    {
      auto & info = HPInfo<TIndex>();
      return info.Range().Contains(i) ? info[i] : HPElementInfo();
    }
    template <typename TIndex>
    void SetHPInfo (TIndex i, HPElementInfo val);
  public:
    template <typename TIndex>
    int GetOrder (TIndex i) const { return GetHPInfo(i).orderx; }
    void GetOrder (ElementIndex i, int & ox, int & oy, int & oz) const
    { auto h = GetHPInfo(i); ox = h.orderx; oy = h.ordery; oz = h.orderz; }
    void GetOrder (SurfaceElementIndex i, int & ox, int & oy, int & oz) const
    { auto h = GetHPInfo(i); ox = h.orderx; oy = h.ordery; oz = 0; }
    void GetOrder (SurfaceElementIndex i, int & ox, int & oy) const
    { auto h = GetHPInfo(i); ox = h.orderx; oy = h.ordery; }
    template <typename TIndex>
    void SetOrder (TIndex i, int order) { SetOrder (i, order, order, order); }
    void SetOrder (ElementIndex i, int ox, int oy, int oz)
    { auto h = GetHPInfo(i); h.orderx = ox; h.ordery = oy; h.orderz = oz; SetHPInfo (i, h); }
    void SetOrder (SurfaceElementIndex i, int ox, int oy, int /* oz */ = 0)
    { auto h = GetHPInfo(i); h.orderx = ox; h.ordery = oy; h.orderz = 1; SetHPInfo (i, h); }
    /// number of the HPRefElement an element was made from (-1 if none)
    template <typename TIndex>
    int GetHpElnr (TIndex i) const { return GetHPInfo(i).hp_elnr; }
    template <typename TIndex>
    void SetHpElnr (TIndex i, int nr) { auto h = GetHPInfo(i); h.hp_elnr = nr; SetHPInfo (i, h); }
    template <typename TIndex>
    void AllocateHPInfo ();
    template <typename TIndex>
    bool HasHPInfo () const { return HPInfo<TIndex>().Size() > 0; }
  
  
    /// number of refinement levels
    // int mglevels;
    // number of vertices on each refinement level:
    Array<size_t> level_nv;
    /// refinement hierarchy
    Array<PointIndices<2>,PointIndex> mlbetweennodes;
    /// parent element of volume element
    Array<ElementIndex, ElementIndex> mlparentelement;
    /// parent element of surface element
    Array<SurfaceElementIndex, SurfaceElementIndex> mlparentsurfaceelement;



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
      hp_volinfo.SetSize(0);
      timestamp = NextTimeStamp();
    }

    ///
    DLL_HEADER void ClearSegments()
    { 
      segments.SetSize(0); 
      hp_seginfo.SetSize(0);
      timestamp = NextTimeStamp();
    }
    
    ///
    bool TestOk () const;

    void SetAllocSize(int nnodes, int nsegs, int nsel, int nel);
    

    DLL_HEADER PointIndex AddPoint (const netgen::Point<3> & p, int layer = 1);
    DLL_HEADER PointIndex AddPoint (const netgen::Point<3> & p, int layer, POINTTYPE type);

    auto GetNP () const { return points.Size(); }

    MeshPoint & Point(PointIndex pi) { return points[pi]; }
    const MeshPoint & Point(PointIndex pi) const { return points[pi]; }

    const MeshPoint & operator[] (PointIndex pi) const { return points[pi]; }
    MeshPoint & operator[] (PointIndex pi) { return points[pi]; }

    const T_POINTS & Points() const { return points; }
    T_POINTS & Points() { return points; }


    DLL_HEADER SegmentIndex AddSegment (const Segment & s);
    void DeleteSegment (SegmentIndex si)
    {
      segments[si][0].Invalidate();
      segments[si][1].Invalidate();
    }

    int GetNSeg () const { return segments.Size(); }
    Segment & LineSegment(SegmentIndex si) { return segments[si]; }
    const Segment & LineSegment(SegmentIndex si) const { return segments[si]; }
    const Segment & operator[] (SegmentIndex si) const { return segments[si]; }
    Segment & operator[] (SegmentIndex si) { return segments[si]; }

    const auto & LineSegments() const { return segments; }
    auto & LineSegments() { return segments; }
    
    Array<Element0d> pointelements;  // only via python interface

    DLL_HEADER SurfaceElementIndex AddSurfaceElement (const Element2dRef & el);
    // write to pre-allocated container, thread-safe
    DLL_HEADER void SetSurfaceElement (SurfaceElementIndex sei, const Element2dRef & el);
    
    void Delete (SurfaceElementIndex eli)
    {
      // for (auto & p : surfelements[eli].PNums()) p.Invalidate();
      surfelements[eli].Delete();
      timestamp = NextTimeStamp();
    }

    auto GetNSE () const { return surfelements.Size(); }

    // [[deprecated("Use mesh[](SurfaceElementIndex) instead !")]]
    Element2dRef SurfaceElement(SurfaceElementIndex i) { return surfelements[i]; }
    // [[deprecated("Use mesh[](SurfaceElementIndex) instead !")]]
    const Element2dRef SurfaceElement(SurfaceElementIndex i) const { return surfelements[i]; }

    const Element2dRef operator[] (SurfaceElementIndex ei) const { return surfelements[ei]; }
    Element2dRef operator[] (SurfaceElementIndex ei) { return surfelements[ei]; }

    const auto & SurfaceElements() const { return surfelements; }
    auto & SurfaceElements() { return surfelements; }

    
    DLL_HEADER void RebuildSurfaceElementLists ();
    /// surface elements of face fi; INVALID: all surface elements
    DLL_HEADER void GetSurfaceElementsOfFace (FaceRegionIndex fi, Array<SurfaceElementIndex> & sei) const;
    void GetSurfaceElementsOfFace (int facenr, Array<SurfaceElementIndex> & sei) const
    { GetSurfaceElementsOfFace (FaceRegionIndex::FromNr1(facenr), sei); }

    DLL_HEADER ElementIndex AddVolumeElement (const ElementRef & el);
    // write to pre-allocated container, thread-safe
    DLL_HEADER void SetVolumeElement (ElementIndex sei, const ElementRef & el);

    auto GetNE () const { return volelements.Size(); }

    // [[deprecated("Use mesh[](VolumeElementIndex) instead !")]]
    ElementRef VolumeElement(ElementIndex i) { return volelements[i]; }
    // [[deprecated("Use mesh[](VolumeElementIndex) instead !")]]
    const ElementRef VolumeElement(ElementIndex i) const { return volelements[i]; }

    const ElementRef operator[] (ElementIndex ei) const { return volelements[ei]; }
    ElementRef operator[] (ElementIndex ei) { return volelements[ei]; }

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
    void FixPoints (const TBitArray<PointIndex> & fixpoints);

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
    const Segment & GetOpenSegment (int nr) { return opensegments[nr-1]; }
    /// face descriptor index for open segment nr (1-based)
    int GetOpenSegmentFace (int nr) { return opensegment_faces[nr-1]; }
  
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
    DLL_HEADER void RestrictLocalH (const netgen::Point<3> & p, double hloc, int layer=1);
    ///
    DLL_HEADER void RestrictLocalHLine (const netgen::Point<3> & p1, const netgen::Point<3> & p2, 
                             double hloc, int layer=1);
    /// number of elements per radius
    DLL_HEADER void CalcLocalHFromSurfaceCurvature(double grading, double elperr, int layer=1);
    ///
    DLL_HEADER void CalcLocalHFromPointDistances(double grading, int layer=1);
    ///
    DLL_HEADER void RestrictLocalH (resthtype rht, int nr, double loch);
    DLL_HEADER void RestrictLocalH (const Element2dRef & sel, double loch);
    DLL_HEADER void RestrictLocalH (const Segment & seg, double loch);
    ///
    DLL_HEADER void LoadLocalMeshSize (const filesystem::path & meshsizefilename);
    ///
    DLL_HEADER void SetGlobalH (double h);
    ///
       DLL_HEADER void SetMinimalH (double h);
    ///
        DLL_HEADER double MaxHDomain (int dom) const;
    ///
        DLL_HEADER void SetMaxHDomain (const Array<double> & mhd);
    ///
    DLL_HEADER double GetH (const netgen::Point<3> & p, int layer=1) const;
    DLL_HEADER double GetH (PointIndex pi) const { return GetH(points[pi], points[pi].GetLayer()); }
    ///
    double GetMinH (const netgen::Point<3> & pmin, const netgen::Point<3> & pmax, int layer=1);
    ///
    bool HasLocalHFunction (int layer=1) { return lochfunc[layer-1] != nullptr; }
    ///
    LocalH & LocalHFunction (int layer=1) { return * lochfunc[layer-1]; }

    shared_ptr<LocalH> & GetLocalH(int layer=1) const
    {
      if(lochfunc.Size() == 1)
        return lochfunc[0];
      return lochfunc[layer-1];
    }
    DLL_HEADER void SetLocalH(shared_ptr<LocalH> loch, int layer=1);

    ///
    bool LocalHFunctionGenerated(int layer=1) const { return (lochfunc[layer-1] != NULL); }

    /// Find bounding box
    DLL_HEADER void GetBox (netgen::Point<3> & pmin, netgen::Point<3> & pmax, int dom = -1) const;

    /// Find bounding box of points of typ ptyp or less
    DLL_HEADER void GetBox (netgen::Point<3> & pmin, netgen::Point<3> & pmax, POINTTYPE ptyp ) const;

    ///
    int GetNOpenElements() const
    { return openelements.Size(); }
    ///
    const Element2dRef & OpenElement(int i) const
    { return openelements[i-1]; }

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
      
      return boundaryedges->Used ({pi1, pi2});
    }

    void DeleteBoundaryEdges ()
    {
        boundaryedges = nullptr;
    }

    bool IsSegment (PointIndex pi1, PointIndex pi2) const
    {
      return segmentht->Used ({pi1, pi2});
    }

    SegmentIndex SegmentNr (PointIndex pi1, PointIndex pi2) const
    {
      return segmentht->Get ({pi1, pi2});
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
    void ImproveMeshJacobian (const MeshingParameters & mp, OPTIMIZEGOAL goal = OPT_QUALITY,
                              const TBitArray<PointIndex> * usepoint = NULL);
    ///
    void ImproveMeshJacobianOnSurface (const MeshingParameters & mp,
                                       const TBitArray<PointIndex> & usepoint, 
                                       const Array< Vec<3>* > & nv,
                                       OPTIMIZEGOAL goal = OPT_QUALITY,
                                       const Array< idmap_type* > * idmaps = NULL);
    /**
       free nodes in environment of openelements 
       for optimiztion
    */
    void FreeOpenElementsEnvironment (int layers);


    DLL_HEADER double CalcTotalBad (const MeshingParameters & mp);
    FlatArray<int> GetQualityHistogram() { return tets_in_qualclass; }

    ///
    bool LegalTet (ElementRef el) const
    {
      if (el.IllegalValid())
        return !el.Illegal();
      return LegalTet2 (el);
    }
    ///
    bool LegalTet2 (ElementRef el) const;


    ///
    // Find trigs with same vertices
    // return: number of illegal trigs
    int FindIllegalTrigs ();

    bool LegalTrig (const Element2dRef & el) const;
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
        DLL_HEADER int MarkIllegalElements (int domain=0);

    /// orient surface mesh, for one sub-domain only
        DLL_HEADER void SurfaceMeshOrientation ();

    /// convert mixed element mesh to tet-mesh
        DLL_HEADER void Split2Tets();


    /// build box-search tree
    DLL_HEADER void BuildElementSearchTree (int dim);
    BoxTree<3, ElementIndex>* GetElementSearchTree () const
    {
        return elementsearchtree_vol.get();
    }

    BoxTree<3, SurfaceElementIndex>* GetSurfaceElementSearchTree () const
    {
      return elementsearchtree_surf.get();
    }

    void SetPointSearchStartElement(const int el) const {ps_startelement = el;}

    /// gives element of point, barycentric coordinates
    DLL_HEADER ElementIndex
    GetElementOfPoint (const netgen::Point<3> & p,
                       double * lami,
                       bool build_searchtree = false,
                       int index = -1,
                       bool allowindex = true,
                       double tol=1e-4) const;
    DLL_HEADER ElementIndex
    GetElementOfPoint (const netgen::Point<3> & p,
                       double * lami,
                       std::optional<FlatArray<int>> indices,
                       bool build_searchtree = 0,
                       bool allowindex = true,
                       double tol=1e-4) const;
    DLL_HEADER SurfaceElementIndex
    GetSurfaceElementOfPoint (const netgen::Point<3> & p,
                              double * lami,
                              bool build_searchtree = false,
                              int index = -1,
                              bool allowindex = true) const;
    DLL_HEADER SurfaceElementIndex
    GetSurfaceElementOfPoint (const netgen::Point<3> & p,
                              double * lami,
                              std::optional<FlatArray<int>> indices,
                              bool build_searchtree = false,
                              bool allowindex = true) const;

    /// give list of vol elements which are int the box(p1,p2)
    void GetIntersectingVolEls(const netgen::Point<3>& p1, const netgen::Point<3>& p2, 
                               Array<ElementIndex> & locels) const;

    /// regions of entity dimension D, indexed by RegionIndex<D>
    template <int D> RegionArray<D> & Regions () { return std::get<D>(regions); }
    template <int D> const RegionArray<D> & Regions () const { return std::get<D>(regions); }

    template <int D> Region<D> & GetRegion (RegionIndex<D> i) { return Regions<D>()[i]; }
    template <int D> const Region<D> & GetRegion (RegionIndex<D> i) const { return Regions<D>()[i]; }

    template <int D> RegionIndex<D> AddRegion (const Region<D> & reg) { return Regions<D>().Append(reg); }

    ///
    FaceRegionIndex AddFaceDescriptor(const FaceRegion& fd)
    { return Regions<2>().Append(fd); }

    EdgeRegionIndex AddEdgeDescriptor(const EdgeRegion & fd)
    { return Regions<1>().Append(fd); }

    auto & GetCommunicator() const { return this->comm; }
    void SetCommunicator(NgMPI_Comm acomm);
    
    DLL_HEADER void SplitFacesByAdjacentDomains();
    DLL_HEADER shared_ptr<Mesh> GetSubMesh(string domains="", string faces="") const;

    /// name of domain domnr (1-based): materials in 3D, face descriptor in 2D, edge descriptor in 1D
    DLL_HEADER void SetMaterial (int domnr, const string & mat);
    DLL_HEADER const string & GetMaterial (int domnr) const;
    DLL_HEADER static string defaultmat;
    /// 3D domain name
    const string * GetMaterialPtr (int domnr) const // 1-based
    {
      return (domnr >= 1 && domnr <= Regions<3>().Size()) ? &Regions<3>()[VolumeRegionIndex::FromNr1(domnr)].GetName() : &defaultmat;
    }
    
    /// 1D meshes only (vertex names); a no-op otherwise
    DLL_HEADER void SetNBCNames ( int nbcn );

    /// name of boundary region nr (0-based): face descriptor nr in 3D, edge descriptor nr in 2D, vertex nr in 1D.
    /// Missing descriptors are created.
    DLL_HEADER void SetBCName ( int bcnr, const string & abcname );
    DLL_HEADER const string & GetBCName ( int bcnr ) const;
    /// boundary names keyed by bc number, the layout of files and archives (empty if nothing is named)
    DLL_HEADER Array<string> BCNamesByNumber () const;
    /// name of the boundary described by face descriptor fdi
    const string & GetBCName (FaceRegionIndex fdi) const { return Regions<2>()[fdi].GetBCName(); }

    /// vertex names of 2D meshes (cd2nr 1-based); edge names of 3D meshes live in the edge descriptors
    DLL_HEADER void SetCD2Name (int cd2nr, const string & abcname);
    DLL_HEADER const string & GetCD2Name (int cd2nr ) const;
    DLL_HEADER static string cd2_default_name;
    size_t GetNCD2Names() const { return dimension == 2 ? Regions<0>().Size() : 0; }

    /// edge descriptor nr (1-based); missing descriptors up to nr are created
    DLL_HEADER EdgeRegion & EnsureEdgeDescriptor (int nr);
  private:
    void SetCD2NameCompat (int cd2nr, const string & name);
    /// names of dimension D in the archive layout Array<optional<string>>
    template <int D> void ArchiveRegionNames (Archive & ar)
    {
      auto names = RegionNames<D>();
      ar & names;
      if (ar.Input()) SetRegionNames<D>(names);
    }
  public:

    DLL_HEADER void SetNCD3Names (int ncd3n);
    DLL_HEADER void SetCD3Name (int cd3nr, const string & abcname);
    DLL_HEADER int AddCD3Name (const string & aname);
    DLL_HEADER const string & GetCD3Name (int cd3nr ) const;
    DLL_HEADER static string cd3_default_name;
    const string * GetCD3NamePtr (int cd3nr ) const
    {
      if (cd3nr >= 0 && cd3nr < Regions<0>().Size()) return &Regions<0>()[VertexRegionIndex::FromNr0(cd3nr)].GetName();
      return &cd3_default_name;
    }
    size_t GetNCD3Names() const { return dimension == 3 ? Regions<0>().Size() : 0; }

    DLL_HEADER static string default_bc;
    const string * GetBCNamePtr (int bcnr) const
    {
      if (dimension == 3)
        return (bcnr >= 0 && bcnr < Regions<2>().Size()) ? &Regions<2>()[FaceRegionIndex::FromNr0(bcnr)].GetBCName() : &default_bc;
      if (dimension == 2)
        return (bcnr >= 0 && bcnr < Regions<1>().Size()) ? &Regions<1>()[EdgeRegionIndex::FromNr0(bcnr)].GetName() : &default_bc;
      return (bcnr >= 0 && bcnr < Regions<0>().Size()) ? &Regions<0>()[VertexRegionIndex::FromNr0(bcnr)].GetName() : &default_bc;
    }

    DLL_HEADER std::string_view GetRegionName(const Segment & el) const;
    DLL_HEADER std::string_view GetRegionName(const Element2dRef & el) const;
    DLL_HEADER std::string_view GetRegionName(const ElementRef & el) const;

    std::string_view GetRegionName(SegmentIndex ei) const { return GetRegionName((*this)[ei]); }
    std::string_view GetRegionName(SurfaceElementIndex ei) const { return GetRegionName((*this)[ei]); }
    std::string_view GetRegionName(ElementIndex ei) const { return GetRegionName((*this)[ei]); }

    DLL_HEADER static string_view defaultmat_sv;
    /// number of regions of entity dimension dim (3D domains, faces, edges, vertices)
    size_t GetNRegions (int dim) const
    {
      switch (dim)
        {
        case 3: return Regions<3>().Size();
        case 2: return Regions<2>().Size();
        case 1: return Regions<1>().Size();
        default: return Regions<0>().Size();
        }
    }
    /// name of region nr (1-based) of entity dimension dim
    std::string_view GetRegionName (int dim, int nr) const
    {
      switch (dim)
        {
        case 3: return GetRegionName<3>(nr);
        case 2: return GetRegionName<2>(nr);
        case 1: return GetRegionName<1>(nr);
        default: return GetRegionName<0>(nr);
        }
    }
    template <int D>
    std::string_view GetRegionName (int nr) const  // 1-based
    {
      auto & regs = Regions<D>();
      return (nr >= 1 && nr <= regs.Size()) ? string_view(regs[RegionIndex<D>::FromNr1(nr)].GetName()) : defaultmat_sv;
    }
    /// region names of dimension D, nullopt where not set
    template <int D>
    Array<optional<string>> RegionNames () const
    {
      Array<optional<string>> names(Regions<D>().Size());
      for (auto i : Regions<D>().Range())
        names[i.Nr0()] = Regions<D>()[i].OptName();
      return names;
    }
    template <int D>
    void SetRegionNames (FlatArray<optional<string>> names)
    {
      Regions<D>().SetSize(names.Size());
      for (auto i : Regions<D>().Range())
        Regions<D>()[i].SetName(names[i.Nr0()]);
    }
    /// domain names (dimension of the mesh) as array, the layout of files, archives and MPI messages
    DLL_HEADER Array<optional<string>> DomainNames () const;
    DLL_HEADER void SetDomainNames (Array<optional<string>> names);
    
    ///
    void ClearFaceDescriptors()
    { Regions<2>().SetSize(0); }

    void FreeFaceDescriptors()
    { Regions<2>() = RegionArray<2>(); }

    ///
    int GetNFD () const
    { return Regions<2>().Size(); }

    const FaceRegion & GetFaceDescriptor (const Element2dRef & el) const
    { return Regions<2>()[el.GetIndex()]; }
    FaceRegion & GetFaceDescriptor (const Element2dRef & el)
    { return Regions<2>()[el.GetIndex()]; }

    /// surface element refers to an existing face descriptor
    bool HasFaceDescriptor (const Element2dRef & el) const
    { return Regions<2>().Range().Contains(el.GetIndex()); }
    
    const FaceRegion & GetFaceDescriptor (FaceRegionIndex i) const
    { return Regions<2>()[i]; }

    auto & FaceDescriptors () const { return Regions<2>(); }

    const EdgeRegion & GetEdgeDescriptor (EdgeRegionIndex i) const
    { return Regions<1>()[i]; }
    EdgeRegion & GetEdgeDescriptor (EdgeRegionIndex i)
    { return Regions<1>()[i]; }
    /// 1-based
    const EdgeRegion & GetEdgeDescriptor (int i) const
    { return Regions<1>()[EdgeRegionIndex::FromNr1(i)]; }
    EdgeRegion & GetEdgeDescriptor (int i)
    { return Regions<1>()[EdgeRegionIndex::FromNr1(i)]; }

    const EdgeRegion & GetEdgeDescriptor (const Segment & seg) const
    { return Regions<1>()[seg.GetIndex()]; }
    EdgeRegion & GetEdgeDescriptor (const Segment & seg)
    { return Regions<1>()[seg.GetIndex()]; }

    /// segment refers to an existing edge descriptor
    bool HasEdgeDescriptor (const Segment & seg) const
    { return Regions<1>().Range().Contains(seg.GetIndex()); }

    int GetNED () const
    { return Regions<1>().Size(); }

    auto & EdgeDescriptors () const { return Regions<1>(); }
    auto & EdgeDescriptors () { return Regions<1>(); }

    void ClearEdgeDescriptors()
    { Regions<1>().SetSize(0); }

    void ReconstructEdgeDescriptors(const Array<std::pair<int,int>, SegmentIndex> * seg_surfnrs = nullptr,
                                    const Array<int, SegmentIndex> * seg_edgenrs = nullptr);

    /// Recompute EdgeRegion::fdindex from segment si values or FD lookup
    void RebuildFDIndices();



    ///
    FaceRegion & GetFaceDescriptor (FaceRegionIndex i)
    { return Regions<2>()[i]; }

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
    //   void GetIdentificationMap (int identnr, Array<int> & identmap) const;
    //   ///
    //   void GetIdentificationPairs (int identnr, Array<IVec<2>> & identpairs) const;
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
    DLL_HEADER void AddPointCurvePoint(const netgen::Point<3> & pt) const;
    DLL_HEADER int GetNumPointCurves(void) const;
    DLL_HEADER int GetNumPointsOfPointCurve(int curve) const;
    DLL_HEADER netgen::Point<3> & GetPointCurvePoint(int curve, int n) const;
    DLL_HEADER void GetPointCurveColor(int curve, double & red, double & green, double & blue) const;




    /// find number of vertices
    DLL_HEADER void ComputeNVertices ();
    /// number of vertices (no edge-midpoints)
    DLL_HEADER int GetNV () const;
    /// remove edge points
    DLL_HEADER void SetNP (int np);

  

    DLL_HEADER Table<ElementIndex, PointIndex> CreatePoint2ElementTable(std::optional<TBitArray<PointIndex>> points = std::nullopt, int domain = 0) const;
    // DLL_HEADER Table<SurfaceElementIndex, PointIndex> CreatePoint2SurfaceElementTable( int faceindex=0 ) const;
    DLL_HEADER Table<SurfaceElementIndex, PointIndex> CreatePoint2SurfaceElementTable( int faceindex=0 ) const;
    DLL_HEADER CompressedTable<SurfaceElementIndex, PointIndex> CreateCompressedPoint2SurfaceElementTable( FaceRegionIndex fi = FaceRegionIndex::INVALID ) const;
    CompressedTable<SurfaceElementIndex, PointIndex> CreateCompressedPoint2SurfaceElementTable( int faceindex ) const
    { return CreateCompressedPoint2SurfaceElementTable (FaceRegionIndex::FromNr1(faceindex)); }

    DLL_HEADER bool PureTrigMesh (int faceindex = 0) const;
    DLL_HEADER bool PureTetMesh () const;


    const MeshTopology & GetTopology () const { return topology; }
    MeshTopology & GetTopology () { return topology; }

    DLL_HEADER void UpdateTopology ();
  
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

      void Add (const Element2dRef & sel)
      {
        if (sel.GetNP() == 3)
          area += Cross ( mesh[sel[1]]-mesh[sel[0]],
                          mesh[sel[2]]-mesh[sel[0]] ).Length() / 2;
        else
          area += Cross (Vec<3> (mesh[sel[0]], mesh[sel[2]]),
                         Vec<3> (mesh[sel[0]], mesh[sel[3]])).Length() / 2;;
      }
      void ReCalc ()
      {
        area = 0;
        /*
        for (auto & el : mesh.SurfaceElements())
          Add (el);
        */
        for (auto el : mesh.SurfaceElements())
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
    std::mutex & Mutex ()   { return mutex; }
    std::mutex & MajorMutex ()   { return majormutex; }


    DLL_HEADER shared_ptr<NetgenGeometry> GetGeometry() const;
    void SetGeometry (shared_ptr<NetgenGeometry> geom) 
    {
      geometry = geom;
    }

    ///
    void SetUserData(const char * id, Array<int> & data);
    ///
    bool GetUserData(const char * id, Array<int> & data, int shift = 0) const;
    ///
    void SetUserData(const char * id, Array<double> & data);
    ///
    bool GetUserData(const char * id, Array<double> & data, int shift = 0) const;

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
    DLL_HEADER void Distribute (Array<int> & volume_weights, Array<int> & surface_weights,
                                Array<int> & segment_weights);


    /// find connection to parallel meshes
    //   void FindExchangePoints () ;

    //   void FindExchangeEdges ();
    //   void FindExchangeFaces ();

    /// use metis to decompose master mesh 
    DLL_HEADER void ParallelMetis (int nproc); 
    DLL_HEADER void ParallelMetis (Array<int> & volume_weights, Array<int> & surface_weights,
                                   Array<int> & segment_weights); 

    void PartHybridMesh (); 
    void PartDualHybridMesh (); 
    void PartDualHybridMesh2D ();

    /// send mesh from master to local procs
    void SendRecvMesh ();

    /// send mesh to parallel machine, keep global mesh at master 
    void SendMesh ( ) const;   
    /// loads a mesh sent from master processor
    void ReceiveParallelMesh ();

    
#else
    void ParallelMetis (int /* nproc */) {}
    void Distribute () {}
    void SendRecvMesh () {}
    void Distribute (Array<int> & volume_weights, Array<int> & surface_weights, 
      Array<int> & segment_weights){ }
#endif

    Array<int, ElementIndex> vol_partition;
    Array<int, SurfaceElementIndex> surf_partition;
    Array<int, SegmentIndex> seg_partition;

    shared_ptr<Mesh> Mirror( netgen::Point<3> p, Vec<3> n );
  };

  inline ostream& operator<<(ostream& ost, const Mesh& mesh)
  {
    ost << "mesh: " << endl;
    mesh.Save(ost);
    return ost;
  }



  // the topology tables are DynStrideArrays with a run-time width >= the element's nedges/nfaces,
  // the first GetNEdges/GetNFaces entries belong to the element
  FlatArray<const EdgeIndex> MeshTopology :: GetEdges (SurfaceElementIndex elnr) const
  {
    return FlatArray<const EdgeIndex>(GetNEdges ( (*mesh)[elnr].GetType()), surfedges[elnr].TailPtr());
  }

  FlatArray<const EdgeIndex> MeshTopology :: GetEdges (ElementIndex elnr) const
  {
    return FlatArray<const EdgeIndex>(GetNEdges ( (*mesh)[elnr].GetType()), edges[elnr].TailPtr());
  }
  
  FlatArray<const FaceIndex> MeshTopology :: GetFaces (ElementIndex elnr) const
  {
    return FlatArray<const FaceIndex>(GetNFaces ( (*mesh)[elnr].GetType()), faces[elnr].TailPtr());
  }

  /// a surface element has one face, a segment one edge
  FlatArray<FaceIndex> MeshTopology :: GetFaces (SurfaceElementIndex elnr) const
  {
    return FlatArray<FaceIndex>(1, &surffaces[elnr]);
  }

  FlatArray<EdgeIndex> MeshTopology :: GetEdges (SegmentIndex segnr) const
  {
    return FlatArray<EdgeIndex>(1, &segedges[segnr]);
  }

  DLL_HEADER void AddFacesBetweenDomains(Mesh & mesh);
}

#endif // NETGEN_MESHCLASS_HPP

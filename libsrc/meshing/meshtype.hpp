#ifndef MESHTYPE
#define MESHTYPE


/**************************************************************************/
/* File:   meshtype.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include <mydefs.hpp>
#include <general/template.hpp>
#include <core/mpi_wrapper.hpp>
#include <gprim/geom3d.hpp>
#include <linalg.hpp>

#include "core/exception.hpp"
#include "msghandler.hpp"

namespace netgen
{

  /*
    Classes for NETGEN
  */



  enum ELEMENT_TYPE : unsigned char { 
    SEGMENT = 1, SEGMENT3 = 2,
    TRIG = 10, QUAD=11, TRIG6 = 12, QUAD6 = 13, QUAD8 = 14,
    TET = 20, TET10 = 21, 
    PYRAMID = 22, PRISM = 23, PRISM12 = 24, PRISM15 = 27, PYRAMID13 = 28,
    HEX = 25, HEX20 = 26, HEX7 = 29
  };

  /*
  typedef int ELEMENT_EDGE[2];      // initial point, end point
  typedef int ELEMENT_FACE[4];      // points, last one is -1 for trig
  */

  struct ELEMENT_EDGE
  {
    int vals[2];
    int & operator[] (size_t i) { return vals[i]; }
    int operator[] (size_t i) const { return vals[i]; }
  };
  
  struct ELEMENT_FACE
  {
    int vals[4];
    int & operator[] (size_t i) { return vals[i]; }
    int operator[] (size_t i) const { return vals[i]; }
  };
  

#define ELEMENT_MAXPOINTS 20
#define ELEMENT2D_MAXPOINTS 8


  enum POINTTYPE : unsigned char { FIXEDPOINT = 1, EDGEPOINT = 2, SURFACEPOINT = 3, INNERPOINT = 4 };
  enum ELEMENTTYPE { FREEELEMENT, FIXEDELEMENT };
  enum OPTIMIZEGOAL { OPT_QUALITY, OPT_CONFORM, OPT_REST, OPT_WORSTCASE, OPT_LEGAL };


  extern DLL_HEADER size_t timestamp;
  inline size_t GetTimeStamp() 
  { 
    return timestamp; 
  }

  inline size_t NextTimeStamp()
  {
    timestamp++;
    return timestamp;
  }
  
  /*
  extern DLL_HEADER int GetTimeStamp();
  extern DLL_HEADER int NextTimeStamp();
  */
  class PointGeomInfo
  {
  public:
    int trignum;   // for STL Meshing
    double u, v;   // for OCC Meshing

    PointGeomInfo () = default;
    // : trignum(-1), u(0), v(0) { ; }
    PointGeomInfo (const PointGeomInfo&) = default;
    PointGeomInfo (PointGeomInfo &&) = default;
    PointGeomInfo & operator= (const PointGeomInfo&) = default;
    PointGeomInfo & operator= (PointGeomInfo&&) = default;
    
  };

  inline ostream & operator<< (ostream & ost, const PointGeomInfo & gi)
  {
    return (ost << gi.trignum << " " << gi.u << " " << gi.v);
  }

  inline istream & operator>> (istream & ist, PointGeomInfo & gi)
  {
    return (ist >> gi.trignum >> gi.u >> gi.v);
  }



  class MultiPointGeomInfo
  {
    ArrayMem<PointGeomInfo, 100> mgi;
  public:
    int AddPointGeomInfo (const PointGeomInfo & gi);
    void Init () { mgi.SetSize(0); }
    void DeleteAll () { mgi.SetSize(0); }

    int GetNPGI () const { return mgi.Size(); }
    const PointGeomInfo & GetPGI (int i) const { return mgi[i-1]; }

    MultiPointGeomInfo () = default;
    MultiPointGeomInfo (const MultiPointGeomInfo&) = default;
    MultiPointGeomInfo (MultiPointGeomInfo &&) = default;
    MultiPointGeomInfo & operator= (const MultiPointGeomInfo&) = delete;
    MultiPointGeomInfo & operator= (MultiPointGeomInfo&&) = default;
  };


  class EdgePointGeomInfo
  {
  public:
    int edgenr;
    int body;    // for ACIS
    double dist; // for 2d meshing
    double u, v; // for OCC Meshing

  public:
    EdgePointGeomInfo ()
      : edgenr(-1), body(0), dist(0.0), u(0.0), v(0.0) { ; }

    EdgePointGeomInfo & operator= (const EdgePointGeomInfo & gi2) = default;
  };

  inline ostream & operator<< (ostream & ost, const EdgePointGeomInfo & gi)
  {
    ost << "epgi: edgnr=" << gi.edgenr << ", dist=" << gi.dist;
    return ost;
  }





  class PointIndex
  {
    int i;
  public:
    class t_invalid { public: constexpr t_invalid() = default; };
    static constexpr t_invalid INVALID{};
    
    PointIndex () = default;
    PointIndex (const PointIndex&) = default;
    PointIndex (PointIndex &&) = default;
    PointIndex & operator= (const PointIndex&) = default;
    PointIndex & operator= (PointIndex&&) = default;
     
    constexpr PointIndex (int ai) : i(ai)
    {
#ifdef DEBUG
      if (ai < PointIndex::BASE)
        cout << "illegal PointIndex, use PointIndex::INVALID instead" << endl;
        // throw Exception("illegal PointIndex, use PointIndex::INVALID instead");
#endif
    }
    constexpr PointIndex (t_invalid inv) : i(PointIndex::BASE-1) { ; }
    // PointIndex & operator= (const PointIndex &ai) { i = ai.i; return *this; }
    constexpr operator int () const { return i; }
    PointIndex operator++ (int) { PointIndex hi(*this); i++; return hi; }
    PointIndex operator-- (int) { PointIndex hi(*this); i--; return hi; }
    PointIndex & operator++ () { i++; return *this; }
    PointIndex operator-- () { i--; return *this; }
    PointIndex operator+= (int add) { i += add; return *this; }
    void Invalidate() { i = PointIndex::BASE-1; }
    bool IsValid() const { return i != PointIndex::BASE-1; }
#ifdef BASE0
    static constexpr size_t BASE = 0;
#else
    static constexpr size_t BASE = 1;
#endif  

    void DoArchive (Archive & ar) { ar & i; }
  };

}

namespace ngcore
{
  template<> 
  constexpr netgen::PointIndex IndexBASE<netgen::PointIndex> () { return netgen::PointIndex(netgen::PointIndex::BASE); }
}

namespace netgen
{

  
  inline istream & operator>> (istream & ist, PointIndex & pi)
  {
    int i; ist >> i; pi = PointIndex(i); return ist;
  }

  inline ostream & operator<< (ostream & ost, const PointIndex & pi)
  {
    return (ost << int(pi));
  }

  template <int N> class PointIndices;
  template <> class PointIndices<2> : public INDEX_2
  {
  public:
    PointIndices () = default;
    PointIndices (INDEX_2 i2) : INDEX_2(i2) { ; }
    PointIndices (PointIndex i1, PointIndex i2) : INDEX_2(i1,i2) { ; } 
    PointIndex operator[] (int i) const { return PointIndex(INDEX_2::operator[](i)); }
    PointIndex & operator[] (int i) { return reinterpret_cast<PointIndex&>(INDEX_2::operator[](i)); }
    static PointIndices Sort(PointIndex i1, PointIndex i2) { return INDEX_2::Sort(i1, i2); }
    template <size_t J>
    PointIndex get() const { return PointIndex(INDEX_2::operator[](J)); }    
  };
}

namespace std
{
  // structured binding support
  template <auto N>
  struct tuple_size<netgen::PointIndices<N>> : std::integral_constant<std::size_t, N> {};
  template<size_t N, auto M> struct tuple_element<N,netgen::PointIndices<M>> { using type = netgen::PointIndex; };
}

namespace netgen
{

  class ElementIndex
  {
    int i;
  public:
    ElementIndex () = default;
    constexpr ElementIndex (int ai) : i(ai) { ; }
    ElementIndex & operator= (const ElementIndex & ai) { i = ai.i; return *this; }
    ElementIndex & operator= (int ai) { i = ai; return *this; }
    constexpr operator int () const { return i; }
    ElementIndex operator++ (int) { return ElementIndex(i++); }    
    ElementIndex operator-- (int) { return ElementIndex(i--); }
    ElementIndex & operator++ () { ++i; return *this; }
    ElementIndex & operator-- () { --i; return *this; }
  };

  inline istream & operator>> (istream & ist, ElementIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const ElementIndex & pi)
  {
    return (ost << int(pi));
  }


  class SurfaceElementIndex
  {
    int i;
  public:
    SurfaceElementIndex () = default;
    constexpr SurfaceElementIndex (int ai) : i(ai) { ; }
    /*
    SurfaceElementIndex & operator= (const SurfaceElementIndex & ai) 
    { i = ai.i; return *this; }
    */
    SurfaceElementIndex & operator= (const SurfaceElementIndex & ai) = default;
    SurfaceElementIndex & operator= (int ai) { i = ai; return *this; }
    constexpr operator int () const { return i; }
    SurfaceElementIndex operator++ (int) { SurfaceElementIndex hi(*this); i++; return hi; }
    SurfaceElementIndex operator-- (int) { SurfaceElementIndex hi(*this); i--; return hi; }
    SurfaceElementIndex & operator++ () { ++i; return *this; }
    SurfaceElementIndex & operator-- () { --i; return *this; }
    SurfaceElementIndex & operator+= (int inc) { i+=inc; return *this; }
    void DoArchive (Archive & ar) { ar & i; }
  };
  
  inline void SetInvalid (SurfaceElementIndex & id) { id = -1; }
  inline bool IsInvalid (SurfaceElementIndex & id) { return id == -1; }

  inline istream & operator>> (istream & ist, SurfaceElementIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const SurfaceElementIndex & pi)
  {
    return (ost << int(pi));
  }

  class SegmentIndex
  {
    int i;
  public:
    SegmentIndex () = default;
    constexpr SegmentIndex (int ai) : i(ai) { ; }
    SegmentIndex & operator= (const SegmentIndex & ai) 
    { i = ai.i; return *this; }
    SegmentIndex & operator= (int ai) { i = ai; return *this; }
    constexpr operator int () const { return i; }
    SegmentIndex& operator++ () { ++i; return *this; }
    SegmentIndex& operator-- () { --i; return *this; }
    SegmentIndex operator++ (int) { return i++; }
    SegmentIndex operator-- (int) { return i--; }
  };

  inline void SetInvalid (SegmentIndex & id) { id = -1; }
  inline bool IsInvalid (SegmentIndex & id) { return id == -1; }


  inline istream & operator>> (istream & ist, SegmentIndex & pi)
  {
    int i; ist >> i; pi = i; return ist;
  }

  inline ostream & operator<< (ostream & ost, const SegmentIndex & pi)
  {
    return (ost << int(pi));
  }




  /**
     Point in the mesh.
     Contains layer (a new feature in 4.3 for overlapping meshes.
  */
  class MeshPoint : public Point<3>
  {
    double singular; // singular factor for hp-refinement
    int layer;
    POINTTYPE type;


  public:
    MeshPoint () 
    { 
      ;
    }

    MeshPoint (const Point<3> & ap, int alayer = 1, POINTTYPE apt = INNERPOINT)
      : Point<3> (ap), singular(0.), layer(alayer), type(apt) 
    { 
      ;
    }
  
    void SetPoint (const Point<3> & ap)
    { 
      Point<3>::operator= (ap); 
      layer = 0; 
      singular = 0; 
    }

    void Scale(double factor) { *testout << "before: " << x[0] << endl; x[0] *= factor; x[1] *= factor; x[2] *= factor; *testout << "after: " << x[0] << endl;}

    int GetLayer() const { return layer; }

    POINTTYPE Type() const { return type; }
    void SetType(POINTTYPE at) { type = at; }
 
    double Singularity() const { return singular; }
    void Singularity(double s) { singular = s; }
    bool IsSingular() const { return (singular != 0.0); }

#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType ( );
#endif

    void DoArchive (Archive & ar)
    {
      // ar & x[0] & x[1] & x[2] & layer & singular;
      // ar.Do(&x[0], 3);
      // ar & layer & singular;
      // ar & (unsigned char&)(type);
      ar.DoPacked (x[0], x[1], x[2], layer, singular, (unsigned char&)(type));
    }
  };

  inline ostream & operator<<(ostream  & s, const MeshPoint & pt)
  { 
    return (s << Point<3> (pt)); 
  }




  // typedef NgArray<MeshPoint, PointIndex::BASE, PointIndex> T_POINTS;
  typedef Array<MeshPoint, PointIndex> T_POINTS;



  /**
     Triangle element for surface mesh generation.
  */
  class Element2d
  { 
    /// point numbers
    PointIndex pnum[ELEMENT2D_MAXPOINTS];
    /// geom info of points
    PointGeomInfo geominfo[ELEMENT2D_MAXPOINTS];

    /// surface nr
    int index;
    ///
    ELEMENT_TYPE typ;
    /// number of points
    int8_t np;
    bool refflag;  // marked for refinement
    bool badel:1;
    bool strongrefflag:1;
    bool deleted:1;  // element is deleted

    // Philippose - 08 August 2010
    // Set a new property for each element, to 
    // control whether it is visible or not
    bool visible:1;  // element visible
    bool is_curved;   // element is (high order) curved
    int8_t newest_vertex = -1; // from refinement via bisection
    /// order for hp-FEM
    unsigned int orderx:6;
    unsigned int ordery:6;

    /// a linked list for all segments in the same face
    SurfaceElementIndex next;
    ///
    int hp_elnr;

  public:
    static auto GetDataLayout()
    {
      return std::map<string, int>({
          { "pnum", offsetof(Element2d, pnum)},
          { "index", offsetof(Element2d, index) },
          { "np", offsetof(Element2d, np) },
          { "refine", offsetof(Element2d, refflag) },
          { "curved", offsetof(Element2d, is_curved)}
        });
    }

    ///
    DLL_HEADER Element2d ();
    Element2d (const Element2d &) = default;
    Element2d (Element2d &&) = default;
    Element2d & operator= (const Element2d &) = default;
    Element2d & operator= (Element2d &&) = default;
    Element2d & operator= (initializer_list<PointIndex> list)
    {
      size_t cnt = 0;
      for (auto val : list)
        pnum[cnt++] = val;
      return *this;
    }
    Element2d & operator= (initializer_list<std::tuple<PointIndex,PointGeomInfo>> list)
    {
      size_t cnt = 0;
      for (auto val : list)
        {
          pnum[cnt] = get<0>(val);
          geominfo[cnt++] = get<1>(val);
        }
      return *this;
    }
    ///
    DLL_HEADER Element2d (int anp);
    ///
    DLL_HEADER Element2d (ELEMENT_TYPE type);
    ///
    DLL_HEADER Element2d (int pi1, int pi2, int pi3);
    ///
    DLL_HEADER Element2d (int pi1, int pi2, int pi3, int pi4);
    ///
    ELEMENT_TYPE GetType () const { return typ; }
    /// 
    void SetType (ELEMENT_TYPE atyp)
    {
      typ = atyp;
      switch (typ)
	{
	case TRIG: np = 3; break;
	case QUAD: np = 4; break;
	case TRIG6: np = 6; break;
	case QUAD6: np = 6; break;
	case QUAD8: np = 8; break;
	default:
	  PrintSysError ("Element2d::SetType, illegal type ", int(typ));
	}
      is_curved = (np >= 4); 
    }
    ///
    int GetNP() const { return np; }
    ///
    int GetNV() const
    {
      if (typ == TRIG || typ == TRIG6)
        return 3;
      else
        {
#ifdef DEBUG
          if (typ != QUAD && typ != QUAD6 && typ != QUAD8)
            PrintSysError ("element2d::GetNV not implemented for typ", int(typ));
#endif
          return 4;
        }
      /*
      switch (typ)
	{
	case TRIG:
	case TRIG6: return 3;
          
	case QUAD:
	case QUAD8:
	case QUAD6: return 4;
	default:
#ifdef DEBUG
	  PrintSysError ("element2d::GetNV not implemented for typ", typ)
#endif
	    ;
	}
      return np;
      */
    }

    ///
    PointIndex & operator[] (int i) { return pnum[i]; }
    ///
    const PointIndex & operator[] (int i) const { return pnum[i]; }

    auto PNums () const { return FlatArray<const PointIndex> (np, &pnum[0]); }
    auto PNums ()  { return FlatArray<PointIndex> (np, &pnum[0]); }
    template <int NP>
    auto PNums() const { return FlatArray<const PointIndex> (NP, &pnum[0]); }
    auto Vertices() const { return FlatArray<const PointIndex> (GetNV(), &pnum[0]); }

    auto GeomInfo() const { return FlatArray<const PointGeomInfo> (np, &geominfo[0]); }
    auto GeomInfo() { return FlatArray<PointGeomInfo> (np, &geominfo[0]); }
    
    ///
    PointIndex & PNum (int i) { return pnum[i-1]; }
    ///
    const PointIndex & PNum (int i) const { return pnum[i-1]; }
    ///
    PointIndex & PNumMod (int i) { return pnum[(i-1) % np]; }
    ///
    const PointIndex & PNumMod (int i) const { return pnum[(i-1) % np]; }
    ///

    ///
    PointGeomInfo & GeomInfoPi (int i) { return geominfo[i-1]; }
    ///
    const PointGeomInfo & GeomInfoPi (int i) const { return geominfo[i-1]; }
    ///
    PointGeomInfo & GeomInfoPiMod (int i) { return geominfo[(i-1) % np]; }
    ///
    const PointGeomInfo & GeomInfoPiMod (int i) const { return geominfo[(i-1) % np]; }

    auto & NewestVertex() { return newest_vertex; }
    auto NewestVertex() const { return newest_vertex; }

    void DoArchive (Archive & ar)
    {
      short _np, _typ;
      bool _curved, _vis, _deleted;
      if (ar.Output())
        { _np = np; _typ = typ; _curved = is_curved;
          _vis = visible; _deleted = deleted; }
      // ar & _np & _typ & index & _curved & _vis & _deleted;
      ar.DoPacked (_np, _typ, index, _curved, _vis, _deleted);
      // ar & next; don't need 
      if (ar.Input())
        { np = _np; typ = ELEMENT_TYPE(_typ); is_curved = _curved;
          visible = _vis; deleted = _deleted; }
      /*
      for (size_t i = 0; i < np; i++)
        ar & pnum[i];
      */
      static_assert(sizeof(int) == sizeof (PointIndex));
      ar.Do( (int*)&pnum[0], np);
    }

#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif
    

    void SetIndex (int si) { index = si; }
    ///
    int GetIndex () const { return index; }

    int GetOrder () const { return orderx; }
    void SetOrder (int aorder) { orderx = ordery = aorder; }


    void GetOrder (int & ox, int & oy) const { ox = orderx, oy =ordery;};
    void GetOrder (int & ox, int & oy, int & oz) const { ox = orderx; oy = ordery; oz=0; }
    void SetOrder (int ox, int oy, int  /* oz */) { orderx = ox; ordery = oy;}
    void SetOrder (int ox, int oy) { orderx = ox; ordery = oy;}

    int GetHpElnr() const { return hp_elnr; }
    void SetHpElnr(int _hp_elnr) { hp_elnr = _hp_elnr; }

    ///
    void GetBox (const T_POINTS & points, Box3d & box) const;
    /// invert orientation
    inline void Invert ();
    ///
    DLL_HEADER void Invert2 ();
    /// first point number is smallest
    inline void NormalizeNumbering ();
    ///
    void NormalizeNumbering2 ();

    bool BadElement() const { return badel; }

    // friend ostream & operator<<(ostream  & s, const Element2d & el);
    friend class Mesh;


    /// get number of 'integration points'
    int GetNIP () const;
    void GetIntegrationPoint (int ip, Point<2> & p, double & weight) const;

    void GetTransformation (int ip, const NgArray<Point<2>> & points,
			    class DenseMatrix & trans) const;
    void GetTransformation (int ip, class DenseMatrix & pmat,
			    class DenseMatrix & trans) const;

    void GetShape (const Point<2> & p, class Vector & shape) const;
    DLL_HEADER void GetShapeNew (const Point<2> & p, class FlatVector & shape) const;
    template <typename T>
    DLL_HEADER void GetShapeNew (const Point<2,T> & p, TFlatVector<T> shape) const;
    /// matrix 2 * np
    DLL_HEADER void GetDShape (const Point<2> & p, class DenseMatrix & dshape) const;
    template <typename T>
    DLL_HEADER void GetDShapeNew (const Point<2,T> & p, class MatrixFixWidth<2,T> & dshape) const;
    
    /// matrix 2 * np
    void GetPointMatrix (const NgArray<Point<2>> & points,
			 class DenseMatrix & pmat) const; 

    void ComputeIntegrationPointData () const;
  

    double CalcJacobianBadness (const NgArray<Point<2>> & points) const;
    double CalcJacobianBadness (const T_POINTS & points, 
				const Vec<3> & n) const;
    double CalcJacobianBadnessDirDeriv (const NgArray<Point<2>> & points,
					int pi, Vec<2> & dir, double & dd) const;


    
    void Delete ()
    {
      deleted = 1;
      // for (PointIndex & p : pnum) p.Invalidate(); 
    }
    
    bool IsDeleted () const 
    {
#ifdef DEBUG
      if (pnum[0] < PointIndex::BASE && !deleted)
	cerr << "Surfelement has illegal pnum, but not marked as deleted" << endl;
#endif    
      return deleted; 
    }

    // Philippose - 08 August 2010
    // Access functions for the new property: visible
    void Visible(bool vis = true) 
    { visible = vis; }
    bool IsVisible () const 
    { return visible; }
   
    void SetRefinementFlag (bool rflag = true) 
    { refflag = rflag; }
    bool TestRefinementFlag () const
    { return refflag; }

    void SetStrongRefinementFlag (bool rflag = true) 
    { strongrefflag = rflag; }
    bool TestStrongRefinementFlag () const
    { return strongrefflag; }


    bool IsCurved () const { return is_curved; }
    void SetCurved (bool acurved) { is_curved = acurved; }
  
    SurfaceElementIndex NextElement() { return next; }

    bool operator==(const Element2d & el2) const;

    int HasFace(const Element2d& el) const;
  };

  ostream & operator<<(ostream  & s, const Element2d & el);





  class IntegrationPointData
  {
  public:
    Point<3> p;
    double weight;
    Vector shape;
    DenseMatrix dshape;
  };








  /**
     Volume element
  */
  class Element
  {
  private:
    /// point numbers
    PointIndex pnum[ELEMENT_MAXPOINTS];
    ///
    ELEMENT_TYPE typ;
    /// number of points (4..tet, 5..pyramid, 6..prism, 8..hex, 10..quad tet, 12..quad prism)
    int8_t np;
    int8_t newest_vertex = -1; // from refinement via bisection
    
    /// sub-domain index
    int index;
    /// order for hp-FEM
    unsigned int orderx:6;
    unsigned int ordery:6;
    unsigned int orderz:6;
    /* unsigned int levelx:6;
       unsigned int levely:6;
       unsigned int levelz:6; */ 
    /// stored shape-badness of element
    float badness;
    bool is_curved;   // element is (high order) curved

    class flagstruct {
    public:
      bool refflag;     // mark element for refinement
      bool marked:1;  // marked for refinement
      bool badel:1;   // angles worse then limit
      bool reverse:1; // for refinement a la Bey
      bool illegal:1; // illegal, will be split or swapped
      bool illegal_valid:1; // is illegal-flag valid ?
      bool badness_valid:1; // is badness valid ?
      bool strongrefflag:1;
      bool deleted:1;   // element is deleted, will be removed from array
      bool fixed:1;     // don't change element in optimization
    };

    flagstruct flags;
    int hp_elnr;
  public:

    static auto GetDataLayout()
    {
      return std::map<string, int>({
          { "pnum", offsetof(Element, pnum)},
          { "index", offsetof(Element, index) },
          { "np", offsetof(Element, np) },
          { "refine", offsetof(Element, flags.refflag) },
          { "curved", offsetof(Element, is_curved)}
        });
    }

    ///
    DLL_HEADER Element () = default;
    Element (const Element &) = default;
    Element (Element &&) = default;
    Element & operator= (const Element &) = default;
    Element & operator= (Element &&) = default;

    ///
    DLL_HEADER Element (int anp);
    ///
    DLL_HEADER Element (ELEMENT_TYPE type);
    ///
    // Element & operator= (const Element & el2);

    const flagstruct& Flags() const { return flags; }
    flagstruct& Flags() { return flags; }
  
    ///
    DLL_HEADER void SetNP (int anp);
    ///
    DLL_HEADER void SetType (ELEMENT_TYPE atyp);
    ///
    int GetNP () const { return np; }
    ///
    uint8_t GetNV() const
    {
      // __assume(typ >= TET && typ <= PYRAMID13);
      switch (typ)
	{
        case TET: 
        case TET10: 
          return 4;
        case PRISM12:
        case PRISM15:
        case PRISM:
	  return 6; 
	case PYRAMID:
        case PYRAMID13:
	  return 5;
	case HEX7:
	  return 7;
	case HEX:
	case HEX20:
	  return 8;
        default: // not a 3D element
#ifdef DEBUG
          PrintSysError ("Element3d::GetNV not implemented for typ ", int(typ));
#endif
          __assume(false);
          return -1;
        }
    }

    DLL_HEADER bool operator==(const Element & el2) const;

    // old style:
    int NP () const { return np; }

    ///
    ELEMENT_TYPE GetType () const { return typ; }

    ///
    PointIndex & operator[] (int i) { return pnum[i]; }
    ///
    const PointIndex & operator[] (int i) const { return pnum[i]; }

    auto PNums () const { return FlatArray<const PointIndex> (np, &pnum[0]); }
    auto PNums () { return FlatArray<PointIndex> (np, &pnum[0]); }    
    template <int NP>
    auto PNums() const { return FlatArray<const PointIndex> (NP, &pnum[0]); }

    FlatArray<const PointIndex> Vertices() const { return { GetNV(), &pnum[0] }; }

    ///
    PointIndex & PNum (int i) { return pnum[i-1]; }
    ///
    const PointIndex & PNum (int i) const { return pnum[i-1]; }
    ///
    PointIndex & PNumMod (int i) { return pnum[(i-1) % np]; }
    ///
    const PointIndex & PNumMod (int i) const { return pnum[(i-1) % np]; }

    auto & NewestVertex() { return newest_vertex; }
    auto NewestVertex() const { return newest_vertex; }

    void DoArchive (Archive & ar)
    {
      short _np, _typ;
      bool _curved;
      if (ar.Output())
        { _np = np; _typ = typ; _curved = is_curved; }
      // ar & _np & _typ & index & _curved;
      ar.DoPacked (_np, _typ, index, _curved);                

      if (ar.Input())
        {
          np = _np;
          typ = ELEMENT_TYPE(_typ);
          is_curved = _curved;
          flags.marked = 1;
          flags.badel = 0;
          flags.reverse = 0;
          flags.illegal = 0;
          flags.illegal_valid = 0;
          flags.badness_valid = 0;
          flags.refflag = 1;
          flags.strongrefflag = false;
          flags.deleted = 0;
          flags.fixed = 0;
        }

      static_assert(sizeof(int) == sizeof (PointIndex));
      ar.Do( (int*)&pnum[0], np);
    }
    
#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif

    ///
    void SetIndex (int si) { index = si; }
    ///
    int GetIndex () const { return index; }

    int GetOrder () const { return orderx; }
    void SetOrder (const int aorder) ; 

    void GetOrder (int & ox, int & oy, int & oz) const { ox = orderx; oy = ordery; oz = orderz; }
    void SetOrder (const int ox, const int oy, const int oz);
    // void GetLevel (int & ox, int & oy, int & oz) const { ox = levelx; oy = levely; oz = levelz; }
    // void SetLevel (int ox, int oy, int oz) { levelx = ox; levely = oy; levelz = oz; }


    ///
    void GetBox (const T_POINTS & points, Box3d & box) const;
    /// Calculates Volume of element
    double Volume (const T_POINTS & points) const;
    ///
    DLL_HEADER void Print (ostream & ost) const;
    ///
    int GetNFaces () const
    {
      switch (typ)
	{
	case TET: 
	case TET10: return 4;
	case PYRAMID: case PYRAMID13: return 5;
	case PRISM:
        case PRISM15:
	case PRISM12: return 5;
        case HEX7: return 6;
        case HEX: case HEX20:
          return 6;
	default:
#ifdef DEBUG
	  PrintSysError ("element3d::GetNFaces not implemented for typ", int(typ))
#endif
	    ;
	}
      return 0;
    }
    ///
    inline void GetFace (int i, Element2d & face) const;
    ///
    DLL_HEADER void GetFace2 (int i, Element2d & face) const;
    ///
    DLL_HEADER void Invert ();

    int GetHpElnr() const { return hp_elnr; }
    void SetHpElnr(int _hp_elnr) { hp_elnr = _hp_elnr; }

    /// split into 4 node tets
    void GetTets (NgArray<Element> & locels) const;
    /// split into 4 node tets, local point nrs
    void GetTetsLocal (NgArray<Element> & locels) const;
    /// returns coordinates of nodes
    // void GetNodesLocal (NgArray<Point<3> > & points) const;
    void GetNodesLocalNew (NgArray<Point<3> > & points) const;

    /// split surface into 3 node trigs
    DLL_HEADER void GetSurfaceTriangles (NgArray<Element2d> & surftrigs) const;


    /// get number of 'integration points'
    int GetNIP () const;
    void GetIntegrationPoint (int ip, Point<3> & p, double & weight) const;

    void GetTransformation (int ip, const T_POINTS & points,
			    class DenseMatrix & trans) const;
    void GetTransformation (int ip, class DenseMatrix & pmat,
			    class DenseMatrix & trans) const;

    void GetShape (const Point<3> & p, class Vector & shape) const;
    // void GetShapeNew (const Point<3> & p, class FlatVector & shape) const;
    template <typename T>
    DLL_HEADER void GetShapeNew (const Point<3,T> & p, TFlatVector<T> shape) const;
    /// matrix 2 * np
    void GetDShape (const Point<3> & p, class DenseMatrix & dshape) const;
    template <typename T>
    void GetDShapeNew (const Point<3,T> & p, class MatrixFixWidth<3,T> & dshape) const;
    /// matrix 3 * np
    void GetPointMatrix (const T_POINTS & points,
			 class DenseMatrix & pmat) const; 

    void ComputeIntegrationPointData () const;
  

    double CalcJacobianBadness (const T_POINTS & points) const;
    double CalcJacobianBadnessDirDeriv (const T_POINTS & points,
					int pi, Vec<3> & dir, double & dd) const;
    double CalcJacobianBadnessGradient (const T_POINTS & points,
					int pi, Vec<3> & grad) const;

    ///
    // friend ostream & operator<<(ostream  & s, const Element & el);

    void SetRefinementFlag (bool rflag = 1) 
    { flags.refflag = rflag; }
    int TestRefinementFlag () const
    { return flags.refflag; }

    void SetStrongRefinementFlag (bool rflag = 1) 
    { flags.strongrefflag = rflag; }
    int TestStrongRefinementFlag () const
    { return flags.strongrefflag; }

    int Illegal () const
    {
      NETGEN_CHECK_SAME(flags.illegal_valid, true);
      return flags.illegal;
    }
    int IllegalValid () const
    { return flags.illegal_valid; }
    void SetIllegal (int aillegal)
    {
      flags.illegal = aillegal ? 1 : 0;
      flags.illegal_valid = 1;
    }
    void SetLegal (int alegal)
    {
      flags.illegal = alegal ? 0 : 1;
      flags.illegal_valid = 1;
    }

    bool BadnessValid()
    { return flags.badness_valid; }

    float GetBadness()
    {
      NETGEN_CHECK_SAME(flags.badness_valid, true);
      return badness;
    }

    void SetBadness(float value)
    {
      badness = value;
      flags.badness_valid = 1;
    }

    void Touch() {
      flags.illegal_valid = 0;
      flags.badness_valid = 0;
    }
  
    void Delete () { flags.deleted = 1; }
    bool IsDeleted () const 
    { 
#ifdef DEBUG
      if (pnum[0] < PointIndex::BASE && !flags.deleted)
	cerr << "Volelement has illegal pnum, but not marked as deleted" << endl;
#endif    

      return flags.deleted; 
    }

    bool IsCurved () const { return is_curved; }
    void SetCurved (bool acurved) { is_curved = acurved; }

  };

  ostream & operator<<(ostream  & s, const Element & el);






  /**
     Edge segment.
  */
  class Segment
  {
  public:
    ///
    DLL_HEADER Segment();
    Segment (const Segment& other) = default;

    // friend ostream & operator<<(ostream  & s, const Segment & seg);

    PointIndex pnums[3];  // p1, p2, pmid

    int edgenr;
    ///
    double singedge_left;
    double singedge_right;

    /// 0.. not first segment of segs, 1..first of class, 2..first of class, inverse
    unsigned int seginfo:2;

    /// surface decoding index
    int si;
    /// co dim 2 decoding index
    int cd2i;
    /// domain number inner side
    int domin;
    /// domain number outer side
    int domout;  
    /// top-level object number of surface
    int tlosurf;
    ///
    PointGeomInfo geominfo[2];

    /// surfaces describing edge
    int surfnr1, surfnr2;
    ///
    EdgePointGeomInfo epgeominfo[2];
    ///
    // int pmid; // for second order
    ///
    int meshdocval;

    bool is_curved;
    int hp_elnr;
    /*
      PointIndex operator[] (int i) const
      { return (i == 0) ? p1 : p2; }

      PointIndex & operator[] (int i) 
      { return (i == 0) ? p1 : p2; }
    */

    Segment& operator=(const Segment & other) = default;

  

    int GetNP() const
    {
      return pnums[2].IsValid() ? 3 : 2;
    }

    auto PNums() const { return FlatArray<const PointIndex> (GetNP(), &pnums[0]); }
    auto PNums() { return FlatArray<PointIndex> (GetNP(), &pnums[0]); }
    
    
    ELEMENT_TYPE GetType() const
    {
      return pnums[2].IsValid() ? SEGMENT3 : SEGMENT;
    }
  
    PointIndex & operator[] (int i) { return pnums[i]; }
    const PointIndex & operator[] (int i) const { return pnums[i]; }


    bool IsCurved () const { return is_curved; }
    void SetCurved (bool acurved) { is_curved = acurved; }
    
    void DoArchive (Archive & ar);
#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif
    
  };

  ostream & operator<<(ostream  & s, const Segment & seg);


  class Element0d
  {
  public:
    PointIndex pnum;
    string name;
    int index;
    Element0d () = default;
    Element0d (PointIndex _pnum, int _index)
      : pnum(_pnum), index(_index) { ; }

#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif
    
    void DoArchive (Archive & ar);
  };

  ostream & operator<<(ostream  & s, const Element0d & el);

  // class Surface;  
  // class FaceDescriptor;

  ///
  class FaceDescriptor
  {
    /// which surface, 0 if not available
    int surfnr;
    /// domain nr inside
    int domin;
    /// domain nr outside
    int domout;
    /// top level object number of surface
    int tlosurf;
    /// boundary condition property
    int bcprop;
    // Philippose - 06/07/2009
    // Add capability to store surface colours along with 
    // other face data
    /// surface colour (Default: R=0.0 ; G=1.0 ; B=0.0)
    Vec<4> surfcolour;
    
    ///
    // static string default_bcname;
    // string * bcname = &default_bcname;
    string bcname = "default";
    /// root of linked list 
    SurfaceElementIndex firstelement;
  
    double domin_singular;
    double domout_singular;

  public:
    DLL_HEADER FaceDescriptor();
    DLL_HEADER FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi);
    DLL_HEADER FaceDescriptor(const Segment & seg);
    DLL_HEADER FaceDescriptor(const FaceDescriptor& other);
    DLL_HEADER ~FaceDescriptor()  { ; }

    DLL_HEADER int SegmentFits (const Segment & seg);

    int SurfNr () const { return surfnr; }
    int DomainIn () const { return domin; }
    int DomainOut () const { return domout; }
    int TLOSurface () const { return tlosurf; }
    int BCProperty () const { return bcprop; }


    double DomainInSingular() const { return domin_singular; }
    double DomainOutSingular() const { return domout_singular; }

    // Philippose - 06/07/2009
    // Get Surface colour
    Vec<4> SurfColour () const { return surfcolour; }
    /* DLL_HEADER */ const string & GetBCName () const { return bcname; }
    // string * BCNamePtr () { return bcname; }
    // const string * BCNamePtr () const  { return bcname; }
    void SetSurfNr (int sn) { surfnr = sn; }
    void SetDomainIn (int di) { domin = di; }
    void SetDomainOut (int dom) { domout = dom; }
    void SetBCProperty (int bc) { bcprop = bc; }
    DLL_HEADER void SetBCName (string * bcn); //  { bcname = bcn; }
    void SetBCName (const string & bcn) { bcname = bcn; }    
    // Philippose - 06/07/2009
    // Set the surface colour
    void SetSurfColour (Vec<4> colour) { surfcolour = colour; }

    void SetDomainInSingular (double v) { domin_singular = v; }
    void SetDomainOutSingular (double v) { domout_singular = v; }

    SurfaceElementIndex FirstElement() { return firstelement; }
    // friend ostream & operator<<(ostream  & s, const FaceDescriptor & fd);
    friend class Mesh;

    void DoArchive (Archive & ar);
  };

  ostream & operator<< (ostream  & s, const FaceDescriptor & fd);

  
 
  class EdgeDescriptor
  {
    int tlosurf;
    int surfnr[2];
  public:
    EdgeDescriptor ()
      : tlosurf(-1)
    { surfnr[0] = surfnr[1] = -1; }

    int SurfNr (int i) const { return surfnr[i]; }
    void SetSurfNr (int i, int nr) { surfnr[i] = nr; }

    int TLOSurface() const { return tlosurf; }
    void SetTLOSurface (int nr) { tlosurf = nr; }
  };



  class DLL_HEADER MeshingParameters
  {
  public:
    /**
       3d optimization strategy:
       // m .. move nodes
       // M .. move nodes, cheap functional
       // s .. swap faces
       // c .. combine elements
       // d .. divide elements
       // D .. divide and join opposite edges, remove element
       // p .. plot, no pause
       // P .. plot, Pause
       // h .. Histogramm, no pause
       // H .. Histogramm, pause
       */
    string optimize3d = "cmdDmustm";
    /// number of 3d optimization steps
    int optsteps3d = 3;
    /**
       2d optimization strategy:
       // s .. swap, opt 6 lines/node
       // S .. swap, optimal elements
       // m .. move nodes
       // p .. plot, no pause
       // P .. plot, pause
       // c .. combine
       **/
    string optimize2d = "smcmSmcmSmcm";
    /// number of 2d optimization steps
    int optsteps2d = 3;
    /// power of error (to approximate max err optimization)
    double opterrpow = 2;
    /// do block filling ?  
    bool blockfill = true;
    /// block filling up to distance
    double filldist = 0.1;
    /// radius of local environment (times h)
    double safety = 5;
    /// radius of active environment (times h)
    double relinnersafety = 3;
    /// use local h ?
    bool uselocalh = true;
    /// grading for local h
    double grading = 0.3;
    /// use delaunay for 3d meshing
    bool delaunay = true;
    /// use delaunay for 2d meshing
    bool delaunay2d = false;
    /// maximal mesh size
    double maxh = 1e10;
    /// minimal mesh size
    double minh = 0.0;
    /// file for meshsize
    string meshsizefilename = "";
    /// restrict h based on close edges
    optional<double> closeedgefac = nullopt;
    /// start surfacemeshing from everywhere in surface
    bool startinsurface = false;
    /// check overlapping surfaces (debug)
    bool checkoverlap = true;
    /// check overlapping surface mesh before volume meshing
    bool checkoverlappingboundary = true;
    /// check chart boundary (sometimes too restrictive)
    bool checkchartboundary = true;
    /// safety factor for curvatures (elements per radius)
    double curvaturesafety = 2;
    /// minimal number of segments per edge
    double segmentsperedge = 1;
    /// use parallel threads
    bool parthread = 0;
    /// weight of element size w.r.t element shape
    double elsizeweight = 0.2;
    /// init with default values

    /// start at step
    int perfstepsstart = 0;
    /// end at step
    int perfstepsend = 6;


    /// from mp3:
    /// give up quality class, 2d meshing
    int giveuptol2d = 200;
    /// give up quality class, 3d meshing
    int giveuptol = 10;
    /// give up quality class for closing open quads, > 100 for
    /// free pyramids
    int giveuptolopenquads = 15;
    /// maximal outer steps
    int maxoutersteps = 10;
    /// class starting star-shape filling
    int starshapeclass = 5;
    /// if non-zero, baseelement must have baseelnp points
    int baseelnp = 0;        
    /// quality tolerances are handled less careful
    int sloppy = 1;
  
    /// limit for max element angle (150-180)
    double badellimit = 175;

    bool check_impossible = false;

    int only3D_domain_nr = 0;
  
    ///
    bool secondorder = false;
    /// high order element curvature
    int elementorder = 1;
    /// quad-dominated surface meshing
    bool quad = false;
    ///
    bool try_hexes = false;
    ///
    bool inverttets = false;
    ///
    bool inverttrigs = false;
    ///
    bool autozrefine = false;

    bool parallel_meshing = true;
    int nthreads = 4;

    Flags geometrySpecificParameters;
    ///
    MeshingParameters ();
    ///
    MeshingParameters (const MeshingParameters & mp2) = default;
    MeshingParameters (MeshingParameters && mp2) = default;
    MeshingParameters & operator= (const MeshingParameters & mp2) = default;
    MeshingParameters & operator= (MeshingParameters && mp2) = default;
    ///
    void Print (ostream & ost) const;
    /// 
    // void CopyFrom(const MeshingParameters & other);

    class MeshSizePoint
    {
    public:
      Point<3> pnt;
      double h;
      int layer = 1;
      MeshSizePoint (Point<3> _pnt, double _h) : pnt(_pnt), h(_h) { ; }
      MeshSizePoint () = default;
      MeshSizePoint (const MeshSizePoint &) = default;
      MeshSizePoint (MeshSizePoint &&) = default;
      MeshSizePoint & operator= (const MeshSizePoint &) = default;
      MeshSizePoint & operator= (MeshSizePoint &&) = default;      
    };
    NgArray<MeshSizePoint> meshsize_points;
    
    void (*render_function)(bool) = NULL;
    void Render(bool blocking = false) const
    {
      if (render_function) 
        (*render_function)(blocking);
    }
  };

  inline ostream & operator<< (ostream & ost, const MeshingParameters & mp)
  {
    mp.Print (ost);
    return ost;
  }

  class DebugParameters 
  {
  public:
    ///
    int debugoutput;
    /// use slow checks
    int slowchecks;
    ///
    int haltsuccess;
    ///
    int haltnosuccess;
    ///
    int haltlargequalclass;
    ///
    int haltsegment;
    ///
    int haltnode;
    ///
    int haltsegmentp1;
    ///
    int haltsegmentp2;
    ///
    int haltexistingline;
    ///
    int haltoverlap;
    ///
    int haltface;
    ///
    int haltfacenr;
    ///
    bool write_mesh_on_error;
    ///
    DebugParameters ();
  };




  inline void Element2d :: Invert()
  {
    if (typ == TRIG)
      Swap (PNum(2), PNum(3));
    else
      Invert2();
  }




  inline void Element2d :: NormalizeNumbering ()
  {
    if (GetNP() == 3)
      {
	if (PNum(1) < PNum(2) && PNum(1) < PNum(3))
	  return;
	else
	  {
	    if (PNum(2) < PNum(3))
	      {
		PointIndex pi1 = PNum(2);
		PNum(2) = PNum(3);
		PNum(3) = PNum(1);
		PNum(1) = pi1;
	      }
	    else
	      {
		PointIndex pi1 = PNum(3);
		PNum(3) = PNum(2);
		PNum(2) = PNum(1);
		PNum(1) = pi1;
	      }
	  }
      }
    else
      NormalizeNumbering2();
  }



  static const int gftetfacesa[4][3] = 
    { { 1, 2, 3 },
      { 2, 0, 3 },
      { 0, 1, 3 },
      { 1, 0, 2 } };

  inline void Element :: GetFace (int i, Element2d & face) const
  {
    if (typ == TET)
      {
	face.SetType(TRIG);
	face[0] = pnum[gftetfacesa[i-1][0]];
	face[1] = pnum[gftetfacesa[i-1][1]];
	face[2] = pnum[gftetfacesa[i-1][2]];
      }
    else
      GetFace2 (i, face);
  }







  /**
     Identification of periodic surfaces, close surfaces, etc. 
  */
  class Identifications
  {
  public:
    enum ID_TYPE : unsigned char { UNDEFINED = 1, PERIODIC = 2, CLOSESURFACES = 3, CLOSEEDGES = 4};
  

  private:
    class Mesh & mesh;

    /// identify points (thin layers, periodic b.c.)  
    INDEX_2_HASHTABLE<int> identifiedpoints;
  
    /// the same, with info about the id-nr
    INDEX_3_HASHTABLE<int> identifiedpoints_nr;

    /// sorted by identification nr
    TABLE<INDEX_2> idpoints_table;

    NgArray<ID_TYPE> type;

    /// number of identifications (or, actually used identifications ?)
    int maxidentnr;
    Array<string> names;

  public:
    ///
    DLL_HEADER Identifications (class Mesh & amesh);
    ///
    DLL_HEADER ~Identifications ();

    DLL_HEADER void Delete ();

    /*
      Identify points pi1 and pi2, due to
      identification nr identnr
    */
    DLL_HEADER void Add (PointIndex pi1, PointIndex pi2, int identnr);
    void Add (PointIndex pi1, PointIndex pi2, string name, ID_TYPE type)
    {
        auto nr = GetNr(name);
        Add(pi1, pi2, nr);
        SetType(nr, type);
    }

    int Get (PointIndex pi1, PointIndex pi2) const;
    int GetSymmetric (PointIndex pi1, PointIndex pi2) const;

    bool Get (PointIndex pi1, PointIndex pi2, int identnr) const;
    bool GetSymmetric (PointIndex pi1, PointIndex pi2, int identnr) const;

    // bool HasIdentifiedPoints() const { return identifiedpoints != nullptr; } 
    ///
    INDEX_2_HASHTABLE<int> & GetIdentifiedPoints () 
    { 
      return identifiedpoints; 
    }

    bool Used (PointIndex pi1, PointIndex pi2)
    {
      return identifiedpoints.Used (INDEX_2 (pi1, pi2));
    }

    bool UsedSymmetric (PointIndex pi1, PointIndex pi2)
    {
      return 
	identifiedpoints.Used (INDEX_2 (pi1, pi2)) ||
	identifiedpoints.Used (INDEX_2 (pi2, pi1));
    }

    ///
    void GetMap (int identnr, NgArray<int,PointIndex::BASE> & identmap, bool symmetric = false) const;
    ///
    ID_TYPE GetType(int identnr) const
    {
      if(identnr <= type.Size())
	return type[identnr-1];
      else
	return UNDEFINED;
    }
    void SetType(int identnr, ID_TYPE t)
    {
      while(type.Size() < identnr)
	type.Append(UNDEFINED);
      type[identnr-1] = t;
    }
    
    ///
    DLL_HEADER void GetPairs (int identnr, NgArray<INDEX_2> & identpairs) const;
    ///
    int GetMaxNr () const { return maxidentnr; }  

    int GetNr(string name)
    {
      if(!names.Contains(name))
         names.Append(name);
      return names.Pos(name)+1;
    }

    /// remove secondorder
    void SetMaxPointNr (int maxpnum);

    DLL_HEADER void Print (ostream & ost) const;

    void DoArchive (Archive & ar);
  };
}


#ifdef PARALLEL
namespace ngcore
{
  template <> struct MPI_typetrait<netgen::PointIndex> {
    static NG_MPI_Datatype MPIType ()  { return NG_MPI_INT; }
  };

  template <> struct MPI_typetrait<netgen::ELEMENT_TYPE> {
    static NG_MPI_Datatype MPIType ()  { return NG_MPI_CHAR; }
  };

  template <> struct MPI_typetrait<netgen::MeshPoint> {
    static NG_MPI_Datatype MPIType ()  { return netgen::MeshPoint::MyGetMPIType(); }
  };

  template <> struct MPI_typetrait<netgen::Element> {
    static NG_MPI_Datatype MPIType ()  { return netgen::Element::MyGetMPIType(); }
  };
  template <> struct MPI_typetrait<netgen::Element2d> {
    static NG_MPI_Datatype MPIType ()  { return netgen::Element2d::MyGetMPIType(); }
  };
  template <> struct MPI_typetrait<netgen::Segment> {
    static NG_MPI_Datatype MPIType ()  { return netgen::Segment::MyGetMPIType(); }
  };
  template <> struct MPI_typetrait<netgen::Element0d> {
    static NG_MPI_Datatype MPIType ()  { return netgen::Element0d::MyGetMPIType(); }
  };

}
#endif


#endif


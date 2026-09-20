#ifndef MESHTYPE
#define MESHTYPE


/**************************************************************************/
/* File:   meshtype.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include <variant>

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


  using ELEMENT_EDGE = std::array<int,2>;
  using ELEMENT_FACE = std::array<int,4>;

  /// number of points / vertices / edges / faces of a 3D element type, 0 for other types
  namespace element3d_info
  {
    constexpr size_t SIZE = 32;   // > max ELEMENT_TYPE value
    constexpr std::array<int8_t,SIZE> np = [] {
      std::array<int8_t,SIZE> t{};
      t[TET] = 4;  t[PYRAMID] = 5;  t[PRISM] = 6;  t[HEX7] = 7;  t[HEX] = 8;
      t[TET10] = 10;  t[PRISM12] = 12;  t[PYRAMID13] = 13;  t[PRISM15] = 15;  t[HEX20] = 20;
      return t; } ();
    constexpr std::array<int8_t,SIZE> nv = [] {
      std::array<int8_t,SIZE> t{};
      t[TET] = t[TET10] = 4;  t[PYRAMID] = t[PYRAMID13] = 5;  t[PRISM] = t[PRISM12] = t[PRISM15] = 6;
      t[HEX7] = 7;  t[HEX] = t[HEX20] = 8;
      return t; } ();
    constexpr std::array<int8_t,SIZE> nedges = [] {
      std::array<int8_t,SIZE> t{};
      t[TET] = t[TET10] = 6;  t[PYRAMID] = t[PYRAMID13] = 8;  t[PRISM] = t[PRISM12] = t[PRISM15] = 9;
      t[HEX7] = 11;  t[HEX] = t[HEX20] = 12;
      return t; } ();
    constexpr std::array<int8_t,SIZE> nfaces = [] {
      std::array<int8_t,SIZE> t{};
      t[TET] = t[TET10] = 4;  t[PYRAMID] = t[PYRAMID13] = 5;  t[PRISM] = t[PRISM12] = t[PRISM15] = 5;
      t[HEX7] = 6;  t[HEX] = t[HEX20] = 6;
      return t; } ();
  }

  /// number of points / vertices of a 2D element type, 0 for other types
  namespace element2d_info
  {
    constexpr size_t SIZE = 32;
    constexpr std::array<int8_t,SIZE> np = [] {
      std::array<int8_t,SIZE> t{};
      t[TRIG] = 3;  t[QUAD] = 4;  t[TRIG6] = 6;  t[QUAD6] = 6;  t[QUAD8] = 8;
      return t; } ();
    constexpr std::array<int8_t,SIZE> nv = [] {
      std::array<int8_t,SIZE> t{};
      t[TRIG] = t[TRIG6] = 3;  t[QUAD] = t[QUAD6] = t[QUAD8] = 4;
      return t; } ();
  }

  /// the 2D element type with np points; 6 points give a TRIG6 (a QUAD6 is only made by type)
  constexpr ELEMENT_TYPE Element2dTypeFromNP (int np)
  {
    switch (np)
      {
      case 3: return TRIG;  case 4: return QUAD;  case 6: return TRIG6;  case 8: return QUAD8;
      }
    return ELEMENT_TYPE(0);
  }

  /// the 3D element type with np points (the mapping is one-to-one)
  constexpr ELEMENT_TYPE Element3dTypeFromNP (int np)
  {
    switch (np)
      {
      case 4: return TET;  case 5: return PYRAMID;  case 6: return PRISM;  case 7: return HEX7;  case 8: return HEX;
      case 10: return TET10;  case 12: return PRISM12;  case 13: return PYRAMID13;  case 15: return PRISM15;  case 20: return HEX20;
      }
    return ELEMENT_TYPE(0);
  }


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
  
  class PointGeomInfo
  {
  public:
    int trignum;   // for STL Meshing
    double u, v;   // for OCC Meshing

    PointGeomInfo () = default;
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
    PointGeomInfo gi;
    double dist; // parameter along edge curve

  public:
    EdgePointGeomInfo ()
      : dist(0.0) { gi.trignum = -1; gi.u = 0.0; gi.v = 0.0; }
    EdgePointGeomInfo (const EdgePointGeomInfo&) = default;
    EdgePointGeomInfo & operator= (const EdgePointGeomInfo & gi2) = default;
  };
  inline ostream & operator<< (ostream & ost, const EdgePointGeomInfo & gi)
  {
    ost << "epgi: dist=" << gi.dist;
    return ost;
  }


  template <typename T, typename TIndex, int BASE_>
  class Index
  {
  public:
    T i;

    static constexpr int BASE = BASE_;
    static constexpr TIndex Base() { return TIndex(BASE_); } 

    class t_invalid { public: constexpr t_invalid() = default; };
    static constexpr t_invalid INVALID{};

    typedef decltype( declval<T>()-declval<T>() ) T_diff;
    
  public:
    constexpr Index () = default;
    constexpr Index (const Index& i2) = default;
    constexpr Index (Index &&) = default;
    Index & operator= (const Index&) = default;
    Index & operator= (Index&&) = default;

    // private:
    constexpr Index (T ai) : i(ai)
    {
#ifdef DEBUG
      if (ai < BASE_)
        cout << "illegal Index, use Index::INVALID instead" << endl;
#endif
    }

  public:
    constexpr Index (t_invalid inv) : i(long(BASE)-1) { ; }
    // protected:
    constexpr operator T () const { return i; }
    explicit constexpr operator T& () { return i; }
  public:
    TIndex operator++ (int) { TIndex hi{*this}; i++; return hi; }
    TIndex operator-- (int) { TIndex hi(*this); i--; return hi; }
    TIndex & operator++ () { i++; return static_cast<TIndex&>(*this); }
    TIndex & operator-- () { i--; return static_cast<TIndex&>(*this); }
    constexpr TIndex operator+= (T_diff add) { i += add; return TIndex{*this}; }
    constexpr TIndex operator-= (T_diff add) { i -= add; return TIndex{*this}; }
    
    constexpr auto operator- (Index i2) const { return i-i2.i; }

    /// 0-based number of this index
    constexpr T Nr0 () const { return i - BASE_; }
    /// index for a 0-based number
    static constexpr TIndex FromNr0 (T i0) { return TIndex(T(BASE_) + i0); }
    /// 1-based number of this index (file formats, external interfaces)
    constexpr T Nr1 () const { return i - BASE_ + 1; }
    /// index for a 1-based number
    static constexpr TIndex FromNr1 (T i1) { return TIndex(T(BASE_) + i1 - 1); }
    /// the underlying integer, for diagnostics and low-level interfaces only
    constexpr T GetRawInteger () const { return i; }

    void Invalidate() { i = long(TIndex::BASE)-1; }
    bool IsValid() const { return i+1 != TIndex::BASE; }
    // operator bool() const { return IsValid(); }

    // archives store 1-based numbers (the historic base-1 raw value), independent of BASE
    void DoArchive (Archive & ar)
    {
      T nr1 = i - T(BASE) + 1;
      ar & nr1;
      if (ar.Input()) i = nr1 - 1 + T(BASE);
    }
  };


  template <typename T, typename TIndex, int Base>  
  constexpr auto operator+ (Index<T,TIndex,Base> ind, int i) { Index<T,TIndex,Base> res(ind); return res += i; }
  template <typename T, typename TIndex, int Base>
  constexpr auto operator+ (Index<T,TIndex,Base> ind, size_t i) { Index<T,TIndex,Base> res(ind); return res += i; }
  template <typename T, typename TIndex, int Base>    
  constexpr TIndex operator+ (int i, Index<T,TIndex,Base> ind) { return ind+i; } // Indexx<T,TIndex,Base> res(ind); return res += i; 
  template <typename T, typename TIndex, int Base>    
  inline TIndex operator+ (size_t i, Index<T,TIndex,Base> ind) { return ind+i; } //  TIndex res(ind); res += i; return res; }
  
  template <typename T, typename TIndex, int Base>    
  constexpr inline auto operator- (Index<T,TIndex,Base> ind, int i) { Index<T,TIndex,Base> res(ind); return res -= i; }  
  
  template <typename T, typename TIndex, int Base>      
  inline bool operator< (Index<T,TIndex,Base> a, Index<T,TIndex,Base> b) { return a-b < 0; }
  template <typename T, typename TIndex, int Base>      
  inline bool operator> (Index<T,TIndex,Base> a, Index<T,TIndex,Base> b) { return a-b > 0; }
  template <typename T, typename TIndex, int Base>      
  inline bool operator>= (Index<T,TIndex,Base> a, Index<T,TIndex,Base> b) { return a-b >= 0; }
  template <typename T, typename TIndex, int Base>      
  inline bool operator<= (Index<T,TIndex,Base> a, Index<T,TIndex,Base> b) { return a-b <= 0; }

  template <typename T, typename TIndex, int Base>      
  constexpr bool operator== (Index<T,TIndex,Base> a, Index<T,TIndex,Base> b) { return a.i == b.i; }
  template <typename T, typename TIndex, int Base>      
  constexpr bool operator!= (Index<T,TIndex,Base> a, Index<T,TIndex,Base> b) { return a.i != b.i; }


  template <typename T, typename TIndex, int Base>      
  inline void SetInvalid (Index<T,TIndex,Base> & id) { id.Invalidate(); }
  template <typename T, typename TIndex, int Base>        
  inline bool IsInvalid (const Index<T,TIndex,Base> & id) { return !id.IsValid(); }
  template <typename T, typename TIndex, int Base>
  inline size_t HashValue (Index<T,TIndex,Base> id, size_t size) { return (113*size_t(id)) % size; }


  
  class PointIndex : public Index<int,PointIndex,NETGEN_POINTINDEX_BASE>
  {
    friend class Index<int,PointIndex,NETGEN_POINTINDEX_BASE>;
    constexpr PointIndex (int ai) : Index(ai) { }   // use IndexBASE<PointIndex>()+nr
  public:
    using Index::Index;
    operator int () const = delete;    // a PointIndex stays a PointIndex
    operator int & () = delete;
  };

  // input-output is 1-based
  inline istream & operator>> (istream & ist, PointIndex & pi)
  {
    int i; ist >> i;
    pi = PointIndex::FromNr1(i);
    return ist;
  }

  inline ostream & operator<< (ostream & ost, const PointIndex & pi)
  {
    // return (ost << int(pi));
    int intpi = pi.Nr1();
    return (ost << intpi);    
  }



  /*
    PointIndices<N> is just IVec<N,PointIndex> - use operator[], not I1()/I2()/...
    To get a sorted one, write PointIndices<N>(...).Sort() or use SortedPointIndices<N>.
   */
  
  template <int N>
  using PointIndices = IVec<N,PointIndex>;

  template <int N>
  class SortedPointIndices : public PointIndices<N>
  {
    using PointIndices<N>::Sort;
    struct ordered_t { };
    constexpr SortedPointIndices (ordered_t, PointIndices<N> pnts)
      : PointIndices<N>(pnts) { }
  public:
    /// for values that are known to be ordered already (skips the sort)
    static constexpr SortedPointIndices Ordered (PointIndices<N> pnts)
    { return SortedPointIndices(ordered_t{}, pnts); }
    constexpr SortedPointIndices (PointIndices<N> pnts)
      : PointIndices<N>(pnts.Sort()) { } 
    
      template <typename ...Pnts>
    constexpr SortedPointIndices (Pnts ...pnts)
      : PointIndices<N>(pnts...)
    { Sort(); }
  };
  
}



namespace ngcore
{

  template <>
  struct CHT_trait<netgen::PointIndex>
  {
    constexpr static inline netgen::PointIndex Invalid() { return netgen::PointIndex::INVALID; }
    constexpr static inline size_t HashValue (const netgen::PointIndex & hash, size_t mask)
    { return (hash-IndexBASE<netgen::PointIndex>()) & mask; }
  };


  template <>
  struct CHT_trait<netgen::PointIndices<2>>
  {
    constexpr static inline netgen::PointIndices<2> Invalid() { return { netgen::PointIndex::INVALID, netgen::PointIndex::INVALID} ; }
    constexpr static inline size_t HashValue (const netgen::PointIndices<2> & hash, size_t mask)
    { return HashValue2(IVec<2>(hash[0]-IndexBASE<netgen::PointIndex>(),
                                hash[1]-IndexBASE<netgen::PointIndex>()), mask); }
  };
  

  template <>
  struct CHT_trait<netgen::SortedPointIndices<2>>
  {
    constexpr static inline netgen::SortedPointIndices<2> Invalid()
    { return netgen::SortedPointIndices<2>::Ordered ({ netgen::PointIndex::INVALID, netgen::PointIndex::INVALID }); }
    constexpr static inline size_t HashValue (const netgen::SortedPointIndices<2> & hash, size_t mask)
    // { return HashValue2(IVec<2,netgen::int>(hash[0], hash[1]), mask); }
    { return CHT_trait<netgen::PointIndices<2>>::HashValue (hash, mask); }
  };
  

  template <>
  struct CHT_trait<netgen::PointIndices<3>>
  {
    constexpr static inline netgen::PointIndices<3> Invalid()
    { return { netgen::PointIndex::INVALID, netgen::PointIndex::INVALID, netgen::PointIndex::INVALID }; }
    constexpr static inline size_t HashValue (const netgen::PointIndices<3> & hash, size_t mask)
    { return HashValue2(IVec<3>(hash[0]-IndexBASE<netgen::PointIndex>(),
                                hash[1]-IndexBASE<netgen::PointIndex>(),
                                hash[2]-IndexBASE<netgen::PointIndex>()), mask); }
  };


  template <>
  struct CHT_trait<netgen::SortedPointIndices<3>>
  {
    constexpr static inline netgen::SortedPointIndices<3> Invalid()
    { return netgen::SortedPointIndices<3>::Ordered ({ netgen::PointIndex::INVALID, netgen::PointIndex::INVALID, netgen::PointIndex::INVALID }); }
    constexpr static inline size_t HashValue (const netgen::SortedPointIndices<3> & hash, size_t mask)
    { return CHT_trait<netgen::PointIndices<3>>::HashValue (hash, mask); }
  };


  template <>
  constexpr inline netgen::PointIndices<3> InvalidHash<netgen::PointIndices<3>> ()
  { return netgen::PointIndices<3>{netgen::PointIndex::INVALID, netgen::PointIndex::INVALID, netgen::PointIndex::INVALID}; }

  /*
  template <>
  constexpr inline netgen::SortedPointIndices<2> InvalidHash<netgen::SortedPointIndices<2>> ()
  //   { return InvalidHash<netgen::PointIndices<2>>(); }
  { return CHT_trait<netgen::PointIndices<2>>::Invalid(); }    
  */
}


namespace std
{
  // structured binding support
  template <auto N>
  struct tuple_size<netgen::PointIndices<N>> : std::integral_constant<std::size_t, N> {};
  template<size_t N, auto M> struct tuple_element<N,netgen::PointIndices<M>> { using type = netgen::PointIndex; };
  template <auto N>
  struct tuple_size<netgen::SortedPointIndices<N>> : std::integral_constant<std::size_t, N> {};
  template<size_t N, auto M> struct tuple_element<N,netgen::SortedPointIndices<M>> { using type = netgen::PointIndex; };
}

namespace netgen
{

  class AnyElementIndex;

  /**
     Element number of a D-dimensional mesh element:
     ElIndex<3> is a volume element, ElIndex<2> a surface element, ElIndex<1> a segment.
  */
  template <int D>
  class ElIndex : public Index<int,ElIndex<D>,0>
  {
    typedef Index<int,ElIndex<D>,0> TBase;
    friend class Index<int,ElIndex<D>,0>;
    constexpr ElIndex (int ai) : TBase(ai) { }   // use IndexBASE<ElIndex<D>>()+nr, or FromNr0/FromNr1
  public:
    using TBase::TBase;
    operator int () const = delete;    // an ElIndex stays an ElIndex
    operator int & () = delete;
    /// narrowing from AnyElementIndex is explicit - name the kind you mean
    explicit constexpr ElIndex (AnyElementIndex bi);
  };

  using ElementIndex = ElIndex<3>;
  using SurfaceElementIndex = ElIndex<2>;
  using SegmentIndex = ElIndex<1>;
}

namespace ngcore
{
  template <int D>
  struct CHT_trait<netgen::ElIndex<D>>
  {
    constexpr static inline netgen::ElIndex<D> Invalid() { return netgen::ElIndex<D>::INVALID; }
    constexpr static inline size_t HashValue (const netgen::ElIndex<D> & hash, size_t mask)
    { return (hash-IndexBASE<netgen::ElIndex<D>>()) & mask; }
  };
}

namespace netgen
{

  template <int D>
  inline istream & operator>> (istream & ist, ElIndex<D> & ei)
  {
    int i; ist >> i; ei = ElIndex<D>::FromNr0(i); return ist;
  }

  template <int D>
  inline ostream & operator<< (ostream & ost, const ElIndex<D> & ei)
  {
    return ost << ei.Nr0();
  }


  /**
     An element number whose kind (volume element, surface element or segment)
     is fixed by the context, not by the value - e.g. HPRefElement::coarse_elnr.
     Widening from a concrete index is implicit, narrowing back is explicit.
  */
  class AnyElementIndex : public Index<int,AnyElementIndex,0>
  {
  public:
    using Index::Index;
    template <int D>
    constexpr AnyElementIndex (ElIndex<D> ei) : Index(ei.Nr0()) { }
  };

  template <int D>
  constexpr ElIndex<D>::ElIndex (AnyElementIndex bi)
    : Index<int,ElIndex<D>,0>(bi.Nr0()) { }


  /**
     Point number in the advancing front (AdFront3), not a mesh point number.
     AdFront3::GetGlobalIndex maps it to the PointIndex of the mesh.
  */
  class Front3PointIndex : public Index<int,Front3PointIndex,0>
  {
    friend class Index<int,Front3PointIndex,0>;
    constexpr Front3PointIndex (int ai) : Index(ai) { }   // use IndexBASE<Front3PointIndex>()+nr
  public:
    using Index::Index;
    operator int () const = delete;
    operator int & () = delete;
  };

  /**
     Point number in the 2D advancing front (AdFront2), not a mesh point number.
     AdFront2::GetGlobalIndex maps it to the PointIndex of the mesh.
  */
  class Front2PointIndex : public Index<int,Front2PointIndex,0>
  {
    friend class Index<int,Front2PointIndex,0>;
    constexpr Front2PointIndex (int ai) : Index(ai) { }   // use IndexBASE<Front2PointIndex>()+nr
  public:
    using Index::Index;
    operator int () const = delete;
    operator int & () = delete;
  };

  /**
     Point number in the local numbering of AdFront3::GetLocals,
     i.e. an index into the local point array handed to the meshing rules.
  */
  class LocalPointIndex : public Index<int,LocalPointIndex,0>
  {
    friend class Index<int,LocalPointIndex,0>;
    constexpr LocalPointIndex (int ai) : Index(ai) { }   // use IndexBASE<LocalPointIndex>()+nr
  public:
    using Index::Index;
    operator int () const = delete;
    operator int & () = delete;
  };

  /**
     Point number within a meshing rule (as parsed from the .rls files),
     mapped to a LocalPointIndex by the rule application.
  */
  class RulePointIndex : public Index<int,RulePointIndex,1>
  {
    friend class Index<int,RulePointIndex,1>;
    constexpr RulePointIndex (int ai) : Index(ai) { }   // use IndexBASE<RulePointIndex>()+nr
  public:
    using Index::Index;
    operator int () const = delete;
    operator int & () = delete;
  };

  /**
     Vertex number within a single element, 1 .. GetNP().
  */
  class ElementVertexIndex : public Index<int,ElementVertexIndex,1>
  {
    friend class Index<int,ElementVertexIndex,1>;
    constexpr ElementVertexIndex (int ai) : Index(ai) { }   // use IndexBASE<ElementVertexIndex>()+nr
  public:
    using Index::Index;
    operator int () const = delete;
    operator int & () = delete;
  };

  /**
     1-based index into the regions of entity dimension D, Mesh::Regions<D>().
     INVALID (= 0) means unset. Only constructed via FromNr0/FromNr1.
  */
  class AnyRegionIndex;

  template <int D>
  class RegionIndex : public Index<int,RegionIndex<D>,1>
  {
    typedef Index<int,RegionIndex<D>,1> TBase;
    friend class Index<int,RegionIndex<D>,1>;
    constexpr RegionIndex (int ai) : TBase(ai) { }
  public:
    using TBase::TBase;
    /// narrowing from AnyRegionIndex - the dimension is fixed by the context
    constexpr RegionIndex (AnyRegionIndex ai);
    operator int () const = delete;    // use Nr1() / Nr0() / IsValid()
    operator int & () = delete;
  };

  using VertexRegionIndex = RegionIndex<0>;
  using EdgeRegionIndex = RegionIndex<1>;
  using FaceRegionIndex = RegionIndex<2>;
  using VolumeRegionIndex = RegionIndex<3>;

  // old names
  using FaceDescriptorIndex = FaceRegionIndex;
  using EdgeDescriptorIndex = EdgeRegionIndex;

  /**
     A region index whose dimension (vertex, edge, face or volume region) is
     fixed by the context, not by the value - e.g. HPRefElement::index.
     Converts implicitly from and to the four RegionIndex<D>.
  */
  class AnyRegionIndex : public Index<int,AnyRegionIndex,1>
  {
  public:
    using Index::Index;
    template <int D>
    constexpr AnyRegionIndex (RegionIndex<D> ri) : Index(ri.Nr1()) { }
    operator int () const = delete;
    operator int & () = delete;
  };

  template <int D>
  constexpr RegionIndex<D>::RegionIndex (AnyRegionIndex ai)
    : TBase(ai.Nr1()) { }

  template <int D>
  inline ostream & operator<< (ostream & ost, const RegionIndex<D> & i) { return ost << i.Nr1(); }
  template <int D>
  inline istream & operator>> (istream & ist, RegionIndex<D> & i)
  {
    int nr; ist >> nr; i = RegionIndex<D>::FromNr1(nr); return ist;
  }

  inline ostream & operator<< (ostream & ost, const Front2PointIndex & fpi)
  {
    return ost << (fpi - IndexBASE<Front2PointIndex>());
  }

  inline ostream & operator<< (ostream & ost, const Front3PointIndex & fpi)
  {
    return ost << (fpi - IndexBASE<Front3PointIndex>());
  }

  inline ostream & operator<< (ostream & ost, const LocalPointIndex & lpi)
  {
    return ost << (lpi - IndexBASE<LocalPointIndex>());
  }

  // rule files number their points 1-based
  inline constexpr RulePointIndex RuleP (int nr) { return RulePointIndex::FromNr1(nr); }

  inline istream & operator>> (istream & ist, RulePointIndex & rpi)
  {
    int i; ist >> i;
    rpi = RulePointIndex::FromNr1(i);
    return ist;
  }

  inline ostream & operator<< (ostream & ost, const RulePointIndex & rpi)
  {
    return ost << (rpi - IndexBASE<RulePointIndex>());
  }

  inline ostream & operator<< (ostream & ost, const ElementVertexIndex & evi)
  {
    return ost << (evi - IndexBASE<ElementVertexIndex>());
  }


template <typename TINDEX>
class MiniElement2dT
{
protected:
  int np;
  TINDEX pnum[4];
  bool deleted;
public:
  MiniElement2dT ()
  { np = 3; deleted = 0; }
  MiniElement2dT (int anp)
  { np = anp; deleted = 0; }

  int GetNP() const { return np; }
  void SetNP (int anp) { np = anp; }
  TINDEX & operator[] (int i) { return pnum[i]; }
  const TINDEX operator[] (int i) const { return pnum[i]; }

  const TINDEX PNum (int i) const { return pnum[i-1]; }
  TINDEX & PNum (int i) { return pnum[i-1]; }
  const TINDEX PNumMod (int i) const { return pnum[(i-1)%np]; }
  auto PNums() { return FlatArray<TINDEX> (np, &pnum[0]); }
  auto PNums() const { return FlatArray<const TINDEX> (np, &pnum[0]); }
  void Delete () { deleted = true; for (TINDEX & p : pnum) p.Invalidate(); }
  bool IsDeleted () const { return deleted; }
};

/// face in local (GetLocals) numbering
using MiniElement2d = MiniElement2dT<LocalPointIndex>;
/// face in advancing-front numbering
using FrontElement2d = MiniElement2dT<Front3PointIndex>;
/// 2d element / face in rule numbering
using RuleElement2d = MiniElement2dT<RulePointIndex>;
/// face of a volume element, in element-vertex numbering
using ElementFace = MiniElement2dT<ElementVertexIndex>;


/// volume element in a non-mesh numbering (local or rule)
template <typename TINDEX>
class MiniElementT
{
  ELEMENT_TYPE typ;
  int np;
  TINDEX pnum[8];
public:
  MiniElementT () : typ(TET), np(4) { }
  MiniElementT (int anp) { SetNP(anp); }

  ELEMENT_TYPE GetType () const { return typ; }
  int GetNP () const { return np; }

  void SetNP (int anp)
  {
    np = anp;
    switch (np)
      {
      case 4: typ = TET; break;
      case 5: typ = PYRAMID; break;
      case 6: typ = PRISM; break;
      case 7: typ = HEX7; break;
      case 8: typ = HEX; break;
      default: cerr << "MiniElementT::SetNP unknown element with " << np << " points" << endl;
      }
  }

  void SetType (ELEMENT_TYPE atyp)
  {
    typ = atyp;
    switch (typ)
      {
      case TET: np = 4; break;
      case PYRAMID: np = 5; break;
      case PRISM: np = 6; break;
      case HEX7: np = 7; break;
      case HEX: np = 8; break;
      default: cerr << "MiniElementT::SetType unknown type " << int(typ) << endl;
      }
  }

  TINDEX & operator[] (int i) { return pnum[i]; }
  TINDEX operator[] (int i) const { return pnum[i]; }
  TINDEX & PNum (int i) { return pnum[i-1]; }
  TINDEX PNum (int i) const { return pnum[i-1]; }
  auto PNums() { return FlatArray<TINDEX> (np, &pnum[0]); }
  auto PNums() const { return FlatArray<const TINDEX> (np, &pnum[0]); }
};

/// volume element in local (GetLocals) numbering, as produced by the meshing rules
using LocalElement = MiniElementT<LocalPointIndex>;
/// volume element in rule numbering
using RuleElement = MiniElementT<RulePointIndex>;
/// sub-tet of a volume element, in element-vertex numbering
using ElementTet = MiniElementT<ElementVertexIndex>;

template <typename TINDEX>
inline ostream & operator<< (ostream & ost, const MiniElementT<TINDEX> & el)
{
  ost << "np = " << el.GetNP();
  for (int j = 0; j < el.GetNP(); j++)
    ost << " " << el[j];
  return ost;
}


template <typename TINDEX>
inline ostream & operator<<(ostream  & s, const MiniElement2dT<TINDEX> & el)
{
  s << "np = " << el.GetNP();
  for (int j = 0; j < el.GetNP(); j++)
    s << " " << el[j];
  return s;
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
    MeshPoint () = default;

    MeshPoint (const Point<3> & ap, int alayer = 1, POINTTYPE apt = INNERPOINT)
      : Point<3> (ap), singular(0.), layer(alayer), type(apt) { }
  
    void SetPoint (const Point<3> & ap)
    { 
      Point<3>::operator= (ap); 
      layer = 0; 
      singular = 0; 
    }

    void Scale(double factor) { x[0] *= factor; x[1] *= factor; x[2] *= factor; }

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
      ar.DoPacked (x[0], x[1], x[2], layer, singular, (unsigned char&)(type));
    }
  };

  inline ostream & operator<<(ostream  & s, const MeshPoint & pt)
  { 
    return (s << Point<3> (pt)); 
  }



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

    /// face descriptor index (1-based)
    FaceRegionIndex index = FaceRegionIndex::INVALID;
    ///
    ELEMENT_TYPE typ;   // number of points and vertices follow from the type: element2d_info
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

    /// a linked list for all elements in the same face
    SurfaceElementIndex next;

  public:
    static auto GetDataLayout()
    {
      return std::map<string, int>({
          { "pnum", offsetof(Element2d, pnum)},
          { "index", offsetof(Element2d, index) },
          { "type", offsetof(Element2d, typ) },
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
    DLL_HEADER Element2d (PointIndex pi1, PointIndex pi2, PointIndex pi3);
    ///
    DLL_HEADER Element2d (PointIndex pi1, PointIndex pi2, PointIndex pi3, PointIndex pi4);
    ///
    ELEMENT_TYPE GetType () const { return typ; }
    /// 
    void SetType (ELEMENT_TYPE atyp)
    {
      typ = atyp;
      if (element2d_info::np[atyp] == 0)
        PrintSysError ("Element2d::SetType, illegal type ", int(atyp));
      is_curved = (GetNP() >= 4); 
    }
    ///
    int GetNP() const { return element2d_info::np[typ]; }
    ///
    int GetNV() const { return element2d_info::nv[typ]; }

    ///
    PointIndex & operator[] (int i) { return pnum[i]; }
    ///
    const PointIndex & operator[] (int i) const { return pnum[i]; }

    auto PNums () const { return FlatArray<const PointIndex> (GetNP(), &pnum[0]); }
    auto PNums ()  { return FlatArray<PointIndex> (GetNP(), &pnum[0]); }
    template <int NP>
    auto PNums() const { return FlatArray<const PointIndex> (NP, &pnum[0]); }
    auto Vertices() const { return FlatArray<const PointIndex> (GetNV(), &pnum[0]); }

    auto GeomInfo() const { return FlatArray<const PointGeomInfo> (GetNP(), &geominfo[0]); }
    auto GeomInfo() { return FlatArray<PointGeomInfo> (GetNP(), &geominfo[0]); }
    
    ///
    PointIndex & PNum (int i) { return pnum[i-1]; }
    /// vertex i of this element
    PointIndex & PNum (ElementVertexIndex i) { return pnum[i-IndexBASE<ElementVertexIndex>()]; }
    ///
    const PointIndex & PNum (int i) const { return pnum[i-1]; }
    const PointIndex & PNum (ElementVertexIndex i) const { return pnum[i-IndexBASE<ElementVertexIndex>()]; }
    ///
    PointIndex & PNumMod (int i) { return pnum[(i-1) % GetNP()]; }
    ///
    const PointIndex & PNumMod (int i) const { return pnum[(i-1) % GetNP()]; }
    ///

    ///
    PointGeomInfo & GeomInfoPi (int i) { return geominfo[i-1]; }
    ///
    const PointGeomInfo & GeomInfoPi (int i) const { return geominfo[i-1]; }
    ///
    PointGeomInfo & GeomInfoPiMod (int i) { return geominfo[(i-1) % GetNP()]; }
    ///
    const PointGeomInfo & GeomInfoPiMod (int i) const { return geominfo[(i-1) % GetNP()]; }

    auto & NewestVertex() { return newest_vertex; }
    auto NewestVertex() const { return newest_vertex; }

    void DoArchive (Archive & ar)
    {
      short _np, _typ;
      bool _curved, _vis, _deleted;
      if (ar.Output())
        { _np = GetNP(); _typ = typ; _curved = is_curved;
          _vis = visible; _deleted = deleted; }
      // ar & _np & _typ & index & _curved & _vis & _deleted;
      ar.DoPacked (_np, _typ, index, _curved, _vis, _deleted);
      // ar & next; don't need 
      if (ar.Input())
        { typ = ELEMENT_TYPE(_typ); is_curved = _curved;
          visible = _vis; deleted = _deleted; }
      /*
      for (size_t i = 0; i < GetNP(); i++)
        ar & pnum[i];
      */
      // archive stores 1-based point numbers, independent of BASE
      int nr1[ELEMENT2D_MAXPOINTS];
      if (ar.Output())
        for (int k = 0; k < GetNP(); k++) nr1[k] = pnum[k].Nr1();
      ar.Do (nr1, GetNP());
      if (ar.Input())
        for (int k = 0; k < GetNP(); k++) pnum[k] = PointIndex::FromNr1(nr1[k]);
    }

#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif
    

    void SetIndex (FaceRegionIndex si) { index = si; }
    ///
    FaceRegionIndex GetIndex () const { return index; }





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

    void GetTransformation (int ip, FlatArray<Point<2>, PointIndex> points,
                            class DenseMatrix & trans) const;
    void GetTransformation (int ip, class DenseMatrix & pmat,
                            class DenseMatrix & trans) const;

    void GetShape (const Point<2> & p, class Vector & shape) const;
    DLL_HEADER void GetShapeNew (const Point<2> & p, class FlatVector & shape) const;
    template <typename T>
    DLL_HEADER void GetShapeNew (const Point<2,T> & p, TFlatVector<T> shape) const;
    /// matrix 2 * GetNP()
    DLL_HEADER void GetDShape (const Point<2> & p, class DenseMatrix & dshape) const;
    template <typename T>
    DLL_HEADER void GetDShapeNew (const Point<2,T> & p, class MatrixFixWidth<2,T> & dshape) const;
    
    /// matrix 2 * GetNP()
    void GetPointMatrix (FlatArray<Point<2>, PointIndex> points,
                         class DenseMatrix & pmat) const;

    void ComputeIntegrationPointData () const;
  

    double CalcJacobianBadness (FlatArray<Point<2>, PointIndex> points) const;
    double CalcJacobianBadness (const T_POINTS & points, 
                                const Vec<3> & n) const;
    double CalcJacobianBadnessDirDeriv (FlatArray<Point<2>, PointIndex> points,
                                        int pi, Vec<2> & dir, double & dd) const;


    
    void Delete ()
    {
      deleted = true;
      // for (PointIndex & p : pnum) p.Invalidate(); 
    }
    
    bool IsDeleted () const 
    {
#ifdef DEBUG
      if ((pnum[0]-IndexBASE<PointIndex>() < 0) && !deleted)
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

  DLL_HEADER ostream & operator<<(ostream  & s, const Element2d & el);
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
  struct ElementHeader
  {
    ELEMENT_TYPE typ;   // number of points, vertices, ... follow from the type: element3d_info
    int8_t newest_vertex = -1; // from refinement via bisection
    /// sub-domain index
    VolumeRegionIndex index;
    bool is_curved;   // element is (high order) curved
    bool refflag;     // mark element for refinement

    class flagstruct {
    public:
      bool marked:1;  // marked for refinement
      bool badel:1;   // angles worse then limit
      bool reverse:1; // for refinement a la Bey
      bool illegal:1; // illegal, will be split or swapped
      bool illegal_valid:1; // is illegal-flag valid ?
      bool strongrefflag:1;
      bool deleted:1;   // element is deleted, will be removed from array
      bool fixed:1;     // don't change element in optimization
    };

    flagstruct flags;
  };

  class Element;

  /*
    Handle to a volume element: header, point numbers and the number of point slots.
    Refers to a slot of the mesh's element array or to the storage of an Element value.
    Copying the handle rebinds it, assigning through it copies the element contents.
    Read access is const, modification is not.
  */
  class ElementRef
  {
  protected:
    ElementHeader * h;
    PointIndex * pn;
    int maxnp;

  public:
    typedef ElementHeader::flagstruct flagstruct;

    ElementRef (ElementHeader * ah, PointIndex * apn, int amaxnp) : h(ah), pn(apn), maxnp(amaxnp) { }
    ElementRef (DynStrideView<ElementHeader, false, PointIndex> v)
      : h(&v.Head()), pn(v.TailPtr<0>()), maxnp(int(v.Width())) { }
    ElementRef (DynStrideView<ElementHeader, true, PointIndex> v)
      : h(const_cast<ElementHeader*>(&v.Head())), pn(const_cast<PointIndex*>(v.TailPtr<0>())), maxnp(int(v.Width())) { }
    ElementRef (const ElementRef &) = default;

    /// copies header and point numbers, throws if the element does not fit
    DLL_HEADER ElementRef & operator= (const ElementRef & el2);

    ElementHeader & Header () const { return *h; }
    int MaxNP () const { return maxnp; }

    const flagstruct& Flags() const { return h->flags; }
    flagstruct& Flags() { return h->flags; }

    DLL_HEADER void SetNP (int anp);
    DLL_HEADER void SetType (ELEMENT_TYPE atyp);
    int GetNP () const { return element3d_info::np[h->typ]; }
    // old style:
    int NP () const { return element3d_info::np[h->typ]; }

    uint8_t GetNV() const { return element3d_info::nv[h->typ]; }

    ELEMENT_TYPE GetType () const { return h->typ; }

    PointIndex & operator[] (int i) { NETGEN_CHECK_RANGE(i, 0, maxnp); return pn[i]; }
    const PointIndex & operator[] (int i) const { NETGEN_CHECK_RANGE(i, 0, maxnp); return pn[i]; }

    auto PNums () const { return FlatArray<const PointIndex> (GetNP(), pn); }
    auto PNums () { return FlatArray<PointIndex> (GetNP(), pn); }
    template <int NP>
    auto PNums() const { return FlatArray<const PointIndex> (NP, pn); }
    FlatArray<const PointIndex> Vertices() const { return { GetNV(), pn }; }

    PointIndex & PNum (int i) { NETGEN_CHECK_RANGE(i, 1, maxnp+1); return pn[i-1]; }
    /// vertex i of this element
    PointIndex & PNum (ElementVertexIndex i) { return pn[i-IndexBASE<ElementVertexIndex>()]; }
    const PointIndex & PNum (int i) const { NETGEN_CHECK_RANGE(i, 1, maxnp+1); return pn[i-1]; }
    const PointIndex & PNum (ElementVertexIndex i) const { return pn[i-IndexBASE<ElementVertexIndex>()]; }
    PointIndex & PNumMod (int i) { return pn[(i-1) % GetNP()]; }
    const PointIndex & PNumMod (int i) const { return pn[(i-1) % GetNP()]; }

    auto & NewestVertex() { return h->newest_vertex; }
    auto NewestVertex() const { return h->newest_vertex; }

    DLL_HEADER void DoArchive (Archive & ar);

    void SetIndex (VolumeRegionIndex si) { h->index = si; }
    VolumeRegionIndex GetIndex () const { return h->index; }


    DLL_HEADER void GetBox (const T_POINTS & points, Box3d & box) const;
    /// Calculates Volume of element
    DLL_HEADER double Volume (const T_POINTS & points) const;
    DLL_HEADER void Print (ostream & ost) const;
    int GetNFaces () const { return element3d_info::nfaces[h->typ]; }
    inline void GetFace (int i, Element2d & face) const;
    DLL_HEADER void GetFace2 (int i, Element2d & face) const;
    DLL_HEADER void Invert ();


    /// split into 4 node tets
    DLL_HEADER void GetTets (Array<Element> & locels) const;
    /// split into 4 node tets, local point nrs
    DLL_HEADER void GetTetsLocal (Array<ElementTet> & locels) const;
    /// returns coordinates of nodes
    DLL_HEADER void GetNodesLocalNew (Array<Point<3> > & points) const;
    /// split surface into 3 node trigs
    DLL_HEADER void GetSurfaceTriangles (Array<ElementFace> & surftrigs) const;

    /// get number of 'integration points'
    DLL_HEADER int GetNIP () const;
    DLL_HEADER void GetIntegrationPoint (int ip, Point<3> & p, double & weight) const;

    DLL_HEADER void GetTransformation (int ip, const T_POINTS & points,
                                       class DenseMatrix & trans) const;
    DLL_HEADER void GetTransformation (int ip, class DenseMatrix & pmat,
                                       class DenseMatrix & trans) const;

    DLL_HEADER void GetShape (const Point<3> & p, class Vector & shape) const;
    template <typename T>
    DLL_HEADER void GetShapeNew (const Point<3,T> & p, TFlatVector<T> shape) const;
    /// matrix 2 * np
    DLL_HEADER void GetDShape (const Point<3> & p, class DenseMatrix & dshape) const;
    template <typename T>
    DLL_HEADER void GetDShapeNew (const Point<3,T> & p, class MatrixFixWidth<3,T> & dshape) const;
    /// matrix 3 * np
    DLL_HEADER void GetPointMatrix (const T_POINTS & points,
                                    class DenseMatrix & pmat) const;

    DLL_HEADER void ComputeIntegrationPointData () const;

    DLL_HEADER double CalcJacobianBadness (const T_POINTS & points) const;
    DLL_HEADER double CalcJacobianBadnessDirDeriv (const T_POINTS & points,
                                                   int pi, Vec<3> & dir, double & dd) const;
    DLL_HEADER double CalcJacobianBadnessGradient (const T_POINTS & points,
                                                   int pi, Vec<3> & grad) const;

    void SetRefinementFlag (bool rflag = 1) { h->refflag = rflag; }
    int TestRefinementFlag () const { return h->refflag; }

    void SetStrongRefinementFlag (bool rflag = 1) { h->flags.strongrefflag = rflag; }
    int TestStrongRefinementFlag () const { return h->flags.strongrefflag; }

    int Illegal () const
    {
      NETGEN_CHECK_SAME(h->flags.illegal_valid, true);
      return h->flags.illegal;
    }
    int IllegalValid () const { return h->flags.illegal_valid; }
    void SetIllegal (int aillegal)
    {
      h->flags.illegal = aillegal ? 1 : 0;
      h->flags.illegal_valid = 1;
    }
    void SetLegal (int alegal)
    {
      h->flags.illegal = alegal ? 0 : 1;
      h->flags.illegal_valid = 1;
    }

    void Touch() { h->flags.illegal_valid = 0; }

    void Delete () { h->flags.deleted = 1; }
    bool IsDeleted () const
    {
#ifdef DEBUG
      if (pn[0]-IndexBASE<PointIndex>() < 0 && !h->flags.deleted)
        cerr << "Volelement has illegal pnum, but not marked as deleted" << endl;
#endif
      return h->flags.deleted;
    }

    bool IsCurved () const { return h->is_curved; }
    void SetCurved (bool acurved) { h->is_curved = acurved; }

    DLL_HEADER bool operator== (const ElementRef & el2) const;
  };


  /// volume element value type: own storage for ELEMENT_MAXPOINTS points, handled through ElementRef
  class Element : public ElementRef
  {
    ElementHeader hstore;
    PointIndex pnstore[ELEMENT_MAXPOINTS];
  public:
    DLL_HEADER Element ();
    DLL_HEADER Element (int anp);
    DLL_HEADER Element (ELEMENT_TYPE type);
    /// copy of a stored element; explicit, a handle is what you usually want
    explicit Element (const ElementRef & el) : Element() { ElementRef::operator= (el); }
    Element (const Element & e2) : ElementRef(&hstore, pnstore, ELEMENT_MAXPOINTS), hstore(e2.hstore)
    { for (int i = 0; i < ELEMENT_MAXPOINTS; i++) pnstore[i] = e2.pnstore[i]; }
    Element & operator= (const Element & e2)
    { hstore = e2.hstore; for (int i = 0; i < ELEMENT_MAXPOINTS; i++) pnstore[i] = e2.pnstore[i]; return *this; }
    Element & operator= (const ElementRef & el) { ElementRef::operator= (el); return *this; }

    static auto GetDataLayout()
    {
      Element hel;
      auto off = [&hel] (const void * p) { return int((const char*)p - (const char*)&hel); };
      return std::map<string, int>({
          { "pnum", off(&hel.pnstore[0]) },
          { "index", off(&hel.hstore.index) },
          { "type", off(&hel.hstore.typ) },
          { "refine", off(&hel.hstore.refflag) },
          { "curved", off(&hel.hstore.is_curved) }
        });
    }

#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif
  };

  /// array of volume elements with run-time number of point slots
  typedef DynStrideArray<ElementHeader, TailList<PointIndex>, ElementIndex> T_VOLELEMENTS_BASE;

  /// adds ElementRef access and iteration to a (flat or owning) strided element array
  // the template parameter must not be called BASE: the strided array has a static member of that name
  template <class TBASE>
  class ElementRefArray : public TBASE
  {
  public:
    using TBASE::TBASE;
    ElementRefArray (const TBASE & b) : TBASE(b) { }

    ElementRef operator[] (typename TBASE::index_type i) { return ElementRef (TBASE::operator[] (i)); }
    const ElementRef operator[] (typename TBASE::index_type i) const { return ElementRef (TBASE::operator[] (i)); }
    ElementRef First () { return ElementRef (TBASE::First()); }
    const ElementRef First () const { return ElementRef (TBASE::First()); }
    ElementRef Last () { return ElementRef (TBASE::Last()); }
    const ElementRef Last () const { return ElementRef (TBASE::Last()); }

    template <class IT>
    class Iterator
    {
      IT it;
    public:
      Iterator (IT ait) : it(ait) { }
      Iterator & operator++ () { ++it; return *this; }
      ElementRef operator* () const { return ElementRef (*it); }
      bool operator!= (const Iterator & it2) const { return it != it2.it; }
      bool operator== (const Iterator & it2) const { return it == it2.it; }
    };
    auto begin () { return Iterator<decltype(TBASE::begin())> (TBASE::begin()); }
    auto end () { return Iterator<decltype(TBASE::end())> (TBASE::end()); }
    auto begin () const { return Iterator<decltype(TBASE::begin())> (TBASE::begin()); }
    auto end () const { return Iterator<decltype(TBASE::end())> (TBASE::end()); }

    auto Range () const { return TBASE::Range(); }
    template <typename... ARGS>
    auto Range (ARGS... args) const
    { return ElementRefArray<FlatDynStrideArray<ElementHeader, TailList<PointIndex>, size_t>> (TBASE::Range (args...)); }
  };

  class VolumeElementArray : public ElementRefArray<T_VOLELEMENTS_BASE>
  {
    typedef ElementRefArray<T_VOLELEMENTS_BASE> TBASE;
  public:
    using TBASE::TBASE;
    using TBASE::Append;

    // grows the width to the element's number of points if needed
    ElementIndex Append (const ElementRef & el)
    {
      if (size_t(el.GetNP()) > Width()) SetWidth (el.GetNP());
      ElementIndex ei = T_VOLELEMENTS_BASE::Append();
      ElementRef v = (*this)[ei];
      v = el;
      for (int k = el.GetNP(); k < int(Width()); k++) v[k].Invalidate();
      return ei;
    }

    void DoArchive (Archive & ar)
    {
      constexpr const char * width_version = "v6.2.2607-105";
      if (ar.Input() && ar.GetVersion("netgen") < width_version)
        {
          size_t s;
          ar & s;
          SetSize (0);
          for (size_t i = 0; i < s; i++)
            {
              Element el;
              el.DoArchive (ar);
              Append (el);
            }
          return;
        }
      ar.NeedsVersion ("netgen", width_version);
      size_t s = Size(), w = Width();
      ar & s & w;
      if (ar.Input()) { SetWidth (w); SetSize (s); }
      for (auto el : *this) el.DoArchive (ar);
    }
  };
  typedef VolumeElementArray T_VOLELEMENTS;

  /// explicit value copy of a stored element (auto x = mesh[i] yields a handle)
  inline Element Copy (const ElementRef & el) { return Element(el); }

  DLL_HEADER ostream & operator<<(ostream  & s, const ElementRef & el);






  /**
     Edge segment.

     How indices are used up to 2026-03-29

     edgenr:
     OCC: the geometry edge, 1-based
     Spline2D:  edge-nr, 1-based
     CSG:  edge counter 
     
     si:
     CSG: one surface of the edge
     Spline2d: bc-number
     OCC: edgenr, 1-based
     
     cd2i:
     ????
     
     epgeominfo:
     OCC: geometry  edgenr, 0-based
     Spline2D: edgenr, 1-based
     

     NGSoleve Interface:
     GetIndex
     mesh.dim == 3 ->  edgenr
     mesh.dim == 2 ->  si

     
     Python interface:
     edgenr -> segmnr
     index -> si
     
     Python ctor:
     index -> si, segmnr
     edgenr -> epgeominfo DECREMENT 1



     NEW from 2026-03-29
     GetEdgeNr()  -> the geometry edge
     GetIndex()  -> the index for boundary conditions
  */
  
  class Segment
  {
  PointIndex pnums[3];
  EdgePointGeomInfo epgeominfo[2]; // combines PointGeomInfo + dist
  /// 1-based edge descriptor index into mesh.Regions<1>() (INVALID = 0)
  EdgeRegionIndex index = EdgeRegionIndex::INVALID;

  public:
    ///
    DLL_HEADER Segment();
    Segment (const Segment& other) = default;
    Segment& operator=(const Segment & other) = default;

    PointIndex & operator[] (int i) { return pnums[i]; }
    const PointIndex & operator[] (int i) const { return pnums[i]; }

    int GetNP() const { return pnums[2].IsValid() ? 3 : 2; }
    auto PNums() const { return FlatArray<const PointIndex> (GetNP(), &pnums[0]); }
    auto PNums() { return FlatArray<PointIndex> (GetNP(), &pnums[0]); }
    auto Vertices() const { return FlatArray<const PointIndex> (2, &pnums[0]); }
    ELEMENT_TYPE GetType() const { return pnums[2].IsValid() ? SEGMENT3 : SEGMENT; }

    PointGeomInfo & GeomInfo (int i) { return epgeominfo[i].gi; }
    const PointGeomInfo & GeomInfo (int i) const { return epgeominfo[i].gi; }

    EdgePointGeomInfo & EPGeomInfo (int i) { return epgeominfo[i]; }
    const EdgePointGeomInfo & EPGeomInfo (int i) const { return epgeominfo[i]; }


    EdgeRegionIndex GetIndex() const { return index; }
    void SetIndex (EdgeRegionIndex i) { index = i; }

    void DoArchive (Archive & ar);
#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif

    static size_t OffsetPnums() { return offsetof(Segment, pnums); }
    static size_t OffsetIndex() { return offsetof(Segment, index); }
  };

  ostream & operator<<(ostream  & s, const Segment & seg);


  class Element0d
  {
  public:
    PointIndex pnum;
    string name;
    VertexRegionIndex index = VertexRegionIndex::INVALID;
    Element0d () = default;
    Element0d (PointIndex _pnum, VertexRegionIndex _index)
      : pnum(_pnum), index(_index) { ; }

    VertexRegionIndex GetIndex () const { return index; }
    void SetIndex (VertexRegionIndex i) { index = i; }

#ifdef PARALLEL
    static NG_MPI_Datatype MyGetMPIType();
#endif
    
    void DoArchive (Archive & ar);
  };

  ostream & operator<<(ostream  & s, const Element0d & el);

  /// common part of the regions of all dimensions
  class RegionBase
  {
    optional<string> name;   // nullopt: not set, reported as "default"
  public:
    DLL_HEADER static const string default_name;

    const string & GetName () const { return name ? *name : default_name; }
    bool HasName () const { return name.has_value(); }
    void SetName (optional<string> aname) { name = std::move(aname); }
    void ResetName () { name = nullopt; }
    const optional<string> & OptName () const { return name; }
  };

  /**
     Geometric entity of dimension D the mesh elements of dimension D belong to:
     Region<3> volume (material), Region<2> face, Region<1> edge, Region<0> vertex.
     Indexed by RegionIndex<D>.
  */
  template <int D> class Region;

  template <>
  class Region<3> : public RegionBase
  {
  public:
    Region () = default;
    explicit Region (optional<string> aname) { SetName(std::move(aname)); }
  };

  template <>
  class Region<0> : public RegionBase
  {
  public:
    Region () = default;
    explicit Region (optional<string> aname) { SetName(std::move(aname)); }
  };

  ///
  template <>
  class Region<2> : public RegionBase
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
    
    /// root of linked list 
    SurfaceElementIndex firstelement;
  
    double domin_singular;
    double domout_singular;

  public:
    DLL_HEADER Region();
    DLL_HEADER Region(int surfnri, int domini, int domouti, int tlosurfi);
    DLL_HEADER Region(const Region& other);
    Region & operator= (const Region & other) = default;

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
    const string & GetBCName () const { return GetName(); }
    void SetSurfNr (int sn) { surfnr = sn; }
    void SetDomainIn (int di) { domin = di; }
    void SetDomainOut (int dom) { domout = dom; }
    void SetBCProperty (int bc) { bcprop = bc; }
    void SetBCName (const string & bcn) { SetName(bcn); }
    // Philippose - 06/07/2009
    // Set the surface colour
    void SetSurfColour (Vec<4> colour) { surfcolour = colour; }

    void SetDomainInSingular (double v) { domin_singular = v; }
    void SetDomainOutSingular (double v) { domout_singular = v; }

    SurfaceElementIndex FirstElement() { return firstelement; }
    friend class Mesh;

    void DoArchive (Archive & ar);
  };

  using FaceRegion = Region<2>;

  ostream & operator<< (ostream  & s, const FaceRegion & fd);

  
 
  template <>
  class Region<1> : public RegionBase
  {
    int edgenr = -1;
    int surfnr[2] = {-1, -1};
    double singedge_left = 0;
    double singedge_right = 0;

    // legacy, kept for CSG compatibility
    int tlosurf = -1;
    int domin = -1, domout = -1;

    /// transient index: face descriptor index (1-based) in 3D, not serialized - recomputed by RebuildFDIndices()
    FaceRegionIndex index_ = FaceRegionIndex::INVALID;

  public:
    Region () = default;
    Region (int edgenri, int surfnr1 = -1, int surfnr2 = -1,
                    int domini = -1, int domouti = -1, int tlosurfi = -1)
      : edgenr(edgenri), surfnr{surfnr1, surfnr2}, tlosurf(tlosurfi), domin(domini), domout(domouti)
    { ; }

    int EdgeNr () const { return edgenr; }
    void SetEdgeNr (int nr) { edgenr = nr; }

    int SurfNr (int i) const { return surfnr[i]; }
    void SetSurfNr (int i, int nr) { surfnr[i] = nr; }

    double SingEdgeLeft () const { return singedge_left; }
    void SetSingEdgeLeft (double s) { singedge_left = s; }

    double SingEdgeRight () const { return singedge_right; }
    void SetSingEdgeRight (double s) { singedge_right = s; }

    // legacy CSG support
    int TLOSurface () const { return tlosurf; }
    void SetTLOSurface (int nr) { tlosurf = nr; }

    int DomainIn () const { return domin; }
    void SetDomainIn (int nr) { domin = nr; }

    int DomainOut () const { return domout; }
    void SetDomainOut (int nr) { domout = nr; }

    /// face descriptor index (1-based), INVALID if not yet set. Transient, recomputed by RebuildFDIndices().
    FaceRegionIndex GetIndex () const { return index_; }
    void SetIndex (FaceRegionIndex i) { index_ = i; }

    // deprecated aliases
    [[deprecated("use GetIndex()")]] int FDIndex () const { return index_.Nr1(); }

    void DoArchive (Archive & ar)
    {
      string aname = GetName();
      ar & edgenr & surfnr[0] & surfnr[1] & aname
        & singedge_left & singedge_right & tlosurf
        & domin & domout;
      if (ar.Input()) SetName(aname);
    }

    friend inline ostream & operator<< (ostream & ost, const Region & ed)
    {
      ost << "EdgeDescriptor(edgenr=" << ed.edgenr << ", surfnr=(" << ed.surfnr[0] << "," << ed.surfnr[1] << "), domin=" << ed.domin << ", domout=" << ed.domout << ", name=" << ed.GetName() << ")";
      return ost;
    }
  };

  using EdgeRegion = Region<1>;
  using VolumeRegion = Region<3>;
  using VertexRegion = Region<0>;

  // old names
  using FaceDescriptor = FaceRegion;
  using EdgeDescriptor = EdgeRegion;

  template <int D>
  using RegionArray = Array<Region<D>, RegionIndex<D>>;

  struct BoundaryLayerParameters
  {
    std::variant<double, std::vector<double>> thickness;
    std::variant<string, int, std::vector<int>> domain;
    std::variant<string, int, std::vector<int>> boundary = ".*";
    std::optional<std::variant<string, std::map<string, string>>> new_material = nullopt;
    std::optional<std::variant<string, std::vector<int>>> project_boundaries = nullopt;
    bool outside = false;
    bool grow_edges = true;
    bool limit_growth_vectors = false; // automatic reduction of layer thickness to avoid intersections
    std::optional<bool> sides_keep_surfaceindex = nullopt; // !outside by default
    bool disable_curving = true; // disable curving affected boundaries/edges (could lead to self-intersecting volume elements)
  };


  ostream & operator<< (ostream & ost, const BoundaryLayerParameters & mp);

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

    Array<BoundaryLayerParameters> boundary_layers;
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
      MeshSizePoint (Point<3> pnt_, double h_, int layer_ = 1) : pnt(pnt_), h(h_), layer(layer_) { ; }
      MeshSizePoint () = default;
      MeshSizePoint (const MeshSizePoint &) = default;
      MeshSizePoint (MeshSizePoint &&) = default;
      MeshSizePoint & operator= (const MeshSizePoint &) = default;
      MeshSizePoint & operator= (MeshSizePoint &&) = default;      
    };
    Array<MeshSizePoint> meshsize_points;
    
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
    PointIndex haltsegmentp1;
    ///
    PointIndex haltsegmentp2;
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

  inline void ElementRef :: GetFace (int i, Element2d & face) const
  {
    if (h->typ == TET)
      {
        face.SetType(TRIG);
        face[0] = pn[gftetfacesa[i-1][0]];
        face[1] = pn[gftetfacesa[i-1][1]];
        face[2] = pn[gftetfacesa[i-1][2]];
      }
    else
      GetFace2 (i, face);
  }



  // typedef Array<PointIndex,PointIndex::BASE> idmap_type;
  typedef Array<PointIndex,PointIndex> idmap_type;
  


  /**
     Identification of periodic surfaces, close surfaces, etc. 
  */
  class Identifications
  {
  public:
    enum ID_TYPE : unsigned char { UNDEFINED = 1, PERIODIC = 2, CLOSESURFACES = 3, CLOSEEDGES = 4, OFFSET_POINT = 5};
  

  private:
    class Mesh & mesh;

    /// identify points (thin layers, periodic b.c.)  
    // INDEX_2_HASHTABLE<int> identifiedpoints;
    ClosedHashTable<PointIndices<2>, int> identifiedpoints;
  
    /// the same, with info about the id-nr
    // INDEX_3_HASHTABLE<int> identifiedpoints_nr;
    ClosedHashTable<std::tuple<PointIndices<2>, int>, int> identifiedpoints_nr;

    /// sorted by identification nr
    TABLE<PointIndices<2>> idpoints_table;

    Array<ID_TYPE> type;

    /// number of identifications (or, actually used identifications ?)
    int maxidentnr;
    Array<string> names;

  public:
    ///
    DLL_HEADER Identifications (class Mesh & amesh);
    ///
    DLL_HEADER ~Identifications ();

    DLL_HEADER void Delete ();

    // Removes identifications if one point is an INNERPOINT
    DLL_HEADER void DeleteInnerPointIdentifications ();

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
    auto & GetIdentifiedPoints ()
    { 
      return identifiedpoints_nr;
    }

    bool Used (PointIndex pi1, PointIndex pi2)
    {
      // return identifiedpoints.Used (IVec<2> (pi1, pi2));
      return identifiedpoints.Used (PointIndices<2>(pi1, pi2));
    }

    bool UsedSymmetric (PointIndex pi1, PointIndex pi2)
    {
      return 
        identifiedpoints.Used (PointIndices<2>(pi1, pi2)) ||
        identifiedpoints.Used (PointIndices<2>(pi2, pi1));
    }

    ///
    void GetMap (int identnr, idmap_type & identmap, bool symmetric = false) const;
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
    DLL_HEADER void GetPairs (int identnr, Array<PointIndices<2>> & identpairs) const;
    DLL_HEADER Array<std::tuple<PointIndices<2>, int>> GetPairs () const;
    ///
    int GetMaxNr () const { return maxidentnr; }  

    int GetNr(string name)
    {
      if(!names.Contains(name))
         names.Append(name);
      return names.Pos(name)+1;
    }
    string GetName(int nr) const
    {
      if (nr <= names.Size())
        return names[nr - 1];
      else
        return "";
    }
    void SetName(int nr, string name)
    {
      while(names.Size() < nr)
        names.Append("");
      names[nr-1] = name;
    }

    /// remove secondorder
    void SetMaxPointNr (int maxpnum);

    void MapPoints(FlatArray<PointIndex, PointIndex> op2np);

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
  template <int D> struct MPI_typetrait<netgen::RegionIndex<D>> {
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


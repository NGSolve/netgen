#ifndef FILE_ADFRONT3
#define FILE_ADFRONT3

/**************************************************************************/
/* File:   adfront3.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
  Advancing front class for volume meshing
*/

#include <gprim/geomobjects.hpp>
#include <gprim/adtree.hpp>
#include "meshtype.hpp"
#include "geomsearch.hpp"

namespace netgen
{

/// Point in advancing front
class FrontPoint3
{
  /// coordinates
  Point<3> p;           
  /// global node index
  PointIndex globalindex;   
  /// number of faces connected to point 
  int nfacetopoint;    
  /// distance to original boundary
  int frontnr;
  /// 
  Front3PointIndex cluster;
public:
  ///
  FrontPoint3 ();
  ///
  FrontPoint3 (const Point<3> & ap, PointIndex agi);
  
  ///
  const Point<3> & P () const
  { return p; }
  ///
  PointIndex GlobalIndex () const
  { return globalindex; }
  
  ///
  void AddFace ()
  { nfacetopoint++; }

  /// if last face is removed, then point is invalidated
  void RemoveFace()
  { 
    nfacetopoint--;
    if (nfacetopoint == 0) nfacetopoint = -1;
  }
  
  ///
  bool Valid () const
  { return nfacetopoint >= 0; }

  ///
  void DecFrontNr (int afrontnr)
  {
    if (frontnr > afrontnr) frontnr = afrontnr;
  }
  
  ///
  int FrontNr () const
  { return frontnr; }

  ///
  friend class AdFront3;
};



/// Face in advancing front
class FrontFace
{
private:
  ///
  FrontElement2d f;
  ///
  int qualclass;
  ///
  char oldfront;
  ///
  int hashvalue;
  ///
  Front3PointIndex cluster;
  
public:
  ///
  FrontFace ();
  ///
  FrontFace (const FrontElement2d & af);
  ///
  const FrontElement2d & Face () const
  { return f; }
  
  ///
  int QualClass () const
  { return qualclass; }

  ///
  void IncrementQualClass ()
  { qualclass++; }

  ///
  void ResetQualClass ()
  {
    if (qualclass > 1)
      {
        qualclass = 1;
        oldfront = 0;
      }
  }
  
  ///
  bool Valid () const
  { return !f.IsDeleted(); }

  ///
  void Invalidate ();

  ///
  int HashValue() const 
  { return hashvalue; }

  ///
  void SetHashValue(int hv) 
  { hashvalue = hv; }

  ///
  friend class AdFront3;

  Front3PointIndex Cluster () const { return cluster; }
};  




/// Advancing front, 3D.
class AdFront3
{
  ///
  // Array<FrontPoint3, PointIndex::BASE, PointIndex> points;
  Array<FrontPoint3, Front3PointIndex> points;
  ///
  Array<FrontFace> faces;
  ///
  Array<Front3PointIndex> delpointl;
  
  /// which points are connected to pi ?
  // TABLE<PointIndex, PointIndex::BASE> * connectedpairs;
  unique_ptr<DynamicTable<Front3PointIndex, Front3PointIndex>> connectedpairs;
  
  /// number of total front faces;
  int nff;
  /// number of quads in front
  int nff4; 
  
  ///
  double vol;
  
  ///
  GeomSearch3d hashtable;
  
  /// 
  int hashon;

  ///
  int hashcreated;
  
  /// counter for rebuilding internal tables
  int rebuildcounter;
  /// last base element
  int lasti;
  /// minimal selection-value of baseelements
  int minval;
  Array<LocalPointIndex, Front3PointIndex> invpindex;   // front -> local
  Array<char, Front3PointIndex> pingroup;
  
  ///
  class BoxTree<3> * facetree;
public:
  
  ///
  AdFront3 ();
  ///
  ~AdFront3 ();
  ///
  void GetPoints (Array<Point<3> > & apoints) const;
  ///
  int GetNP() const 
  { return points.Size(); }
  ///
  const Point<3> & GetPoint (Front3PointIndex pi) const
  { return points[pi].P(); }
  ///
  int GetNF() const
  { return nff; }
  /// 1-based
  const FrontElement2d & GetFace (int i) const
  { return faces[i-1].Face(); }
  const auto & Faces() const { return faces; }
  ///
  void Print () const;
  ///
  bool Empty () const
  { return nff == 0; }
  ///
  bool Empty (int elnp) const
  {
    if (elnp == 4)
      return (nff4 == 0);
    return (nff - nff4 == 0);
  }
  ///
  int SelectBaseElement ();

  ///
  void CreateTrees ();

  ///
  void GetIntersectingFaces (const Point<3> & pmin, const Point<3> & pmax, 
                             Array<int> & ifaces) const;

  bool PointInsideGroup(const Array<Front3PointIndex, LocalPointIndex> &grouppindex,
                        const Array<MiniElement2d>& groupfaces) const;

  ///
  void GetFaceBoundingBox (int i, Box3d & box) const;

  ///
  int GetLocals (int baseelement,
                 Array<Point<3>, LocalPointIndex> & locpoints,
                 Array<MiniElement2d> & locfaces,   // local index
                 Array<Front3PointIndex, LocalPointIndex> & pindex,   // local -> front
                 Array<int> & findex,
                 ClosedHashTable<IVec<2>,int> & connectedpairs,
                 float xh,
                 float relh,
                 int& facesplit);
  
  ///
  void GetGroup (int fi,
                 Array<MeshPoint, LocalPointIndex> & grouppoints,
                 Array<MiniElement2d> & groupelements,
                 Array<Front3PointIndex, LocalPointIndex> & pindex,
                 Array<int> & findex);

  ///
  void DeleteFace (int fi);
  ///
  Front3PointIndex AddPoint (const Point<3> & p, PointIndex globind);
  ///
  int AddFace (const FrontElement2d & e);
  ///
  int AddConnectedPair (IVec<2,Front3PointIndex> pair);
  ///
  void IncrementClass (int fi)
  { faces[fi-1].IncrementQualClass(); }

  ///
  void ResetClass (int fi)
  { faces[fi-1].ResetQualClass(); }

  ///
  void SetStartFront (int baseelnp = 0);

  /// is Point p inside Surface ?
  bool Inside (const Point<3> & p) const;
  /// both points on same side ?
  int SameSide (const Point<3> & lp1, const Point<3> & lp2, 
                const Array<int> * testfaces = NULL) const;


  ///
  PointIndex GetGlobalIndex (Front3PointIndex pi) const
  { return points[pi].GlobalIndex(); }
  ///
  double Volume () const
  { return vol; }


private:
  void RebuildInternalTables();
};

} // namespace netgen
#endif

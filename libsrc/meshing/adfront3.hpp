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
  PointIndex cluster;
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



class MiniElement2d
{
protected:
  int np;
  PointIndex pnum[4]; // can be global or local nums
  bool deleted;
public:
  MiniElement2d ()
  { np = 3; deleted = 0; }
  MiniElement2d (int anp)
  { np = anp; deleted = 0; }

  int GetNP() const { return np; }
  PointIndex & operator[] (int i) { return pnum[i]; }
  const PointIndex operator[] (int i) const { return pnum[i]; }

  const PointIndex PNum (int i) const { return pnum[i-1]; }
  PointIndex & PNum (int i) { return pnum[i-1]; }
  const PointIndex PNumMod (int i) const { return pnum[(i-1)%np]; }
  auto PNums() { return FlatArray<PointIndex> (np, &pnum[0]); }
  auto PNums() const { return FlatArray<const PointIndex> (np, &pnum[0]); }
  void Delete () { deleted = true; for (PointIndex & p : pnum) p.Invalidate(); }
  bool IsDeleted () const { return deleted; }
};


inline ostream & operator<<(ostream  & s, const MiniElement2d & el)
{
  s << "np = " << el.GetNP();
  for (int j = 0; j < el.GetNP(); j++)
    s << " " << el[j];
  return s;
}




/// Face in advancing front
class FrontFace
{
private:
  ///
  MiniElement2d f;
  ///
  int qualclass;
  ///
  char oldfront;
  ///
  int hashvalue;
  ///
  PointIndex cluster;
  
public:
  ///
  FrontFace ();
  ///
  FrontFace (const MiniElement2d & af);
  ///
  const MiniElement2d & Face () const
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

  PointIndex Cluster () const { return cluster; }
};  




/// Advancing front, 3D.
class AdFront3
{
  ///
  // NgArray<FrontPoint3, PointIndex::BASE, PointIndex> points;
  Array<FrontPoint3, PointIndex> points
;
  ///
  NgArray<FrontFace> faces;
  ///
  Array<PointIndex> delpointl;
  
  /// which points are connected to pi ?
  // TABLE<PointIndex, PointIndex::BASE> * connectedpairs;
  unique_ptr<DynamicTable<PointIndex, PointIndex>> connectedpairs;
  
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
  Array<PointIndex, PointIndex> invpindex;
  Array<char, PointIndex> pingroup;
  
  ///
  class BoxTree<3> * facetree;
public:
  
  ///
  AdFront3 ();
  ///
  ~AdFront3 ();
  ///
  void GetPoints (NgArray<Point<3> > & apoints) const;
  ///
  int GetNP() const 
  { return points.Size(); }
  ///
  const Point<3> & GetPoint (PointIndex pi) const
  { return points[pi].P(); }
  ///
  int GetNF() const
  { return nff; }
  /// 1-based
  const MiniElement2d & GetFace (int i) const
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
			     NgArray<int> & ifaces) const;

  bool PointInsideGroup(const Array<PointIndex, PointIndex> &grouppindex,
                        const Array<MiniElement2d>& groupfaces) const;

  ///
  void GetFaceBoundingBox (int i, Box3d & box) const;

  ///
  int GetLocals (int baseelement,
		 Array<Point3d, PointIndex> & locpoints,
                 Array<MiniElement2d> & locfaces,   // local index
                 Array<PointIndex, PointIndex> & pindex,
                 Array<INDEX> & findex,
		 INDEX_2_HASHTABLE<int> & connectedpairs,
                 float xh,
		 float relh,
		 INDEX& facesplit);
  
  ///
  void GetGroup (int fi,
                 Array<MeshPoint, PointIndex> & grouppoints,
                 Array<MiniElement2d> & groupelements,
                 Array<PointIndex, PointIndex> & pindex,
                 Array<INDEX> & findex);

  ///
  void DeleteFace (INDEX fi);
  ///
  PointIndex AddPoint (const Point<3> & p, PointIndex globind);
  ///
  INDEX AddFace (const MiniElement2d & e);
  ///
  INDEX AddConnectedPair (PointIndices<2> pair);
  ///
  void IncrementClass (INDEX fi)
  { faces[fi-1].IncrementQualClass(); }

  ///
  void ResetClass (INDEX fi)
  { faces[fi-1].ResetQualClass(); }

  ///
  void SetStartFront (int baseelnp = 0);

  /// is Point p inside Surface ?
  bool Inside (const Point<3> & p) const;
  /// both points on same side ?
  int SameSide (const Point<3> & lp1, const Point<3> & lp2, 
		const NgArray<int> * testfaces = NULL) const;


  ///
  PointIndex GetGlobalIndex (PointIndex pi) const
  { return points[pi].GlobalIndex(); }
  ///
  double Volume () const
  { return vol; }


private:
  void RebuildInternalTables();
};

} // namespace netgen
#endif

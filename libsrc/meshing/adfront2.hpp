#ifndef NETGEN_ADFRONT2_HPP
#define NETGEN_ADFRONT2_HPP

/**************************************************************************/
/* File:   adfront2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/


/**

    Advancing front class for surfaces

*/

#include <gprim/geomobjects.hpp>
#include <gprim/adtree.hpp>
#include "meshtype.hpp"

namespace netgen
{
  ///
  class FrontPoint2
  {
    /// coordinates
    Point<3> p;            
    /// global node index
    PointIndex globalindex;   
    /// number of front lines connected to point 
    int nlinetopoint;    
    /// distance to original boundary
    int frontnr;          

    bool onsurface;

  public:
    ///
    MultiPointGeomInfo * mgi;

    ///
    FrontPoint2 ()
    {
      globalindex.Invalidate(); //  = -1;
      nlinetopoint = 0;
      frontnr = INT_MAX-10;    // attention: overflow on calculating  INT_MAX + 1
      mgi = NULL;
      onsurface = true;
    }

    ///
    FrontPoint2 (const Point<3> & ap, PointIndex agi,
		 MultiPointGeomInfo * amgi, bool aonsurface = true);
    ///
    ~FrontPoint2 () { ; }

    ///
    const Point<3> & P () const { return p; }
    ///
    operator const Point<3> & () const { return p; }
    ///
    PointIndex GlobalIndex () const { return globalindex; }

    ///
    void AddLine () { nlinetopoint++; }
    ///
    void RemoveLine ()
    {
      nlinetopoint--;
      if (nlinetopoint == 0)
	nlinetopoint = -1;
    }

    ///
    bool Valid () const
    { return nlinetopoint >= 0; }

    ///
    bool OnSurface() const
    { return onsurface; }

    ///
    void DecFrontNr (int afrontnr)
    {
      if (frontnr > afrontnr) frontnr = afrontnr;
    }
    
    ///
    int FrontNr () const { return frontnr; }
  };

  
  ///
  class FrontLine
  {
  private:
    /// front point indizes
    IVec<2,Front2PointIndex> l;
    /// quality class 
    int lineclass;      
    /// geometry specific data
    PointGeomInfo geominfo[2];
  public:

    FrontLine ()
    {
      lineclass = 1;
    }

    ///
    FrontLine (const IVec<2,Front2PointIndex> & al)
      : l(al), lineclass(1) { } 

    ///
    const auto & L () const { return l; }
    ///
    int LineClass() const { return lineclass; }

    ///
    void IncrementClass ()
    {
      lineclass++;
    }
    ///
    void ResetClass ()
    {
      lineclass = 1;
    }

    ///
    bool Valid () const
    {
      return l[0].IsValid();
    }
    ///
    void Invalidate ()
    {
      l[0] = Front2PointIndex::INVALID;
      l[1] = Front2PointIndex::INVALID;
      lineclass = 1000;
    }

    void SetGeomInfo (const PointGeomInfo & gi1, const PointGeomInfo & gi2)
      {
	geominfo[0] = gi1;
	geominfo[1] = gi2;
      }

    const PointGeomInfo * GetGeomInfo () const
    { return geominfo; }
    
    const PointGeomInfo & GetGeomInfo (int endp) const
    { return geominfo[endp-1]; }

    friend class AdFront2;
  };


class AdFront2
{

  ///
  Array<FrontPoint2, Front2PointIndex> points;  /// front points
  Array<FrontLine> lines;     /// front lines

  Box3d boundingbox;
  BoxTree<3> linesearchtree;       /// search tree for lines
  Point3dTree pointsearchtree;    /// search tree for points
  Point3dTree cpointsearchtree;   /// search tree for cone points (not used ???)

  Array<Front2PointIndex> delpointl;     /// list of deleted front points
  Array<int> dellinel;      /// list of deleted front lines

  int nfl;                  /// number of front lines;
  INDEX_2_HASHTABLE<int> * allflines; /// all front lines ever have been

  Array<int, Front2PointIndex> invpindex;   // front -> local number (0/-1 markers)

  int minval;
  int starti;

public:
  ///
  //  AdFront2 ();
  AdFront2 (const Box3d & aboundingbox);
  ///
  ~AdFront2 ();

  ///
  // void GetPoints (NgArray<Point<3> > & apoints) const;
  ///
  void Print (ostream & ost) const;

  ///
  bool Empty () const
  {
    return nfl == 0;
  }
  ///
  int GetNFL () const { return nfl; }

  const FrontLine & GetLine (int nr) const { return lines[nr]; }
  const FrontPoint2 & GetPoint (Front2PointIndex nr) const { return points[nr]; }
  const auto & GetLines () const { return lines; }

  ///
  int SelectBaseLine (Point<3> & p1, Point<3> & p2, 
		      const PointGeomInfo *& geominfo1,
		      const PointGeomInfo *& geominfo2,
		      int & qualclass);

  ///
  int GetLocals (int baseline, 
		 Array<Point<3>, LocalPointIndex> & locpoints,
		 Array<MultiPointGeomInfo, LocalPointIndex> & pgeominfo,
                 NgArray<IVec<2,LocalPointIndex>> & loclines,
                 Array<Front2PointIndex, LocalPointIndex> & pindex,   // local -> front
                 NgArray<int> & lindex,
                 double xh);

  ///
  void DeleteLine (int li);
  ///
  Front2PointIndex AddPoint (const Point<3> & p, PointIndex globind, 
                MultiPointGeomInfo * mgi = NULL,
                bool pointonsurface = true);
  ///
  int AddLine (Front2PointIndex pi1, Front2PointIndex pi2, 
               const PointGeomInfo & gi1, const PointGeomInfo & gi2);
  ///
  /// gpi are mesh point numbers
  int ExistsLine (PointIndex gpi1, PointIndex gpi2);

  ///
  void IncrementClass (int li)
  {
    lines[li].IncrementClass();
  }

  ///
  void ResetClass (int li)
  {
    lines[li].ResetClass();
  }

  ///
  const PointGeomInfo & GetLineGeomInfo (int li, int lend) const
    { return lines[li].GetGeomInfo (lend); }
  ///

  PointIndex GetGlobalIndex (Front2PointIndex pi) const
  {
    return points[pi].GlobalIndex();
  }


  /// is Point p inside Surface (flat geometry only)
  bool Inside (const Point<2> & p) const;

  bool SameSide (const Point<2> & lp1, const Point<2> & lp2, 
                 const FlatArray<int> * /* testfaces */ = NULL) const;
  /*
  {
    return Inside (lp1) == Inside (lp2);
  }
  */

  ///
  void SetStartFront ();
  ///
  void PrintOpenSegments (ostream & ost) const;
};

} // namespace netgen
#endif // NETGEN_ADFRONT2_HPP

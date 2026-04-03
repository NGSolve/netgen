#ifndef FILE_STLTOPOLOGY
#define FILE_STLTOPOLOGY

/**************************************************************************/
/* File:   stltopology.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   26. Jul. 99                                                    */
/**************************************************************************/

/*
  The STLTopology contains topologic information as
  triangle->point, point->triangles, triangle->edge, 2-points->edge,...
*/


namespace netgen {

class STLGeometry;

  // #define STLBASE 1

class STLPointId
{
  int i;
public:
  STLPointId () { ; }
  constexpr STLPointId (int ai) : i(ai) { ; }
  STLPointId & operator= (const STLPointId & ai) { i = ai.i; return *this; }
  STLPointId & operator= (int ai) { i = ai; return *this; }
  constexpr operator int () const { return i; }
  STLPointId operator++ (int) { return STLPointId(i++); }    
  STLPointId operator-- (int) { return STLPointId(i--); }
  STLPointId & operator++ () { ++i; return *this; }
  STLPointId & operator-- () { --i; return *this; }
  
  void DoArchive(Archive& ar) { ar & i; }
};



class STLTrigId
{
  int i;
public:
  STLTrigId () { ; }
  constexpr STLTrigId (int ai) : i(ai) { ; }
  STLTrigId & operator= (const STLTrigId & ai) { i = ai.i; return *this; }
  STLTrigId & operator= (int ai) { i = ai; return *this; }
  constexpr operator int () const { return i; }

  STLTrigId operator++ (int) { return STLTrigId(i++); }    
  STLTrigId operator-- (int) { return STLTrigId(i--); }
  STLTrigId & operator++ () { ++i; return *this; }
  STLTrigId & operator-- () { --i; return *this; }

  int operator- (STLTrigId i2) const { return i-i2.i; }
};

  inline void SetInvalid (STLTrigId & id) { id = 0; }
  inline bool IsInvalid (STLTrigId & id) { return id == 0; }

}


namespace ngcore
{
  template<> 
  constexpr netgen::STLPointId IndexBASE<netgen::STLPointId> () { return netgen::STLPointId(1); }
  template<> 
  constexpr netgen::STLTrigId IndexBASE<netgen::STLTrigId> () { return netgen::STLTrigId(1); }
}



namespace netgen {


// triangle structure for loading stl files
class STLReadTriangle
{
  Vec<3> normal;
  Point<3> pts[3];
public:
  STLReadTriangle (const Point<3> * apts, const Vec<3> & anormal);
  STLReadTriangle () {};
  const Point<3> & operator[] (int i) const { return pts[i]; }
  const Vec<3> & Normal() const { return normal; }
};



class STLTriangle
{
  // topology edges of triangle, edge[i] opposite to point[i]
  int topedges[3];
  // neighbour triangles, trig[i] opposite to point[i]
  int nbtrigs[2][3]; 
  // normalized stored normal vector ??
  Vec<3> normal;
  // point numbers of triangle
  STLPointId pts[3];
  // front-side and back-side domains
  int domains[2];


public:

  Box<3> box;
  Point<3> center;
  double rad;
  int facenum;

  struct 
  {
    unsigned int toperror : 1;
  } flags;




  STLTriangle (const STLPointId * apts);
  STLTriangle ()
  {
    pts[0]=0;pts[1]=0;pts[2]=0;
    nbtrigs[0][0] = nbtrigs[0][1] = nbtrigs[0][2] = 0.;
    nbtrigs[1][0] = nbtrigs[1][1] = nbtrigs[1][2] = 0.;
  }

  void DoArchive(Archive& ar)
  {
    ar.Do(&topedges[0],3);
    ar.Do(&nbtrigs[0][0], 6);
    // ar.Do(&pts[0],3);
    ar & pts[0] & pts[1] & pts[2];
    ar.Do(&domains[0],2);
    size_t i = flags.toperror;
    ar & normal & box & center & rad & facenum & i;
    flags.toperror = i;
  }

  STLPointId operator[] (int i) const { return pts[i]; }
  STLPointId & operator[] (int i) { return pts[i]; }

  int EdgeNum(int i) const { return topedges[(i-1)]; }
  int & EdgeNum(int i) { return topedges[(i-1)]; }

  int NBTrig (bool side, int i) const { return nbtrigs[side][i]; }
  int & NBTrig (bool side, int i) { return nbtrigs[side][i]; }

  
  int Domain (bool side) const { return domains[side]; }
  int & Domain (bool side) { return domains[side]; }



  // obsolete:
  STLPointId PNum(int i) const { return pts[(i-1)]; }
  STLPointId & PNum(int i) { return pts[(i-1)]; }
  STLPointId PNumMod(int i) const { return pts[(i-1)%3]; }
  STLPointId & PNumMod(int i)  { return pts[(i-1)%3]; }
  FlatArray<const STLPointId> PNums() const { return { 3, pts }; }

  
  int EdgeNumMod(int i) const { return topedges[(i-1)%3]; }
  int & EdgeNumMod(int i)  { return topedges[(i-1)%3]; }

  int NBTrigNum(int i) const { return nbtrigs[0][(i-1)]; }
  int & NBTrigNum(int i) { return nbtrigs[0][(i-1)]; }
  int NBTrigNumMod(int i) const { return nbtrigs[0][(i-1)%3]; }
  int & NBTrigNumMod(int i)  { return nbtrigs[0][(i-1)%3]; }
  

  // consistently oriented neighbour:
  int IsNeighbourFrom(const STLTriangle& t) const;
  // opposite to consistently oriented neighbour:
  int IsWrongNeighbourFrom(const STLTriangle& t) const;

  ///Get the two points of neighbour-Triangles in orientation of this-Triangle
  void GetNeighbourPoints(const STLTriangle& t, STLPointId & p1, STLPointId & p2) const;
  int GetNeighbourPointsAndOpposite(const STLTriangle& t, STLPointId & p1, STLPointId & p2, STLPointId & po) const;



  // NON-normalized geometry - normal vector
  Vec<3> GeomNormal(const Array<Point<3>,STLPointId>& ap) const;
  
  // Stored normal vector, normalized
  void SetNormal (const Vec<3> & n);
  const Vec<3> & Normal () const { return normal; }


  void ChangeOrientation(); 

  //project with a certain normal vector in plane
  void ProjectInPlain(const Array<Point<3>, STLPointId>& ap, 
		      const Vec<3> & n, Point<3> & pp) const;
  //project with the triangle's normal vector in plane
  void ProjectInPlain(const Array<Point<3>, STLPointId> & ap, Point<3> & pp) const;


  /*
    Project the point pp along the nproj into the plane of
    the triangle. The triangle normal is given by ntrig to 
    avoid numerical instabilities.
    The local coordinates lam are defined by

    pp(input) = P1 + lam1 v1 + lam2 v2 + lam3 n

    the result is
    
    pp(output) = P1 + lam1 v1 + lam2 v2
  */
  int ProjectInPlain (const Array<Point<3>,STLPointId>& ap, 
		      const Vec<3> & nproj, 
		      Point<3> & pp, Vec<3> & lam) const;

  bool PointInside(const Array<Point<3>,STLPointId>& ap, const Point<3> & pp) const;

  //get nearest point on triangle and distance to it
  double GetNearestPoint(const Array<Point<3>,STLPointId>& ap, 
			 Point<3> & p3d) const;

  double Area(const Array<Point<3>,STLPointId>& ap) const;

  double MinHeight(const Array<Point<3>,STLPointId>& ap) const;
  double MaxLength(const Array<Point<3>,STLPointId>& ap) const; 
  //max length of a side of triangle

  int GetFaceNum() {return facenum;}
  void SetFaceNum(int i) {facenum = i;}

  bool HasEdge(STLPointId p1, STLPointId p2) const;
};


/**
   Topology Edge:
   Useful unside a face.
   A edges sharing more than 2 faces: trigs are undefined 
 */
class STLTopEdge 
{
  STLPointId pts[2];  
  int trigs[2];  
  double cosangle;
  int status;  // excluded, confirmed, candidate, undefined
public:
  STLTopEdge ();
  STLTopEdge (STLPointId p1, STLPointId p2, int trig1, int trig2);

  STLPointId operator[] (int i) const { return pts[i]; }
  STLPointId & operator[] (int i) { return pts[i]; }


  STLPointId PNum(int i) const { return pts[(i-1)]; }
  STLPointId & PNum(int i) { return pts[(i-1)]; }
  STLPointId PNumMod(int i) const { return pts[(i-1)%2]; }
  STLPointId & PNumMod(int i)  { return pts[(i-1)%2]; }

  int TrigNum(int i) const { return trigs[(i-1)]; }
  int & TrigNum(int i) { return trigs[(i-1)]; }
  int TrigNumMod(int i) const { return trigs[(i-1)%2]; }
  int & TrigNumMod(int i)  { return trigs[(i-1)%2]; }

  void SetCosAngle (double ca) { cosangle = ca; }
  double CosAngle () const { return cosangle; }
  double Angle () const { return acos (cosangle); }

  void SetStatus (int stat) { status = stat; }
  int GetStatus () const { return status; }
};



ostream& operator<<(ostream& os, const STLTriangle& t);







class STLTopology
{
protected:
  Array<STLTriangle, STLTrigId> trias;
  NgArray<STLTopEdge> topedges;
  Array<Point<3>, STLPointId> points;
  bool surface = false;

  // mapping of sorted pair of points to topedge
  INDEX_2_HASHTABLE<int> * ht_topedges;
  // mapping of node to trigs
  TABLE<int, IndexBASE<STLPointId>()> trigsperpoint; 
  // mapping of node to edges
  TABLE<int> topedgesperpoint; 
  
  // searchtree for trigs and points

  BoxTree<3> * searchtree; // ADT
  Point3dTree * pointtree;

  Box<3> boundingbox;
  double pointtol;

public:
  enum STL_GEOM_STATUS { STL_GOOD, STL_WARNING, STL_ERROR };

protected:
  STL_GEOM_STATUS status;
  string statustext;
  
  bool topology_ok;
  bool orientation_ok;

public:
  STLTopology();
  virtual ~STLTopology();

  static STLGeometry * LoadNaomi (istream & ist);
  DLL_HEADER static STLGeometry * Load (istream & ist, bool surface=false);
  static STLGeometry * LoadBinary (istream & ist);

  void Save (const filesystem::path & filename) const;
  void SaveBinary (const filesystem::path & filename, const char* aname) const;
  void SaveSTLE (const filesystem::path & filename) const; // stores trigs and edges

  bool IsSurfaceSTL() const { return surface; }
  void SetSurfaceSTL( bool surface_ ) { surface = surface_; }

  virtual void DoArchive(Archive& ar)
  {
    ar & trias & points & boundingbox & pointtol;
    if(ar.Input())
      FindNeighbourTrigs();
  }
  
  virtual void InitSTLGeometry (const NgArray<STLReadTriangle> & readtrigs);

  virtual void TopologyChanged() {}; //do some things, if topology changed!

  /// Generate topology tables
  void FindNeighbourTrigs();

  
  void GetTrianglesInBox (const Box<3> & box,
			  NgArray<int> & trias) const;


  int GetNP() const { return points.Size(); }
  int AddPoint(const Point<3> & p) { points.Append(p); return points.Size(); }
  const Point<3> & GetPoint(STLPointId nr) const { return points[nr]; } // .Get(nr); }
  int GetPointNum (const Point<3> & p);
  void SetPoint(STLPointId nr, const Point<3> & p) { points[nr] = p; } // { points.Elem(nr) = p; }
  auto & GetPoints() const { return points; }

  const Point<3> & operator[] (STLPointId i) const { return points[i]; }
  Point<3> & operator[] (STLPointId i) { return points[i]; }




  int GetNT() const { return trias.Size(); }
  void AddTriangle(const STLTriangle& t);
  const STLTriangle & GetTriangle (STLTrigId nr) const { return trias[nr]; } // .Get(nr); }
  STLTriangle & GetTriangle (STLTrigId nr) { return trias[nr]; } // .Elem(nr); }
  
  const STLTriangle & operator[] (STLTrigId i) const { return trias[i]; }
  STLTriangle & operator[] (STLTrigId i) { return trias[i]; }


  int GetNTE() const { return topedges.Size(); }
  const STLTopEdge & GetTopEdge (int nr) const { return topedges.Get(nr); }
  STLTopEdge & GetTopEdge (int nr)  { return topedges.Elem(nr); }
  DLL_HEADER int GetTopEdgeNum (int pi1, int pi2) const;


  int NOTrigsPerPoint(int pn) { return trigsperpoint.EntrySize(pn); }
  int TrigPerPoint(int pn, int i) { return trigsperpoint.Get(pn, i); }


  int NTopEdgesPerPoint (int pn) const { return topedgesperpoint.EntrySize(pn); }
  int TopEdgePerPoint (int pn, int ei) const { return topedgesperpoint.Get(pn, ei); }

  
  bool Topology_Ok() const { return topology_ok; }
  bool Orientation_Ok() const { return orientation_ok; }

  STL_GEOM_STATUS GetStatus () const { return status; }
  const string & GetStatusText () const { return statustext; }

  DLL_HEADER void InvertTrig (int trig);
  DLL_HEADER void DeleteTrig (int trig);
  DLL_HEADER void OrientAfterTrig (int trig);


  // Table will be constructed, if topology is not ok
  /// neighbourtrigs for surfacetrigs
  TABLE<STLTrigId> neighbourtrigs;

  /// get nr-th neighbour Triangle for triangle trig
  int NONeighbourTrigs(STLTrigId trig) const { return neighbourtrigs.EntrySize(int(trig)); }
  STLTrigId NeighbourTrig(STLTrigId trig, int nr) const { return neighbourtrigs.Get(int(trig),nr); }
  int NeighbourTrigSorted(int trig, int nr) const;
  void AddNeighbourTrig(STLTrigId i, STLTrigId nt) { neighbourtrigs.Add1(int(i), nt); }




  int GetLeftTrig (int p1, int p2) const;
  int GetRightTrig (int p1, int p2) const;

  const Box<3> & GetBoundingBox () const { return boundingbox; }
};

} // namespace netgen

#endif

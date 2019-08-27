#ifndef FILE_STLTOOL
#define FILE_STLTOOL


//#include "gprim/gprim.hh"

/**************************************************************************/
/* File:   stlgeom.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/



// use one normal vector for whole chart
extern int usechartnormal;
extern int chartdebug;

extern int geomsearchtreeon;
extern int AddPointIfNotExists(NgArray<Point3d>& ap, const Point3d& p, double eps = 1e-8);
//get distance from line lp1-lp2 to point p
extern double GetDistFromLine(const Point<3>& lp1, const Point<3>& lp2, Point<3>& p);
extern double GetDistFromInfiniteLine(const Point<3>& lp1, const Point<3>& lp2, const Point<3>& p);


extern void FIOReadInt(istream& ios, int& i);
extern void FIOWriteInt(ostream& ios, const int& i);
extern void FIOReadDouble(istream& ios, double& i);
extern void FIOWriteDouble(ostream& ios, const double& i);
extern void FIOReadFloat(istream& ios, float& i);
extern void FIOWriteFloat(ostream& ios, const float& i);
extern void FIOReadString(istream& ios, char* str, int len);
extern void FIOReadStringE(istream& ios, char* str, int len);
extern void FIOWriteString(ostream& ios, char* str, int len);


typedef NgArray <int> * ArrayINTPTR;

class STLGeometry;
class STLParameters;

class STLChart
{
private:
  STLGeometry * geometry;
  NgArray<int> charttrigs; // trigs which only belong to this chart
  NgArray<int> outertrigs; // trigs which belong to other charts
  BoxTree<3> * searchtree; // ADT containing outer trigs

  NgArray<twoint> olimit; //outer limit of outer chart
  NgArray<twoint> ilimit; //outer limit of inner chart
  const STLParameters& stlparam;


public:
  
  STLChart(STLGeometry * ageometry, const STLParameters& astlparam);
  ~STLChart();
  void AddChartTrig(int i);
  void AddOuterTrig(int i);
  
  int IsInWholeChart(int nr) const;

  int GetChartTrig(int i) const {return charttrigs.Get(i);}
  int GetOuterTrig(int i) const {return outertrigs.Get(i);}
  //get all trigs:
  int GetTrig(int i) const
    {
      if (i <= charttrigs.Size()) {return charttrigs.Get(i);}
      else {return outertrigs.Get(i-charttrigs.Size());}
    }
  
  int GetNChartT() const {return charttrigs.Size();}
  int GetNOuterT() const {return outertrigs.Size();}
  int GetNT() const {return charttrigs.Size()+outertrigs.Size(); }

  void GetTrianglesInBox (const Point3d & pmin,
			  const Point3d & pmax,
			  NgArray<int> & trias) const;
  void AddOLimit(twoint l) {olimit.Append(l);}
  void AddILimit(twoint l) {ilimit.Append(l);}

  void ClearOLimit() {olimit.SetSize(0);}
  void ClearILimit() {ilimit.SetSize(0);}

  int GetNOLimit() const {return olimit.Size();}
  int GetNILimit() const {return ilimit.Size();}

  twoint GetOLimit(int i) const {return olimit.Get(i);}
  twoint GetILimit(int i) const {return ilimit.Get(i);}

  //move triangles trigs (local chart-trig numbers) to outer chart
  void MoveToOuterChart(const NgArray<int>& trigs);
  void DelChartTrigs(const NgArray<int>& trigs);


  // define local coordinate system, JS:
private:
  Vec<3> normal;
  Point<3> pref;
  Vec<3> t1, t2;
public:
  void SetNormal (const Point<3> & apref, const Vec<3> & anormal);
  const Vec<3> & GetNormal () const { return normal; }
  Point<2> Project2d (const Point<3> & p3d) const
  {
    Vec<3> v = p3d-pref;
    return Point<2> (t1 * v, t2 * v);
  }
};

class STLBoundarySeg
{
  Point<3> p1, p2, center;
  Point<2> p2d1, p2d2;
  Box<2> boundingbox;
  //  Point<2> p2dmin, p2dmax;

  double rad;
  int i1, i2;
  int smoothedge;
public:
  STLBoundarySeg () { ; }
  STLBoundarySeg (int ai1, int ai2, const NgArray<Point<3> > & points,
		  const STLChart * chart)
    : p1(points.Get(ai1)), p2(points.Get(ai2)),
      i1(ai1), i2(ai2)
  {
    center = ::netgen::Center (p1, p2);
    rad = Dist (p1, center);
    
    p2d1 = chart->Project2d (p1);
    p2d2 = chart->Project2d (p2);
    
    boundingbox.Set (p2d1);
    boundingbox.Add (p2d2);
  }

  int operator== (const STLBoundarySeg & s2) const
    { return i1 == s2.i1 && i2 == s2.i2; }
  void Swap ();
  int I1() const { return i1; }
  int I2() const { return i2; }
  const Point<3> & P1() const { return p1; }
  const Point<3> & P2() const { return p2; }
  const Point<2> & P2D1() const { return p2d1; }
  const Point<2> & P2D2() const { return p2d2; }
  const Point<2> & P2DMin() const { return boundingbox.PMin(); }
  const Point<2> & P2DMax() const { return boundingbox.PMax(); }
  const Point<3> & Center() const { return center; }
  const Box<2> & BoundingBox() const { return boundingbox; }
  double Radius () const { return rad; }

  void SetSmoothEdge (int se) { smoothedge = se; }
  int IsSmoothEdge () const { return smoothedge; }
  friend class STLBoundary;
};

class STLBoundary
{
private:
  STLGeometry * geometry;
  const STLChart * chart;
  NgArray<STLBoundarySeg> boundary;
  ClosedHashTable<INDEX_2, STLBoundarySeg> boundary_ht;  
  BoxTree<2,INDEX_2> * searchtree = nullptr;
public:
  STLBoundary(STLGeometry * ageometry);
  ~STLBoundary() { delete searchtree; }

  void Clear() {boundary.SetSize(0); boundary_ht = ClosedHashTable<INDEX_2,STLBoundarySeg>(); }
  void SetChart (const STLChart * achart) { chart = achart; }
  //don't check, if already exists!
  void AddNewSegment(const STLBoundarySeg & seg) {boundary.Append(seg);};
  //check if segment exists
  void AddOrDelSegment(const STLBoundarySeg & seg);
  //addordelsegment for all 3 triangle segments!
  void AddTriangle(const STLTriangle & t);
  int NOSegments() {return boundary.Size();};
  const STLBoundarySeg & GetSegment(int i) {return boundary.Get(i);}

  void BuildSearchTree();
  void DeleteSearchTree();
  int TestSeg(const Point<3> & p1, const Point<3> & p2, const Vec<3> & sn, 
	      double sinchartangle, int divisions, NgArray<Point<3> >& points,
	      double eps);

  int TestSegChartNV(const Point3d& p1, const Point3d& p2, const Vec3d& sn);
};


class STLDoctorParams
{
public:
  int drawmeshededges;
  double geom_tol_fact;

  double longlinefact;
  int showexcluded;

  int selectmode; //0==trig, 1==edge, 2==point, 3==multiedge, 4==line cluster
  int edgeselectmode;

  int useexternaledges;
  int showfaces;
  int showedgecornerpoints;
  int showtouchedtrigchart;
  int conecheck;
  int spiralcheck;
  int selecttrig;
  int nodeofseltrig;
  int selectwithmouse;
  int showmarkedtrigs;
  double dirtytrigfact;
  double smoothangle;

  double smoothnormalsweight;

  int showvicinity;
  int vicinity;
  ///
  STLDoctorParams();
  ///
  void Print (ostream & ost) const;
};

DLL_HEADER extern STLDoctorParams stldoctor;



// TODO change enable flag to optional parameters
class DLL_HEADER STLParameters
{
public:
  /// angle for edge detection
  double yangle = 30.;
  double contyangle = 20.; //edges continued with contyangle
  /// angle of geometry edge at which the mesher should set a point
  double edgecornerangle = 60.;
  /// angle inside on chart
  double chartangle = 15.;
  /// angle for overlapping parts of char
  double outerchartangle = 70.;
  /// 0 .. no, 1 .. local, (2 .. global)
  int usesearchtree = 0;
  ///
  double resthatlasfac = 2.; 
  bool resthatlasenable = true;
  double atlasminh = 0.1;

  double resthsurfcurvfac = 2.; 
  bool resthsurfcurvenable = false;

  double resthchartdistfac = 1.2;
  bool resthchartdistenable = true;

  double resthcloseedgefac = 1.;
  bool resthcloseedgeenable = true;
  
  double resthedgeanglefac = 1.;
  bool resthedgeangleenable = false;
  
  double resthsurfmeshcurvfac = 1.;
  bool resthsurfmeshcurvenable = false;
  
  double resthlinelengthfac = 0.5;
  bool resthlinelengthenable = true;

  ///
  bool recalc_h_opt = true;
  ///
  STLParameters();
  ///
  void Print (ostream & ost) const;
};


void STLMeshing (STLGeometry & geom,
		 Mesh & mesh,
                 const MeshingParameters& mparam,
                 const STLParameters& stlpar);


int STLSurfaceMeshing (STLGeometry & geom,
                       Mesh & mesh,
                       const MeshingParameters& mparam,
                       const STLParameters& stlpar);

void STLSurfaceOptimization (STLGeometry & geom,
			     Mesh & mesh,
			     const MeshingParameters & mparam);




#endif

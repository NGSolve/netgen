#ifndef FILE_EDGEFLW
#define FILE_EDGEFLW

/**************************************************************************/
/* File:   edgeflw.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

namespace netgen
{

  /// Lightweight struct for CSG refedges (replaces Segment as refedge type)
  struct RefEdge
  {
    int si = 0;           // surface index (class representant)
    int surfnr1 = -1;    // original surface 1
    int surfnr2 = -1;    // original surface 2
    int domin = -1;       // domain in (TLO index, or -1)
    int domout = -1;      // domain out (TLO index, or -1)
    int tlosurf = -1;     // TLO surface index, or -1
    int edgenr = 0;       // edge number (1-based)
    int index_ = 0;      // edge descriptor index (1-based, 0 = invalid)

    int GetIndex() const { return index_; }
    void SetIndex(int i) { index_ = i; }

    friend ostream & operator<< (ostream & s, const RefEdge & re)
    {
      s << "GitRefEdge(si=" << re.si
        << ", surfnr=" << re.surfnr1 << "/" << re.surfnr2
        << ", domin=" << re.domin << ", domout=" << re.domout
        << ", tlosurf=" << re.tlosurf << ", edgenr=" << re.edgenr
        << ", index=" << re.index_ << ")";
      return s;
    }
  };


  /*
  
  Edge - following function and
  Projection to edge of implicitly given edge

  */
 

  /**
     Calculates edges.
     The edges of a solid geometry are computed. Special
     points have to be given.
  */
  extern void CalcEdges (const CSGeometry & geometry,
			 const Array<SpecialPoint> & specpoints,
			 double h, Mesh & mesh);





  class EdgeCalculation
  {
    const CSGeometry & geometry;
    Array<SpecialPoint> & specpoints;
    Point3dTree<> * searchtree;
    Point3dTree<PointIndex> * meshpoint_tree;
    int cntedge;

  public:
    Array<char, SegmentIndex> seg_seginfo;

  private:
    double ideps;
    MeshingParameters & mparam;

  public:
    EdgeCalculation (const CSGeometry & ageometry,
		     Array<SpecialPoint> & aspecpoints,
                     MeshingParameters & amparam);

    ~EdgeCalculation();

    void SetIdEps(const double epsin) {ideps = epsin;}

    void Calc(double h, Mesh & mesh);


  private:
    void CalcEdges1 (double h, Mesh & mesh);
  

    void FollowEdge (int pi1, int & ep, int & pos,
		     // const Array<SpecialPoint> & hsp,
		     const Array<int> & hsp,
		     double h, const Mesh & mesh,
		     Array<Point<3> > & edgepoints,
		     Array<double> & curvelength);
		   

    void AnalyzeEdge (int s1, int s2, int s1_rep, int s2_rep, int pos, int layer,
		      const Array<Point<3> > & edgepoints,
		      Array<RefEdge> & refedges,
		      Array<bool> & refedgesinv);

    void StoreEdge (const Array<RefEdge> & refedges,
		    const Array<bool> & refedgesinv,
		    const Array<Point<3> > & edgepoints,
		    const Array<double> & curvelength,
		    int layer,
		    Mesh & mesh);

    void StoreShortEdge (const Array<RefEdge> & refedges,
			 const Array<bool> & refedgesinv,
			 const Array<Point<3> > & edgepoints,
			 const Array<double> & curvelength,
			 int layer,
			 Mesh & mesh);

    void CopyEdge (const Array<RefEdge> & refedges,
		   const Array<bool> & refedgesinv,
		   int copyfromedge, 
		   const Point<3> & fromstart, const Point<3> & fromend,
		   const Point<3> & tostart, const Point<3> & toend,
		   int copyedgeidentification,
		   int layer,
		   Mesh & mesh);

  
    void SplitEqualOneSegEdges (Mesh & mesh);
    void FindClosedSurfaces (double h, Mesh & mesh);


  public:
    bool point_on_edge_problem;

  };

}


#endif

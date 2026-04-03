#ifndef FILE_EDGEFLW
#define FILE_EDGEFLW

/**************************************************************************/
/* File:   edgeflw.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

namespace netgen
{



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
			 const NgArray<SpecialPoint> & specpoints,
			 double h, Mesh & mesh);





  class EdgeCalculation
  {
    const CSGeometry & geometry;
    NgArray<SpecialPoint> & specpoints;
    Point3dTree * searchtree;
    Point3dTree * meshpoint_tree;
    int cntedge;

    double ideps;
    MeshingParameters & mparam;

  public:
    EdgeCalculation (const CSGeometry & ageometry,
		     NgArray<SpecialPoint> & aspecpoints,
                     MeshingParameters & amparam);

    ~EdgeCalculation();

    void SetIdEps(const double epsin) {ideps = epsin;}

    void Calc(double h, Mesh & mesh);


  private:
    void CalcEdges1 (double h, Mesh & mesh);
  

    void FollowEdge (int pi1, int & ep, int & pos,
		     // const NgArray<SpecialPoint> & hsp,
		     const NgArray<int> & hsp,
		     double h, const Mesh & mesh,
		     NgArray<Point<3> > & edgepoints,
		     NgArray<double> & curvelength);
		   

    void AnalyzeEdge (int s1, int s2, int s1_rep, int s2_rep, int pos, int layer,
		      const NgArray<Point<3> > & edgepoints,
		      NgArray<Segment> & refedges,
		      NgArray<bool> & refedgesinv);

    void StoreEdge (const NgArray<Segment> & refedges,
		    const NgArray<bool> & refedgesinv,
		    const NgArray<Point<3> > & edgepoints,
		    const NgArray<double> & curvelength,
		    int layer,
		    Mesh & mesh);

    void StoreShortEdge (const NgArray<Segment> & refedges,
			 const NgArray<bool> & refedgesinv,
			 const NgArray<Point<3> > & edgepoints,
			 const NgArray<double> & curvelength,
			 int layer,
			 Mesh & mesh);

    void CopyEdge (const NgArray<Segment> & refedges,
		   const NgArray<bool> & refedgesinv,
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

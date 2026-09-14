#ifndef NETGEN_MESHING2_HPP
#define NETGEN_MESHING2_HPP

/**************************************************************************/
/* File:   meshing2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/


#include "adfront2.hpp"
#include "ruler2.hpp"
#include "basegeom.hpp"

namespace netgen
{


enum MESHING2_RESULT
{
  MESHING2_OK = 0,
  MESHING2_GIVEUP = 1
};


/*
   
The basis class for 2D mesh generation. 
Has the method GenerateMesh

For surface mesh generation, or non-Euklidean meshing,
derive from Meshing2, and replace transformation.

*/

class Meshing2
{
  /// the current advancing front
  AdFront2 adfront;
  /// mesh point number -> front point number
  Array<Front2PointIndex, PointIndex> glob2front;
  /// rules for mesh generation
  Array<unique_ptr<netrule>> rules;
  /// statistics
  Array<int> ruleused, canuse, foundmap;
  /// 
  Box<3> boundingbox;
  ///
  double starttime;
  ///
  double maxarea;

  Vec<3> ex, ey, ez;
  Point<3> p1, p2;

  const NetgenGeometry& geo;

public:
  ///
  DLL_HEADER Meshing2 (const NetgenGeometry& geo,
                       const MeshingParameters & mp,
                       const Box<3> & aboundingbox);

  ///
  DLL_HEADER virtual ~Meshing2 ();

  /// Load rules, either from file, or compiled rules
  void LoadRules (const char * filename, bool quad);

  /// 
  DLL_HEADER MESHING2_RESULT GenerateMesh (Mesh & mesh, const MeshingParameters & mp, double gh, int facenr, int layer=1);

  DLL_HEADER void Delaunay (Mesh & mesh, int domainnr, const MeshingParameters & mp);
  DLL_HEADER void BlockFillLocalH (Mesh & mesh, const MeshingParameters & mp);


  ///
  DLL_HEADER Front2PointIndex AddPoint (const Point<3> & p, PointIndex globind, MultiPointGeomInfo * mgi = NULL,
		 bool pointonsurface = true);
  DLL_HEADER PointIndex GetGlobalIndex(Front2PointIndex pi) const;

  ///
  /// pi1, pi2 are mesh point numbers
  DLL_HEADER void AddBoundaryElement (PointIndex pi1, PointIndex pi2,
			   const PointGeomInfo & gi1, const PointGeomInfo & gi2);
  /// fpi1, fpi2 are front point numbers (as returned by AddPoint)
  DLL_HEADER void AddBoundaryElement (Front2PointIndex fpi1, Front2PointIndex fpi2,
			   const PointGeomInfo & gi1, const PointGeomInfo & gi2);
  
  ///
  void SetStartTime (double astarttime);

  ///
  void SetMaxArea (double amaxarea);

protected:
  ///
  virtual void StartMesh ();
  ///
  virtual void EndMesh ();
  ///
  virtual double CalcLocalH (const Point<3> & p, double gh) const;

  ///
  virtual void DefineTransformation (const Point<3> & p1, const Point<3> & p2,
				     const PointGeomInfo * geominfo1,
				     const PointGeomInfo * geominfo2);
  ///
  virtual void TransformToPlain (const Point<3> & locpoint, const MultiPointGeomInfo &  geominfo,
				 Point<2> & plainpoint, double h, int & zone);
  /// return 0 .. ok
  /// return >0 .. cannot transform point to true surface
  virtual int TransformFromPlain (const Point<2>& plainpoint,
				  Point<3> & locpoint, 
				  PointGeomInfo & geominfo, 
				  double h);
  
  /// projects to surface
  /// return 0 .. ok
  virtual int BelongsToActiveChart (const Point<3> & p, 
				    const PointGeomInfo & gi);

  /// computes geoinfo data for line with respect to
  /// selected chart
  virtual int ComputePointGeomInfo (const Point<3> & p, 
				    PointGeomInfo & gi);

  /// Tries to select unique geominfo on active chart
  /// return 0: success
  /// return 1: failed
  virtual int ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
					PointGeomInfo & pgi);



  /*
    tests, whether endpoint (= 1 or 2) of line segment p1-p2
    is inside of the selected chart. The endpoint must be on the
    chart
   */
  virtual int IsLineVertexOnChart (const Point<3> & p1, const Point<3> & p2,
				   int endpoint, const PointGeomInfo & geominfo);

  /*
    get (projected) boundary of current chart
   */
  virtual void GetChartBoundary (Array<Point<2>> & points, 
				 Array<Point<3>> & points3d,
				 Array<INDEX_2> & lines, double p) const;

  virtual double Area () const;


/** Applies 2D rules.
 Tests all 2D rules */
  int ApplyRules (Array<Point<2>, LocalPointIndex> & lpoints, 
		  Array<int, LocalPointIndex> & legalpoints,
		  int maxlegalpoint,
		  Array<IVec<2,LocalPointIndex>> & llines,
		  int maxlegelline,
		  Array<MiniElement2d> & elements, Array<INDEX> & dellines,
		  int tolerance,
		  const MeshingParameters & mp);
  

};
} // namespace netgen

#endif // NETGEN_MESHING2_HPP

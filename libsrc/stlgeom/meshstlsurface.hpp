#ifndef FILE_MESHSTLSURF
#define FILE_MESHSTLSURF

/* *************************************************************************/
/* File:   meshstlsurf.hpp                                                 */
/* Author: Johannes Gerstmayr, Joachim Schoeberl                           */
/* Date:   01. Aug. 99                                                     */
/* *************************************************************************/

/*

The interface between mesh generation and stl geometry

*/


/// 
class MeshingSTLSurface : public Meshing2
{
  ///
  STLGeometry & geom;
  ///
  int transformationtrig;
public:
  ///
  MeshingSTLSurface (STLGeometry & ageom, const MeshingParameters & mp);

protected:
  ///
  void DefineTransformation (const Point<3> & p1, const Point<3> & p2,
                             const PointGeomInfo * geominfo1,
                             const PointGeomInfo * geominfo2) override;
  ///
  void TransformToPlain (const Point<3> & locpoint, const MultiPointGeomInfo & geominfo,
                         Point<2> & plainpoint, double h, int & zone) override;
  ///
  int TransformFromPlain (const Point<2>& plainpoint,
                          Point<3> & locpoint, 
                          PointGeomInfo & gi,
                          double h) override;
  ///
  int BelongsToActiveChart (const Point3d & p, 
                            const PointGeomInfo & gi) override;

  ///
  int ComputePointGeomInfo (const Point3d & p, PointGeomInfo & gi) override;
  ///
  int ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
                                PointGeomInfo & pgi) override;

  ///
  int IsLineVertexOnChart (const Point3d & p1, const Point3d & p2,
                           int endpoint, const PointGeomInfo & gi) override;

  void GetChartBoundary (NgArray<Point<2>> & points, 
                         NgArray<Point<3>> & poitns3d,
                         NgArray<INDEX_2> & lines, double h) const override;

  ///
  double CalcLocalH (const Point<3> & p, double gh) const override;

  ///
  double Area () const override;
};

#endif


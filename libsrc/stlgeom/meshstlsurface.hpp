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
  int BelongsToActiveChart (const Point<3> & p, 
                            const PointGeomInfo & gi) override;

  ///
  int ComputePointGeomInfo (const Point<3> & p, PointGeomInfo & gi) override;
  ///
  int ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
                                PointGeomInfo & pgi) override;

  ///
  int IsLineVertexOnChart (const Point<3> & p1, const Point<3> & p2,
                           int endpoint, const PointGeomInfo & gi) override;

  void GetChartBoundary (Array<Point<2>> & points, 
                         Array<Point<3>> & poitns3d,
                         Array<IVec<2>> & lines, double h) const override;

  ///
  double CalcLocalH (const Point<3> & p, double gh) const override;

  ///
  double Area () const override;
};

#endif


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



///
class MeshOptimizeSTLSurface : public MeshOptimize2d
  {
  ///
    STLGeometry & geom;

public:
    ///
    MeshOptimizeSTLSurface (STLGeometry & ageom); 
   
    ///
    virtual void SelectSurfaceOfPoint (const Point<3> & p,
				       const PointGeomInfo & gi);
    ///
    virtual void ProjectPoint (INDEX surfind, Point<3> & p) const;
    ///
    virtual void ProjectPoint2 (INDEX surfind, INDEX surfind2, Point<3> & p) const;
    ///
    virtual int CalcPointGeomInfo(PointGeomInfo& gi, const Point<3> & p3) const;
    ///
    void GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const override;
    void GetNormalVector(INDEX surfind, const Point<3>  & p, PointGeomInfo & gi, Vec<3> & n) const override;
};




class RefinementSTLGeometry : public Refinement
{
  const STLGeometry & geom;

public:
  RefinementSTLGeometry (const STLGeometry & ageom);
  virtual ~RefinementSTLGeometry ();
  
  virtual void PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
			     int surfi, 
			     const PointGeomInfo & gi1, 
			     const PointGeomInfo & gi2,
			     Point<3> & newp, PointGeomInfo & newgi) const override;

  virtual void PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
			     int surfi1, int surfi2, 
			     const EdgePointGeomInfo & ap1, 
			     const EdgePointGeomInfo & ap2,
			     Point<3> & newp, EdgePointGeomInfo & newgi) const override;

  virtual void ProjectToSurface (Point<3> & p, int surfi) const override;
  virtual void ProjectToSurface (Point<3> & p, int surfi, PointGeomInfo & gi) const override;
};



#endif


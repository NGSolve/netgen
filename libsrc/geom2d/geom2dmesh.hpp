#ifndef FILE_GEOM2DMESH
#define FILE_GEOM2DMESH

/**************************************************************************/
/* File:   geom2dmesh.hh                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   22. Jan. 01                                                    */
/**************************************************************************/


namespace netgen
{

  class Refinement2d : public Refinement
  {
    const class SplineGeometry2d & geometry;

  public:
    Refinement2d (const class SplineGeometry2d & ageometry);
    virtual ~Refinement2d ();
  
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


    virtual Vec<3> GetTangent (const Point<3> & p, int surfi1, int surfi2,
			       const EdgePointGeomInfo & ap1) const override;

    virtual Vec<3> GetNormal (const Point<3> & p, int surfi1, 
			      const PointGeomInfo & gi) const override;

    virtual void ProjectToSurface (Point<3> & p, int surfi, PointGeomInfo & /* gi */) const override;

    virtual void ProjectToEdge (Point<3> & p, int surfi1, int surfi2, 
				const EdgePointGeomInfo & egi) const override;
  };


}



#endif

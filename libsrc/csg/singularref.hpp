#ifndef FILE_SINGULARREF
#define FILE_SINGULARREF

/**************************************************************************/
/* File:   singularref.hh                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Sep. 99                                                    */
/**************************************************************************/

namespace netgen
{


  /**
     Control for local refinement
  */



  /**
     Singular Face.
     Causes a boundary layer mesh refinement.
     All elements in subdomain domnr will get a boundary layer
     on faces sharing the solid sol
  */
  class DLL_HEADER SingularFace 
  {
  public:
    int domnr;
    const Solid *sol;
    double factor; 
    // NgArray<Point<3> > points;
    // NgArray<INDEX_2> segms;
  public:
    SingularFace (int adomnr, const Solid * asol, double sf)
      : domnr(adomnr), sol(asol), factor(sf) { ; }
    const Solid * GetSolid() const { return sol; }
    int GetDomainNr () const { return domnr; }
  };


  ///
  class DLL_HEADER SingularEdge
  {
  public:
    double beta;
    int domnr;
    const CSGeometry& geom;
    const Solid *sol1, *sol2;
    NgArray<Point<3> > points;
    NgArray<INDEX_2> segms;
    double factor; 

    double maxhinit;
  public:
    SingularEdge (double abeta, int adomnr, 
		  const CSGeometry & ageom,
		  const Solid * asol1, const Solid * asol2, double sf,
		  const double maxh_at_initialization = -1);
    void FindPointsOnEdge (class Mesh & mesh);
    void SetMeshSize (class Mesh & mesh, double globalh);
  };


  ///
  class DLL_HEADER SingularPoint
  {
  public:
    double beta;
    const Solid *sol1, *sol2, *sol3;
    NgArray<Point<3> > points;
    double factor; 
 
  public:
    SingularPoint (double abeta, const Solid * asol1, const Solid * asol2,
		   const Solid * asol3, double sf);
    void FindPoints (class Mesh & mesh);
    void SetMeshSize (class Mesh & mesh, double globalh);
  };

}


#endif

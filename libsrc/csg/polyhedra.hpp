#ifndef FILE_POLYHEDRA
#define FILE_POLYHEDRA


/**************************************************************************/
/* File:   polyhedra.hh                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Mar. 2000                                                  */
/**************************************************************************/

namespace netgen
{

  /*

  Polyhedral primitive
  
  */

  class Polyhedra : public Primitive
  {
    class Face {
    public:
      int pnums[3];
      int planenr;

      int inputnr;

      Box<3> bbox;
      //    Point<3> center;
      Vec<3> v1, v2;   // edges
      Vec<3> w1, w2;   // pseudo-inverse
      Vec<3> n;        // normal to face
      Vec<3> nn;       // normed normal

      Face () { ; }
      Face (int pi1, int pi2, int pi3, 
	    const NgArray<Point<3> > & points,
	    int ainputnr);
    };

    NgArray<Point<3> > points;
    NgArray<Face> faces;
    NgArray<Plane*> planes;
    Box<3> poly_bbox;

    double eps_base1;

  public:
    Polyhedra ();
    virtual ~Polyhedra () override;
    static Primitive * CreateDefault ();

    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const override;
    virtual INSOLID_TYPE PointInSolid (const Point<3> & p,
				       double eps) const override;
    virtual INSOLID_TYPE VecInSolidNew (const Point<3> & p,
                                        const Vec<3> & v,
                                        double eps, bool printing = false) const;
    virtual INSOLID_TYPE VecInSolidOld (const Point<3> & p,
				     const Vec<3> & v,
				     double eps) const;
    
    virtual INSOLID_TYPE VecInSolid (const Point<3> & p,
				     const Vec<3> & v,
				     double eps) const override;

    virtual INSOLID_TYPE VecInSolid2 (const Point<3> & p,
				      const Vec<3> & v1,
				      const Vec<3> & v2,
				      double eps) const override;
    
    virtual INSOLID_TYPE VecInSolid3 (const Point<3> & p,
				      const Vec<3> & v1,
				      const Vec<3> & v2,
				      double eps) const override;

    virtual INSOLID_TYPE VecInSolid4 (const Point<3> & p,
				      const Vec<3> & v,
				      const Vec<3> & v2,
				      const Vec<3> & m,
				      double eps) const override;
    
    virtual void GetTangentialSurfaceIndices (const Point<3> & p, 
					      NgArray<int> & surfind, double eps) const override;


    virtual void GetTangentialVecSurfaceIndices2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
						  NgArray<int> & surfind, double eps) const override;

    virtual void CalcSpecialPoints (NgArray<Point<3> > & pts) const override;
    virtual void AnalyzeSpecialPoint (const Point<3> & pt, 
				      NgArray<Point<3> > & specpts) const override;
    virtual Vec<3> SpecialPointTangentialVector (const Point<3> & p, int s1, int s2) const override;

    virtual int GetNSurfaces() const override
    { return planes.Size(); }
    virtual Surface & GetSurface (int i) override
    { return *planes[i]; }
    virtual const Surface & GetSurface (int i) const override
    { return *planes[i]; }

    virtual void GetPrimitiveData (const char *& classname, NgArray<double> & coeffs) const override;
    virtual void SetPrimitiveData (NgArray<double> & coeffs) override;

    virtual void Reduce (const BoxSphere<3> & box) override;
    virtual void UnReduce () override;

    int AddPoint (const Point<3> & p);
    int AddFace (int pi1, int pi2, int pi3, int inputnum);

    void GetPolySurfs(NgArray < NgArray<int> * > & polysurfs);
  
  protected:
    int FaceBoxIntersection (int fnr, const BoxSphere<3> & box) const;
    //  void CalcData();
  };

}

#endif

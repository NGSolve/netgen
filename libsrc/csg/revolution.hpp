#ifndef _REVOLUTION_HPP
#define _REVOLUTION_HPP

namespace netgen
{

  class Revolution;

  class RevolutionFace : public Surface
  {
  private:
    bool isfirst, islast;
    const SplineSeg<2> * spline;
    bool deletable;

    Point<3> p0;
    Vec<3> v_axis;

    int id;

    // coefficient for implicizt polynomial
    mutable Vector spline_coefficient;
    mutable Vector spline_coefficient_shifted;


    NgArray < Vec<2>* > checklines_vec;
    NgArray < Point<2>* > checklines_start;
    NgArray < Vec<2>* > checklines_normal;
  
  private:
    void Init (void);

  public:
    void CalcProj(const Point<3> & point3d, Point<2> & point2d) const;
    void CalcProj(const Point<3> & point3d, Point<2> & point2d,
		  const Vec<3> & vector3d, Vec<2> & vector2d) const;
    void CalcProj0(const Vec<3> & point3d_minus_p0, Point<2> & point2d) const;

  public:
    RevolutionFace(const SplineSeg<2> & spline_in,
		   const Point<3> & p,
		   const Vec<3> & vec,
		   bool first = false,
		   bool last = false,
		   const int id_in = 0);

    RevolutionFace(const NgArray<double> & raw_data);
    // default constructor for archive
    RevolutionFace() {}

    ~RevolutionFace();

    virtual void DoArchive(Archive& ar)
    {
      Surface::DoArchive(ar);
      ar & isfirst & islast & spline & deletable & p0 & v_axis & id & spline_coefficient
        & spline_coefficient_shifted & checklines_vec & checklines_start & checklines_normal;
    }

    virtual int IsIdentic (const Surface & s2, int & inv, double eps) const; 
  
    virtual double CalcFunctionValue (const Point<3> & point) const;
    virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
    virtual void CalcHesse (const Point<3> & point, Mat<3> & hesse) const;
    virtual double HesseNorm () const;

    virtual double MaxCurvature () const;
    //virtual double MaxCurvatureLoc (const Point<3> & /* c */ , 
    //				  double /* rad */) const;
    
    Point<3> P0() const { return p0; }
    Vec<3> Axis() const { return v_axis; }
    
    virtual void Project (Point<3> & p) const;

    virtual Point<3> GetSurfacePoint () const;
    virtual void Print (ostream & str) const;
  
    virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					   const Box<3> & boundingbox, 
					   double facets) const;

    bool BoxIntersectsFace (const Box<3> & box) const;
    /*
      bool BoxIntersectsFace (const BoxSphere<2> & box, bool & uncertain) const;
      bool BoxIntersectsFace (const BoxSphere<3> & box, bool & uncertain) const;
    */

    const SplineSeg<2> & GetSpline(void) const {return *spline;}

    /* INSOLID_TYPE */ bool PointInFace (const Point<3> & p, const double eps) const;

    void GetRawData(NgArray<double> & data) const;

  };



  /*

  Primitive of revolution
  
  */


  class Revolution : public Primitive
  {
  private:
    Point<3> p0,p1;
    Vec<3> v_axis;

    // 1 ... torus-like
    // 2 ... sphere-like
    int type;
  

    Array<RevolutionFace*> faces;
    shared_ptr<SplineGeometry<2>> splinegeo;

    mutable int intersecting_face;

  public:
    Revolution(const Point<3> & p0_in,
	       const Point<3> & p1_in,
	       shared_ptr<SplineGeometry<2>> spline_in);
    // default constructor for archive
    Revolution() {}

    ~Revolution();

    virtual void DoArchive(Archive& ar)
    {
      Primitive::DoArchive(ar);
      ar & p0 & p1 & v_axis & type & faces & intersecting_face;
    }
  
    /*
      Check, whether box intersects solid defined by surface.

      return values:
      0 .. box outside solid \\
      1 .. box in solid \\
      2 .. can't decide (allowed, iff box is close to solid)
    */
    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
    virtual INSOLID_TYPE PointInSolid (const Point<3> & p,
				       double eps) const;

    virtual void GetTangentialSurfaceIndices (const Point<3> & p, 
					      NgArray<int> & surfind, double eps) const;
    
    virtual INSOLID_TYPE VecInSolid (const Point<3> & p,
				     const Vec<3> & v,
				     double eps) const;

    // checks if lim s->0 lim t->0  p + t(v1 + s v2) in solid
    virtual INSOLID_TYPE VecInSolid2 (const Point<3> & p,
				      const Vec<3> & v1,
				      const Vec<3> & v2,
				      double eps) const;

    virtual void GetTangentialVecSurfaceIndices2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
						  NgArray<int> & surfind, double eps) const;

    
    virtual int GetNSurfaces() const;
    virtual Surface & GetSurface (int i = 0);
    virtual const Surface & GetSurface (int i = 0) const;


    virtual void Reduce (const BoxSphere<3> & box);
    virtual void UnReduce ();
  
	     
  };

}


#endif

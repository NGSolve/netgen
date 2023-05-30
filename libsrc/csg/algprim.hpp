#ifndef FILE_ALGPRIM
#define FILE_ALGPRIM


/**************************************************************************/
/* File:   algprim.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   1. Dez. 95                                                     */
/**************************************************************************/

namespace netgen
{

  /*

  Quadric Surfaces (Plane, Sphere, Cylinder)
  
  */


  /**
     A quadric surface.
     surface defined by
     cxx x^2 + cyy y^2 + czz z^2 + cxy x y + cxz x z + cyz y z +
     cx x + cy y + cz z + c1 = 0.
  **/
  class QuadraticSurface : public  OneSurfacePrimitive
  {
  protected:
    double cxx, cyy, czz, cxy, cxz, cyz, cx, cy, cz, c1;

  public:
    virtual double CalcFunctionValue (const Point<3> & point) const;
    virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
    virtual void CalcHesse (const Point<3> & point, Mat<3> & hesse) const;
    /*
      virtual int RootInBox (const Box<3> & box) 
      const { return 0; }
      virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) 
      const { return DOES_INTERSECT; }
    */
    virtual double HesseNorm () const { return cxx + cyy + czz; }

    virtual Point<3> GetSurfacePoint () const;


    virtual void Print (ostream & str) const;
    virtual void Read (istream & ist);
    void PrintCoeff (ostream & ost) const;
    virtual void DoArchive(Archive& ar)
    {
      OneSurfacePrimitive::DoArchive(ar);
      ar & cxx & cyy & czz & cxy & cxz & cyz & cx & cy & cz & c1;
    }
  };


  /// A Plane (i.e., the plane and everything behind it).
  class Plane : public QuadraticSurface
  {
  protected:
    /// a point in the plane
    Point<3> p;
    /// outward normal vector 
    Vec<3> n;

    double eps_base;

  public:
    ///
    Plane (const Point<3> & ap, Vec<3> an);
    // default constructor for archive
    Plane() {}

    virtual void DoArchive(Archive& ar)
    {
      QuadraticSurface::DoArchive(ar);
      ar & p & n & eps_base;
    }
    Point<3> P() const { return p; }
    Vec<3> N() const { return n; }
    virtual void GetPrimitiveData (const char *& classname, 
				   NgArray<double> & coeffs) const;
    virtual void SetPrimitiveData (NgArray<double> & coeffs);
    static Primitive * CreateDefault ();

    virtual Primitive * Copy () const;
    virtual void Print (ostream & str) const;

    virtual void Transform (Transformation<3> & trans);


    virtual int IsIdentic (const Surface & s2, int & inv, double eps) const;

    ///
    virtual void DefineTangentialPlane (const Point<3> & ap1, 
					const Point<3> & ap2);
    ///
    virtual void ToPlane (const Point<3> & p3d, 
			  Point<2> & pplane, double h,
			  int & zone) const;
    ///
    virtual void FromPlane (const Point<2> & pplane, 
			    Point<3> & p3d, 
			    double h) const;
    ///
    virtual void Project (Point<3> & p) const;

    ///
    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;

    ///
    inline virtual double CalcFunctionValue (const Point<3> & p3d) const
    {return cx * p3d(0) + cy * p3d(1) + cz * p3d(2) + c1;}
    ///
    virtual void CalcGradient (const Point<3> & point, 
			       Vec<3> & grad) const;
    ///
    virtual void CalcHesse (const Point<3> & point, 
			    Mat<3> & hesse) const;
    ///
    virtual double HesseNorm () const;
    ///
    virtual Point<3> GetSurfacePoint () const;
    ///
    virtual void GetTriangleApproximation 
    (TriangleApproximation & tas, 
     const Box<3> & boundingbox, double facets) const;
  protected:
    void CalcData();
  };

  // typedef Plane Plane;


  ///
  class Sphere : public QuadraticSurface
  {
    ///
    Point<3> c;
    ///
    double r, invr;
  public:
    ///
    Sphere (const Point<3> & ac, double ar);
    // default constructor for archive
    Sphere() {}

    virtual void DoArchive(Archive& ar)
    {
      QuadraticSurface::DoArchive(ar);
      ar & c & r & invr;
    }

    virtual void GetPrimitiveData (const char *& classname, 
				   NgArray<double> & coeffs) const;
    virtual void SetPrimitiveData (NgArray<double> & coeffs);
    static Primitive * CreateDefault ();

    virtual Primitive * Copy () const;
    virtual void Transform (Transformation<3> & trans);

    virtual double CalcFunctionValue (const Point<3> & point) const;


    virtual int IsIdentic (const Surface & s2, int & inv, double eps) const;

    ///
    virtual void DefineTangentialPlane (const Point<3> & ap1, 
					const Point<3> & ap2);
    ///
    virtual void ToPlane (const Point<3> & p3d, 
			  Point<2> & pplane, double h,
			  int & zone) const;
    ///
    virtual void FromPlane (const Point<2> & pplane, 
			    Point<3> & p, double h) const;
    ///
    virtual void Project (Point<3> & p) const;

    ///
    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
    ///
    virtual double HesseNorm () const;
    ///
    virtual Point<3> GetSurfacePoint () const;
    ///
    const Point<3> & Center () const { return c; }
    ///
    double Radius () const { return r; }

    ///
    virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					   const Box<3> & bbox, 
					   double facets) const;
  };


  ///
  class Cylinder : public QuadraticSurface
  {
    ///
    Point<3> a, b;
    ///
    double r;
    ///
    Vec<3> vab;

  public:
    Cylinder (const Point<3> & aa, const Point<3> & ab, double ar);
    Cylinder (NgArray<double> & coeffs);
    // default constructor for archive
    Cylinder() {}

    virtual void DoArchive(Archive& ar)
    {
      QuadraticSurface::DoArchive(ar);
      ar & a & b & r & vab;
    }
    Point<3> A() const { return a; }
    Point<3> B() const { return b; }
    double R() const { return r; }
    virtual void GetPrimitiveData (const char *& classname, NgArray<double> & coeffs) const;
    virtual void SetPrimitiveData (NgArray<double> & coeffs);
    static Primitive * CreateDefault ();

    virtual Primitive * Copy () const;
    virtual void Print (ostream & str) const;

    virtual void Transform (Transformation<3> & trans);

    ///
    virtual int IsIdentic (const Surface & s2, int & inv, double eps) const;
    ///
    virtual void DefineTangentialPlane (const Point<3> & ap1, 
					const Point<3> & ap2);
    ///
    virtual void ToPlane (const Point<3> & p, 
			  Point<2> & pplane, 
			  double h,
			  int & zone) const;
    ///
    virtual void FromPlane (const Point<2> & pplane, 
			    Point<3> & p, 
			    double h) const;
    ///
    virtual void Project (Point<3> & p) const;

    ///
    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
    ///
    virtual double HesseNorm () const;
    ///
    virtual Point<3> GetSurfacePoint () const;
    ///
    virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					   const Box<3> & bbox, 
					   double facets) const;
  };





  ///
  class EllipticCylinder : public QuadraticSurface
  {
  private:
    ///
    Point<3> a;
    ///
    Vec<3> vl, vs;
    ///
    Vec<3> vab, t0vec, t1vec;
    ///
    double vabl, t0, t1;
  public:
    ///
    EllipticCylinder (const Point<3> & aa,
		      const Vec<3> & avl, const Vec<3> & avs);
    EllipticCylinder (NgArray<double> & coeffs);
    // default constructor for archive
    EllipticCylinder() {}

    virtual void DoArchive(Archive& ar)
    {
      QuadraticSurface::DoArchive(ar);
      ar & a & vl & vs & vab & t0vec & t1vec & vabl & t0 & t1;
    }

    // static Primitive * CreateDefault ();
    virtual void GetPrimitiveData (const char *& classname, NgArray<double> & coeffs) const;
    virtual void SetPrimitiveData (NgArray<double> & coeffs);

    ///
    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
    ///
    virtual double HesseNorm () const;
    ///
    virtual Point<3> GetSurfacePoint () const;

    virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					   const Box<3> & bbox, 
					   double facets) const;

    virtual int IsIdentic (const Surface & s2, int & inv, double eps) const;

    virtual double MaxCurvature () const;

    virtual double MaxCurvatureLoc (const Point<3> & /* c */ , 
				    double /* rad */) const;


  private:
    void CalcData();
  };






  ///
  class Ellipsoid : public QuadraticSurface
  {
  private:
    ///
    Point<3> a;
    ///
    Vec<3> v1, v2, v3;
    ///
    double rmin;
  public:
    ///
    Ellipsoid (const Point<3> & aa,
	       const Vec<3> & av1, 
	       const Vec<3> & av2,
	       const Vec<3> & av3);
    // default constructor for archive
    Ellipsoid() {}

    void DoArchive(Archive& ar) override
    {
      QuadraticSurface::DoArchive(ar);
      ar & a & v1 & v2 & v3 & rmin;
    }
    ///
    INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const override;
    ///
    double HesseNorm () const override;
    ///
    double MaxCurvature () const override;
    ///
    Point<3> GetSurfacePoint () const override;

    void GetTriangleApproximation (TriangleApproximation & tas,
					   const Box<3> & bbox, 
					   double facets) const override;

    void GetPrimitiveData (const char *& classname, NgArray<double> & coeffs) const override;
    void SetPrimitiveData (NgArray<double> & coeffs) override;

  private:
    void CalcData();
  };








  ///
  class Cone : public QuadraticSurface
  {
    ///
    Point<3> a, b;
    ///
    double ra, rb, minr;
    ///
    Vec<3> vab, t0vec, t1vec;
    ///
    double vabl, t0, t1;
    double cosphi;
  public:
    ///
    Cone (const Point<3> & aa, const Point<3> & ab, double ara, double arb);
    ///
    // default constructor for archive
    Cone() {}

    virtual void DoArchive(Archive& ar)
    {
      QuadraticSurface::DoArchive(ar);
      ar & a & b & ra & rb & minr & vab & t0vec & t1vec & vabl & t0 & t1 & cosphi;
    }
    static Primitive * CreateDefault ();
    virtual void GetPrimitiveData (const char *& classname, NgArray<double> & coeffs) const;
    virtual void SetPrimitiveData (NgArray<double> & coeffs);

    ///
    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
    ///
    virtual double HesseNorm () const;

    virtual double LocH (const Point<3> & p, double x, 
			 double c,                          
                         const MeshingParameters & mparam,
                         double hmax) const;

    ///
    virtual Point<3> GetSurfacePoint () const;

    virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					   const Box<3> & bbox, 
					   double facets) const;

  private:
    void CalcData();
  };


///

/// Elliptic Cone
/// Josephat Kalezhi (kalezhi@cbu.ac.zm)
/// February 21st, 2018
///

  ///
  class EllipticCone : public QuadraticSurface
  {
	Point<3> a;
	Vec<3> vl, vs;
	double h, vlr;

  public:
  ///
  EllipticCone (const Point<3> & aa, const Vec<3> & avl,
	      const Vec<3> & avs, double ah, double avlr);
    // default constructor for archive
    EllipticCone() {}

    virtual void DoArchive(Archive& ar)
    {
      QuadraticSurface::DoArchive(ar);
      ar & a & vl & vs & h & vlr;
    }
  static Primitive * CreateDefault ();
  virtual void GetPrimitiveData (const char *& classname, NgArray<double> & coeffs) const;
  virtual void SetPrimitiveData (NgArray<double> & coeffs);

  ///
  virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
  ///
  virtual double HesseNorm () const;
  virtual double MaxCurvature () const;
  virtual double MaxCurvatureLoc (const Point<3> & /* c */ ,
				  double /* rad */) const;
  ///
  virtual Point<3> GetSurfacePoint () const;

  virtual void GetTriangleApproximation (TriangleApproximation & tas,
					 const Box<3> & bbox,
					 double facets) const;

  private:
  void CalcData();
  };


  /** Torus 
  /// Lorenzo Codecasa (codecasa@elet.polimi.it)
  /// April 27th, 2005 
  */
  class Torus : public OneSurfacePrimitive
  { 
    /// center of the torus
    Point<3> c;
    /// vector normal to the symmetry plane of the torus
    Vec<3> n;
    /// Large radius of the torus
    double R;
    /// Small radius of the torus
    double r;
  
  public:
    /// OK
    Torus (const Point<3> & ac, const Vec<3> & an, double aR, double ar);
    // default constructor for archive
    Torus() {}

    virtual void DoArchive(Archive& ar)
    {
      OneSurfacePrimitive::DoArchive(ar);
      ar & c & n & R & r;
    }
    /// OK
    const Point<3> & Center () const { return c; }
    /// OK
    const Vec<3> & NormalToPlane () const { return n; }
    /// OK
    double LargeRadius () const { return R; }
    /// OK
    double SmallRadius () const { return r; }
    /// OK
    virtual double CalcFunctionValue (const Point<3> & point) const;
    /// OK
    virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
    /// OK
    virtual void CalcHesse (const Point<3> & point, Mat<3> & hesse) const;
    /// OK
    virtual double HesseNorm () const;
    /// OK
    virtual Point<3> GetSurfacePoint () const;
    /// OK
    virtual void GetPrimitiveData (const char *& classname, 
				   NgArray<double> & coeffs) const;
    /// OK			 
    virtual void SetPrimitiveData (NgArray<double> & coeffs);
    /// OK
    static Primitive * CreateDefault ();
    /// OK
    virtual Primitive * Copy () const;
    /// OK
    virtual void Transform (Transformation<3> & trans);
    /// OK
    virtual int IsIdentic (const Surface & s2, int & inv, double eps) const;
    /// OK
    /// virtual void DefineTangentialPlane (const Point<3> & ap1, 
    //				      const Point<3> & ap2);
    /// OK
    /// virtual void ToPlane (const Point<3> & p3d, 
    ///			Point<2> & pplane, 
    ///			double h, int & zone) const;
    /// OK
    /// virtual void FromPlane (const Point<2> & pplane, 
    //			  Point<3> & p, double h) const;
    /// OK
    /// virtual void Project (Point<3> & p) const;
    /// OK
    virtual INSOLID_TYPE BoxInSolid (const BoxSphere<3> & box) const;
    /// OK
    virtual void GetTriangleApproximation (TriangleApproximation & tas, 
					   const Box<3> & bbox, 
					   double facets) const;
    /// OK		 
    virtual void Print (ostream & ist) const;
    /// OK
    virtual void Read (istream & ist);
  };

  /// ...end


}





#endif

#ifndef FILE_GEOM2D
#define FILE_GEOM2D

/* *************************************************************************/
/* File:   geom2d.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   5. Aug. 95                                                      */
/* *************************************************************************/

#include <mydefs.hpp>

#include <general/template.hpp>
#include "geomobjects.hpp"
#include <meshing/global.hpp>

namespace netgen 
{

  /* Geometric Algorithms */

#define EPSGEOM 1E-5




  class LINE2D;
  class Line2d;
  class PLine2d;
  class TRIANGLE2D;
  class PTRIANGLE2D;


  /// 2d cross product
  inline double Cross (const Vec<2> & v1, const Vec<2> & v2)
  { return v1(0) * v2(1) - v1(1) * v2(0); }

  DLL_HEADER double Angle (const Vec<2> & v);
  DLL_HEADER double FastAngle (const Vec<2> & v);
  DLL_HEADER double Angle (const Vec<2> & v1, const Vec<2> & v2);
  DLL_HEADER double FastAngle (const Vec<2> & v1, const Vec<2> & v2);
  double Dist2(const Line2d & g, const Line2d & h );            // GH
  int Near (const Point<2> & p1, const Point<2> & p2, const double eps);

  int Parallel (const Line2d & l1, const Line2d & l2, double peps = EPSGEOM);
  int IsOnLine (const Line2d & l, const Point<2> & p, double heps = EPSGEOM);
  int IsOnLongLine (const Line2d & l, const Point<2> & p);
  int Hit (const Line2d & l1, const Line2d & l2, double heps = EPSGEOM);
  ostream & operator<<(ostream  & s, const Line2d & l);
  DLL_HEADER Point<2> CrossPoint (const PLine2d & l1, const PLine2d & l2);
  DLL_HEADER Point<2> CrossPoint (const Line2d & l1, const Line2d & l2);
  int Parallel (const PLine2d & l1, const PLine2d & l2, double peps = EPSGEOM);
  int IsOnLine (const PLine2d & l, const Point<2> & p, double heps = EPSGEOM);
  int IsOnLongLine (const PLine2d & l, const Point<2> & p);
  int Hit (const PLine2d & l1, const Line2d & l2, double heps = EPSGEOM);
  ostream & operator<<(ostream  & s, const Line2d & l);
  ostream & operator<<(ostream  & s, const TRIANGLE2D & t); 
  ostream & operator<<(ostream & s, const PTRIANGLE2D & t);




  ///
  class Line2d
  {
  protected:
    ///
    Point<2> p1, p2;

  public:
    ///
    Line2d() : p1(), p2() { };
    ///
    Line2d(const Point<2> & ap1, const Point<2> & ap2)
    { p1 = ap1; p2 = ap2; }

    ///
    Line2d & operator= (const Line2d & l2)
    { p1 = l2.p1; p2 = l2.p2; return *this;}

    ///
    Point<2> & P1() { return p1; }
    ///
    Point<2> & P2() { return p2; }
    ///
    const Point<2> & P1() const { return p1; }
    ///
    const Point<2> & P2() const { return p2; }

    ///
    double XMax() const { return max2 (p1(0), p2(0)); }
    ///
    double YMax() const { return max2 (p1(1), p2(1)); }
    ///
    double XMin() const { return min2 (p1(0), p2(0)); }
    ///
    double YMin() const { return min2 (p1(1), p2(1)); }

    ///
    Vec<2> Delta () const { return Vec<2> (p2(0)-p1(0), p2(1)-p1(1)); }
    ///
    double Length () const { return Delta().Length(); }
    ///
    double Length2 () const
    { return sqr (p1(0) - p2(0)) +
        sqr (p1(1) - p2(1)); }

    void GetNormal (Line2d & n) const;                                  // GH
    Vec<2> NormalDelta () const;                                                // GH

    /// square of the distance between two 2d-lines.
    friend double Dist2(const Line2d & g, const Line2d & h );           // GH

    ///
    friend DLL_HEADER Point<2> CrossPoint (const Line2d & l1, const Line2d & l2);
    /// returns 1 iff parallel
    friend int CrossPointBarycentric (const Line2d & l1, const Line2d & l2,
                                      double & lam1, double & lam2, double eps);
    
    ///
    friend int Parallel (const Line2d & l1, const Line2d & l2, double peps);
    ///
    friend int IsOnLine (const Line2d & l, const Point<2> & p, double heps);
    ///
    friend int IsOnLongLine (const Line2d & l, const Point<2> & p);
    ///
    friend int Hit (const Line2d & l1, const Line2d & l2, double heps);

    ///
    friend ostream & operator<<(ostream  & s, const Line2d & l);
  };


#ifdef NONE
  ///
  class PLine2d
  {
  protected:
    ///
    Point<2> const * p1, *p2;

  public:
    ///
    PLine2d() { };
    ///
    PLine2d(Point<2> const * ap1, Point<2> const * ap2)
    { p1 = ap1; p2 = ap2; }

    ///
    PLine2d & operator= (const PLine2d & l2)
    { p1 = l2.p1; p2 = l2.p2; return *this;}

    ///
    const Point<2> *& P1() { return p1; }
    ///
    const Point<2> *& P2() { return p2; }
    ///
    const Point<2> & P1() const { return *p1; }
    ///
    const Point<2> & P2() const { return *p2; }

    ///
    double XMax() const { return max2 (p1->X(), p2->X()); }
    ///
    double YMax() const { return max2 (p1->Y(), p2->Y()); }
    ///
    double XMin() const { return min2 (p1->X(), p2->X()); }
    ///
    double YMin() const { return min2 (p1->Y(), p2->Y()); }


    ///
    Vec<2> Delta () const { return Vec<2> (p2->X()-p1->X(), p2->Y()-p1->Y()); }
    ///
    double Length () const { return Delta().Length(); }
    ///
    double Length2 () const
    { return sqr (p1->X() - p2->X()) +
        sqr (p1->Y() - p2->Y()); }


    
    ///
    friend Point<2> CrossPoint (const PLine2d & l1, const PLine2d & l2);
    ///
    friend int Parallel (const PLine2d & l1, const PLine2d & l2, double peps);
    ///
    friend int IsOnLine (const PLine2d & l, const Point<2> & p, double heps);
    ///
    friend int IsOnLongLine (const PLine2d & l, const Point<2> & p);
    ///
    friend int Hit (const PLine2d & l1, const Line2d & l2, double heps);

    ///
    friend ostream & operator<<(ostream  & s, const Line2d & l);
  };



  ///
  class ILINE
  {
    ///
    int i[2];

  public:
    ///
    ILINE() {};
    ///
    ILINE(int i1, int i2) { i[0] = i1; i[1] = i2; }
    ///
    ILINE(const ILINE & l) { i[0] = l.i[0]; i[1] = l.i[1]; }

    ///
    ILINE & operator= (const ILINE & l)
    { i[0] = l.i[0]; i[1] = l.i[1]; return *this; }

    ///
    const int & I(int ai) const { return i[ai-1]; }
    ///
    const int & X() const { return i[0]; }
    ///
    const int & Y() const { return i[1]; }
    ///
    const int & I1() const { return i[0]; }
    ///
    const int & I2() const { return i[1]; }

    ///
    int & I(int ai) { return i[ai-1]; }
    ///
    int & X() { return i[0]; }
    ///
    int & Y() { return i[1]; }
    ///
    int & I1() { return i[0]; }
    ///
    int & I2() { return i[1]; }
  };




  ///
  class TRIANGLE2D
  {
  private:
    ///
    Point<2> p1, p2, p3;

  public:
    ///
    TRIANGLE2D() { };
    ///
    TRIANGLE2D (const Point<2> & ap1, const Point<2> & ap2,
                const Point<2> & ap3)
    { p1 = ap1; p2 = ap2; p3 = ap3;}

    ///
    TRIANGLE2D & operator= (const TRIANGLE2D & t2)
    { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }

    ///
    Point<2> & P1() { return p1; }
    ///
    Point<2> & P2() { return p2; }
    ///
    Point<2> & P3() { return p3; }
    ///
    const Point<2> & P1() const { return p1; }
    ///
    const Point<2> & P2() const { return p2; }
    ///
    const Point<2> & P3() const { return p3; }

    ///
    double XMax() const { return max3 (p1.X(), p2.X(), p3.X()); }
    ///
    double YMax() const { return max3 (p1.Y(), p2.Y(), p3.Y()); }
    ///
    double XMin() const { return min3 (p1.X(), p2.X(), p3.X()); }
    ///
    double YMin() const { return min3 (p1.Y(), p2.Y(), p3.Y()); }

    ///
    inline Point<2> Center () const
    { return Point<2>( (p1.X()+p2.X()+p3.X())/3, (p1.Y()+p2.Y()+p3.Y())/3); }

    ///
    int Regular() const;
    /// 
    int CW () const;
    ///
    int CCW () const;

    ///
    int IsOn (const Point<2> & p) const;
    ///
    int IsIn (const Point<2> & p) const;
    ///
    friend ostream & operator<<(ostream  & s, const TRIANGLE2D & t);
  };


  ///
  class PTRIANGLE2D
  {
  private:
    ///
    Point<2> const *p1, *p2, *p3;

  public:
    ///
    PTRIANGLE2D() { };
    ///
    PTRIANGLE2D (const Point<2> * ap1, const Point<2> * ap2,
                 const Point<2> * ap3)
    { p1 = ap1; p2 = ap2; p3 = ap3;}

    ///
    PTRIANGLE2D & operator= (const PTRIANGLE2D & t2)
    { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }

    ///
    const Point<2> *& P1() { return p1; }
    ///
    const Point<2> *& P2() { return p2; }
    ///
    const Point<2> *& P3() { return p3; }
    ///
    const Point<2> * P1() const { return p1; }
    ///
    const Point<2> * P2() const { return p2; }
    ///
    const Point<2> * P3() const { return p3; }

    ///
    double XMax() const { return max3 (p1->X(), p2->X(), p3->X()); }
    ///
    double YMax() const { return max3 (p1->Y(), p2->Y(), p3->Y()); }
    ///
    double XMin() const { return min3 (p1->X(), p2->X(), p3->X()); }
    ///
    double YMin() const { return min3 (p1->Y(), p2->Y(), p3->Y()); }

    ///
    Point<2> Center () const
    { return Point<2>( (p1->X()+p2->X()+p3->X())/3, (p1->Y()+p2->Y()+p3->Y())/3); }


    ///
    int Regular() const;
    ///
    int CW () const;
    ///
    int CCW () const;

    ///
    int IsOn (const Point<2> & p) const;
    ///
    int IsIn (const Point<2> & p) const;
    ///
    friend ostream & operator<<(ostream & s, const PTRIANGLE2D & t);
  };
#endif


  /** Cheap approximation to atan2.
      A monotone function of atan2(x,y) is computed.
  */
  extern double Fastatan2 (double x, double y);







#ifdef none
  inline int TRIANGLE2D :: Regular() const
  {
    return fabs(Cross ( p2 - p1, p3 - p2)) > EPSGEOM;
  }


  inline int TRIANGLE2D :: CW () const
  {
    return Cross ( p2 - p1, p3 - p2) < 0;
  }


  inline int TRIANGLE2D :: CCW () const
  {
    return Cross ( p2 - p1, p3 - p2) > 0;
  }




  inline int PTRIANGLE2D :: Regular() const
  {
    return fabs(Cross ( *p2 - *p1, *p3 - *p2)) > EPSGEOM;
  }


  inline int PTRIANGLE2D :: CW () const
  {
    return Cross ( *p2 - *p1, *p3 - *p2) < 0;
  }


  inline int PTRIANGLE2D :: CCW () const
  {
    return Cross ( *p2 - *p1, *p3 - *p2) > 0;
  }


#endif


  ///
  class Mat2d
  {
  protected:
    ///
    double coeff[4];

  public:
    ///
    Mat2d() { coeff[0] = coeff[1] = coeff[2] = coeff[3] = 0; }
    ///
    Mat2d(double a11, double a12, double a21, double a22)
    { coeff[0] = a11; coeff[1] = a12; coeff[2] = a21; coeff[3] = a22; }
    ///
    Mat2d(const Mat2d & m2)
    { for (int i = 0; i < 4; i++) coeff[i] = m2.Get(i); }

    ///
    double & Elem (int i, int j) { return coeff[2*(i-1)+j-1]; }
    ///
    double & Elem (int i) {return coeff[i]; }
    ///
    double Get (int i, int j) const { return coeff[2*(i-1)+j-1]; }
    ///
    double Get (int i) const {return coeff[i]; }

    ///  
    double Det () const { return coeff[0] * coeff[3] - coeff[1] * coeff[2]; }

    ///
    void Mult (const Vec<2> & v, Vec<2> & prod) const;
    ///
    void MultTrans (const Vec<2> & v , Vec<2> & prod) const;
    ///
    void Solve (const Vec<2> & rhs, Vec<2> & x) const;
    /// Solves mat * x = rhs, but using a positive definite matrix instead of mat
    void SolvePositiveDefinite (const Vec<2> & rhs, Vec<2> & x) const;
    /// add a term \alpha * v * v^T
    void AddDiadicProduct (double alpha, Vec<2> & v);
  };



  inline void Mat2d :: Mult (const Vec<2> & v, Vec<2> & prod) const
  {
    prod(0) = coeff[0] * v(0) + coeff[1] * v(1);
    prod(1) = coeff[2] * v(0) + coeff[3] * v(1);
  }


  inline  void Mat2d :: MultTrans (const Vec<2> & v, Vec<2> & prod) const
  {
    prod(0) = coeff[0] * v(0) + coeff[2] * v(1);
    prod(1) = coeff[1] * v(0) + coeff[3] * v(1);
  }



  inline void Mat2d :: Solve (const Vec<2> & rhs, Vec<2> & x) const
  {
    double det = Det();
  
    if (det == 0)
      throw Exception ("Mat2d::Solve: zero determinant");
    x(0) = (coeff[3] * rhs(0) - coeff[1] * rhs(1)) / det;
    x(1) = (-coeff[2] * rhs(0) + coeff[0] * rhs(1)) / det;
  }


  inline void Mat2d :: SolvePositiveDefinite (const Vec<2> & rhs, Vec<2> & x) const
  {
    double a = max2(coeff[0], 1e-8);
    double b = coeff[1] / a;
    double c = coeff[2] / a;
    double d = max2(coeff[3] - a *b * c, 1e-8);

    x(0) = (rhs(0) - b * rhs(1)) / a;
    x(1) = rhs(1) / d - c * x(0);
  }


  inline void Mat2d :: AddDiadicProduct (double alpha, Vec<2> & v)
  {
    coeff[0] += alpha * v(0) * v(0);
    coeff[1] += alpha * v(0) * v(1);
    coeff[2] += alpha * v(1) * v(0);
    coeff[3] += alpha * v(1) * v(1);
  }

}

#endif

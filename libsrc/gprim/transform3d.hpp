#ifndef FILE_TRANSFORM3D
#define FILE_TRANSFORM3D

/* *************************************************************************/
/* File:   transform3d.hh                                                  */
/* Author: Joachim Schoeberl                                               */
/* Date:   22. Mar. 98                                                     */
/* *************************************************************************/

/*
  Affine - Linear mapping in 3D space
 */

#include "geom3d.hpp"
#include "geomfuncs.hpp"

namespace netgen
{

class Transformation3d;
ostream & operator<< (ostream & ost, Transformation3d & trans);

class Transformation3d
{
  double lin[3][3];
  double offset[3];
public:
  ///
  Transformation3d ();
  /// Unit tet is mapped to tet described by pp
  Transformation3d (const Point<3> ** pp);
  /// Unit tet is mapped to tet described by pp
  Transformation3d (const Point<3> pp[]);
  /// translation
  Transformation3d (const Vec<3> & translate);
  /// rotation with ...
  Transformation3d (const Point<3> & c, double alpha, double beta, double gamma);
  /// 
  void CalcInverse (Transformation3d & inv) const;
  /// this = ta x tb
  void Combine (const Transformation3d & ta, const Transformation3d & tb);
  /// dir = 1..3 (== x..z)
  void SetAxisRotation (int dir, double alpha);
  ///
  void Transform (const Point<3> & from, Point<3> & to) const
    {
      for (int i = 1; i <= 3; i++)
	{
	  to(i-1) = offset[i-1] + lin[i-1][0] * from(0) + 
	    lin[i-1][1] * from(1) + lin[i-1][2] * from(2);
	}
    }

  ///
  void Transform (Point<3> & p) const
  {
    Point<3> hp;
    Transform (p, hp);
    p = hp;
  }

  /// transform vector, apply only linear part, not offset
  void Transform (const Vec<3> & from, Vec<3> & to) const
    {
      for (int i = 1; i <= 3; i++)
	{
	  to(i-1) = lin[i-1][0] * from(0) + 
	    lin[i-1][1] * from(1) + lin[i-1][2] * from(2);
	}
    }
  friend ostream & operator<< (ostream & ost, Transformation3d & trans);
};














template <int D>
class Transformation
{
  Mat<D> m;
  Vec<D> v;
public:
  ///
  Transformation () { m = 0; v = 0; }

  /// Unit tet is mapped to tet described by pp
  Transformation (const Point<D> * pp);

  /// translation
  Transformation (const Vec<D> & translate)
  {
    v = translate;
    m = 0;
    for (int i = 0; i < D; i++)
      m(i,i) = 1;
  }

  Transformation (const Point<D> & c, const Vec<3> & axes, double angle);
  
  // rotation with ...
  Transformation (const Point<D> & c, double alpha, double beta, double gamma)
  {
    // total = T_c x Rot_0 x T_c^{-1}
    // Use Euler angles, see many books from tech mech, e.g. 
    // Shabana "multibody systems"
    
    Vec<D> vc(c);
    Transformation<D> tc(vc);
    Transformation<D> tcinv(-vc);
    // tc.CalcInverse (tcinv);
    
    Transformation<D> r1, r2, r3, ht, ht2;
    r1.SetAxisRotation (3, alpha);
    r2.SetAxisRotation (1, beta);
    r3.SetAxisRotation (3, gamma);
    
    ht.Combine (tc, r3);
    ht2.Combine (ht, r2);
    ht.Combine (ht2, r1);
    Combine (ht, tcinv);
    
    // cout << "Rotation - Transformation:" << (*this) << endl;
    //  (*testout) << "Rotation - Transformation:" << (*this) << endl;
  }

  Mat<D> & GetMatrix() { return m; }
  Vec<D> & GetVector() { return v; }

  void DoArchive(Archive& ar)
  {
      ar & m & v;
  }

  /// 
  Transformation CalcInverse () const
  {
    Transformation inv;
    // inv.m = Inv(m);
    ::netgen::CalcInverse (m, inv.m);
    inv.v = inv.m * (-v);
    return inv;
  }

  /// this = ta x tb
  void Combine (const Transformation & ta, const Transformation & tb)
  {
    v = ta.v + ta.m * tb.v;
    m = ta.m * tb.m;
  }



  /// dir = 1..3 (== x..z)
  void SetAxisRotation (int dir, double alpha)
  {
    double co = cos(alpha);
    double si = sin(alpha);
    dir--;
    int pos1 = (dir+1) % 3;
    int pos2 = (dir+2) % 3;
    
    int i, j;
    for (i = 0; i <= 2; i++)
    {
      v(i) = 0;
      for (j = 0; j <= 2; j++)
	m(i,j) = 0;
    }
    
    m(dir,dir) = 1;
    m(pos1, pos1) = co;
    m(pos2, pos2) = co;
    m(pos1, pos2) = si;
    m(pos2, pos1) = -si;
  }

  ///
  void Transform (const Point<D> & from, Point<D> & to) const
  {
    to = Point<D> (v + m * Vec<D>(from));
  }

  void Transform (Point<D> & p) const
  {
    p = Point<D> (v + m * Vec<D>(p));
  }



  /// transform vector, apply only linear part, not offset
  void Transform (const Vec<D> & from, Vec<D> & to) const
  {
    to = m * from;
  }

  Point<D> operator() (Point<D> from) const { Point<D> to; Transform(from, to); return to; }
  Vec<D> operator() (Vec<D> from) const { Vec<D> to; Transform(from, to); return to; }
};

template <int D>
ostream & operator<< (ostream & ost, Transformation<D> & trans);


}

#endif

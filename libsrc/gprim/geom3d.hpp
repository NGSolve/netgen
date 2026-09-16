#ifndef FILE_GEOM3D
#define FILE_GEOM3D

/* *************************************************************************/
/* File:   geom3d.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   5. Aug. 95                                                      */
/* *************************************************************************/

#include <mydefs.hpp>
#include "geom2d.hpp"

namespace netgen
{



  /// a unit vector normal to v
  DLL_HEADER void GetNormal (const Vec<3> & v, Vec<3> & n);
  inline Vec<3> GetNormal (const Vec<3> & v) { Vec<3> n; GetNormal (v, n); return n; }

  double Angle (const Vec<3> & v);
  double FastAngle (const Vec<3> & v);
  double Angle (const Vec<3> & v1, const Vec<3> & v2);
  double FastAngle (const Vec<3> & v1, const Vec<3> & v2);
  ostream & operator<<(ostream  & s, const Vec<3> & v);
  void Transpose (Vec<3> & v1, Vec<3> & v2, Vec<3> & v3);
  int SolveLinearSystem (const Vec<3> & col1,
                         const Vec<3> & col2,
                         const Vec<3> & col3,
                         const Vec<3> & rhs,
                         Vec<3> & sol);
  int SolveLinearSystemLS (const Vec<3> & col1,
                           const Vec<3> & col2,
                           const Vec<2> & rhs,
                           Vec<3> & sol);
  int SolveLinearSystemLS2 (const Vec<3> & col1,
                            const Vec<3> & col2,
                            const Vec<2> & rhs, 
                            Vec<3> & sol,
                            double & x, double & y);
  int PseudoInverse (const Vec<3> & col1,
                     const Vec<3> & col2,
                     Vec<3> & inv1,
                     Vec<3> & inv2);







  class QuadraticFunction3d
  {
    double c0, cx, cy, cz;
    double cxx, cyy, czz, cxy, cxz, cyz;

  public:
    QuadraticFunction3d (const Point<3> & p, const Vec<3> & v);
    double Eval (const Point<3> & p)
    {
      return 
        c0 
        + p(0) * (cx + cxx * p(0) + cxy * p(1) + cxz * p(2))
        + p(1) * (cy + cyy * p(1) + cyz * p(2))
        + p(2) * (cz + czz * p(2));
    }
  };




  ///
  class Box3d
  {
  protected:
    ///
    double minx[3], maxx[3];

  public:
    ///
    Box3d () { };
    ///
    DLL_HEADER Box3d ( double aminx, double amaxx,
            double aminy, double amaxy,
            double aminz, double amaxz );
    ///
    DLL_HEADER Box3d ( const Box3d & b2 );
    ///
    DLL_HEADER Box3d (const Point<3>& p1, const Point<3>& p2);
    ///
    DLL_HEADER Box3d (const Box<3> & b2);
    ///
    double MinX () const { return minx[0]; }
    ///
    double MaxX () const { return maxx[0]; }
    ///
    double MinY () const { return minx[1]; }
    ///
    double MaxY () const { return maxx[1]; }
    ///
    double MinZ () const { return minx[2]; }
    ///
    double MaxZ () const { return maxx[2]; }

    ///
    double Mini (int i) const { return minx[i-1]; }
    ///
    double Maxi (int i) const { return maxx[i-1]; }

    ///
    Point<3> PMin () const { return Point<3>(minx[0], minx[1], minx[2]); }
    ///
    Point<3> PMax () const { return Point<3>(maxx[0], maxx[1], maxx[2]); }

    ///
    void GetPointNr (int i, Point<3> & point) const;
    /// increase Box at each side with dist 
    void Increase (double dist);
    /// increase Box by factor rel
    void IncreaseRel (double rel);
    /// return 1 if closures are intersecting
    int Intersect (const Box3d & box2) const
    {
      if (minx[0] > box2.maxx[0] || maxx[0] < box2.minx[0] ||
          minx[1] > box2.maxx[1] || maxx[1] < box2.minx[1] ||
          minx[2] > box2.maxx[2] || maxx[2] < box2.minx[2])
        return 0;
      return 1;
    }
    /// return 1 if point p in closure
    int IsIn (const Point<3> & p) const
    {
      if (minx[0] <= p(0) && maxx[0] >= p(0) &&
          minx[1] <= p(1) && maxx[1] >= p(1) &&
          minx[2] <= p(2) && maxx[2] >= p(2))
        return 1;
      return 0;
    }
    ///
    inline void SetPoint (const Point<3> & p)
    {
      minx[0] = maxx[0] = p(0);
      minx[1] = maxx[1] = p(1);
      minx[2] = maxx[2] = p(2);    
    }

    ///
    inline void AddPoint (const Point<3> & p)
    {
      if (p(0) < minx[0]) minx[0] = p(0);
      if (p(0) > maxx[0]) maxx[0] = p(0);
      if (p(1) < minx[1]) minx[1] = p(1);
      if (p(1) > maxx[1]) maxx[1] = p(1);
      if (p(2) < minx[2]) minx[2] = p(2);
      if (p(2) > maxx[2]) maxx[2] = p(2);
    }

    ///
    const Box3d& operator+=(const Box3d& b);

    ///
    Point<3> MaxCoords() const;
    ///
    Point<3> MinCoords() const;

    /// Make a negative sized box;
    //  void CreateNegMinMaxBox();
  
    ///
    Point<3> CalcCenter () const { return Point<3>(0.5*(minx[0] + maxx[0]),
                                                 0.5*(minx[1] + maxx[1]),
                                                 0.5*(minx[2] + maxx[2])); }
    ///
    double CalcDiam () const { return sqrt(sqr(maxx[0]-minx[0])+
                                           sqr(maxx[1]-minx[1])+
                                           sqr(maxx[2]-minx[2])); }

    ///
    void WriteData(ofstream& fout) const;
    ///
    void ReadData(ifstream& fin);
  };


  class Box3dSphere : public Box3d
  {
  protected:
    ///
    double diam, inner;
    ///
    Point<3> c = Point<3>(0,0,0);
  public:
    ///
    Box3dSphere () { };
    ///
    Box3dSphere ( double aminx, double amaxx,
                  double aminy, double amaxy,
                  double aminz, double amaxz);
    ///
    const Point<3> & Center () const { return c; }

    ///
    double Diam () const { return diam; }
    ///
    double Inner () const { return inner; }
    ///
    void GetSubBox (int i, Box3dSphere & sbox) const;

    // private:
    ///
    void CalcDiamCenter ();
  };




  ///
  class referencetransform
  {
    ///
    Vec<3> ex, ey, ez;
    ///
    Vec<3> exh, eyh, ezh;
    ///
    Vec<3> ex_h, ey_h, ez_h;
    ///
    Point<3> rp = Point<3>(0,0,0);
    ///
    double h;

  public:

    ///
    void Set (const Point<3> & p1, const Point<3> & p2,
              const Point<3> & p3, double ah);

    ///
    void ToPlain (const Point<3> & p, Point<3> & pp) const;
    ///
    void ToPlain (const Array<Point<3>> & p, Array<Point<3>> & pp) const;
    ///
    void FromPlain (const Point<3> & pp, Point<3> & p) const;
  };

}


#endif

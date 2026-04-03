#ifndef FILE_OBJECTS
#define FILE_OBJECTS

/* *************************************************************************/
/* File:   geomobjects.hpp                                                 */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

#include <core/array.hpp>

#include <general/ngarray.hpp>

namespace netgen
{
  using namespace ngcore;

  template <int D, typename T = double> class Vec;
  template <int D, typename T = double> class Point;


  template <int D, typename T>
  class Point
  {

  protected:
    T x[D];

  public:
    Point () { ; }
    Point (T ax) { for (int i = 0; i < D; i++) x[i] = ax; }
    Point (T ax, T ay) 
    { 
      // static_assert(D==2, "Point<D> constructor with 2 args called");
      x[0] = ax; x[1] = ay; 
    }
    Point (T ax, T ay, T az) 
    {
      // static_assert(D==3, "Point<D> constructor with 3 args called");
      x[0] = ax; x[1] = ay; x[2] = az; 
    }
    Point (T ax, T ay, T az, T au)
    { x[0] = ax; x[1] = ay; x[2] = az; x[3] = au;}

    template <typename T2>
    Point (const Point<D,T2> & p2)
    { for (int i = 0; i < D; i++) x[i] = p2(i); }

    explicit Point (const Vec<D,T> & v)
    { for (int i = 0; i < D; i++) x[i] = v(i); }


    template <typename T2>
    Point & operator= (const Point<D,T2> & p2)
    {
      for (int i = 0; i < D; i++) x[i] = p2(i);
      return *this;
    }

    Point & operator= (T val)
    {
      for (int i = 0; i < D; i++) x[i] = val;
      return *this;
    }

    T & operator() (int i) { return x[i]; }
    const T & operator() (int i) const { return x[i]; }

    T& operator[] (int i) { return x[i]; }
    const T& operator[] (int i) const { return x[i]; }

    operator const T* () const { return x; }

    template<typename ARCHIVE>
    void DoArchive(ARCHIVE& archive)
    {
      for(int i=0; i<D; i++)
        archive & x[i];
    }
  };

  template <int D, typename T>
  class Vec
  {

  protected:
    T x[D];

  public:
    Vec () { ; } // for (int i = 0; i < D; i++) x[i] = 0; }
    Vec (T ax) { for (int i = 0; i < D; i++) x[i] = ax; }
    Vec (T ax, T ay) 
    { 
      // static_assert(D==2, "Vec<D> constructor with 2 args called");
      x[0] = ax; x[1] = ay; 
    }
    Vec (T ax, T ay, T az)
    { 
      // static_assert(D==3, "Vec<D> constructor with 3 args called");
      x[0] = ax; x[1] = ay; x[2] = az; 
    }
    Vec (T ax, T ay, T az, T au)
    { x[0] = ax; x[1] = ay; x[2] = az; x[3] = au; }

    Vec (const Vec<D> & p2)
    { for (int i = 0; i < D; i++) x[i] = p2.x[i]; }

    explicit Vec (const Point<D,T> & p)
    { for (int i = 0; i < D; i++) x[i] = p(i); }

    explicit Vec(const Point<D,T>& p1, const Point<D,T>& p2)
    { for(int i=0; i<D; i++) x[i] = p2(i)-p1(i); }

    template <typename T2>
    Vec & operator= (const Vec<D,T2> & p2)
    {
      for (int i = 0; i < D; i++) x[i] = p2(i);
      return *this;
    }

    Vec & operator= (T s)
    {
      for (int i = 0; i < D; i++) x[i] = s;
      return *this;
    }

    bool operator== (const Vec<D,T> &a) const
    {
      bool res = true;
      for (auto i : Range(D))
        res &= (x[i]==a.x[i]);
      return res;
    }

    T & operator() (int i) { return x[i]; }
    const T & operator() (int i) const { return x[i]; }

    T& operator[] (int i) { return x[i]; }
    const T& operator[] (int i) const { return x[i]; }

    operator const T* () const { return x; }

    template <typename ARCHIVE>
    void DoArchive(ARCHIVE& archive)
    {
      for(int i=0; i<D; i++)
        archive & x[i];
    }

    T Length () const
    {
      T l = 0;
      for (int i = 0; i < D; i++)
	l += x[i] * x[i];
      return sqrt (l);
    }

    T Length2 () const
    {
      T l = 0;
      for (int i = 0; i < D; i++)
	l += x[i] * x[i];
      return l;
    }

    Vec & Normalize ()
    {
      T l = Length();
      // if (l != 0)
      for (int i = 0; i < D; i++)
        x[i] /= (l+1e-40);
      return *this;
    }

    Vec<D> GetNormal () const;
  };

  template <int D>
  inline ostream & operator<< (ostream & ost, const Vec<D> & a)
  {
    ost << "(";
    for (int i = 0; i < D-1; i++)
      ost << a(i) << ", ";
    ost << a(D-1) << ")";
    return ost;
  }

  template <int D>
  inline ostream & operator<< (ostream & ost, const Point<D> & a)
  {
    ost << "(";
    for (int i = 0; i < D-1; i++)
      ost << a(i) << ", ";
    ost << a(D-1) << ")";
    return ost;
  }

  template<int D>
  inline Vec<D> operator-(const Point<D>& p1, const Point<D>& p2)
  {
    Vec<D> result;
    for(auto i : Range(D))
      result[i] = p1[i] - p2[i];
    return result;
  }

  template<int D>
  inline Vec<D> operator*(const Vec<D>& v, double d)
  {
    Vec<D> result;
    for(auto i : Range(D))
      result[i] = d*v[i];
    return result;
  }

  inline double Cross2(const Vec<2>& v1, const Vec<2>& v2)
  {
    return v1[0] * v2[1] - v1[1] * v2[0];
  }

  // are points clockwise?
  inline bool CW(const Point<2>& p1, const Point<2>& p2,
                  const Point<2>& p3)
  {
    return Cross2(p2-p1, p3-p2) < 0;
  }

  // are points counterclockwise?
  inline bool CCW(const Point<2>& p1, const Point<2>& p2,
                  const Point<2>& p3)
  {
    return Cross2(p2-p1, p3-p2) > 0;
  }

  // are strictly points counterclockwise?
  inline bool CCW(const Point<2>& p1, const Point<2>& p2,
                  const Point<2>& p3, double eps)
  {
    auto v1 = p2-p1;
    auto v2 = p3-p2;
    return Cross2(v1, v2) > eps*eps*max2(v1.Length2(),
                                         v2.Length2());
  }


  template <int H, int W=H, typename T = double>
  class Mat
  {

  protected:
    T x[H*W];

  public:
    Mat () { ; }
    Mat (const Mat & b)
    { for (int i = 0; i < H*W; i++) x[i] = b.x[i]; }
  
    Mat & operator= (T s)
    {
      for (int i = 0; i < H*W; i++) x[i] = s;
      return *this;
    }

    Mat & operator= (const Mat & b)
    {
      for (int i = 0; i < H*W; i++) x[i] = b.x[i]; 
      return *this;
    }

    T & operator() (int i, int j) { return x[i*W+j]; }
    const T & operator() (int i, int j) const { return x[i*W+j]; }
    T & operator() (int i) { return x[i]; }
    const T & operator() (int i) const { return x[i]; }

    Vec<H,T> Col (int i) const
    {
      Vec<H,T> hv; 
      for (int j = 0; j < H; j++)
	hv(j) = x[j*W+i];
      return hv; 
    }

    Vec<W,T> Row (int i) const
    {
      Vec<W,T> hv; 
      for (int j = 0; j < W; j++)
	hv(j) = x[i*W+j];
      return hv; 
    }

    void Solve (const Vec<H,T> & rhs, Vec<W,T> & sol) const
    {
      Mat<W,H,T> inv;
      CalcInverse (*this, inv);
      sol = inv * rhs;
    }

    template <typename ARCHIVE>
    void DoArchive(ARCHIVE & ar)
    {
      ar.Do(x, H*W);
    }
  };




  template <int D>
  class Box
  {
  protected:
    Point<D> pmin, pmax;
  public:
    Box () { ; }

    Box ( const Point<D> & p1)
    {
      for (int i = 0; i < D; i++)
	pmin(i) = pmax(i) = p1(i);
    }


    Box ( const Point<D> & p1, const Point<D> & p2)
    {
      for (int i = 0; i < D; i++)
	{
	  pmin(i) = min2(p1(i), p2(i));
	  pmax(i) = max2(p1(i), p2(i));
	}
    }

    Box (const Point<D> & p1, const Point<D> & p2, const Point<D> & p3)
      : Box(p1,p2)
    {
      Add (p3);
    }

    enum EB_TYPE { EMPTY_BOX = 1 };
    Box ( EB_TYPE et ) 
    {
      for (int i = 0; i < D; i++)
        {
          pmin(i) = 1e99;
          pmax(i) = -1e99;
        }
    }

    const Point<D> & PMin () const { return pmin; }
    const Point<D> & PMax () const { return pmax; }
  
    void Set (const Point<D> & p)
    { pmin = pmax = p; }

    void Add (const Point<D> & p)
    { 
      for (int i = 0; i < D; i++)
	{
	  if (p(i) < pmin(i)) pmin(i) = p(i);
	  /* else */ if (p(i) > pmax(i)) pmax(i) = p(i);
          // optimization invalid for empty-box !
	}
    }

    template <typename T1, typename T2>
    void Set (const NgIndirectArray<T1, T2> & points)
    {
      // Set (points[points.Begin()]);
      Set (points[*points.Range().begin()]);
      // for (int i = points.Begin()+1; i < points.End(); i++)
      for (int i : points.Range().Modify(1,0))
        Add (points[i]);
    }

    template <typename T1, typename T2>
    void Add (const NgIndirectArray<T1, T2> & points)
    {
      // for (int i = points.Begin(); i < points.End(); i++)
      for (int i : points.Range())
        Add (points[i]);
    }


    Point<D> Center () const 
    { 
      Point<D> c;
      for (int i = 0; i < D; i++)
	c(i) = 0.5 * (pmin(i)+pmax(i)); 
      return c;
    }
    double Diam () const { return Abs (pmax-pmin); }

    Point<D> GetPointNr (int nr) const
    {
      Point<D> p;
      for (int i = 0; i < D; i++)
	{
	  p(i) = (nr & 1) ? pmax(i) : pmin(i);
	  nr >>= 1;
	}
      return p;
    }


    bool Intersect (const Box<D> & box2) const
    {
      for (int i = 0; i < D; i++)
	if (pmin(i) > box2.pmax(i) ||
	    pmax(i) < box2.pmin(i)) return 0;
      return 1;
    }


    bool IsIn (const Point<D> & p) const
    {
      for (int i = 0; i < D; i++)
	if (p(i) < pmin(i) || p(i) > pmax(i)) return false;
      return true;
    }

    // is point in eps-increased box
    bool IsIn (const Point<D> & p, double eps) const
    {
      for (int i = 0; i < D; i++)
	if (p(i) < pmin(i)-eps || p(i) > pmax(i)+eps) return false;
      return true;
    }


    void Increase (double dist)
    {
      for (int i = 0; i < D; i++)
	{
	  pmin(i) -= dist;
	  pmax(i) += dist;
	}
    }

    void Scale (double factor)
    {
      auto center = Center();
      pmin = center + factor*(pmin-center);
      pmax = center + factor*(pmax-center);
    }

    template <typename ARCHIVE>
    void DoArchive(ARCHIVE & archive)
    { archive & pmin & pmax; }
  };




  template <int D>
  class BoxSphere : public Box<D>
  {
  protected:
    ///
    Point<D> c;
    ///
    double diam;
    ///
    double inner;
  public:
    ///
    BoxSphere () { };
    ///
    BoxSphere (const Box<D> & box) 
      : Box<D> (box) 
    { 
      CalcDiamCenter();
    };

    ///
    BoxSphere ( Point<D> apmin, Point<D> apmax )
      : Box<D> (apmin, apmax)
    {
      CalcDiamCenter();
    }

    ///
    const Point<D> & Center () const { return c; }
    ///
    double Diam () const { return diam; }
    ///
    double Inner () const { return inner; }


    ///
    void GetSubBox (int nr, BoxSphere & sbox) const
    {
      for (int i = 0; i < D; i++)
	{
	  if (nr & 1)
	    {
	      sbox.pmin(i) = c(i);
	      sbox.pmax(i) = this->pmax(i);
	    }
	  else
	    {
	      sbox.pmin(i) = this->pmin(i);
	      sbox.pmax(i) = c(i);
	    }
	  sbox.c(i) = 0.5 * (sbox.pmin(i) + sbox.pmax(i));
	  nr >>= 1;
	}
      sbox.diam = 0.5 * diam;
      sbox.inner = 0.5 * inner;
    }


    ///
    void CalcDiamCenter ()
    {
      c = Box<D>::Center ();
      diam = Dist (this->pmin, this->pmax);

      inner = this->pmax(0) - this->pmin(0);
      for (int i = 1; i < D; i++)
	if (this->pmax(i) - this->pmin(i) < inner)
	  inner = this->pmax(i) - this->pmin(i);
    }

  };


#ifdef PARALLEL_OLD
  template <> 
  inline MPI_Datatype MyGetMPIType<Vec<3, double> > ()
  {
    static MPI_Datatype MPI_T = 0;
    if (!MPI_T)
      {
	MPI_Type_contiguous ( 3, MPI_DOUBLE, &MPI_T);
	MPI_Type_commit ( &MPI_T );
      }
    return MPI_T;
  };
#endif
  

}


#endif

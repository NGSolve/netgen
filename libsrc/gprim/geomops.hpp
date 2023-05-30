#ifndef FILE_GEOMOPS
#define FILE_GEOMOPS

/* *************************************************************************/
/* File:   geomops.hpp                                                     */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/


namespace netgen
{

  /*

  Point - Vector operations

  */


  template <int D, typename T>
  inline Vec<D,T> operator+ (Vec<D,T> a, Vec<D,T> b)
  {
    Vec<D,T> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) + b(i);
    return res;
  }



  template <int D, typename T>
  inline Point<D,T> operator+ (Point<D,T> a, Vec<D,T> b)
  {
    Point<D,T> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) + b(i);
    return res;
  }



  template <int D, typename T>
  inline Vec<D,T> operator- (Point<D,T> a, Point<D,T> b)
  {
    Vec<D,T> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) - b(i);
    return res;
  }

  template <int D, typename T>
  inline Point<D,T> operator- (Point<D,T> a, Vec<D,T> b)
  {
    Point<D,T> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) - b(i);
    return res;
  }

  template <int D, typename T>
  inline Vec<D,T> operator- (Vec<D,T> a, Vec<D,T> b)
  {
    Vec<D,T> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) - b(i);
    return res;
  }



  template <int D, typename T>
  inline Vec<D,T> operator* (T s, Vec<D,T> b)
  {
    Vec<D,T> res;
    for (int i = 0; i < D; i++)
      res(i) = s * b(i);
    return res;
  }

  template <int D, int SW>
  inline Vec<D, SIMD<double,SW>> operator* (SIMD<double,SW> s, Vec<D,double> b)
  {
    Vec<D,SIMD<double,SW>> res;
    for (int i = 0; i < D; i++)
      res(i) = s * b(i);
    return res;
  }
  

  template <int D>
  inline double operator* (Vec<D> a, Vec<D> b)
  {
    double sum = 0;
    for (int i = 0; i < D; i++)
      sum += a(i) * b(i);
    return sum;
  }



  template <int D, typename T>
  inline Vec<D,T> operator- (Vec<D,T> b)
  {
    Vec<D,T> res;
    for (int i = 0; i < D; i++)
      res(i) = -b(i);
    return res;
  }


  template <int D, typename T>
  inline Point<D,T> & operator+= (Point<D,T> & a, Vec<D,T> b)
  {
    for (int i = 0; i < D; i++)
      a(i) += b(i);
    return a;
  }

  template <int D, typename T>
  inline Vec<D,T> & operator+= (Vec<D,T> & a, Vec<D> b)
  {
    for (int i = 0; i < D; i++)
      a(i) += b(i);
    return a;
  }


  template <int D>
  inline Point<D> & operator-= (Point<D> & a, const Vec<D> & b)
  {
    for (int i = 0; i < D; i++)
      a(i) -= b(i);
    return a;
  }

  template <int D>
  inline Vec<D> & operator-= (Vec<D> & a, const Vec<D> & b)
  {
    for (int i = 0; i < D; i++)
      a(i) -= b(i);
    return a;
  }



  template <int D, typename T1, typename T2>
  inline Vec<D,T1> & operator*= (Vec<D,T1> & a, T2 s)
  {
    for (int i = 0; i < D; i++)
      a(i) *= s;
    return a;
  }


  template <int D>
  inline Vec<D> & operator/= (Vec<D> & a, double s)
  {
    for (int i = 0; i < D; i++)
      a(i) /= s;
    return a;
  }




  // Matrix - Vector operations

  /*
    template <int H, int W>
    inline Vec<H> operator* (const Mat<H,W> & m, const Vec<W> & v)
    {
    Vec<H> res;
    for (int i = 0; i < H; i++)
    {
    res(i) = 0;
    for (int j = 0; j < W; j++)
    res(i) += m(i,j) * v(j);
    }
    return res;
    }
  */

  // thanks to VC60 partial template specialization features !!!

  inline Vec<2> operator* (const Mat<2,2> & m, const Vec<2> & v)
  {
    Vec<2> res;
    for (int i = 0; i < 2; i++)
      {
	res(i) = 0;
	for (int j = 0; j < 2; j++)
	  res(i) += m(i,j) * v(j);
      }
    return res;
  }

  inline Vec<2> operator* (const Mat<2,3> & m, const Vec<3> & v)
  {
    Vec<2> res;
    for (int i = 0; i < 2; i++)
      {
	res(i) = 0;
	for (int j = 0; j < 3; j++)
	  res(i) += m(i,j) * v(j);
      }
    return res;
  }


  inline Vec<3> operator* (const Mat<3,2> & m, const Vec<2> & v)
  {
    Vec<3> res;
    for (int i = 0; i < 3; i++)
      {
	res(i) = 0;
	for (int j = 0; j < 2; j++)
	  res(i) += m(i,j) * v(j);
      }
    return res;
  }


  inline Vec<3> operator* (const Mat<3,3> & m, const Vec<3> & v)
  {
    Vec<3> res;
    for (int i = 0; i < 3; i++)
      {
	res(i) = 0;
	for (int j = 0; j < 3; j++)
	  res(i) += m(i,j) * v(j);
      }
    return res;
  }







  /*
    template <int H1, int W1, int H2, int W2>
    inline Mat<H1,W2> operator* (const Mat<H1,W1> & a, const Mat<H2,W2> & b)
    {
    Mat<H1,W2> m;
    for (int i = 0; i < H1; i++)
    for (int j = 0; j < W2; j++)
    {
    double sum = 0;
    for (int k = 0; k < W1; k++)
    sum += a(i,k) * b(k, j);
    m(i,j) = sum; 
    }
    return m;
    }
  */

  template <typename T>
  inline Mat<2,2,T> operator* (const Mat<2,2,T> & a, const Mat<2,2,T> & b)
  {
    Mat<2,2,T> m;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
	{
          T sum(0);
	  for (int k = 0; k < 2; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }

  inline Mat<2,2> operator* (const Mat<2,3> & a, const Mat<3,2> & b)
  {
    Mat<2,2> m;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
	{
	  double sum = 0;
	  for (int k = 0; k < 3; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }

  template <typename T>
  inline Mat<3,2,T> operator* (const Mat<3,2,T> & a, const Mat<2,2,T> & b)
  {
    Mat<3,2,T> m;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
	{
	  T sum(0.0);
	  for (int k = 0; k < 2; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }



  inline Mat<2,3> operator* (const Mat<2,2> & a, const Mat<2,3> & b)
  {
    Mat<2,3> m;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 3; j++)
	{
	  double sum = 0;
	  for (int k = 0; k < 2; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }

  template <typename T>
  inline Mat<3,3,T> operator* (const Mat<3,3,T> & a, const Mat<3,3,T> & b)
  {
    Mat<3,3,T> m;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	{
	  T sum = T(0);
	  for (int k = 0; k < 3; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }








  template <int H, int W>
  inline Mat<W,H> Trans (const Mat<H,W> & m)
  {
    Mat<W,H> res;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
	res(j,i) = m(i,j);
    return res;
  }











  template <int D, typename T>
  inline ostream & operator<< (ostream & ost, const Vec<D,T> & a)
  {
    ost << "(";
    for (int i = 0; i < D-1; i++)
      ost << a(i) << ", ";
    ost << a(D-1) << ")";
    return ost;
  }

  template <int D, typename T>
  inline ostream & operator<< (ostream & ost, const Point<D,T> & a)
  {
    ost << "(";
    for (int i = 0; i < D-1; i++)
      ost << a(i) << ", ";
    ost << a(D-1) << ")";
    return ost;
  }

  template <int D>
  inline ostream & operator<< (ostream & ost, const Box<D> & b)
  {
    ost << b.PMin() << " - " << b.PMax();
    return ost;
  }

  template <int H, int W, typename T>
  inline ostream & operator<< (ostream & ost, const Mat<H,W,T> & m)
  {
    ost << "(";
    for (int i = 0; i < H; i++)
      {
	for (int j = 0; j < W; j++)
	  ost << m(i,j) << "   ";
	ost << endl;
      }
    return ost;
  }


}

#endif

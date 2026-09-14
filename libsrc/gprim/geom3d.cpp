#include <algorithm>
#include <mystdlib.h>

#include <myadt.hpp>
#include <gprim.hpp>

namespace netgen
{
ostream & operator<<(ostream  & s, const Point<3> & p)
  {
  return s << "(" << p(0) << ", " << p(1) << ", " << p(2) << ")";
  }

ostream & operator<<(ostream  & s, const Vec<3> & v)
  {
  return s << "(" << v(0) << ", " << v(1) << ", " << v(2) << ")";
  }

double Angle (const Vec<3> & v1, const Vec<3> & v2)
{
  double co = (v1 * v2) / (v1.Length() * v2.Length());
  if (co > 1) co = 1;
  if (co < -1) co = -1;
  return acos ( co );
}


void GetNormal (const Vec<3> & v, Vec<3> & n)
  {
  if (fabs (v(0)) > fabs (v(2)))
    {
    n(0) = -v(1);
    n(1) = v(0);
    n(2) = 0;
    }
  else
    {
    n(0) = 0;
    n(1) = v(2);
    n(2) = -v(1);
    }
  double len = n.Length();
  if (len == 0)
    {
    n(0) = 1;
    n(1) = n(2) = 0;
    }
  else
    n /= len;
  }

/*
ostream & operator<<(ostream  & s, const ROTDenseMatrix3D & r)
  {
  return s << "{ (" << r.txx << ", " << r.txy << ", " << r.txz << ") , ("
                    << r.tyx << ", " << r.tyy << ", " << r.tyz << ") , ("
                    << r.tzx << ", " << r.tzy << ", " << r.tzz << ") }";
  }
*/

/*
Vec<3> operator- (const Point<3> & p1, const Point<3> & p2)
  {
  return Vec<3> (p1.X() - p2.X(), p1.Y() - p2.Y(),p1.Z() - p2.Z());
  }

Point<3> operator- (const Point<3> & p1, const Vec<3> & v)
  {
  return Point<3> (p1.X() - v.X(), p1.Y() - v.Y(),p1.Z() - v.Z());
  }

Point<3> operator+ (const Point<3> & p1, const Vec<3> & v)
  {
  return Point<3> (p1.X() + v.X(), p1.Y() + v.Y(),p1.Z() + v.Z());
  }

Vec<3> operator- (const Vec<3> & v1, const Vec<3> & v2)
  {
  return Vec<3> (v1.X() - v2.X(), v1.Y() - v2.Y(),v1.Z() - v2.Z());
  }

Vec<3> operator+ (const Vec<3> & v1, const Vec<3> & v2)
  {
  return Vec<3> (v1.X() + v2.X(), v1.Y() + v2.Y(),v1.Z() + v2.Z());
  }

Vec<3> operator* (double scal, const Vec<3> & v)
  {
  return Vec<3> (scal * v.X(), scal * v.Y(), scal * v.Z());
  }
*/
/*
double operator* (const Vec<3> & v1, const Vec<3> & v2)
  {
  return v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z();
  }

double Cross (const Vec<3> & v1, const Vec<3> & v2)
  {
  return v1.X() * v2.Y() - v1.Y() * v2.X();
  }
*/

/*
void ROTDenseMatrix3D :: CalcRotMat(double ag, double bg, double lg, double size2, Vec<3> r)
  {
  size = size2;
  txx=size * ( cos(bg) * cos(lg) );
  txy=size * ( cos(bg) * sin(lg) );
  txz=size * (-sin(bg)           );

  tyx=size * ( sin(ag) * sin(bg) * cos(lg) - cos(ag) * sin(lg) );
  tyy=size * ( sin(ag) * sin(bg) * sin(lg) + cos(ag) * cos(lg) );
  tyz=size * ( sin(ag) * cos(bg)                               );

  tzx=size * ( cos(ag) * sin(bg) * cos(lg) + sin(ag) * sin(lg) );
  tzy=size * ( cos(ag) * sin(bg) * sin(lg) - sin(ag) * cos(lg) );
  tzz=size * ( cos(ag) * cos(bg)                               );

  deltaR=r;
  }
ROTDenseMatrix3D :: ROTDenseMatrix3D(double ag, double bg, double lg, double size2, Vec<3> r)
  {CalcRotMat(ag, bg, lg, size2, r); }

ROTDenseMatrix3D :: ROTDenseMatrix3D(Vec<3> rot2)
  {
  Vec<3> r2(0,0,0);
  CalcRotMat(rot2.X(), rot2.Y(), rot2.Z(), 1, r2);
  }

ROTDenseMatrix3D ROTDenseMatrix3D :: INV()
  {
  ROTDenseMatrix3D rinv(txx/sqr(size),tyx/sqr(size),tzx/sqr(size),
                   txy/sqr(size),tyy/sqr(size),tzy/sqr(size),
                   txz/sqr(size),tyz/sqr(size),tzz/sqr(size),
                   1/size,deltaR);
  return rinv;
  }

Vec<3> operator* (const ROTDenseMatrix3D & r, const Vec<3> & v)
  {
  return Vec<3> (r.XX() * v.X() + r.XY() * v.Y() + r.XZ() * v.Z(),
                r.YX() * v.X() + r.YY() * v.Y() + r.YZ() * v.Z(),
                r.ZX() * v.X() + r.ZY() * v.Y() + r.ZZ() * v.Z() );
  }

Point<3> operator* (const ROTDenseMatrix3D & r, const Point<3> & p)
  {
  return Point<3> (r.XX() * p.X() + r.XY() * p.Y() + r.XZ() * p.Z(),
                  r.YX() * p.X() + r.YY() * p.Y() + r.YZ() * p.Z(),
                  r.ZX() * p.X() + r.ZY() * p.Y() + r.ZZ() * p.Z() );
  }
*/







Box3d :: Box3d ( double aminx, double amaxx,
                 double aminy, double amaxy,
                 double aminz, double amaxz )
{
  minx[0] = aminx; maxx[0] = amaxx;
  minx[1] = aminy; maxx[1] = amaxy;
  minx[2] = aminz; maxx[2] = amaxz;
}

Box3d :: Box3d ( const Box3d & b2 )
{
  for (int i = 0; i < 3; i++)
    {
      minx[i] = b2.minx[i];
      maxx[i] = b2.maxx[i];
    }
}

Box3d :: Box3d ( const Box<3> & b2 )
{
  for (int i = 0; i < 3; i++)
    {
      minx[i] = b2.PMin()(i);
      maxx[i] = b2.PMax()(i);
    }
}


/*
int Box3d :: Intersect (const Box3d & box2) const
{
  int i;
  for (i = 0; i <= 2; i++)
    if (minx[i] > box2.maxx[i] || maxx[i] < box2.minx[i])
      return 0;
  return 1;
}
*/

/*
void Box3d :: SetPoint (const Point<3> & p)
{
  minx[0] = maxx[0] = p.X();
  minx[1] = maxx[1] = p.Y();
  minx[2] = maxx[2] = p.Z();
}

void Box3d :: AddPoint (const Point<3> & p)
{
  if (p.X() < minx[0]) minx[0] = p.X();
  if (p.X() > maxx[0]) maxx[0] = p.X();
  if (p.Y() < minx[1]) minx[1] = p.Y();
  if (p.Y() > maxx[1]) maxx[1] = p.Y();
  if (p.Z() < minx[2]) minx[2] = p.Z();
  if (p.Z() > maxx[2]) maxx[2] = p.Z();
}
*/

void Box3d :: GetPointNr (int i, Point<3> & point) const
{
  i--;
  point(0) = (i & 1) ? maxx[0] : minx[0];
  point(1) = (i & 2) ? maxx[1] : minx[1];
  point(2) = (i & 4) ? maxx[2] : minx[2];
}


void Box3d :: Increase (double d)
{
  for (int i = 0; i <= 2; i++)
    {
      minx[i] -= d;
      maxx[i] += d;
    }
}

void Box3d :: IncreaseRel (double /* rel */)
{
  for (int i = 0; i <= 2; i++)
    {
      double d = 0.5 * (maxx[i] - minx[i]);
      minx[i] -= d;
      maxx[i] += d;
    }
}


Box3d :: Box3d (const Point<3>& p1, const Point<3>& p2)
{
  minx[0] = min2 (p1(0), p2(0));
  minx[1] = min2 (p1(1), p2(1));
  minx[2] = min2 (p1(2), p2(2));
  maxx[0] = max2 (p1(0), p2(0));
  maxx[1] = max2 (p1(1), p2(1));
  maxx[2] = max2 (p1(2), p2(2));
}

const Box3d& Box3d :: operator+=(const Box3d& b)
{
  minx[0] = min2 (minx[0], b.minx[0]);
  minx[1] = min2 (minx[1], b.minx[1]);
  minx[2] = min2 (minx[2], b.minx[2]);
  maxx[0] = max2 (maxx[0], b.maxx[0]);
  maxx[1] = max2 (maxx[1], b.maxx[1]);
  maxx[2] = max2 (maxx[2], b.maxx[2]);

  return *this;
}

Point<3> Box3d :: MaxCoords() const
{
  return Point<3>(maxx[0], maxx[1], maxx[2]);
}

Point<3> Box3d :: MinCoords() const
{
  return Point<3>(minx[0], minx[1], minx[2]);
}

/*
void Box3d :: CreateNegMinMaxBox()
{
  minx[0] = MAXDOUBLE;
  minx[1] = MAXDOUBLE;
  minx[2] = MAXDOUBLE;
  maxx[0] = MINDOUBLE;
  maxx[1] = MINDOUBLE;
  maxx[2] = MINDOUBLE;

}
*/

void Box3d :: WriteData(ofstream& fout) const
{
  for(int i = 0; i < 3; i++)
    {
      fout << minx[i] << " " << maxx[i] << " ";
    }
  fout << "\n";
}

void Box3d :: ReadData(ifstream& fin)
{
  for(int i = 0; i < 3; i++)
    {
      fin >> minx[i];
      fin >> maxx[i];
    }
}




Box3dSphere :: Box3dSphere ( double aminx, double amaxx,
			     double aminy, double amaxy,
			     double aminz, double amaxz )
  : Box3d (aminx, amaxx, aminy, amaxy, aminz, amaxz)
{
  CalcDiamCenter ();
}


void Box3dSphere :: CalcDiamCenter ()
{
  diam = sqrt( sqr (maxx[0] - minx[0]) +
	       sqr (maxx[1] - minx[1]) + 
	       sqr (maxx[2] - minx[2]));
  
  c(0) = 0.5 * (minx[0] + maxx[0]);
  c(1) = 0.5 * (minx[1] + maxx[1]);
  c(2) = 0.5 * (minx[2] + maxx[2]);
  
  inner = min2 ( min2 (maxx[0] - minx[0], maxx[1] - minx[1]), maxx[2] - minx[2]) / 2;
}


void Box3dSphere :: GetSubBox (int i, Box3dSphere & sbox) const
{
  i--;
  if (i & 1)
    {
      sbox.minx[0] = c(0);
      sbox.maxx[0] = maxx[0];
    }
  else
    {
      sbox.minx[0] = minx[0];
      sbox.maxx[0] = c(0);
    }
  if (i & 2)
    {
      sbox.minx[1] = c(1);
      sbox.maxx[1] = maxx[1];
    }
  else
    {
      sbox.minx[1] = minx[1];
      sbox.maxx[1] = c(1);
    }
  if (i & 4)
    {
      sbox.minx[2] = c(2);
      sbox.maxx[2] = maxx[2];
    }
  else
    {
      sbox.minx[2] = minx[2];
      sbox.maxx[2] = c(2);
    }
  
  //  sbox.CalcDiamCenter ();

  sbox.c(0) = 0.5 * (sbox.minx[0] + sbox.maxx[0]);
  sbox.c(1) = 0.5 * (sbox.minx[1] + sbox.maxx[1]);
  sbox.c(2) = 0.5 * (sbox.minx[2] + sbox.maxx[2]);
  sbox.diam = 0.5 * diam;
  sbox.inner = 0.5 * inner;
}




/*
double Determinant (const Vec<3> & col1,
		    const Vec<3> & col2,
		    const Vec<3> & col3)
{
  return
    col1(0) * ( col2(1) * col3(2) - col2(2) * col3(1)) +
    col1(1) * ( col2(2) * col3(0) - col2(0) * col3(2)) +
    col1(2) * ( col2(0) * col3(1) - col2(1) * col3(0));
}
*/

void Transpose (Vec<3> & v1, Vec<3> & v2, Vec<3> & v3)
{
  Swap (v1(1), v2(0));
  Swap (v1(2), v3(0));
  Swap (v2(2), v3(1));
}


/*
  gcc4.8.3  warning: array subscript is above array bounds [-Warray-bounds]
 */
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif

int SolveLinearSystem (const Vec<3> & col1, const Vec<3> & col2,
		       const Vec<3> & col3, const Vec<3> & rhs,
		       Vec<3> & sol)
{
  // changed by MW
  double matrix[3][3];
  double locrhs[3];
  int retval = 0;

  for(int i=0; i<3; i++)
    {
      matrix[i][0] = col1(i);
      matrix[i][1] = col2(i);
      matrix[i][2] = col3(i);
      locrhs[i] = rhs(i);
    }

  for(int i=0; i<2; i++)
    {
      int pivot = i;
      double maxv = fabs(matrix[i][i]);
      for(int j=i+1; j<3; j++)
	if(fabs(matrix[j][i]) > maxv)
	  {
	    maxv = fabs(matrix[j][i]);
	    pivot = j;
	  }

      if(fabs(maxv) > 1e-40)
	{
	  if(pivot != i)
	    {
	      swap(matrix[i][0],matrix[pivot][0]);
	      swap(matrix[i][1],matrix[pivot][1]);
	      swap(matrix[i][2],matrix[pivot][2]);
	      swap(locrhs[i],locrhs[pivot]);
	    }
	  for(int j=i+1; j<3; j++)
	    {
	      double fac = matrix[j][i] / matrix[i][i];
	      
	      for(int k=i+1; k<3; k++)
		matrix[j][k] -= fac*matrix[i][k];
	      locrhs[j] -= fac*locrhs[i];
	    }
	}
      else
	retval = 1;
    }

  if(fabs(matrix[2][2]) < 1e-40)
    retval = 1;

  if(retval != 0)
    return retval;
  

  for(int i=2; i>=0; i--)
    {
      double sum = locrhs[i];
      for(int j=2; j>i; j--)
	sum -= matrix[i][j]*sol(j);

      sol(i) = sum/matrix[i][i];
    }

  return 0;
  
  
  


  /*
  double det = Determinant (col1, col2, col3);
  if (fabs (det) < 1e-40)
    return 1;
  
  sol.X() = Determinant (rhs, col2, col3) / det;
  sol.Y() = Determinant (col1, rhs, col3) / det;
  sol.Z() = Determinant (col1, col2, rhs) / det;

  return 0;
  */
  /*
  Vec<3> cr;
  Cross (col1, col2, cr);
  double det = cr * col3;

  if (fabs (det) < 1e-40)
    return 1;

  if (fabs(cr.Z()) > 1e-12)
    {
      // solve for 3. component
      sol.Z() = (cr * rhs) / det;
      
      // 2x2 system for 1. and 2. component
      double res1 = rhs.X() - sol.Z() * col3.X();
      double res2 = rhs.Y() - sol.Z() * col3.Y();
      
      sol.X() = (col2.Y() * res1 - col2.X() * res2) / cr.Z();
      sol.Y() = (col1.X() * res2 - col1.Y() * res1) / cr.Z();
  
    }
  else
    {
      det = Determinant (col1, col2, col3);
      if (fabs (det) < 1e-40)
	return 1;
      
      sol.X() = Determinant (rhs, col2, col3) / det;
      sol.Y() = Determinant (col1, rhs, col3) / det;
      sol.Z() = Determinant (col1, col2, rhs) / det;
    }

  return 0;
  */
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif



int SolveLinearSystemLS (const Vec<3> & col1,
			 const Vec<3> & col2,
			 const Vec<2> & rhs,
			 Vec<3> & sol)
{
  double a11 = col1 * col1;
  double a12 = col1 * col2;
  double a22 = col2 * col2;
  
  double det = a11 * a22 - a12 * a12;

  if (det*det <= 1e-24 * a11 * a22)
    {
      sol = Vec<3> (0, 0, 0);
      return 1;
    }
  
  Vec<2> invrhs;
  invrhs(0) = ( a22 * rhs(0) - a12 * rhs(1)) / det;
  invrhs(1) = (-a12 * rhs(0) + a11 * rhs(1)) / det;

  sol(0) = invrhs(0) * col1(0) + invrhs(1) * col2(0);
  sol(1) = invrhs(0) * col1(1) + invrhs(1) * col2(1);
  sol(2) = invrhs(0) * col1(2) + invrhs(1) * col2(2);

  return 0;

  /*
  Vec<3> inv1, inv2;
  int err = 
    PseudoInverse (col1, col2, inv1, inv2);

   sol = rhs.X() * inv1 + rhs.Y() * inv2;
   return err;
  */
}

int SolveLinearSystemLS2 (const Vec<3> & col1,
			 const Vec<3> & col2,
			 const Vec<2> & rhs,
			 Vec<3> & sol, double & x, double & y)
{
  double a11 = col1 * col1;
  double a12 = col1 * col2;
  double a22 = col2 * col2;
  
  double det = a11 * a22 - a12 * a12;

  if (fabs (det) <= 1e-12 * col1.Length() * col2.Length() || 
      col1.Length2() == 0 || col2.Length2() == 0)
    {
      sol = Vec<3> (0, 0, 0);
      x = 0; y = 0;
      return 1;
    }
  
  Vec<2> invrhs;
  invrhs(0) = ( a22 * rhs(0) - a12 * rhs(1)) / det;
  invrhs(1) = (-a12 * rhs(0) + a11 * rhs(1)) / det;

  sol(0) = invrhs(0) * col1(0) + invrhs(1) * col2(0);
  sol(1) = invrhs(0) * col1(1) + invrhs(1) * col2(1);
  sol(2) = invrhs(0) * col1(2) + invrhs(1) * col2(2);

  x = invrhs(0);
  y = invrhs(1);

  return 0;

  /*
  Vec<3> inv1, inv2;
  int err = 
    PseudoInverse (col1, col2, inv1, inv2);

   sol = rhs.X() * inv1 + rhs.Y() * inv2;
   return err;
  */
}

int PseudoInverse (const Vec<3> & col1,
		   const Vec<3> & col2,
		   Vec<3> & inv1,
		   Vec<3> & inv2)
{
  double a11 = col1 * col1;
  double a12 = col1 * col2;
  double a22 = col2 * col2;
  
  double det = a11 * a22 - a12 * a12;

  if (fabs (det) < 1e-12 * col1.Length() * col2.Length())
    {
      inv1 = Vec<3> (0, 0, 0);
      inv2 = Vec<3> (0, 0, 0);
      return 1;
    }

  double ia11 = a22 / det;
  double ia12 = -a12 / det;
  double ia22 = a11 / det;

  inv1 = ia11 * col1 + ia12 * col2;
  inv2 = ia12 * col1 + ia22 * col2;

  return 0;
}




QuadraticFunction3d :: 
QuadraticFunction3d (const Point<3> & p, const Vec<3> & v)
{
  Vec<3> hv(v);
  hv /= (hv.Length() + 1e-12);
  Vec<3> t1, t2;
  GetNormal (hv, t1);
  Cross (hv, t1, t2);
  
  double t1p = t1(0) * p(0) + t1(1) * p(1) + t1(2) * p(2);
  double t2p = t2(0) * p(0) + t2(1) * p(1) + t2(2) * p(2);
  c0 = sqr (t1p) + sqr (t2p);
  cx = -2 * (t1p * t1(0) + t2p * t2(0));
  cy = -2 * (t1p * t1(1) + t2p * t2(1));
  cz = -2 * (t1p * t1(2) + t2p * t2(2));

  cxx = t1(0) * t1(0) + t2(0) * t2(0);
  cyy = t1(1) * t1(1) + t2(1) * t2(1);
  czz = t1(2) * t1(2) + t2(2) * t2(2);

  cxy = 2 * t1(0) * t1(1) + 2 * t2(0) * t2(1);
  cxz = 2 * t1(0) * t1(2) + 2 * t2(0) * t2(2);
  cyz = 2 * t1(1) * t1(2) + 2 * t2(1) * t2(2);

  /*
  (*testout) << "c0 = " << c0
	     << " clin = " << cx << " " << cy << " " << cz 
	     << " cq = " << cxx << " " << cyy << " " << czz
	     << cxy << " " << cyz << " " << cyz << endl;
  */
}

// QuadraticFunction3d gqf (Point<3> (0,0,0), Vec<3> (1, 0, 0));





void referencetransform :: Set (const Point<3> & p1, const Point<3> & p2,
                                const Point<3> & p3, double ah)
{
  ex = p2 - p1;
  ex /= ex.Length();
  ey = p3 - p1;
  ey -= (ex * ey) * ex;
  ey /= ey.Length();
  ez = Cross (ex, ey);
  rp = p1;
  h = ah;

  exh = ah * ex;
  eyh = ah * ey;
  ezh = ah * ez;
  ah = 1 / ah;
  ex_h = ah * ex;
  ey_h = ah * ey;
  ez_h = ah * ez;
}

void referencetransform :: ToPlain (const Point<3> & p, Point<3> & pp) const
{
  Vec<3> v;
  v = p - rp;
  pp(0) = (ex_h * v);
  pp(1) = (ey_h * v);
  pp(2) = (ez_h * v);
}

void referencetransform :: ToPlain (const Array<Point<3>> & p,
                                    Array<Point<3>> & pp) const
{
  Vec<3> v;
  int i;

  pp.SetSize (p.Size());
  for (i = 1; i <= p.Size(); i++)
    {
      v = p[i-1] - rp;
      pp[i-1](0) = (ex_h * v);
      pp[i-1](1) = (ey_h * v);
      pp[i-1](2) = (ez_h * v);
    }
}

void referencetransform :: FromPlain (const Point<3> & pp, Point<3> & p) const
{
  Vec<3> v;
  //  v = (h * pp.X()) * ex + (h * pp.Y()) * ey + (h * pp.Z()) * ez;
  //  p = rp + v;
  v(0) = pp(0) * exh(0) + pp(1) * eyh(0) + pp(2) * ezh(0);
  v(1) = pp(0) * exh(1) + pp(1) * eyh(1) + pp(2) * ezh(1);
  v(2) = pp(0) * exh(2) + pp(1) * eyh(2) + pp(2) * ezh(2);
  p(0) = rp(0) + v(0);
  p(1) = rp(1) + v(1);
  p(2) = rp(2) + v(2);
}


}

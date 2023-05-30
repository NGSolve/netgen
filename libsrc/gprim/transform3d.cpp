#include <mystdlib.h>

#include <myadt.hpp>
#include <gprim.hpp>
#include <linalg.hpp>

namespace netgen
{

Transformation3d :: Transformation3d ()
{
  for (int i = 0; i < 3; i++)
    {
      offset[i] = 0;
      for (int j = 0; j < 3; j++)
	lin[i][j] = 0;
    }
}

Transformation3d :: Transformation3d (const Vec3d & translate)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lin[i][j] = 0;
  for (int i = 0; i < 3; i++)
    {
      offset[i] = translate.X(i+1);
      lin[i][i] = 1;
    }
}


Transformation3d :: 
Transformation3d (const Point3d & c, double alpha, 
		  double beta, double gamma)
{
  // total = T_c x Rot_0 x T_c^{-1}
  // Use Euler angles, see many books from tech mech, e.g. 
  // Shabana "multibody systems"

  Transformation3d tc(c);
  Transformation3d tcinv;
  tc.CalcInverse (tcinv);

  Transformation3d r1, r2, r3, ht, ht2;
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




Transformation3d :: Transformation3d (const Point3d ** pp)
{
  for (int i = 1; i <= 3; i++)
    {
      offset[i-1] = (*pp[0]).X(i);
      for (int j = 1; j <= 3; j++)
	lin[i-1][j-1] = (*pp[j]).X(i) - (*pp[0]).X(i);
    }
}

Transformation3d :: Transformation3d (const Point3d pp[])
{
  for (int i = 1; i <= 3; i++)
    {
      offset[i-1] = pp[0].X(i);
      for (int j = 1; j <= 3; j++)
	lin[i-1][j-1] = pp[j].X(i) - pp[0].X(i);
    }
}


void Transformation3d :: CalcInverse (Transformation3d & inv) const
{
  static DenseMatrix a(3), inva(3);
  static Vector b(3), sol(3);
  
  for (int i = 0; i < 3; i++)
    {
      b(i) = offset[i];
      for (int j = 0; j < 3; j++)
	a(i, j) = lin[i][j];
    }

  ::netgen::CalcInverse (a, inva);
  inva.Mult (b, sol);

  for (int i = 0; i < 3; i++)
    {
      inv.offset[i] = -sol(i);
      for (int j = 0; j < 3; j++)
	inv.lin[i][j] = inva(i, j);
    }
}


void  Transformation3d:: 
Combine (const Transformation3d & ta, const Transformation3d & tb)
{
  // o = o_a+ m_a o_b
  // m = m_a m_b

  for (int i = 0; i <= 2; i++)
    {
      offset[i] = ta.offset[i];
      for (int j = 0; j <= 2; j++)
	offset[i] += ta.lin[i][j] * tb.offset[j];
    }
  
  for (int i = 0; i <= 2; i++)
    for (int j = 0; j <= 2; j++)
      {
	lin[i][j] = 0;
	for (int k = 0; k <= 2; k++)
	  lin[i][j] += ta.lin[i][k] * tb.lin[k][j];
      }
}
void Transformation3d :: SetAxisRotation (int dir, double alpha)
{
  double co = cos(alpha);
  double si = sin(alpha);
  dir--;
  int pos1 = (dir+1) % 3;
  int pos2 = (dir+2) % 3;

  int i, j;
  for (i = 0; i <= 2; i++)
    {
      offset[i] = 0;
      for (j = 0; j <= 2; j++)
	lin[i][j] = 0;
    }

  lin[dir][dir] = 1;
  lin[pos1][pos1] = co;
  lin[pos2][pos2] = co;
  lin[pos1][pos2] = si;
  lin[pos2][pos1] = -si;
}

ostream & operator<< (ostream & ost, Transformation3d & trans)
{
  ost << "offset = ";
  for (int i = 0; i <= 2; i++)
    ost << trans.offset[i] << " ";
  ost << endl << "linear = " << endl;
  for (int i = 0; i <= 2; i++)
    {
      for (int j = 0; j <= 2; j++)
	ost << trans.lin[i][j] << " ";
      ost << endl;
    }
  return ost;
}

  template <>
  Transformation<3> :: Transformation (const Point<3> & c, const Vec<3> & axes, double angle)
  {
    Vec<3> vc(c);
    Transformation<3> tc(vc);
    Transformation<3> tcinv(-vc);
    Transformation<3> r, ht, ht2;

    // r.SetAxisRotation (3, alpha);
    Vec<3> naxes = axes;
    naxes.Normalize();
    Vec<3> n1 = naxes.GetNormal();
    Vec<3> n2 = Cross(naxes, n1);
    r.v = Vec<3>(0,0,0);
    double co = cos(angle);
    double si = sin(angle);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        r.m(i,j) = naxes(i)*naxes(j) + co*(n1(i)*n1(j)+n2(i)*n2(j)) + si*( (n2(i)*n1(j)-n2(j)*n1(i)) );

    ht.Combine (tc, r);
    Combine (ht, tcinv);
  }

  
  template <int D>
  Transformation<D> :: Transformation (const Point<D> * pp)
  {
    v = Vec<D> (pp[0]);
    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++)
        m(j,i) = pp[i+1](j)-pp[0](j);
  }
  
  template class Transformation<3>;
}

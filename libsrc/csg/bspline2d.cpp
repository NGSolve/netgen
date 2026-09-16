#include <mystdlib.h>

#include <csg.hpp>

namespace netgen
{

BSplineCurve2d :: BSplineCurve2d ()
{
  redlevel = 0;
}


void BSplineCurve2d :: AddPoint (const Point<2> & apoint)
{
  points.Append (apoint);
  intervallused.Append (0);
}

bool BSplineCurve2d :: Inside (const Point<2> & p, double & dist) const
{
  Point<2> hp = p;
  double t = ProjectParam (p);
  hp = Eval(t);
  Vec<2> v = EvalPrime (t);

  Vec<2> n (v(0), -v(1));
  
  cout << "p = " << p << ", hp = " << hp << endl;
  dist = Dist (p, hp);
  double scal = (hp-p) * n;
  cout << "scal = " << scal << endl;

  return scal >= 0;
}
  
double BSplineCurve2d :: ProjectParam (const Point<2> & p) const
{
  double t, dt, mindist, mint = 0.0;
  int n1;
  
  mindist = 1e10;
  dt = 0.2;
  for (n1 = 1; n1 <= points.Size(); n1++)
    if (intervallused[n1-1] == 0)
      for (t = n1; t <= n1+1; t += dt)
        if (Dist (Eval(t), p) < mindist)
          {
            mint = t;
            mindist = Dist (Eval(t), p);
          }
    
  if (mindist > 1e9) 
    {
      for (t = 0; t <= points.Size(); t += dt)
        if (Dist (Eval(t), p) < mindist)
          {
            mint = t;
            mindist = Dist (Eval(t), p);
          }   
    }

  while (Dist (Eval (mint-dt), p) < mindist)
    {
      mindist = Dist (Eval (mint-dt), p);
      mint -= dt;
    }
  while (Dist (Eval (mint+dt), p) < mindist)
    {
      mindist = Dist (Eval (mint+dt), p);
      mint += dt;
    }


  return NumericalProjectParam (p, mint-dt, mint+dt);  
}
  
  
// t \in (n1, n2)  
  
Point<2> BSplineCurve2d :: Eval (double t) const
{
  int n, n1, n2, n3, n4;
  double loct, b1, b2, b3, b4;
  Point<2> hp;

  static int cnt = 0;
  cnt++;
  if (cnt % 100000 == 0) cout << "cnt = " << cnt << endl;
  
  n = int(t);   
  loct = t - n;
  
  b1 = 0.25 * (1 - loct) * (1 - loct);
  b4 = 0.25 * loct * loct;
  b2 = 0.5 - b4;
  b3 = 0.5 - b1;
  
  n1 = (n + 10 * points.Size() -1) % points.Size() + 1;
  n2 = n1+1;
  if (n2 > points.Size()) n2 = 1;
  n3 = n2+1;
  if (n3 > points.Size()) n3 = 1;
  n4 = n3+1;
  if (n4 > points.Size()) n4 = 1;

  //  cout << "t = " << t << " n = " << n << " loct = " << loct 
  //      << " n1 = " << n1 << endl;

  
  hp(0) = b1 * points[n1-1](0) + b2 * points[n2-1](0) +
    b3 * points[n3-1](0) + b4 * points[n4-1](0);
  hp(1) = b1 * points[n1-1](1) + b2 * points[n2-1](1) +
    b3 * points[n3-1](1) + b4 * points[n4-1](1);
  return hp;
}
  
Vec<2> BSplineCurve2d :: EvalPrime (double t) const
{
  int n, n1, n2, n3, n4;
  double loct, db1, db2, db3, db4;
  Vec<2> hv;
  
  n = int(t);   
  loct = t - n;
  
  db1 = 0.5 * (loct - 1);
  db4 = 0.5 * loct;
  db2 = -db4;
  db3 = -db1;
  
  n1 = (n + 10 * points.Size() -1) % points.Size() + 1;
  n2 = n1+1;
  if (n2 > points.Size()) n2 = 1;
  n3 = n2+1;
  if (n3 > points.Size()) n3 = 1;
  n4 = n3+1;
  if (n4 > points.Size()) n4 = 1;
  
  hv(0) = db1 * points[n1-1](0) + db2 * points[n2-1](0) +
    db3 * points[n3-1](0) + db4 * points[n4-1](0);
  hv(1) = db1 * points[n1-1](1) + db2 * points[n2-1](1) +
    db3 * points[n3-1](1) + db4 * points[n4-1](1);
  return hv;
}

Vec<2> BSplineCurve2d :: EvalPrimePrime (double t) const
{
  int n, n1, n2, n3, n4;
  double ddb1, ddb2, ddb3, ddb4;
  Vec<2> hv;
  
  n = int(t);   
  //  double loct = t - n;
  
  ddb1 = 0.5;
  ddb4 = 0.5;
  ddb2 = -0.5;
  ddb3 = -0.5;
  
  n1 = (n + 10 * points.Size() -1) % points.Size() + 1;
  n2 = n1+1;
  if (n2 > points.Size()) n2 = 1;
  n3 = n2+1;
  if (n3 > points.Size()) n3 = 1;
  n4 = n3+1;
  if (n4 > points.Size()) n4 = 1;
  
  hv(0) = ddb1 * points[n1-1](0) + ddb2 * points[n2-1](0) +
    ddb3 * points[n3-1](0) + ddb4 * points[n4-1](0);
  hv(1) = ddb1 * points[n1-1](1) + ddb2 * points[n2-1](1) +
    ddb3 * points[n3-1](1) + ddb4 * points[n4-1](1);
  return hv;
}
  

int BSplineCurve2d :: SectionUsed (double t) const
{
  int n1 = int(t);   
  n1 = (n1 + 10 * points.Size() - 1) % points.Size() + 1;
  return (intervallused[n1-1] == 0);
}

void BSplineCurve2d :: Reduce (const Point<2> & p, double rad)
{
  int n1, n;
  int j; 
  double minx, miny, maxx, maxy;
  
  //  (*testout) << "Reduce: " << p << "," << rad << endl;
  
  redlevel++;
  
  for (n1 = 1; n1 <= points.Size(); n1++)
    {
      if (intervallused[n1-1] != 0) continue;
    
      minx = maxx = points[n1-1](0);
      miny = maxy = points[n1-1](1);
    
      n = n1;
      for (j = 1; j <= 3; j++)
        {
          n++;
          if (n > points.Size()) n = 1;
          if (points[n-1](0) < minx) minx = points[n-1](0);
          if (points[n-1](1) < miny) miny = points[n-1](1);
          if (points[n-1](0) > maxx) maxx = points[n-1](0);
          if (points[n-1](1) > maxy) maxy = points[n-1](1);
        }
      
      if (minx > p(0) + rad || maxx < p(0) - rad ||
          miny > p(1) + rad || maxy < p(1) - rad)
        {
          intervallused[n1-1] = redlevel;
          //      (*testout) << 0;
        }
      else
        {
          //      (*testout) << 1;
          intervallused[n1-1] = 0;
        }
    }
  //  (*testout) << endl;
}

void BSplineCurve2d :: UnReduce () 
{
  int i;
  for (i = 1; i <= intervallused.Size(); i++)
    if (intervallused[i-1] == redlevel)
      intervallused[i-1] = 0;
  redlevel--;
}
  
void BSplineCurve2d :: Print (ostream & ost) const
{
  ost << "SplineCurve: " << points.Size() << " points." << endl;
  for (int i = 1; i <= points.Size(); i++)
    ost << "P" << i << " = " << points[i-1] << endl;
}
}

#include <mystdlib.h>
#include "meshing.hpp"



namespace netgen
{


vnetrule :: vnetrule ()
{
  name = new char[1];
  name[0] = char(0);
  quality = 0;
}

vnetrule :: ~vnetrule ()
{
  // if (strlen(name)) 
  delete [] name;
  for (int i = 1; i <= freefaces.Size(); i++)
    delete freefaces[i-1];
  for (int i = 1; i <= freesets.Size(); i++)
    delete freesets[i-1];
  for (int i = 1; i <= freeedges.Size(); i++)
    delete freeedges[i-1];
  for (int i = 1; i <= freefaceinequ.Size(); i++)
    delete freefaceinequ[i-1];
  delete oldutofreezone;
  delete oldutofreezonelimit;
}

int vnetrule :: TestFlag (char flag) const
{
  for (int i = 1; i <= flags.Size(); i++)
    if (flags[i-1] == flag) return 1;
  return 0;
}


void vnetrule :: SetFreeZoneTransformation (const Vector & allp, int tolclass)
{
  int i, j;
  // double nx, ny, nz, v1x, v1y, v1z, v2x, v2y, v2z;
  double nl;
  const threeint * ti;
  int fs;

  double lam1 = 1.0/(2 * tolclass - 1);
  double lam2 = 1-lam1;

  transfreezone.SetSize (freezone.Size());
  
  int np = points.Size();
  int nfp = freezone.Size();
  Vector vp(np), vfp1(nfp), vfp2(nfp);


  for (i = 1; i <= 3; i++)
    {
      for (j = 1; j <= np; j++)
        vp(j-1) = allp(i+3*j-3-1);

      oldutofreezone->Mult (vp, vfp1);
      oldutofreezonelimit->Mult (vp, vfp2);

      vfp1 *= lam1;
      vfp1.Add (lam2, vfp2);

      for (j = 1; j <= nfp; j++)
        transfreezone[j-1](i-1) = vfp1(j-1);
    }

  // MARK(setfz2);


  fzbox.SetPoint (transfreezone[0]);
  for (i = 2; i <= freezone.Size(); i++)
    fzbox.AddPoint (transfreezone[i-1]);
  fzbox.IncreaseRel(1e-8);
  
  
  // MARK(setfz3);


  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      Array<threeint> & freesetfaces = *freefaces[fs-1];
      DenseMatrix & freesetinequ = *freefaceinequ[fs-1];
      
      for (i = 1; i <= freesetfaces.Size(); i++)
        {
          ti = &freesetfaces[i-1];
          const Point<3> & p1 = transfreezone[(ti->i1)-1];
          const Point<3> & p2 = transfreezone[(ti->i2)-1];
          const Point<3> & p3 = transfreezone[(ti->i3)-1];

          Vec<3> v1(p1, p2);   
          Vec<3> v2(p1, p3);   
          Vec<3> n;
          Cross (v1, v2, n);

          nl = n.Length();

          if (nl < 1e-10)
            {
              freesetinequ.Set(1, 1, 0);
              freesetinequ.Set(1, 2, 0);
              freesetinequ.Set(1, 3, 0);
              freesetinequ.Set(1, 4, -1);
            }
          else
            {
              //              n /= nl;
              
              freesetinequ.Set(i, 1, n(0)/nl);
              freesetinequ.Set(i, 2, n(1)/nl);
              freesetinequ.Set(i, 3, n(2)/nl);
              freesetinequ.Set(i, 4,
                               -(p1(0) * n(0) + p1(1) * n(1) + p1(2) * n(2)) / nl);
            }
        }
    }

  /*
  (*testout) << "Transformed freezone: " << endl;
  for (i = 1; i <= transfreezone.Size(); i++)
    (*testout) << transfreezone.Get(i) << " ";
  (*testout) << endl;
  */
}

int vnetrule :: ConvexFreeZone () const
{
  int i, j, k, fs;

  // (*mycout) << "Convex free zone...\n";
  
  int ret1=1;
  // int ret2=1;

  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      const DenseMatrix & freesetinequ = *freefaceinequ[fs-1];

      // const Array<int> & freeset = *freesets.Get(fs);
      const Array<twoint> & freesetedges = *freeedges[fs-1];
      // const Array<threeint> & freesetfaces = *freefaces.Get(fs);
      
      for (i = 1; i <= freesetedges.Size(); i++)
        {
          j = freesetedges[i-1].i1;    //triangle j with opposite point k
          k = freesetedges[i-1].i2;
          
          if ( freesetinequ.Get(j, 1) * transfreezone[k-1](0) +
               freesetinequ.Get(j, 2) * transfreezone[k-1](1) +
               freesetinequ.Get(j, 3) * transfreezone[k-1](2) +
               freesetinequ.Get(j, 4) > 0 )
            {
              ret1=0;
            }
        }
      
    }

  return ret1;
}


int vnetrule :: IsInFreeZone (const Point<3> & p) const
{
  int i, fs;
  char inthis;
  
  
  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      inthis = 1;
      Array<threeint> & freesetfaces = *freefaces[fs-1];
      DenseMatrix & freesetinequ = *freefaceinequ[fs-1];
      
      for (i = 1; i <= freesetfaces.Size() && inthis; i++)
        {
          if (freesetinequ.Get(i, 1) * p(0) + freesetinequ.Get(i, 2) * p(1) +
              freesetinequ.Get(i, 3) * p(2) + freesetinequ.Get(i, 4) > 0)
            inthis = 0;
        }
      
      if (inthis) return 1;
    }
  
  return 0;
}


int vnetrule :: IsTriangleInFreeZone (const Point<3> & p1, 
                                      const Point<3> & p2,
                                      const Point<3> & p3, 
                                      const Array<int> & pi, int newone)
{
  int fs;
  int infreeset, cannot = 0;


  ArrayMem<int,3> pfi(3), pfi2(3);

  // convert from local index to freeset index
  int i, j;
  for (i = 1; i <= 3; i++)
    {
      pfi[i-1] = 0;
      if (pi[i-1])
        {
          for (j = 1; j <= freezonepi.Size(); j++)
            if (freezonepi[j-1] == pi[i-1])
              pfi[i-1] = j;
        }
    }

  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      const Array<int> & freeseti = *freesets[fs-1];
      for (i = 1; i <= 3; i++)
        {
          pfi2[i-1] = 0;
          for (j = 1; j <= freeseti.Size(); j++)
            if (pfi[i-1] == freeseti[j-1])
              pfi2[i-1] = pfi[i-1];
        }

      infreeset = IsTriangleInFreeSet(p1, p2, p3, fs, pfi2, newone);
      if (infreeset == 1) return 1;
      if (infreeset == -1) cannot = -1;
    }
  
  return cannot;
}



int vnetrule :: IsTriangleInFreeSet (const Point<3> & p1, const Point<3> & p2,
                                     const Point<3> & p3, int fs,
                                     const Array<int> & pi, int newone)
{
  int i, ii;
  Vec<3> n;
  int allleft, allright;
  int hos1, hos2, hos3, os1, os2, os3;
  double hf, lam1, lam2, f, c1, c2, alpha;
  double v1n, v2n, h11, h12, h22, dflam1, dflam2;
  double lam1old, lam2old, fold;
  double hpx, hpy, hpz, v1x, v1y, v1z, v2x, v2y, v2z;
  int act1, act2, act3, it;
  int cntout;
  Array<int> activefaces;
  int isin;
  

  // MARK(triinfz);
  
  Array<threeint> & freesetfaces = *freefaces[fs-1];
  DenseMatrix & freesetinequ = *freefaceinequ[fs-1];
  

  int cnt = 0;
  for (i = 1; i <= 3; i++)
    if (pi[i-1]) cnt++;

  /*
  (*testout) << "trig in free set : " << p1 << " - " << p2 << " - " << p3 << endl;
  (*testout) << "common points: " << cnt << endl;
  */
  if (!newone)
    cnt = 0;

  if (cnt == 1)
    {
      // MARK(triinfz1);

      int upi = 0, lpiu = 0;
      for (i = 1; i <= 3; i++)
        if (pi[i-1])
          {
            upi = i;
            lpiu = pi[i-1];
          }

      Vec<3> v1, v2;
      switch (upi)
        {
        case 1:
          {
            v1 = p2 - p1;
            v2 = p3 - p1;
            break;
          }
        case 2:
          {
            v1 = p3 - p2;
            v2 = p1 - p2;
            break;
          }
        case 3:
          {
            v1 = p1 - p3;
            v2 = p2 - p3;
            break;
          }
        }

      v1 /= v1.Length();
      v2 /= v2.Length();
      Cross (v1, v2, n);
      n /= n.Length();

      //      (*testout) << "Test new: " << endl;
      for (i = 1; i <= freesetfaces.Size(); i++)
        {
          if ( (freesetfaces[i-1].i1 == lpiu) || 
               (freesetfaces[i-1].i2 == lpiu) ||
               (freesetfaces[i-1].i3 == lpiu) )
            {
              // freeface has point


              Vec<3> a (freesetinequ.Get(i, 1),
                       freesetinequ.Get(i, 2),
                       freesetinequ.Get(i, 3));
              
              //              if (1 - fabs (a * n) < 1e-8 ) 
              //                continue;

              Vec<3> an;
              Cross (a, n, an);
              double lan = an.Length();
              if (lan < 1e-10)
                continue;

              an /= lan;
              
              int out1 = (a * v1) > 0;
              int out2 = (a * v2) > 0;
              //              (*testout) << "out1, out2 = " << out1 << ", " << out2 << endl;
              if (out1 && out2)
                return 0;

              if (!out1 && !out2) 
                continue;


              //              if ( ( (an * v1) < 0) &&  ( (an * v2) < 0) )   // falsch !!!!
              //                an *= -1;

              // solve  an = lam1 v1 + lam2 v2
              double vii11 = v1 * v1;
              double vii12 = v1 * v2;
              double vii22 = v2 * v2;
              double det = vii11 * vii22 - vii12 * vii12;
              if ( fabs (det) < 1e-10 )
                continue;
              double rs1 = an * v1;
              double rs2 = an * v2;
              
              double lambda1 = rs1 * vii22 - rs2 * vii12;
              double lambda2 = rs2 * vii11 - rs1 * vii12;

              if (fabs (lambda1) > fabs (lambda2))
                {
                  if (lambda1 < 0)
                    an *= -1;
                }
              else
                {
                  if (lambda2 < 0)
                    an *= -1;
                }


              if (lambda1 * lambda2 < 0 && 0)
                {
                  if (fabs (lambda1) > 1e-14 && fabs (lambda2) > 1e-14)
                    {
                      //                      (*mycout) << "lambda1 lambda2 < 0" << endl;
                      (*testout) << "lambdai different" << endl;
                      (*testout) << "v1 = " << v1 << endl;
                      (*testout) << "v2 = " << v2 << endl;
                      (*testout) << "n = " << n << endl;
                      (*testout) << "a = " << a << endl;
                      (*testout) << "an = " << an << endl;
                      (*testout) << "a * v1 = " << (a * v1) << endl;
                      (*testout) << "a * v2 = " << (a * v2) << endl;
                      (*testout) << "an * v1 = " << (an * v1) << endl;
                      (*testout) << "an * v2 = " << (an * v2) << endl;
                      
                      (*testout) << "vii = " << vii11 << ", " << vii12 << ", " << vii22 << endl;
                      (*testout) << "lambdai = " << lambda1 << ", " << lambda2 << endl;
                      (*testout) << "rs = " << rs1 << ", " << rs2 << endl;
                      continue;
                    }
                }

              if (out1)
                v1 = an;
              else
                v2 = an;
            }
        }
      
      return 1;

      /*
      (*testout) << "overlap trig " << p1 << p2 << p3 << endl;
      (*testout) << "upi = " << upi << endl;
      (*testout) << "v1 = " << v1 << " v2 = " << v2 << endl;
      */

      switch (upi)
        {
        case 1:
          {
            v1 = p2 - p1;
            v2 = p3 - p1;
            break;
          }
        case 2:
          {
            v1 = p3 - p2;
            v2 = p1 - p2;
            break;
          }
        case 3:
          {
            v1 = p1 - p3;
            v2 = p2 - p3;
            break;
          }
        }

      v1 /= v1.Length();
      v2 /= v2.Length();
      Cross (v1, v2, n);
      n /= n.Length();

      //      (*testout) << "orig v1, v2 = " << v1 << ", " << v2 << endl;

      
      for (i = 1; i <= freesetfaces.Size(); i++)
        {
          if ( (freesetfaces[i-1].i1 == lpiu) || 
               (freesetfaces[i-1].i2 == lpiu) ||
               (freesetfaces[i-1].i3 == lpiu) )
            {
              /*
              (*testout) << "v1, v2, now = " << v1 << ", " << v2 << endl;

              // freeface has point
              (*testout) << "freesetface: "
                         << freesetfaces.Get(i).i1 << " "
                         << freesetfaces.Get(i).i2 << " "
                         << freesetfaces.Get(i).i3 << " ";
              */

              Vec<3> a (freesetinequ.Get(i, 1),
                       freesetinequ.Get(i, 2),
                       freesetinequ.Get(i, 3));
              //              (*testout) << "a = " <<  a << endl;


              Vec<3> an;
              Cross (a, n, an);
              double lan = an.Length();
              
              //              (*testout) << "an = " << an << endl;

              if (lan < 1e-10)
                continue;

              an /= lan;

              //              (*testout) << "a*v1 = " << (a*v1) << " a*v2 = " << (a*v2) << endl;
              
              int out1 = (a * v1) > 0;
              // int out2 = (a * v2) > 0;


              //              (*testout) << "out1, 2 = " << out1 << ", " << out2 << endl;

              
              double vii11 = v1 * v1;
              double vii12 = v1 * v2;
              double vii22 = v2 * v2;
              double det = vii11 * vii22 - vii12 * vii12;
              if ( fabs (det) < 1e-10 )
                continue;
              double rs1 = an * v1;
              double rs2 = an * v2;
              
              double lambda1 = rs1 * vii22 - rs2 * vii12;
              double lambda2 = rs2 * vii11 - rs1 * vii12;

              //              (*testout) << "lambda1, lambda2 = " << lambda1 << ", " << lambda2 << endl;


              if (fabs (lambda1) > fabs (lambda2))
                {
                  if (lambda1 < 0)
                    an *= -1;
                }
              else
                {
                  if (lambda2 < 0)
                    an *= -1;
                }


              if (lambda1 * lambda2 < 0)
                {
                  if (fabs (lambda1) > 1e-14 && fabs (lambda2) > 1e-14)
                    {
                      //                      (*mycout) << "lambda1 lambda2 < 0" << endl;
                      (*testout) << "lambdai different" << endl;
                      (*testout) << "v1 = " << v1 << endl;
                      (*testout) << "v2 = " << v2 << endl;
                      (*testout) << "n = " << n << endl;
                      (*testout) << "a = " << a << endl;
                      (*testout) << "an = " << an << endl;
                      (*testout) << "a * v1 = " << (a * v1) << endl;
                      (*testout) << "a * v2 = " << (a * v2) << endl;
                      (*testout) << "an * v1 = " << (an * v1) << endl;
                      (*testout) << "an * v2 = " << (an * v2) << endl;
                      
                      (*testout) << "vii = " << vii11 << ", " << vii12 << ", " << vii22 << endl;
                      (*testout) << "lambdai = " << lambda1 << ", " << lambda2 << endl;
                      (*testout) << "rs = " << rs1 << ", " << rs2 << endl;
                      continue;
                    }
                }

              if (out1)
                v1 = an;
              else
                v2 = an;



            }
        }

      return 1;
    }



  if (cnt == 2)
    {
      //      (*testout) << "tripoitns: " << p1 << " " << p2 << " " << p3 << endl;

      // MARK(triinfz2);

      int pi1 = 0, pi2 = 0, pi3 = 0;
      Vec<3> a1, a2;  // outer normals
      Vec<3> trivec;  // vector from common edge to third point of triangle
      for (i = 1; i <= 3; i++)
        if (pi[i-1])
          {
            pi2 = pi1;
            pi1 = pi[i-1];
          }
        else
          pi3 = i;

      switch (pi3)
        {
        case 1: trivec = (p1 - p2); break;
        case 2: trivec = (p2 - p3); break;
        case 3: trivec = (p3 - p2); break;
        }

      Array<int> lpi(freezonepi.Size());
      for (i = 1; i <= lpi.Size(); i++)
        lpi[i-1] = 0;
      lpi[pi1-1] = 1;
      lpi[pi2-1] = 1;
      
      int ff1 = 0, ff2 = 0;
      for (i = 1; i <= freesetfaces.Size(); i++)
        {
          if (lpi[freesetfaces[i-1].i1-1] + 
              lpi[freesetfaces[i-1].i2-1] + 
              lpi[freesetfaces[i-1].i3-1] == 2)
            {
              ff2 = ff1;
              ff1 = i;
            }
        }

      if (ff2 == 0)
        return 1;

      a1 = Vec<3> (freesetinequ.Get(ff1, 1),
                  freesetinequ.Get(ff1, 2),
                  freesetinequ.Get(ff1, 3));
      a2 = Vec<3> (freesetinequ.Get(ff2, 1),
                  freesetinequ.Get(ff2, 2),
                  freesetinequ.Get(ff2, 3));

      if ( ( (a1 * trivec) > 0) || ( (a2 * trivec) > 0))
        return 0;

      return 1;
    }


  if (cnt == 3)
    {
      // MARK(triinfz3);  

      Array<int> lpi(freezonepi.Size());
      for (i = 1; i <= lpi.Size(); i++)
        lpi[i-1] = 0;

      for (i = 1; i <= 3; i++)
        lpi[pi[i-1]-1] = 1;
      
      for (i = 1; i <= freesetfaces.Size(); i++)
        {
          if (lpi[freesetfaces[i-1].i1-1] + 
              lpi[freesetfaces[i-1].i2-1] + 
              lpi[freesetfaces[i-1].i3-1] == 3)
            {
              return 0;
            }
        }
      return 1;
    }

  // MARK(triinfz0);  

  
  os1 = os2 = os3 = 0;
  activefaces.SetSize(0);

  // is point inside ?

  for (i = 1; i <= freesetfaces.Size(); i++)
    {
      hos1 = freesetinequ.Get(i, 1) * p1(0) +
        freesetinequ.Get(i, 2) * p1(1) +
        freesetinequ.Get(i, 3) * p1(2) +
        freesetinequ.Get(i, 4) > -1E-5;
      
      hos2 = freesetinequ.Get(i, 1) * p2(0) +
        freesetinequ.Get(i, 2) * p2(1) +
        freesetinequ.Get(i, 3) * p2(2) +
        freesetinequ.Get(i, 4) > -1E-5;
      
      hos3 = freesetinequ.Get(i, 1) * p3(0) +
        freesetinequ.Get(i, 2) * p3(1) +
        freesetinequ.Get(i, 3) * p3(2) +
        freesetinequ.Get(i, 4) > -1E-5;
      
      if (hos1 && hos2 && hos3) return 0;
      
      if (hos1) os1 = 1;
      if (hos2) os2 = 1;
      if (hos3) os3 = 1;
      
      if (hos1 || hos2 || hos3) activefaces.Append (i);
    }
  
  if (!os1 || !os2 || !os3) return 1;

  v1x = p2(0) - p1(0);
  v1y = p2(1) - p1(1);
  v1z = p2(2) - p1(2);

  v2x = p3(0) - p1(0);
  v2y = p3(1) - p1(1);
  v2z = p3(2) - p1(2);

  n(0) = v1y * v2z - v1z * v2y;
  n(1) = v1z * v2x - v1x * v2z;
  n(2) = v1x * v2y - v1y * v2x;
  n /= n.Length();

  allleft = allright = 1;
  for (i = 1; i <= transfreezone.Size() && (allleft || allright); i++)
    {
      const Point<3> & p = transfreezone[i-1];
      float scal = (p(0) - p1(0)) * n(0) +
        (p(1) - p1(1)) * n(1) +
        (p(2) - p1(2)) * n(2);

      if ( scal >  1E-8 ) allleft = 0;
      if ( scal < -1E-8 ) allright = 0;
    }

  if (allleft || allright) return 0;


  lam1old = lam2old = lam1 = lam2 = 1.0 / 3.0;


  //  testout << endl << endl << "Start minimizing" << endl;

  it = 0;
  int minit;
  minit = 1000;
  fold = 1E10;

  

  while (1)
    {
      it++;

      if (it > 1000) return -1;

      if (lam1 < 0) lam1 = 0;
      if (lam2 < 0) lam2 = 0;
      if (lam1 + lam2 > 1) lam1 = 1 - lam2;

      if (it > minit)
        {
          (*testout) << "it = " << it << endl;
          (*testout) << "lam1/2 = " << lam1 << "  " << lam2 << endl;
        }

      hpx = p1(0) + lam1 * v1x + lam2 * v2x;
      hpy = p1(1) + lam1 * v1y + lam2 * v2y;
      hpz = p1(2) + lam1 * v1z + lam2 * v2z;

      f = 0;

      h11 = h12 = h22 = dflam1 = dflam2 = 0;
      cntout = 0;

      isin = 1;

      for (i = 1; i <= activefaces.Size(); i++)
        {
          ii = activefaces[i-1];

          hf = freesetinequ.Get(ii, 1) * hpx +
            freesetinequ.Get(ii, 2) * hpy +
            freesetinequ.Get(ii, 3) * hpz +
            freesetinequ.Get(ii, 4);

          if (hf > -1E-7) isin = 0;

          hf += 1E-4;
          if (hf > 0)
            {
              f += hf * hf;

              v1n = freesetinequ.Get(ii, 1) * v1x +
                freesetinequ.Get(ii, 2) * v1y +
                freesetinequ.Get(ii, 3) * v1z;
              v2n = freesetinequ.Get(ii, 1) * v2x +
                freesetinequ.Get(ii, 2) * v2y +
                freesetinequ.Get(ii, 3) * v2z;

              h11 += 2 * v1n * v1n;
              h12 += 2 * v1n * v2n;
              h22 += 2 * v2n * v2n;
              dflam1 += 2 * hf * v1n;
              dflam2 += 2 * hf * v2n;
              cntout++;
            }
        }

      if (isin) return 1;

      if (it > minit)
        {
          (*testout) << "f = " << f
                     << "  dfdlam = " << dflam1 << "  " << dflam2 << endl;
          (*testout) << "h = " << h11 << "  " << h12 << "  " << h22 << endl;
          (*testout) << "active: " << cntout << endl;
          (*testout) << "lam1-lam1old = " << (lam1 - lam1old) << endl;
          (*testout) << "lam2-lam2old = " << (lam2 - lam2old) << endl;
        }


      if (f >= fold)
        {
          lam1 = 0.100000000000000 * lam1 + 0.9000000000000000 * lam1old;
          lam2 = 0.100000000000000 * lam2 + 0.9000000000000000 * lam2old;
        }
      else
        {
          lam1old = lam1;
          lam2old = lam2;
          fold = f;


          if (f < 1E-9) return 1;

          h11 += 1E-10;
          h22 += 1E-10;
          c1 = - ( h22 * dflam1 - h12 * dflam2) / (h11 * h22 - h12 * h12);
          c2 = - (-h12 * dflam1 + h11 * dflam2) / (h11 * h22 - h12 * h12);
          alpha = 1;


          if (it > minit)
            (*testout) << "c1/2 = " << c1 << "  " << c2 << endl;

          act1 = lam1 <= 1E-6 && c1 <= 0;
          act2 = lam2 <= 1E-6 && c2 <= 0;
          act3 = lam1 + lam2 >= 1 - 1E-6 && c1 + c2 >= 0;

          if (it > minit)
            (*testout) << "act1,2,3 = " << act1 << act2 << act3 << endl;

          if ( (act1 && act2) || (act1 && act3) || (act2 && act3) ) return 0;

          if (act1)
            {
              c1 = 0;
              c2 = - dflam2 / h22;
            }

          if (act2)
            {
              c1 = - dflam1 / h11;
              c2 = 0;
            }

          if (act3)
            {
              c1 = - (dflam1 - dflam2) / (h11 + h22 - 2 * h12);
              c2 = -c1;
            }

          if (it > minit)
            (*testout) << "c1/2 now = " << c1 << "  " << c2 << endl;


          if (f > 100 * sqrt (sqr (c1) + sqr (c2))) return 0;


          if (lam1 + alpha * c1 < 0 && !act1)
            alpha = -lam1 / c1;
          if (lam2 + alpha * c2 < 0 && !act2)
            alpha = -lam2 / c2;
          if (lam1 + lam2 + alpha * (c1 + c2) > 1 && !act3)
            alpha = (1 - lam1 - lam2) / (c1 + c2);

          if (it > minit)
            (*testout) << "alpha = " << alpha << endl;

          lam1 += alpha * c1;
          lam2 += alpha * c2;
        }
    }
}




int vnetrule :: IsQuadInFreeZone (const Point<3> & p1, 
                                  const Point<3> & p2,
                                  const Point<3> & p3, 
                                  const Point<3> & p4, 
                                  const Array<int> & pi, int newone)
{
  int fs;
  int infreeset, cannot = 0;


  ArrayMem<int,4> pfi(4), pfi2(4);

  // convert from local index to freeset index
  int i, j;
  for (i = 1; i <= 4; i++)
    {
      pfi[i-1] = 0;
      if (pi[i-1])
        {
          for (j = 1; j <= freezonepi.Size(); j++)
            if (freezonepi[j-1] == pi[i-1])
              pfi[i-1] = j;
        }
    }

  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      const Array<int> & freeseti = *freesets[fs-1];
      for (i = 1; i <= 4; i++)
        {
          pfi2[i-1] = 0;
          for (j = 1; j <= freeseti.Size(); j++)
            if (pfi[i-1] == freeseti[j-1])
              pfi2[i-1] = pfi[i-1];
        }

      infreeset = IsQuadInFreeSet(p1, p2, p3, p4, fs, pfi2, newone);
      if (infreeset == 1) return 1;
      if (infreeset == -1) cannot = -1;
    }
  
  return cannot;
}


int vnetrule :: IsQuadInFreeSet (const Point<3> & p1, const Point<3> & p2,
                                 const Point<3> & p3, const Point<3> & p4, 
                                 int fs, const Array<int> & pi, int newone)
{
  int i;
  
  int cnt = 0;
  for (i = 1; i <= 4; i++)
    if (pi[i-1]) cnt++;
  
  /*
  (*testout) << "test quad in freeset: " << p1 << " - " << p2 << " - " << p3 << " - " << p4 << endl;
  (*testout) << "pi = ";
  for (i = 1; i <= pi.Size(); i++)
    (*testout) << pi.Get(i) << " ";
  (*testout) << endl;
  (*testout) << "cnt = " << cnt  << endl;
  */
  if (cnt == 4)
    {
      return 1;
    }

  if (cnt == 3)
    {
      return 1;
    }

  ArrayMem<int,3> pi3(3);
  int res;

  pi3[0] = pi[0];
  pi3[1] = pi[1];
  pi3[2] = pi[2];
  res = IsTriangleInFreeSet (p1, p2, p3, fs, pi3, newone);
  if (res) return res;


  pi3[0] = pi[1];
  pi3[1] = pi[2];
  pi3[2] = pi[3];
  res = IsTriangleInFreeSet (p2, p3, p4, fs, pi3, newone);
  if (res) return res;

  pi3[0] = pi[2];
  pi3[1] = pi[3];
  pi3[2] = pi[0];
  res = IsTriangleInFreeSet (p3, p4, p1, fs, pi3, newone);
  if (res) return res;

  pi3[0] = pi[3];
  pi3[1] = pi[0];
  pi3[2] = pi[1];
  res = IsTriangleInFreeSet (p4, p1, p2, fs, pi3, newone);
  return res;
}












float vnetrule :: CalcPointDist (RulePointIndex pi, const Point<3> & p) const
{
  float dx = p(0) - points[pi](0);
  float dy = p(1) - points[pi](1);
  float dz = p(2) - points[pi](2);
  
  return tolerances[pi] * (dx * dx + dy * dy + dz * dz);
}


int vnetrule :: TestOk () const
{
  Array<int, RulePointIndex> cntpused(points.Size());
  Array<RulePointIndex> edge1, edge2;
  Array<int> delf(faces.Size());
  int i, j, k;
  RulePointIndex pi1, pi2;
  int found;

  cntpused = 0;
  for (i = 1; i <= faces.Size(); i++)
    delf[i-1] = 0;
  for (i = 1; i <= delfaces.Size(); i++)
    delf[delfaces[i-1]-1] = 1;


  for (i = 1; i <= faces.Size(); i++)
    if (delf[i-1] || i > noldf)
      for (j = 1; j <= faces[i-1].GetNP(); j++)
        cntpused[faces[i-1].PNum(j)]++;

  for (auto pi : cntpused.Range())
    if (cntpused[pi] > 0 && cntpused[pi] < 2)
      {
        return 0;
      }


  //  (*testout) << endl;
  for (i = 1; i <= faces.Size(); i++)
    {
      //      (*testout) << "face " << i << endl;
      for (j = 1; j <= faces[i-1].GetNP(); j++)
        {
          pi1.Invalidate(); pi2.Invalidate();
          if (delf[i-1])
            {
              pi1 = faces[i-1].PNumMod(j);
              pi2 = faces[i-1].PNumMod(j+1);
            }
          if (i > noldf)
            {
              pi1 = faces[i-1].PNumMod(j+1);
              pi2 = faces[i-1].PNumMod(j);
            }

          found = 0;
          if (pi1.IsValid())
            {
              for (k = 1; k <= edge1.Size(); k++)
                if (edge1[k-1] == pi1 && edge2[k-1] == pi2)
                  {
                    found = 1;
                    edge1.DeleteElement(k-1);
                    edge2.DeleteElement(k-1);
                    k--;
                    //              (*testout) << "Del edge " << pi1 << "-" << pi2 << endl;
                  }
              if (!found)
                {
                  edge1.Append (pi2);
                  edge2.Append (pi1);
                  //              (*testout) << "Add edge " << pi1 << "-" << pi2 << endl;
                }
            }
        }
    }


  if (edge1.Size() > 0)
    {
      return 0;
    }

  /*
    cntpused.SetSize(freezone.Size());
    for (i = 1; i <= cntpused.Size(); i++)
    cntpused[i] = 0;

    for (i = 1; i <= freefaces.Size(); i++)
    {
    cntpused[freefaces[i].i1]++;
    cntpused[freefaces[i].i2]++;
    cntpused[freefaces[i].i3]++;
    }

    for (i = 1; i <= cntpused.Size(); i++)
    if (cntpused[i] < 3)
    {
    (*mycout) << "Fall 3" << endl;
    return 0;
    }



    for (i = 1; i <= freefaces.Size(); i++)
    {
    for (j = 1; j <= 3; j++)
    {
    if (j == 1)
    {
    pi1 = freefaces[i].i1;
    pi2 = freefaces[i].i2;
    }
    if (j == 2)
    {
    pi1 = freefaces[i].i2;
    pi2 = freefaces[i].i3;
    }
    if (j == 3)
    {
    pi1 = freefaces[i].i3;
    pi2 = freefaces[i].i1;
    }

    found = 0;
    for (k = 1; k <= edge1.Size(); k++)
    if (edge1[k] == pi1 && edge2[k] == pi2)
    {
    found = 1;
    edge1.DeleteElement(k);
    edge2.DeleteElement(k);
    k--;
    }

    if (!found)
    {
    edge1.Append (pi2);
    edge2.Append (pi1);
    }
    }
    }

    if (edge1.Size() > 0)
    {
    (*mycout) << "Fall 4" << endl;
    return 0;
    }
    */
  return 1;
}


int vnetrule :: IsDelFace (int fn) const
{
  int i;
  for (i = 1; i <= GetNDelF(); i++)
    if (GetDelFace(i) == fn) return 1;
  return 0;
}

}

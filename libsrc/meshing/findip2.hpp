#ifndef NETGEN_FINDIP2_HPP
#define NETGEN_FINDIP2_HPP

// find inner point

namespace netgen
{

template <typename POINTArray, typename FACEArray>
inline int FindInnerPoint2 (POINTArray & points,
                            FACEArray & faces,
                            Point<3> & p)
{
  static Timer timer("FindInnerPoint2");
  RegionTimer reg (timer);

  Array<Vec<3>> a;
  Array<double> c;
  Mat<3> m, inv;
  Vec<3> rs, x, pmin;

  int nf = faces.Size();

  a.SetSize (nf);
  c.SetSize (nf);

  for (int i = 0; i < nf; i++)
    {
      Point<3> p1 = points.Get(faces[i][0]);
      a[i] = Cross (points.Get(faces[i][1]) - p1,
                    points.Get(faces[i][2]) - p1);
      a[i] /= a[i].Length();
      c[i] = - (a[i](0) * p1(0) + a[i](1) * p1(1) + a[i](2) * p1(2));
    }


  x = 0;
  
  
  double hmax = 0;
  for (int i = 0; i < nf; i++)
    {
      const Element2dRef & el = faces[i];
      for (int j = 1; j <= 3; j++)
        {
          double hi = Dist (points.Get(el.PNumMod(j)),
                            points.Get(el.PNumMod(j+1)));
          if (hi > hmax) hmax = hi;
        }
    }

  double fmin = 0;

  for (int i1 = 1; i1 <= nf; i1++)
    for (int i2 = i1+1; i2 <= nf; i2++)
      for (int i3 = i2+1; i3 <= nf; i3++)
        for (int i4 = i3+1; i4 <= nf; i4++)
          {
            m(0, 0) = a[i1-1](0) - a[i2-1](0);
            m(0, 1) = a[i1-1](1) - a[i2-1](1);
            m(0, 2) = a[i1-1](2) - a[i2-1](2);
            rs(0) = c[i2-1] - c[i1-1];

            m(1, 0) = a[i1-1](0) - a[i3-1](0);
            m(1, 1) = a[i1-1](1) - a[i3-1](1);
            m(1, 2) = a[i1-1](2) - a[i3-1](2);
            rs(1) = c[i3-1] - c[i1-1];

            m(2, 0) = a[i1-1](0) - a[i4-1](0);
            m(2, 1) = a[i1-1](1) - a[i4-1](1);
            m(2, 2) = a[i1-1](2) - a[i4-1](2);
            rs(2) = c[i4-1] - c[i1-1];


            if (fabs (Det (m)) > 1e-10)
              {
                CalcInverse (m, inv);
                x = inv * rs;

                double f = -1e10;
                for (int i = 0; i < nf; i++)
                  {
                    double hd = 
                      x(0) * a[i](0) + x(1) * a[i](1) + x(2) * a[i](2) + c[i];
                    if (hd > f) f = hd;
                    if (hd > fmin) break;
                  }

                if (f < fmin)
                  {
                    fmin = f;
                    pmin = x;
                  }
              }
          }

  p = Point<3> (pmin(0), pmin(1), pmin(2));
  (*testout) << "fmin = " << fmin << endl;
  return (fmin < -1e-3 * hmax);
}

} // namespace netgen
#endif // NETGEN_FINDIP2_HPP

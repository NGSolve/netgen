#include <mystdlib.h>

#include "meshing.hpp"
#ifdef SOLIDGEOM
#include <csg.hpp>
#endif
#include <opti.hpp>
#include <core/array.hpp>
#include <core/taskmanager.hpp>


namespace netgen
{
  using namespace ngcore;
  

  double MinFunctionSum :: Func (const Vector & x) const
  {
    double retval = 0;
    for(int i=0; i<functions.Size(); i++)
      retval += functions[i]->Func(x);
      
    return retval;
  }

  void MinFunctionSum :: Grad (const Vector & x, Vector & g) const
  {
    g = 0.;
    VectorMem<3> gi;
    for(int i=0; i<functions.Size(); i++)
      {
        functions[i]->Grad(x,gi);
        for(int j=0; j<g.Size(); j++)
          g[j] += gi[j];
      }
  }
      

  double MinFunctionSum :: FuncGrad (const Vector & x, Vector & g) const
  {
    double retval = 0;
    g = 0.;
    VectorMem<3> gi;
    for(int i=0; i<functions.Size(); i++)
      {
        retval += functions[i]->FuncGrad(x,gi);
        for(int j=0; j<g.Size(); j++)
          g[j] += gi[j];
      }
    return retval;
  }

  double MinFunctionSum :: FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    double retval = 0;
    deriv = 0.;
    double derivi;
    for(int i=0; i<functions.Size(); i++)
      {
        retval += functions[i]->FuncDeriv(x,dir,derivi);
        deriv += derivi;
      }
    return retval;
  }

  double MinFunctionSum :: GradStopping (const Vector & x) const
  {
    double minfs(0), mini;
    for(int i=0; i<functions.Size(); i++)
      {
        mini = functions[i]->GradStopping(x);
        if(i==0 || mini < minfs)
          minfs = mini;
      }
    return minfs;
  }


  void MinFunctionSum :: AddFunction(MinFunction & fun)
  {
    functions.Append(&fun);
  }
  
  const MinFunction & MinFunctionSum :: Function(int i) const
  {
    return *functions[i];
  }
  MinFunction & MinFunctionSum :: Function(int i)
  {
    return *functions[i];
  }

  PointFunction1 :: PointFunction1 (Mesh::T_POINTS & apoints, 
                                    const Array<PointIndices<3>> & afaces,
                                    const MeshingParameters & amp,
                                    double ah)
    : points(apoints), faces(afaces), mp(amp)
  {
    h = ah;
  }
  

  double PointFunction1 :: Func (const Vector & vp) const
  {
    double badness = 0;
    Point<3> pp(vp(0), vp(1), vp(2));

    for (int j = 0; j < faces.Size(); j++)
      {
        const PointIndices<3> & el = faces[j];

        double bad = CalcTetBadness (points[el[0]],
                                     points[el[2]],
                                     points[el[1]], 
                                     pp, 0, mp);
        badness += bad;
      }
 
    return badness;
  }


  double PointFunction1 :: 
  FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    VectorMem<3> hx;
    const double eps = 1e-6;

    double dirlen = dir.L2Norm();
    if (dirlen < 1e-14)
      {
        deriv = 0;
        return Func(x);
      }

    hx.Set(1, x);
    hx.Add(eps * h / dirlen, dir);
    double fr = Func (hx);
    hx.Set(1, x);
    hx.Add(-eps * h / dirlen, dir);
    double fl = Func (hx);

    deriv = (fr - fl) / (2 * eps * h) * dirlen;

    return Func(x);
  }


  double PointFunction1 :: FuncGrad (const Vector & x, Vector & g) const
  {
    VectorMem<3> hx;
    double eps = 1e-6;

    hx = x;
    for (int i = 0; i < 3; i++)
      {
        hx(i) = x(i) + eps * h;
        double fr = Func (hx);
        hx(i) = x(i) - eps * h;
        double fl = Func (hx);
        hx(i) = x(i);

        g(i) = (fr - fl) / (2 * eps * h);
      }

    return Func(x);
  }

  double PointFunction1 :: GradStopping (const Vector & x) const
  {
    double f = Func(x);
    return 1e-8 * f * f;
  }


  /* Cheap Functional depending of inner point inside triangular surface */

  // is it used ????
  class CheapPointFunction1 : public MinFunction
  {
    Mesh::T_POINTS & points;
    const Array<PointIndices<3>> & faces;
    DenseMatrix m;
    double h;
  public:
    CheapPointFunction1 (Mesh::T_POINTS & apoints, 
                         const Array<PointIndices<3>> & afaces,
                         double ah);
  
    virtual double Func (const Vector & x) const;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
  };

  CheapPointFunction1 :: CheapPointFunction1 (Mesh::T_POINTS & apoints, 
                                              const Array<PointIndices<3>> & afaces,
                                              double ah)
    : points(apoints), faces(afaces)
  {
    h = ah;
  

    int nf = faces.Size();

    m.SetSize (nf, 4);
  
    for (int i = 1; i <= nf; i++)
      {
        const Point<3> & p1 = points[faces[i-1][0]];
        const Point<3> & p2 = points[faces[i-1][1]];
        const Point<3> & p3 = points[faces[i-1][2]];
        Vec<3> v1 (p1, p2);
        Vec<3> v2 (p1, p3);
        Vec<3> n;
        Cross (v1, v2, n);
        n /= n.Length();

        m.Elem(i, 1) = n(0);
        m.Elem(i, 2) = n(1);
        m.Elem(i, 3) = n(2);
        m.Elem(i, 4) = - (n(0) * p1(0) + n(1) * p1(1) + n(2) * p1(2));
      } 
  }
  
  double CheapPointFunction1 :: Func (const Vector & vp) const
  {

    /*
      int j;
      double badness = 0;
      Point<3> pp(vp.Get(1), vp.Get(2), vp.Get(3));

      for (j = 1; j <= faces.Size(); j++)
      {
      const IVec<3> & el = faces.Get(j);

      double bad = CalcTetBadness (points.Get(el[0]), 
      points.Get(el[2]), 
      points.Get(el[1]), 
      pp, 0);
      badness += bad;
      }
    */

    double badness = 0;
    VectorMem<4> hv;
    Vector res(m.Height());

    for (int i = 0;i < 3; i++)
      hv(i) = vp(i);
    hv(3) = 1;
    m.Mult (hv, res);

    for (int i = 0; i < res.Size(); i++)
      {
        if (res(i) < 1e-10)
          badness += 1e24;
        else
          badness += 1 / res(i);
      }
 
    return badness;
  }


  double CheapPointFunction1 :: FuncGrad (const Vector & x, Vector & g) const
  {
    VectorMem<3> hx;
    double eps = 1e-6;

    hx = x;
    for (int i = 0; i < 3; i++)
      {
        hx(i) = x(i) + eps * h;
        double fr = Func (hx);
        hx(i) = x(i) - eps * h;
        double fl = Func (hx);
        hx(i) = x(i);

        g(i) = (fr - fl) / (2 * eps * h);
      }

    return Func(x);
  }














  /* ************* PointFunction **************************** */


  class PointFunction 
  {
  public:
    Mesh::T_POINTS & points;
    const Array<Element, ElementIndex> & elements;
    Table<ElementIndex, PointIndex> &elementsonpoint;
    bool own_elementsonpoint;
    const MeshingParameters & mp;
    PointIndex actpind;
    double h;
  
  public:
    PointFunction (Mesh & mesh, const MeshingParameters & amp);
    PointFunction (const PointFunction & pf);
    virtual ~PointFunction () { if(own_elementsonpoint) delete &elementsonpoint; }
    virtual void SetPointIndex (PointIndex aactpind);
    void SetLocalH (double ah) { h = ah; }
    double GetLocalH () const { return h; }
    const Table<ElementIndex, PointIndex> & GetPointToElementTable() { return elementsonpoint; };
    virtual double PointFunctionValue (const Point<3> & pp) const;
    virtual double PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const;
    virtual double PointFunctionValueDeriv (const Point<3> & pp, const Vec<3> & dir, double & deriv) const;

    int MovePointToInner ();
  };


  PointFunction :: PointFunction (const PointFunction & pf)
    : points(pf.points), elements(pf.elements), elementsonpoint(pf.elementsonpoint), own_elementsonpoint(false), mp(pf.mp)
  { }

  PointFunction :: PointFunction (Mesh & mesh, const MeshingParameters & amp)
    : points(mesh.Points()), elements(mesh.VolumeElements()), elementsonpoint(* new Table<ElementIndex,PointIndex>()), own_elementsonpoint(true), mp(amp)
  {
    static Timer tim("PointFunction - build elementsonpoint table"); RegionTimer reg(tim);

    Array<bool, PointIndex> non_tet_points(points.Size());
    non_tet_points = false;
    // Don't optimize if point is adjacent to a non-tet element
    ParallelForRange(elements.Range(), [&] (auto myrange)
        {
          for(auto ei : myrange)
            {
              const auto & el = elements[ei];
              if(el.NP()!=4)
                for(auto pi : el.PNums())
                  non_tet_points[pi] = true;
            }
       });

    elementsonpoint = ngcore::CreateSortedTable<ElementIndex, PointIndex>( elements.Range(),
               [&](auto & table, ElementIndex ei)
               {
                 const auto & el = elements[ei];

                 if(el.NP()!=4 || (mp.only3D_domain_nr && mp.only3D_domain_nr != el.GetIndex()) )
                   return;

                 for (PointIndex pi : el.PNums())
                   if(!non_tet_points[pi])
                     table.Add (pi, ei);
               }, points.Size());
  }

  void PointFunction :: SetPointIndex (PointIndex aactpind)
  {
    actpind = aactpind; 
  }  

  double PointFunction :: PointFunctionValue (const Point<3> & pp) const
  {
    double badness;
    Point<3> hp;

    badness = 0;

    hp = points[actpind];
    points[actpind] = Point<3> (pp);

    for (auto ei : elementsonpoint[actpind])
      {
        const Element & el = elements[ei];
        badness += CalcTetBadness (points[el[0]], points[el[1]], 
                                   points[el[2]], points[el[3]], -1, mp);
      }
  
    points[actpind] = Point<3> (hp); 
    return badness;
  }


  double PointFunction :: PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const
  {
    double f = 0;

    Point<3> hp = points[actpind];
    Vec<3> vgradi, vgrad(0,0,0);
    points[actpind] = Point<3> (pp);

    for (auto ei : elementsonpoint[actpind])
      {
        const Element & el = elements[ei];
        for (int k = 0; k < 4; k++)
          if (el[k] == actpind)
            {
              f += CalcTetBadnessGrad (points[el[0]], points[el[1]], 
                                       points[el[2]], points[el[3]], 
                                       -1, k+1, vgradi, mp);

              vgrad += vgradi;
            }
      }

    points[actpind] = Point<3> (hp); 

    grad = vgrad;
    return f;
  }


  double PointFunction :: PointFunctionValueDeriv (const Point<3> & pp, const Vec<3> & dir,
                                                   double & deriv) const
  {
    Vec<3> vgradi, vgrad(0,0,0);

    Point<3> hp = points[actpind];
    points[actpind] = pp;
    double f = 0;

    for (auto ei : elementsonpoint[actpind])
      {
        const Element & el = elements[ei];

        for (int k = 1; k <= 4; k++)
          if (el.PNum(k) == actpind)
            {
              f += CalcTetBadnessGrad (points[el[0]], 
                                       points[el[1]], 
                                       points[el[2]], 
                                       points[el[3]], -1, k, vgradi, mp);

              vgrad += vgradi;
            }
      }

    points[actpind] = Point<3> (hp); 
    deriv = dir * vgrad;
    return f;
  }

  int PointFunction :: MovePointToInner ()
  {
    // try point movement 
    Array<Element2d> faces;
  
    for (auto ei : elementsonpoint[actpind])
      {
        const Element & el = elements[ei];
      
        for (int k = 1; k <= 4; k++)
          if (el.PNum(k) == actpind)
            {
              Element2d face(TRIG);
              el.GetFace (k, face);
              Swap (face[1], face[2]);
              faces.Append (face);
            }
      }
  
    Point<3> hp;
    int hi = FindInnerPoint (points, faces, hp);
    if (hi)
      {
        // cout << "inner point found" << endl;
        points[actpind] = Point<3> (hp);
      }
    else
      ;
    //      cout << "no inner point found" << endl;

    /*
    Point<3> hp2;
    int hi2 = FindInnerPoint (points, faces, hp2);
    if (hi2)
      {
        cout << "new: inner point found" << endl;
      }
    else
      cout << "new: no inner point found" << endl;
  
    (*testout) << "hi(orig) = " << hi << ", hi(new) = " << hi2;
    if (hi != hi2) (*testout) << "hi different" << endl;
    */

    return hi;
  }






  class CheapPointFunction : public PointFunction
  {
    DenseMatrix m;
  public:
    CheapPointFunction (Mesh & mesh, const MeshingParameters & amp);
    virtual void SetPointIndex (PointIndex aactpind);
    virtual double PointFunctionValue (const Point<3> & pp) const;
    virtual double PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const;
  };


  CheapPointFunction :: CheapPointFunction (Mesh & mesh, const MeshingParameters & amp)
    : PointFunction (mesh, amp)
  {
    ;
  }


  void CheapPointFunction :: SetPointIndex (PointIndex aactpind)
  {
    actpind = aactpind; 

    int ne = elementsonpoint[actpind].Size();
    PointIndex pi1, pi2, pi3;

    m.SetSize (ne, 4);

    for (int i = 0; i < ne; i++)
      {
        pi1 = PointIndex::INVALID;
        pi2 = PointIndex::INVALID;
        pi3 = PointIndex::INVALID;

        const Element & el = elements[elementsonpoint[actpind][i]];
        for (int j = 1; j <= 4; j++)
          if (el.PNum(j) != actpind)
            {
              pi3 = pi2;
              pi2 = pi1;
              pi1 = el.PNum(j);
            }

        const Point<3> & p1 = points[pi1];
        Vec<3> v1 (p1, points[pi2]);
        Vec<3> v2 (p1, points[pi3]);
        Vec<3> n;
        Cross (v1, v2, n);
        n /= n.Length();

        Vec<3> v (p1, points[actpind]);
        double c = v * n;
      
        if (c < 0)
          n *= -1;    
      
        // n is inner normal

        m.Elem(i+1, 1) = n(0);
        m.Elem(i+1, 2) = n(1);
        m.Elem(i+1, 3) = n(2);
        m.Elem(i+1, 4) = - (n(0) * p1(0) + n(1) * p1(1) + n(2) * p1(2));
      }
  }

  double CheapPointFunction :: PointFunctionValue (const Point<3> & pp) const
  {
    VectorMem<4> p4;
    Vector di;
    int n = m.Height();

    p4(0) = pp(0);
    p4(1) = pp(1);
    p4(2) = pp(2);
    p4(3) = 1;

    di.SetSize (n);
    m.Mult (p4, di);
  
    double sum = 0;
    for (int i = 0; i < n; i++)
      {
        if (di(i) > 0)
          sum += 1 / di(i);
        else
          return 1e16;
      }
    return sum;
  }




  double CheapPointFunction :: PointFunctionValueGrad (const Point<3> & pp, Vec<3> & grad) const
  {
    VectorMem<4> p4;
    Vector di;

    int n = m.Height();

    p4(0) = pp(0);
    p4(1) = pp(1);
    p4(2) = pp(2);
    p4(3) = 1;

    di.SetSize (n);
    m.Mult (p4, di);
  
    double sum = 0;
    grad = 0;
    for (int i = 0; i < n; i++)
      {
        if (di(i) > 0)
          {
            double idi = 1 / di(i);
            sum += idi;
            grad(0) -= idi * idi * m(i, 0);
            grad(1) -= idi * idi * m(i, 1);
            grad(2) -= idi * idi * m(i, 2);
          }
        else
          {
            return 1e16;
          }
      }
    return sum;
  }








  class Opti3FreeMinFunction : public MinFunction
  { 
    const PointFunction & pf;
    Point<3> sp1;
  
  public:
    Opti3FreeMinFunction (const PointFunction & apf);
    void SetPoint (const Point<3> & asp1) { sp1 = asp1; }
    virtual double Func (const Vector & x) const;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
    virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;  
    virtual double GradStopping (const Vector & x) const;
    virtual void ApproximateHesse (const Vector & x,
                                   DenseMatrix & hesse) const;
  };

  Opti3FreeMinFunction :: Opti3FreeMinFunction (const PointFunction & apf)
    : pf(apf)
  {
    ;
  }

  double Opti3FreeMinFunction :: Func (const Vector & x) const
  {
    Point<3> pp;
    for (int j = 0; j < 3; j++)
      pp(j) = sp1(j) + x(j);
    return pf.PointFunctionValue (pp);
  }
  
  double Opti3FreeMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
  {
    Vec<3> vgrad;
    Point<3> pp;

    for (int j = 0; j < 3; j++)
      pp(j) = sp1(j) + x(j);

    double val = pf.PointFunctionValueGrad (pp, vgrad);

    for (int j = 0; j < 3; j++)
      grad(j) = vgrad(j);

    return val;
  }

  double Opti3FreeMinFunction :: FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
  {
    Point<3> pp;

    for (int j = 0; j < 3; j++)
      pp(j) = sp1(j) + x(j);

    Vec<3> vdir;
    for (int j = 0; j < 3; j++)
      vdir(j) = dir(j);

    return pf.PointFunctionValueDeriv (pp, vdir, deriv);
  }
  
  double Opti3FreeMinFunction :: GradStopping (const Vector & x) const
  {
    double f = Func(x);
    return 1e-3 * f / pf.GetLocalH();
  }


  void Opti3FreeMinFunction :: ApproximateHesse (const Vector & x,
                                                 DenseMatrix & hesse) const
  {
    int n = x.Size();

    Vector hx;
    hx.SetSize(n);

    double eps = 1e-8;
    double f, f11, f22; //, f12, f21

    f = Func(x);
  
    for (int i = 1; i <= n; i++)
      {
        for (int j = 1; j < i; j++)
          {
            /*
              hx = x;
              hx.Elem(i) = x.Get(i) + eps;
              hx.Elem(j) = x.Get(j) + eps;
              f11 = Func(hx);
              hx.Elem(i) = x.Get(i) + eps;
              hx.Elem(j) = x.Get(j) - eps;
              f12 = Func(hx);
              hx.Elem(i) = x.Get(i) - eps;
              hx.Elem(j) = x.Get(j) + eps;
              f21 = Func(hx);
              hx.Elem(i) = x.Get(i) - eps;
              hx.Elem(j) = x.Get(j) - eps;
              f22 = Func(hx);
            */
            hesse.Elem(i, j) = hesse.Elem(j, i) = 0;
            //      (f11 + f22 - f12 - f21) / (2 * eps * eps);
          }

        hx = x;
        hx(i-1) = x(i-1) + eps;
        f11 = Func(hx);
        hx(i-1) = x(i-1) - eps;
        f22 = Func(hx);

        hesse.Elem(i, i) = (f11 + f22 - 2 * f) / (eps * eps) + 1e-12;
      }
  }






#ifdef SOLIDGEOM
  class Opti3SurfaceMinFunction : public MinFunction
  {
    const PointFunction & pf;
    Point<3> sp1 = Point<3>(0,0,0);
    const Surface * surf;
    Vec<3> t1, t2;
  
  public:
    Opti3SurfaceMinFunction (const PointFunction & apf);
  
    void SetPoint (const Surface * asurf, const Point<3> & asp1);

    void CalcNewPoint (const Vector & x, Point<3> & np) const; 
    virtual double Func (const Vector & x) const;
    virtual double FuncGrad (const Vector & x, Vector & g) const;
  };


  Opti3SurfaceMinFunction :: Opti3SurfaceMinFunction (const PointFunction & apf)
    : MinFunction(), pf(apf)
  {
    ;
  }

  void Opti3SurfaceMinFunction :: SetPoint (const Surface * asurf, const Point<3> & asp1)
  { 
    Vec<3> n = Vec<3>(0,0,0);
    sp1 = asp1; 
    surf = asurf;
  
    Vec<3> hn;
    surf -> GetNormalVector (sp1, hn);
    n = hn;

    n.GetNormal (t1);
    t1 /= t1.Length();
    t2 = Cross (n, t1);
  }

  
  void Opti3SurfaceMinFunction :: CalcNewPoint (const Vector & x, 
                                                Point<3> & np) const
  {
    np.X() = sp1.X() + x.Get(1) * t1.X() + x.Get(2) * t2.X();
    np.Y() = sp1.Y() + x.Get(1) * t1.Y() + x.Get(2) * t2.Y();
    np.Z() = sp1.Z() + x.Get(1) * t1.Z() + x.Get(2) * t2.Z();

    Point<3> hnp = np;
    surf -> Project (hnp);
    np = hnp;
  }


  double Opti3SurfaceMinFunction :: Func (const Vector & x) const
  {
    Point<3> pp1;

    CalcNewPoint (x, pp1);
    return pf.PointFunctionValue (pp1);
  }



  double Opti3SurfaceMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
  {
    Vec<3> n, vgrad;
    Point<3> pp1;
    VectorMem<3> freegrad;

    CalcNewPoint (x, pp1);

    double badness = pf.PointFunctionValueGrad (pp1, freegrad);
    vgrad.X() = freegrad.Get(1);
    vgrad.Y() = freegrad.Get(2);
    vgrad.Z() = freegrad.Get(3);

    Vec<3> hn;
    surf -> GetNormalVector (pp1, hn);
    n = hn;

    vgrad -= (vgrad * n) * n;

    grad.Elem(1) = vgrad * t1;
    grad.Elem(2) = vgrad * t2;
    
    return badness;
  }
#endif
  
  
  


  
  
  
#ifdef SOLIDGEOM
  class Opti3EdgeMinFunction : public MinFunction
  {
    const PointFunction & pf;
    Point<3> sp1 = Point<3>(0,0,0);
    const Surface *surf1, *surf2;
    Vec<3> t1 = Vec<3>(0,0,0);
  
  public:
    Opti3EdgeMinFunction (const PointFunction & apf);
  
    void SetPoint (const Surface * asurf1, const Surface * asurf2,
                   const Point<3> & asp1);
    void CalcNewPoint (const Vector & x, Point<3> & np) const; 
    virtual double FuncGrad (const Vector & x, Vector & g) const;
    virtual double Func (const Vector & x) const;
  };

  Opti3EdgeMinFunction :: Opti3EdgeMinFunction (const PointFunction & apf)
    : MinFunction(), pf(apf)
  {
    ;
  }
  
  void Opti3EdgeMinFunction :: SetPoint (const Surface * asurf1, 
                                         const Surface * asurf2, 
                                         const Point<3> & asp1) 
  { 
    Vec<3> n1, n2;
    sp1 = asp1; 
    surf1 = asurf1;
    surf2 = asurf2;

    Vec<3> hn1, hn2;
    surf1 -> GetNormalVector (sp1, hn1);
    surf2 -> GetNormalVector (sp1, hn2);
    n1 = hn1;
    n2 = hn2;
    t1 = Cross (n1, n2);
  }

  void Opti3EdgeMinFunction :: CalcNewPoint (const Vector & x,
                                             Point<3> & np) const
{
  np.X() = sp1.X() + x.Get(1) * t1.X();
  np.Y() = sp1.Y() + x.Get(1) * t1.Y();
  np.Z() = sp1.Z() + x.Get(1) * t1.Z();
  Point<3> hnp = np;
  ProjectToEdge (surf1, surf2, hnp);
  np = hnp;
}   

double Opti3EdgeMinFunction :: Func (const Vector & x) const
{
  Vector g(x.Size());
  return FuncGrad (x, g);
}


double Opti3EdgeMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
{
  Vec<3> n1, n2, v1, vgrad;
  Point<3> pp1 = Point<3>(0,0,0);
  double badness;
  VectorMem<3> freegrad;

  CalcNewPoint (x, pp1);


  badness = pf.PointFunctionValueGrad (pp1, freegrad);

  vgrad.X() = freegrad.Get(1);
  vgrad.Y() = freegrad.Get(2);
  vgrad.Z() = freegrad.Get(3);

  Vec<3> hn1, hn2;
  surf1 -> GetNormalVector (pp1, hn1);
  surf2 -> GetNormalVector (pp1, hn2);
  n1 = hn1;
  n2 = hn2;

  v1 = Cross (n1, n2);
  v1 /= v1.Length();

  grad.Elem(1) = (vgrad * v1) * (t1 * v1);
  return badness;
}
#endif





int WrongOrientation (const Mesh::T_POINTS & points, const Element & el)
{
  const Point<3> & p1 = points[el[0]];
  const Point<3> & p2 = points[el[1]];
  const Point<3> & p3 = points[el[2]];
  const Point<3> & p4 = points[el[3]];

  Vec<3> v1(p1, p2);
  Vec<3> v2(p1, p3);
  Vec<3> v3(p1, p4);
  Vec<3> n;

  Cross (v1, v2, n);
  double vol = n * v3;

  return (vol > 0);
}











/* ************* JacobianPointFunction **************************** */




// class JacobianPointFunction : public MinFunction
// {
// public:
//   Mesh::T_POINTS & points;
//   const Array<Element> & elements;
//   TABLE<int> elementsonpoint;
//   PointIndex actpind;
  
// public:
//   JacobianPointFunction (Mesh::T_POINTS & apoints, 
//                       const Array<Element> & aelements);
  
//   virtual void SetPointIndex (PointIndex aactpind);
//   virtual double Func (const Vector & x) const;
//   virtual double FuncGrad (const Vector & x, Vector & g) const;
//   virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
// };


JacobianPointFunction :: 
JacobianPointFunction (Mesh::T_POINTS & apoints, 
                       const Array<Element, ElementIndex> & aelements)
  : points(apoints), elements(aelements)
{
  elementsonpoint = ngcore::CreateSortedTable<ElementIndex, PointIndex>
    ( elements.Range(),
      [&](auto & table, ElementIndex ei)
      {
        for (PointIndex pi : elements[ei].PNums())
          table.Add (pi, ei);
      }, apoints.Size());

  onplane = false;
}

void JacobianPointFunction :: SetPointIndex (PointIndex aactpind)
{
  actpind = aactpind; 
}  


double JacobianPointFunction :: Func (const Vector & v) const
{
  double badness = 0;

  Point<3> hp = points[actpind];

  points[actpind] = hp + Vec<3> (v(0), v(1), v(2));

  if(onplane)
    points[actpind] -= (v(0)*nv(0)+v(1)*nv(1)+v(2)*nv(2)) * nv;


  for (auto eli : elementsonpoint[actpind])
      badness += elements[eli].CalcJacobianBadness (points);
  
  points[actpind] = hp; 

  return badness;
}





double JacobianPointFunction :: 
FuncGrad (const Vector & x, Vector & g) const
{
  int lpi;
  double badness = 0;//, hbad;

  Point<3> hp = points[actpind];
  points[actpind] = hp + Vec<3> (x(0), x(1), x(2));

  if(onplane)
    points[actpind] -= (x(0)*nv(0)+x(1)*nv(1)+x(2)*nv(2)) * nv;

  Vec<3> hderiv;
  //Vec<3> vdir;
  g.SetSize(3);
  g = 0;

  for (auto ei : elementsonpoint[actpind])
    {
      const Element & el = elements[ei];

      lpi = 0;
      for (int k = 1; k <= el.GetNP(); k++)
        if (el.PNum(k) == actpind)
          lpi = k;
      if (!lpi) cerr << "loc point not found" << endl;

      badness += elements[ei].
        CalcJacobianBadnessGradient (points, lpi, hderiv);

      for(int k=0; k<3; k++)
        g(k) += hderiv(k);
        
      /*
      for (k = 1; k <= 3; k++)
        {
          vdir = Vec<3>(0,0,0);
          vdir.X(k) = 1;

          hbad = elements.Get(eli).
            CalcJacobianBadnessDirDeriv (points, lpi, vdir, hderiv);
          //(*testout) << "hderiv " << k << ": " << hderiv << endl;
          g.Elem(k) += hderiv;
          if (k == 1)
            badness += hbad;
        }
      */
    }

  if(onplane)
    {
      double scal = nv(0)*g(0) + nv(1)*g(1) + nv(2)*g(2);
      g(0) -= scal*nv(0);
      g(1) -= scal*nv(1);
      g(2) -= scal*nv(2);
    }

  //(*testout) << "g = " << g << endl;

  
  points[actpind] = hp; 

  return badness;
}


double JacobianPointFunction :: 
FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const
{
  int lpi;
  double badness = 0;

  Point<3> hp = points[actpind];
  points[actpind] = Point<3> (hp + Vec<3> (x(0), x(1), x(2)));

  if(onplane)
    points[actpind] -= (Vec<3> (x(0), x(1), x(2))*nv) * nv;

  double hderiv;
  deriv = 0;
  Vec<3> vdir(dir(0), dir(1), dir(2));
 
  if(onplane)
    {
      double scal = vdir * nv;
      vdir -= scal*nv;
    }

  for (auto ei : elementsonpoint[actpind])
    {
      const Element & el = elements[ei];

      lpi = 0;
      for (int k = 1; k <= el.GetNP(); k++)
        if (el.PNum(k) == actpind)
          lpi = k;
      if (!lpi) cerr << "loc point not found" << endl;

      badness += elements[ei].
        CalcJacobianBadnessDirDeriv (points, lpi, vdir, hderiv);
      deriv += hderiv;
    }
  
  points[actpind] = hp; 

  return badness;
  
}













  
void Mesh :: ImproveMesh (const MeshingParameters & mp, OPTIMIZEGOAL goal)
{
  static Timer t("Mesh::ImproveMesh"); RegionTimer reg(t);
  static Timer tcoloring("coloring");
  static Timer tcalcbadmax("Calc badmax");
  static Timer topt("optimize");
  static Timer trange("range");
  static Timer tloch("loch");

  BuildBoundaryEdges(false);

  (*testout) << "Improve Mesh" << "\n";
  PrintMessage (3, "ImproveMesh");

  // int np = GetNP();
  int ne = GetNE();

  PointFunction pf_glob(*this, mp);

  auto & elementsonpoint = pf_glob.GetPointToElementTable();

  const auto & getDofs = [&] (int i)
  {
      return elementsonpoint[PointIndex::FromNr0(i)];
  };

  Array<int> colors(points.Size());

  tcoloring.Start();
  int ncolors = ngcore::ComputeColoring( colors, ne, getDofs );
  auto color_table = CreateTable<PointIndex, int>( points.Size(),
         [&] ( auto & table, int i )
          {
            PointIndex pi = PointIndex::FromNr0(i);
            table.Add(colors[i], pi);
          }, ncolors);

  tcoloring.Stop();

  if (goal == OPT_QUALITY)
    {
      double bad1 = CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad1 << endl;
      PrintMessage (5, "Total badness = ", bad1);
    }


  (*testout) << setprecision(8);

  Array<double, PointIndex> pointh (points.Size());

  if(HasLocalHFunction())
    {
      RegionTimer rt(tloch);
      ParallelForRange(points.Range(), [&] (auto myrange)
         {
           for(auto pi : myrange)
             pointh[pi] = GetH(pi);
         });
    }
  else
    {
      pointh = 0;
      for (Element & el : VolumeElements())
        {
          double h = pow(el.Volume(points),1./3.);
          for (PointIndex pi : el.PNums())
            if (h > pointh[pi])
              pointh[pi] = h;
        }
    }

  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Smooth Mesh";

  topt.Start();
  for (auto icolor : Range(ncolors))
  {
      if (multithread.terminate)
          throw NgException ("Meshing stopped");

      ParallelForRange( color_table[icolor].Range(), [&](auto myrange)
      {
        RegionTracer reg(ngcore::TaskManager::GetThreadId(), trange, myrange.Size());
        Vector x(3);

        PointFunction pf{pf_glob};

        Opti3FreeMinFunction freeminf(pf);

        OptiParameters par;
        par.maxit_linsearch = 20;
        par.maxit_bfgs = 20;

        for (auto i : myrange)
        {
          PointIndex pi = color_table[icolor][i];
          if ( (*this)[pi].Type() == INNERPOINT )
          {
            double lh = pointh[pi];
            pf.SetLocalH (lh);
            par.typx = lh;

            freeminf.SetPoint (points[pi]);
            pf.SetPointIndex (pi);

            x = 0;
            int pok;
            pok = freeminf.Func (x) < 1e10;

            if (!pok)
              {
                pok = pf.MovePointToInner ();

                freeminf.SetPoint (points[pi]);
                pf.SetPointIndex (pi);
              }

            if (pok)
              {
                //*testout << "start BFGS, pok" << endl;
                BFGS (x, freeminf, par);
                //*testout << "BFGS complete, pok" << endl;
                points[pi](0) += x(0);
                points[pi](1) += x(1);
                points[pi](2) += x(2);
              }
          }
        }
      }, 4*ngcore::TaskManager::GetNumThreads());
  }
  topt.Stop();

  multithread.task = savetask;

  if (goal == OPT_QUALITY)
    {
      double bad1 = CalcTotalBad (mp);
      (*testout) << "Total badness = " << bad1 << endl;
      PrintMessage (5, "Total badness = ", bad1);
    }
}



// Improve Condition number of Jacobian, any elements  
void Mesh :: ImproveMeshJacobian (const MeshingParameters & mp,
                                  OPTIMIZEGOAL goal, const TBitArray<PointIndex> * usepoint)
{
  // int i, j;
  
  (*testout) << "Improve Mesh Jacobian" << "\n";
  PrintMessage (3, "ImproveMesh Jacobian");

  int np = GetNP();
  int ne = GetNE();

  
  Vector x(3);
  
  (*testout) << setprecision(8);
  
  JacobianPointFunction pf(points, volelements);
  

  OptiParameters par;
  par.maxit_linsearch = 20;
  par.maxit_bfgs = 20;
  
  TBitArray<PointIndex> badnodes(np);
  badnodes.Clear();

  for (ElementIndex i : T_Range<ElementIndex>(ne))
    {
      const Element & el = (*this)[i];
      double bad = el.CalcJacobianBadness (Points());
      if (bad > 1)
        for (int j = 1; j <= el.GetNP(); j++)
          badnodes.SetBit (el.PNum(j));
    }

  Array<double, PointIndex> pointh (points.Size());

  if(HasLocalHFunction())
    {
      // for(i = 1; i<=points.Size(); i++)
      for (PointIndex pi : points.Range())
        pointh[pi] = GetH(pi);
    }
  else
    {
      pointh = 0;
      for (const Element & el : VolumeElements())
        {
          double h = pow(el.Volume(points),1./3.);
          for(int j=1; j<=el.GetNV(); j++)
            if(h > pointh[el.PNum(j)])
              pointh[el.PNum(j)] = h;
        }
    }
 


  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Smooth Mesh Jacobian";
  
  // for (PointIndex pi = points.Begin(); i < points.End(); pi++)
  for (PointIndex pi : points.Range())
    {
      if ((*this)[pi].Type() != INNERPOINT)
        continue;

      if(usepoint && !usepoint->Test(pi))
        continue;

      //(*testout) << "improvejac, p = " << i << endl;

      if (goal == OPT_WORSTCASE && !badnodes.Test(pi))
        continue;
      //        (*testout) << "smooth p " << i << endl;

      /*
        if (multithread.terminate)
        break;
      */
      if (multithread.terminate)
        throw NgException ("Meshing stopped");

      multithread.percent = 100.0 * (pi-IndexBASE<PointIndex>()) / points.Size();

      if (points.Size() < 1000)
        PrintDot ();
      else
        if ((pi-IndexBASE<PointIndex>()) % 10 == 0)
          PrintDot ('+');

      double lh = pointh[pi];
      par.typx = lh;

      pf.SetPointIndex (pi);

      x = 0;
      int pok = (pf.Func (x) < 1e10); 

      if (pok)
        {
          //*testout << "start BFGS, Jacobian" << endl;
          BFGS (x, pf, par);
          //*testout << "end BFGS, Jacobian" << endl;
          points[pi](0) += x(0);
          points[pi](1) += x(1);
          points[pi](2) += x(2);
        }
      else
        {
          cout << "el not ok" << endl;
        }
    }
  PrintDot ('\n');
  

  multithread.task = savetask;
}




// Improve Condition number of Jacobian, any elements  
void Mesh :: ImproveMeshJacobianOnSurface (const MeshingParameters & mp,
                                           const TBitArray<PointIndex> & usepoint, 
                                           const Array< Vec<3>* > & nv,
                                           OPTIMIZEGOAL goal,
                                           const Array< idmap_type* > * idmaps)
{
  // int i, j;
  
  (*testout) << "Improve Mesh Jacobian" << "\n";
  PrintMessage (3, "ImproveMesh Jacobian");

  int np = GetNP();
  int ne = GetNE();

  
  Vector x(3);
  
  (*testout).precision(8);
  
  JacobianPointFunction pf(points, volelements);

  Array< idmap_type* > locidmaps;
  const Array< idmap_type* > * used_idmaps;

  if(idmaps)
    used_idmaps = idmaps;
  else
    {
      used_idmaps = &locidmaps;
      
      for(int i=1; i<=GetIdentifications().GetMaxNr(); i++)
        {
          if(GetIdentifications().GetType(i) == Identifications::PERIODIC)
            {
              locidmaps.Append(new idmap_type);
              GetIdentifications().GetMap(i,*locidmaps.Last(),true);
            }
        }
    }

  
  bool usesum = (used_idmaps->Size() > 0);
  MinFunctionSum pf_sum;
  
  JacobianPointFunction * pf2ptr = NULL;
  if(usesum)
    {
      pf2ptr = new JacobianPointFunction(points, volelements);
      pf_sum.AddFunction(pf);
      pf_sum.AddFunction(*pf2ptr);
    }
  

  OptiParameters par;
  par.maxit_linsearch = 20;
  par.maxit_bfgs = 20;
  
  TBitArray<PointIndex> badnodes(np);
  badnodes.Clear();

  for (ElementIndex i : T_Range<ElementIndex>(ne))
    {
      const Element & el = (*this)[i];
      double bad = el.CalcJacobianBadness (Points());
      if (bad > 1)
        for (int j = 1; j <= el.GetNP(); j++)
          badnodes.SetBit (el.PNum(j));
    }

  Array<double, PointIndex> pointh (points.Size());
 
  if(HasLocalHFunction())
    {
      // for(i=1; i<=points.Size(); i++)
      for (PointIndex pi : points.Range())
        pointh[pi] = GetH(pi);
    }
  else
    {
      pointh = 0;
      for (const Element & el : VolumeElements())
        {
          double h = pow(el.Volume(points),1./3.);
          for(int j=1; j<=el.GetNV(); j++)
            if(h > pointh[el.PNum(j)])
              pointh[el.PNum(j)] = h;
        }
    }


  const char * savetask = multithread.task;
  multithread.task = "Optimize Volume: Smooth Mesh Jacobian";
  
  // for (PointIndex pi = points.Begin(); pi <= points.End(); pi++)
  for (PointIndex pi : points.Range())
    if ( usepoint.Test(pi) )
      {
        //(*testout) << "improvejac, p = " << i << endl;

        if (goal == OPT_WORSTCASE && !badnodes.Test(pi))
          continue;
        //      (*testout) << "smooth p " << i << endl;

        /*
        if (multithread.terminate)
          break;
        */
        if (multithread.terminate)
          throw NgException ("Meshing stopped");

        multithread.percent = 100.0 * (pi-IndexBASE<PointIndex>()) / points.Size();

        if (points.Size() < 1000)
          PrintDot ();
        else
          if ((pi-IndexBASE<PointIndex>()) % 10 == 0)
            PrintDot ('+');

        double lh = pointh[pi];//GetH(points.Get(i));
        par.typx = lh;

        pf.SetPointIndex (pi);

        constexpr PointIndex state0(PointIndex::INVALID);
        constexpr PointIndex statem1 = state0-1;
        
        PointIndex brother = statem1;
        if(usesum)
          {
            for(int j=0; brother == statem1 && j<used_idmaps->Size(); j++)
              {
                if(pi < (*used_idmaps)[j]->Size() + IndexBASE<PointIndex>())
                  {
                    brother = (*(*used_idmaps)[j])[pi];
                    if(brother == pi || brother == state0)
                      brother = statem1;
                  }
              }
            // if(brother >= pi)
            if(brother-pi >= 0)
              {
                pf2ptr->SetPointIndex(brother);
                pf2ptr->SetNV(*nv[brother-IndexBASE<PointIndex>()]);
              }
          }

        // if(usesum && brother < pi)
        if(usesum && (brother-pi < 0))
          continue;

        //pf.UnSetNV(); x = 0;
        //(*testout) << "before " << pf.Func(x);

        pf.SetNV(*nv[pi-IndexBASE<PointIndex>()]);

        x = 0;
        int pok = (brother == statem1) ? (pf.Func (x) < 1e10) : (pf_sum.Func (x) < 1e10);

        if (pok)
          {
            
            if(brother == statem1)
              BFGS (x, pf, par);
            else
              BFGS (x, pf_sum, par);


            for(int j=0; j<3; j++)
              points[pi](j) += x(j);// - scal*nv[i-1].X(j);

            if(brother != statem1)
              for(int j=0; j<3; j++)
                points[brother](j) += x(j);// - scal*nv[brother-1].X(j);


          }
        else
          {
            cout << "el not ok" << endl;
            (*testout) << "el not ok" << endl
                       << "   func " << ((brother == statem1) ? pf.Func(x) : pf_sum.Func (x)) << endl;
            if(brother != statem1)
              (*testout) << "   func1 " << pf.Func(x) << endl
                         << "   func2 " << pf2ptr->Func(x) << endl;
          }
      }
  
  PrintDot ('\n');

  delete pf2ptr;
  for(int i=0; i<locidmaps.Size(); i++)
    delete locidmaps[i];

  multithread.task = savetask;
}




}

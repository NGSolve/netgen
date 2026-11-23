#include <mystdlib.h>
#include "meshing.hpp"
#include "extend_curving.hpp"

#include <map>
#include <mutex>
#include <array>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

namespace netgen
{

  // ============================================================
  //   Bernstein polynomials on the tetrahedron
  // ============================================================

  template <typename T, int N>
  class LazyTable
  {
    std::atomic<T*> slots[N];
    std::mutex m;
  public:
    LazyTable () { for (auto & s : slots) s.store (nullptr); }
    ~LazyTable () { for (auto & s : slots) delete s.load(); }
    template <typename CREATE>
    const T & Get (int i, CREATE create)
    {
      if (i < 0 || i >= N) throw NgException ("ExtendCurving: table index out of range");
      T * p = slots[i].load (std::memory_order_acquire);
      if (p) return *p;
      std::lock_guard<std::mutex> lock(m);
      p = slots[i].load (std::memory_order_relaxed);
      if (!p)
        {
          p = create().release();
          slots[i].store (p, std::memory_order_release);
        }
      return *p;
    }
  };
  static constexpr int MAXDEG = 40;

  static double Binom (int n, int k)
  {
    if (k < 0 || k > n) return 0;
    double r = 1;
    for (int i = 1; i <= k; i++) r = r * (n-k+i) / i;
    return r;
  }

  static Vec<4> V4 (double a, double b, double c, double d)
  { Vec<4> v; v(0) = a; v(1) = b; v(2) = c; v(3) = d; return v; }

  // multi-indices alpha, |alpha| = d, lattice points and multinomials
  class TetBernstein
  {
  public:
    int d, n;
    Array<std::array<int,4>> alpha;
    Array<int> lookup;       // (a0,a1,a2) -> index
    Array<double> multinom;  // d! / alpha!
    Array<int> vertex;       // indices of (d,0,0,0) ... (0,0,0,d)

    TetBernstein (int ad) : d(ad)
    {
      lookup.SetSize ((d+1)*(d+1)*(d+1));
      lookup = -1;
      for (int a0 = d; a0 >= 0; a0--)
        for (int a1 = d-a0; a1 >= 0; a1--)
          for (int a2 = d-a0-a1; a2 >= 0; a2--)
            {
              lookup[(a0*(d+1)+a1)*(d+1)+a2] = alpha.Size();
              alpha.Append ({a0, a1, a2, d-a0-a1-a2});
            }
      n = alpha.Size();
      multinom.SetSize(n);
      for (int i = 0; i < n; i++)
        {
          auto a = alpha[i];
          multinom[i] = Binom(d, a[0]) * Binom(d-a[0], a[1]) * Binom(d-a[0]-a[1], a[2]);
        }
      vertex.SetSize(4);
      vertex[0] = Index(d,0,0); vertex[1] = Index(0,d,0);
      vertex[2] = Index(0,0,d); vertex[3] = Index(0,0,0);
    }

    int Index (int a0, int a1, int a2) const
    { return lookup[(a0*(d+1)+a1)*(d+1)+a2]; }

    Point<3> Lattice (int i) const
    { return Point<3> (double(alpha[i][0])/d, double(alpha[i][1])/d, double(alpha[i][2])/d); }

    double Basis (int i, const Vec<4> & lam) const
    {
      double v = multinom[i];
      for (int k = 0; k < 4; k++)
        for (int r = 0; r < alpha[i][k]; r++)
          v *= lam(k);
      return v;
    }

    static const TetBernstein & Get (int d)
    {
      static LazyTable<TetBernstein, MAXDEG> cache;
      return cache.Get (d, [d]() { return make_unique<TetBernstein>(d); });
    }
  };

  // inverse of the Bernstein-Vandermonde matrix at the lattice of degree d
  static const DenseMatrix & BernsteinVandermondeInverse (int d)
  {
    static LazyTable<DenseMatrix, MAXDEG> cache;
    return cache.Get (d, [d]()
      {
        const auto & tb = TetBernstein::Get(d);
        DenseMatrix V(tb.n, tb.n);
        for (int i = 0; i < tb.n; i++)
          {
            Point<3> xi = tb.Lattice(i);
            Vec<4> lam = V4 (xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2));
            for (int j = 0; j < tb.n; j++)
              V(i,j) = tb.Basis (j, lam);
          }
        auto p = make_unique<DenseMatrix>(tb.n, tb.n);
        CalcInverse (V, *p);
        return p;
      });
  }

  // B^m_alpha B^n_beta = w B^{m+n}_{alpha+beta}
  struct BernsteinProduct
  {
    Array<int> ia, ib, ic;
    Array<double> w;

    static const BernsteinProduct & Get (int m, int n)
    {
      static LazyTable<BernsteinProduct, MAXDEG*MAXDEG> cache;
      return cache.Get (m*MAXDEG+n, [m,n]()
        {
          auto p = make_unique<BernsteinProduct>();
          const auto & ta = TetBernstein::Get(m);
          const auto & tb = TetBernstein::Get(n);
          const auto & tc = TetBernstein::Get(m+n);
          double scale = 1.0 / Binom(m+n, m);
          for (int i = 0; i < ta.n; i++)
            for (int j = 0; j < tb.n; j++)
              {
                auto a = ta.alpha[i], b = tb.alpha[j];
                double wij = scale;
                for (int k = 0; k < 4; k++)
                  wij *= Binom (a[k]+b[k], a[k]);
                p->ia.Append(i);
                p->ib.Append(j);
                p->ic.Append(tc.Index(a[0]+b[0], a[1]+b[1], a[2]+b[2]));
                p->w.Append(wij);
              }
          return p;
        });
    }

    template <typename T>
    void MultAdd (FlatArray<T> a, FlatArray<T> b, double sign, FlatArray<T> c) const
    {
      for (int k = 0; k < ia.Size(); k++)
        c[ic[k]] += sign * w[k] * a[ia[k]] * b[ib[k]];
    }
  };

  // one de Casteljau step with point u (barycentric), degree e -> e-1
  static void DeCasteljauStep (int e, FlatArray<double> in, const Vec<4> & u, FlatArray<double> out)
  {
    const auto & te = TetBernstein::Get(e);
    const auto & tl = TetBernstein::Get(e-1);
    for (int i = 0; i < tl.n; i++)
      {
        auto b = tl.alpha[i];
        out[i] = u(0) * in[te.Index(b[0]+1, b[1], b[2])]
               + u(1) * in[te.Index(b[0], b[1]+1, b[2])]
               + u(2) * in[te.Index(b[0], b[1], b[2]+1)]
               + u(3) * in[te.Index(b[0], b[1], b[2])];
      }
  }

  // Bernstein coefficients of the restriction to the sub-tet with vertices v
  // (barycentric w.r.t. the parent), via the blossom
  static void RestrictBernstein (int d, FlatArray<double> c, const Vec<4> v[4], FlatArray<double> out)
  {
    const auto & td = TetBernstein::Get(d);
    Array<double> tmp(c.Size()), tmp2(c.Size());
    for (int i = 0; i < td.n; i++)
      {
        tmp = c;
        int e = d;
        for (int k = 0; k < 4; k++)
          for (int r = 0; r < td.alpha[i][k]; r++)
            {
              DeCasteljauStep (e, tmp.Range(0, TetBernstein::Get(e).n), v[k],
                               tmp2.Range(0, TetBernstein::Get(e-1).n));
              e--;
              std::swap (tmp, tmp2);
            }
        out[i] = tmp[0];
      }
  }

  // octasection of the reference tet, children in barycentric coordinates
  static void TetChildren (const Vec<4> v[4], Vec<4> children[8][4]);

  // coefficients of the 8 children as linear maps of the parent coefficients
  static const Array<DenseMatrix> & SubdivisionMatrices (int d)
  {
    static LazyTable<Array<DenseMatrix>, MAXDEG> cache;
    return cache.Get (d, [d]()
      {
        const auto & td = TetBernstein::Get(d);
        Vec<4> ref[4] = { V4(1,0,0,0), V4(0,1,0,0), V4(0,0,1,0), V4(0,0,0,1) };
        Vec<4> children[8][4];
        TetChildren (ref, children);
        auto p = make_unique<Array<DenseMatrix>>();
        Array<double> unit(td.n), col(td.n);
        for (int c = 0; c < 8; c++)
          {
            DenseMatrix M(td.n, td.n);
            for (int j = 0; j < td.n; j++)
              {
                unit = 0.0; unit[j] = 1;
                RestrictBernstein (d, unit, children[c], col);
                for (int i = 0; i < td.n; i++) M(i,j) = col[i];
              }
            p->Append (std::move(M));
          }
        return p;
      });
  }

  static void TetChildren (const Vec<4> v[4], Vec<4> children[8][4])
  {
    auto mid = [&](int i, int j) { return Vec<4> (0.5*(v[i]+v[j])); };
    Vec<4> m01 = mid(0,1), m02 = mid(0,2), m03 = mid(0,3),
      m12 = mid(1,2), m13 = mid(1,3), m23 = mid(2,3);
    Vec<4> c[8][4] = {
      { v[0], m01, m02, m03 }, { m01, v[1], m12, m13 },
      { m02, m12, v[2], m23 }, { m03, m13, m23, v[3] },
      { m02, m13, m01, m03 }, { m02, m13, m03, m23 },
      { m02, m13, m23, m12 }, { m02, m13, m12, m01 } };
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 4; j++)
        children[i][j] = c[i][j];
  }

  // partial derivatives of a degree p polynomial, out[k*n1+i] for d/dxi_k, degree p-1
  static void BernsteinGradient (int p, FlatArray<double> b, FlatArray<double> out)
  {
    const auto & tp = TetBernstein::Get(p);
    const auto & tp1 = TetBernstein::Get(p-1);
    for (int i = 0; i < tp1.n; i++)
      {
        auto be = tp1.alpha[i];
        int i3 = tp.Index(be[0], be[1], be[2]);
        out[i] = p * (b[tp.Index(be[0]+1, be[1], be[2])] - b[i3]);
        out[tp1.n+i] = p * (b[tp.Index(be[0], be[1]+1, be[2])] - b[i3]);
        out[2*tp1.n+i] = p * (b[tp.Index(be[0], be[1], be[2]+1)] - b[i3]);
      }
  }

  // cofactors C(r,k) = d det J / d J(r,k) of the first nrows rows, J(r,k) at
  // (3r+k)*n1 of degree p-1, cofactors at (3r+k)*n2 of degree 2(p-1)
  static void BernsteinCofactors (int p, FlatArray<double> J, int nrows, FlatArray<double> cof)
  {
    const auto & p11 = BernsteinProduct::Get(p-1, p-1);
    int n1 = TetBernstein::Get(p-1).n, n2 = TetBernstein::Get(2*(p-1)).n;
    auto Jrk = [&](int r, int k) { return J.Range((3*r+k)*n1, (3*r+k+1)*n1); };
    cof.Range(0, 3*nrows*n2) = 0.0;
    for (int r = 0; r < nrows; r++)
      for (int k = 0; k < 3; k++)
        {
          int r1 = (r+1)%3, r2 = (r+2)%3, k1 = (k+1)%3, k2 = (k+2)%3;
          auto C = cof.Range((3*r+k)*n2, (3*r+k+1)*n2);
          p11.MultAdd (Jrk(r1,k1), Jrk(r2,k2), 1.0, C);
          p11.MultAdd (Jrk(r1,k2), Jrk(r2,k1), -1.0, C);
        }
  }

  static void BernsteinMinBounds (int d, FlatArray<double> c, int depth, double t,
                                  double & lower, double & upper)
  {
    const auto & td = TetBernstein::Get(d);
    double lo = c[0], cen = 0;
    for (int i = 0; i < td.n; i++)
      {
        lo = min (lo, c[i]);
        cen += c[i] * td.multinom[i];
      }
    cen /= pow (4.0, d);
    double up = cen;
    for (int k = 0; k < 4; k++) up = min (up, c[td.vertex[k]]);

    // t = NaN: no threshold, refine until the gap is small
    bool decide = !std::isnan(t);
    lower = lo; upper = up;
    if ((decide && (lo >= t || up < t)) || depth == 0 || up - lo <= 1e-3 * max (1.0, fabs(up))) return;

    const auto & sub = SubdivisionMatrices (d);
    Array<double> cc(8*td.n);
    double clo[8];
    int order[8];
    for (int i = 0; i < 8; i++)
      {
        auto ci = cc.Range(i*td.n, (i+1)*td.n);
        const DenseMatrix & M = sub[i];
        for (int r = 0; r < td.n; r++)
          {
            double v = 0;
            for (int k = 0; k < td.n; k++) v += M(r,k) * c[k];
            ci[r] = v;
          }
        clo[i] = ci[0];
        for (int j = 1; j < td.n; j++) clo[i] = min (clo[i], ci[j]);
        order[i] = i;
      }
    std::sort (order, order+8, [&](int a, int b) { return clo[a] < clo[b]; });

    lower = 1e300; upper = up;
    for (int k = 0; k < 8; k++)
      {
        int i = order[k];
        if ((decide && clo[i] >= t) || clo[i] >= lower)
          { lower = min (lower, clo[i]); continue; }
        double l, u;
        BernsteinMinBounds (d, cc.Range(i*td.n, (i+1)*td.n), depth-1, t, l, u);
        lower = min (lower, l);
        upper = min (upper, u);
        if (decide && upper < t) break;
      }
  }



  // ============================================================
  //   quadrature on the reference elements
  // ============================================================

  static void GaussRule (int n, Array<double> & xi, Array<double> & wi)
  {
    // copy of the rule in curvedelems.cpp, on (0,1)
    xi.SetSize (n); wi.SetSize (n);
    int m = (n+1)/2;
    for (int i = 1; i <= m; i++)
      {
        double z = cos (M_PI * (i - 0.25) / (n + 0.5)), z1, pp, p1, p2, p3;
        while (true)
          {
            p1 = 1; p2 = 0;
            for (int j = 1; j <= n; j++)
              { p3 = p2; p2 = p1; p1 = ((2*j-1) * z * p2 - (j-1) * p3) / j; }
            pp = n * (z*p1-p2) / (z*z - 1);
            z1 = z; z = z1-p1/pp;
            if (fabs (z - z1) < 1e-14) break;
          }
        xi[i-1] = 0.5 * (1 - z); xi[n-i] = 0.5 * (1 + z);
        wi[i-1] = wi[n-i] = 1.0 / ((1 - z*z) * pp * pp);
      }
  }

  struct IntRule
  {
    Array<Point<3>> pts;
    Array<double> w;
  };

  // collapsed Gauss rules in netgen's reference coordinates
  static const IntRule & GetIntRule (ELEMENT_TYPE et, int n)
  {
    static LazyTable<IntRule, MAXDEG*MAXDEG> cache;
    return cache.Get (int(et)*MAXDEG+n, [et,n]()
      {
    auto p = make_unique<IntRule>();
    if (et == TET)
      {
        // lattice including vertices, edges and faces: the barrier has to
        // see folds there, Gauss points never do
        const auto & tb = TetBernstein::Get(n);
        for (int i = 0; i < tb.n; i++)
          { p->pts.Append (tb.Lattice(i)); p->w.Append (1.0/6/tb.n); }
        return p;
      }
    Array<double> x, w;
    GaussRule (n, x, w);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
          {
            double u = x[i], v = x[j], s = x[k], wgt = w[i]*w[j]*w[k];
            switch (et)
              {
              case TET: case TET10:
                p->pts.Append (Point<3> (u, v*(1-u), s*(1-u)*(1-v)));
                p->w.Append (wgt * (1-u)*(1-u)*(1-v));
                break;
              case PRISM: case PRISM12: case PRISM15:
                p->pts.Append (Point<3> (u, v*(1-u), s));
                p->w.Append (wgt * (1-u));
                break;
              case HEX: case HEX20:
                p->pts.Append (Point<3> (u, v, s));
                p->w.Append (wgt);
                break;
              case PYRAMID: case PYRAMID13:
                p->pts.Append (Point<3> (u*(1-s), v*(1-s), s));
                p->w.Append (wgt * (1-s)*(1-s));
                break;
              default:
                throw NgException ("ExtendCurving: unsupported element type");
              }
          }
    return p;
      });
  }



  // ============================================================
  //   distortion energy density
  // ============================================================

  inline double Val (double x) { return x; }
  template <int D> inline double Val (const AutoDiffDiff<D> & x) { return x.Value(); }

  /*
    f(F) = |F|^2 / (3 h(det F)^{2/3}),  h(s) = (s + sqrt(s^2+4 delta^2))/2
  */
  struct Distortion
  {
    double delta2;

    template <typename T>
    T operator() (const T * F) const
    {
      T N = F[0]*F[0];
      for (int i = 1; i < 9; i++)
        N += F[i]*F[i];
      T s = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[1]*(F[3]*F[8]-F[5]*F[6]) + F[2]*(F[3]*F[7]-F[4]*F[6]);
      T h = 0.5 * (s + sqrt (s*s + 4*delta2));
      if (Val(h) <= 0) return T(1e300);
      return N / (3 * pow (h, T(2.0/3)));
    }

    double Eval (const Mat<3> & F, Mat<3> * G, double * D2) const
    {
      if (!G && !D2)
        {
          double Fv[9];
          for (int i = 0; i < 9; i++) Fv[i] = F(i/3, i%3);
          return (*this)(Fv);
        }
      AutoDiffDiff<9> Fv[9];
      for (int i = 0; i < 9; i++) Fv[i] = AutoDiffDiff<9> (F(i/3, i%3), i);
      AutoDiffDiff<9> f = (*this)(Fv);
      if (f.Value() >= 1e300) return 1e300;
      if (G)
        for (int i = 0; i < 9; i++)
          (*G)(i/3, i%3) = f.DValue(i);
      if (D2)
        for (int i = 0; i < 9; i++)
          for (int j = 0; j < 9; j++)
            D2[9*i+j] = f.DDValue(i,j);
      return f.Value();
    }
  };

  // small dense SPD solve, returns false if not positive definite
  static bool CholeskySolve (int n, FlatArray<double> A, FlatArray<double> b, FlatArray<double> x)
  {
    Array<double> L(n*n);
    L = 0.0;
    for (int j = 0; j < n; j++)
      {
        double s = A[j*n+j];
        for (int k = 0; k < j; k++) s -= L[j*n+k]*L[j*n+k];
        if (s <= 0) return false;
        L[j*n+j] = sqrt(s);
        for (int i = j+1; i < n; i++)
          {
            double t = A[i*n+j];
            for (int k = 0; k < j; k++) t -= L[i*n+k]*L[j*n+k];
            L[i*n+j] = t / L[j*n+j];
          }
      }
    for (int i = 0; i < n; i++)
      {
        double t = b[i];
        for (int k = 0; k < i; k++) t -= L[i*n+k]*x[k];
        x[i] = t / L[i*n+i];
      }
    for (int i = n-1; i >= 0; i--)
      {
        double t = x[i];
        for (int k = i+1; k < n; k++) t -= L[k*n+i]*x[k];
        x[i] = t / L[i*n+i];
      }
    return true;
  }



  // ============================================================
  //   the optimizer
  // ============================================================

  class CurvingOptimizer
  {
    Mesh & mesh;
    CurvedElements & ce;
    const MeshTopology & top;
    const ExtendCurvingParameters & params;
    int order;
    int quad_n;
    double bezier_weight = 1e3;

    Table<ElementIndex, EdgeIndex> edge2el;
    Table<ElementIndex, FaceIndex> face2el;
    Table<SurfaceElementIndex, EdgeIndex> edge2sel;
    // vertices that must not move: boundary and identified points
    Array<bool, PointIndex> vfixed;

    Array<bool, EdgeIndex> efixed, efree, ereleased;
    Array<bool, FaceIndex> ffixed, ffree, freleased;
    Array<int, EdgeIndex> edepth;
    Array<int, FaceIndex> fdepth;
    Array<int, ElementIndex> eldepth;      // -1: not in the bent region
    Array<double, ElementIndex> lower, upper;
    Array<double, ElementIndex> elvol;     // straight volumes

    // reference values of released boundary coefficients
    Array<Vec<3>> edgecoeffs0, facecoeffs0;
    // fidelity weight per released entity, lowered only around elements
    // that are still invalid
    Array<double, EdgeIndex> elambda;
    Array<double, FaceIndex> flambda;

    struct Block { bool face; int nr; };
    Array<Block> blocks;
    Table<int, int> colors;
    // active set: a block is stepped only if an entity of one of its
    // elements moved in the previous sweep (or the energy changed)
    Array<int, EdgeIndex> eblock;
    Array<int, FaceIndex> fblock;
    Array<bool> dirty;
    std::atomic<size_t> cnt_skipped{0};

    std::map<int, Array<Point<3>>> surf_samples;
    std::atomic<size_t> cnt_steps{0}, cnt_nochol{0}, cnt_nodesc{0}, cnt_nols{0};
    // thread-summed seconds of the block step phases
    std::atomic<double> t_collect{0}, t_assemble{0}, t_solve{0};
    static void AtomicAdd (std::atomic<double> & a, double v)
    { double e = a.load(); while (!a.compare_exchange_weak (e, e+v)); }
    Array<bool, ElementIndex> touched;

  public:
    ExtendCurvingInfo info;
    void ElementBounds (ElementIndex ei, double t, double & lo, double & up, int depth = -1) const;

    CurvingOptimizer (Mesh & amesh, const ExtendCurvingParameters & aparams)
      : mesh(amesh), ce(amesh.GetCurvedElements()), top(amesh.GetTopology()), params(aparams)
    {
      order = ce.order;
      quad_n = max (order, 3);   // lattice degree of the energy, needs face and interior points
    }

    // ------------------------------------------------------------
    //   setup
    // ------------------------------------------------------------

    void Setup ()
    {
      size_t ne = mesh.GetNE(), nedges = top.GetNEdges(), nfaces = top.GetNFaces();

      // sorted rows: the block sums must not depend on the thread timing
      edge2el = CreateSortedTable<ElementIndex, EdgeIndex>
        (mesh.VolumeElements().Range(),
         [&](auto & table, ElementIndex ei)
         { for (auto e : top.GetEdges(ei)) table.Add (e, ei); }, nedges);
      face2el = CreateSortedTable<ElementIndex, FaceIndex>
        (mesh.VolumeElements().Range(),
         [&](auto & table, ElementIndex ei)
         { for (auto f : top.GetFaces(ei)) table.Add (f, ei); }, nfaces);
      edge2sel = CreateSortedTable<SurfaceElementIndex, EdgeIndex>
        (mesh.SurfaceElements().Range(),
         [&](auto & table, SurfaceElementIndex sei)
         { for (auto e : top.GetEdges(sei)) table.Add (e, sei); }, nedges);

      efixed.SetSize(nedges); efixed = false;
      ffixed.SetSize(nfaces); ffixed = false;
      for (auto sei : mesh.SurfaceElements().Range())
        {
          for (auto e : top.GetEdges(sei)) efixed[e] = true;
          ffixed[top.GetFace(sei)] = true;
        }
      for (auto si : mesh.LineSegments().Range())
        efixed[top.GetEdge(si)] = true;

      // identified entities stay as they are
      Array<bool, PointIndex> & ident = vfixed;
      ident.SetSize (mesh.GetNP());
      ident = false;
      bool have_ident = false;
      for (int nr = 1; nr <= mesh.GetIdentifications().GetMaxNr(); nr++)
        if (mesh.GetIdentifications().GetType(nr) == Identifications::PERIODIC)
          {
            idmap_type map;
            mesh.GetIdentifications().GetMap (nr, map, true);
            for (auto pi : map.Range())
              if (map[pi].IsValid() && map[pi] != pi)
                { ident[pi] = true; have_ident = true; }
          }
      if (have_ident)
        {
          for (auto e : T_Range<EdgeIndex>(nedges))
            {
              auto [p1, p2] = top.GetEdgeVertices(e);
              if (ident[p1] && ident[p2]) efixed[e] = true;
            }
          for (auto f : T_Range<FaceIndex>(nfaces))
            {
              bool all = true;
              for (auto v : top.GetFaceVertices(f)) if (!ident[v]) all = false;
              if (all) ffixed[f] = true;
            }
        }

      for (auto pi : mesh.Points().Range())
        if (mesh[pi].Type() != INNERPOINT) vfixed[pi] = true;

      efree.SetSize(nedges); efree = false;
      ffree.SetSize(nfaces); ffree = false;
      ereleased.SetSize(nedges); ereleased = false;
      freleased.SetSize(nfaces); freleased = false;
      edepth.SetSize(nedges); edepth = -1;
      fdepth.SetSize(nfaces); fdepth = -1;
      eldepth.SetSize(ne); eldepth = -1;
      lower.SetSize(ne); upper.SetSize(ne);
      elvol.SetSize(ne);

      ParallelForRange (mesh.VolumeElements().Range(), [&](auto r)
        {
          for (auto ei : r)
            {
              ElementBounds (ei, params.quality, lower[ei], upper[ei]);
              elvol[ei] = mesh[ei].Volume (mesh.Points());
            }
        });
    }

    // ------------------------------------------------------------
    //   element evaluation
    // ------------------------------------------------------------

    struct ElemEval
    {
      CurvedElements::ElementInfo info;
      ELEMENT_TYPE type;
      ArrayMem<Vec<3>,100> coefs;
      ArrayMem<double,300> dsmem;
      unique_ptr<MatrixFixWidth<3>> dshapes;

      double DShape (int j, int k) const { return (*dshapes)(j,k); }

      void Jacobians (const CurvedElements & ce, const Point<3> & xi, Mat<3> & jc, Mat<3> & js)
      {
        ce.CalcElementDShapes (info, xi, *dshapes);
        jc = 0.0; js = 0.0;
        for (int j = 0; j < info.ndof; j++)
          for (int c = 0; c < 3; c++)
            for (int k = 0; k < 3; k++)
              jc(c,k) += coefs[j](c) * (*dshapes)(j,k);
        for (int j = 0; j < info.nv; j++)
          for (int c = 0; c < 3; c++)
            for (int k = 0; k < 3; k++)
              js(c,k) += coefs[j](c) * (*dshapes)(j,k);
      }
    };

    void SetupEval (ElementIndex ei, ElemEval & ev) const
    {
      ce.SetupElementInfo (ei, ev.info);
      ev.type = mesh[ei].GetType();
      ev.coefs.SetSize (ev.info.ndof);
      ce.GetCoefficients (ev.info, &ev.coefs[0]);
      ev.dsmem.SetSize (3*ev.info.ndof);
      ev.dshapes = make_unique<MatrixFixWidth<3>> (ev.info.ndof, &ev.dsmem[0]);
    }

    // position of the entity's coefficients in the element coefficient array
    int CoefficientOffset (const CurvedElements::ElementInfo & info, bool face, int nr) const
    {
      int pos = info.nv;
      for (int i = 0; i < info.nedges; i++)
        {
          if (!face && info.edgenrs[i].Nr0() == nr) return pos;
          pos += ce.edgecoeffsindex[info.edgenrs[i]+1] - ce.edgecoeffsindex[info.edgenrs[i]];
        }
      for (int i = 0; i < info.nfaces; i++)
        {
          if (face && info.facenrs[i].Nr0() == nr) return pos;
          pos += ce.facecoeffsindex[info.facenrs[i]+1] - ce.facecoeffsindex[info.facenrs[i]];
        }
      return -1;
    }

    /*
      certified bounds of min det J_c / det J_s over the element, tets via
      Bezier coefficients, other types sampled
    */
    void ElementBoundsImpl (ElementIndex ei, double t, double & lo, double & up, int depth) const
    {
      if (depth < 0) depth = params.bezier_depth;
      ElemEval ev;
      SetupEval (ei, ev);
      auto type = ev.type;
      if (type == TET && ce.order >= 2)
        {
          int p = ce.order;
          const auto & tp = TetBernstein::Get(p);
          const auto & tp1 = TetBernstein::Get(p-1);
          const auto & tdet = TetBernstein::Get(3*(p-1));
          const auto & Vinv = BernsteinVandermondeInverse(p);

          // map at the lattice -> Bernstein coefficients of the map
          ArrayMem<double,200> shapes(ev.info.ndof);
          Array<Vec<3>> xlat(tp.n), b(tp.n);
          for (int i = 0; i < tp.n; i++)
            {
              ce.CalcElementShapes (ev.info, tp.Lattice(i), TFlatVector<double>(ev.info.ndof, &shapes[0]));
              Vec<3> x(0,0,0);
              for (int j = 0; j < ev.info.ndof; j++) x += shapes[j] * ev.coefs[j];
              xlat[i] = x;
            }
          for (int i = 0; i < tp.n; i++)
            {
              Vec<3> s(0,0,0);
              for (int j = 0; j < tp.n; j++) s += Vinv(i,j) * xlat[j];
              b[i] = s;
            }

          // Jacobian entries, degree p-1
          Array<double> J(9*tp1.n), bc(tp.n);
          for (int c = 0; c < 3; c++)
            {
              for (int i = 0; i < tp.n; i++) bc[i] = b[i](c);
              BernsteinGradient (p, bc, J.Range(3*c*tp1.n, 3*(c+1)*tp1.n));
            }

          // cofactors of the first row, degree 2(p-1), and the determinant
          const auto & p12 = BernsteinProduct::Get(p-1, 2*(p-1));
          int n2 = TetBernstein::Get(2*(p-1)).n;
          Array<double> cof(3*n2), det(tdet.n);
          BernsteinCofactors (p, J, 1, cof);
          det = 0.0;
          for (int k = 0; k < 3; k++)
            p12.MultAdd (J.Range(k*tp1.n, (k+1)*tp1.n), cof.Range(k*n2, (k+1)*n2), 1.0, det);

          // straight det
          Mat<3> jc, js;
          ev.Jacobians (ce, Point<3>(0.25,0.25,0.25), jc, js);
          double dets = Det(js);
          if (fabs(dets) < 1e-300) { lo = up = -1e300; return; }
          for (int i = 0; i < tdet.n; i++) det[i] /= dets;

          BernsteinMinBounds (3*(p-1), det, depth, t, lo, up);
          return;
        }

      // sampled
      const auto & ir = GetIntRule (type, 2*order+2);
      Mat<3> jc, js;
      double q = 1e300;
      for (int i = 0; i < ir.pts.Size(); i++)
        {
          ev.Jacobians (ce, ir.pts[i], jc, js);
          q = min (q, Det(jc)/Det(js));
        }
      lo = up = q;
    }

    // ------------------------------------------------------------
    //   block Newton step
    // ------------------------------------------------------------

    // per quadrature point data for a block, map restricted to the block's coefficients
    struct IPData
    {
      Mat<3> Ffix;        // (J_c without block part) J_s^{-1}
      double w;           // weight incl. det J_s and depth factor
      Vec<3> g[8];        // physical gradients of the block shapes
    };

    struct BlockWork
    {
      Array<ElementIndex> els;
      Array<IPData> ips;
      Array<int> ipstart;
      // per tet: Bernstein coefficients of J_c without the block part, of
      // the gradients of the block shapes, penalty weight and 1/det J_s
      Array<double> bzJ, bzdJ, bzw, bzinv;
    };

    int NCoefs (const Block & bl) const
    {
      return bl.face ? ce.facecoeffsindex[FaceIndex::FromNr0(bl.nr)+1] - ce.facecoeffsindex[FaceIndex::FromNr0(bl.nr)]
                     : ce.edgecoeffsindex[EdgeIndex::FromNr0(bl.nr)+1] - ce.edgecoeffsindex[EdgeIndex::FromNr0(bl.nr)];
    }

    Vec<3> * Coefs (const Block & bl)
    {
      return bl.face ? &ce.facecoeffs[ce.facecoeffsindex[FaceIndex::FromNr0(bl.nr)]]
                     : &ce.edgecoeffs[ce.edgecoeffsindex[EdgeIndex::FromNr0(bl.nr)]];
    }

    const Vec<3> * Coefs0 (const Block & bl) const
    {
      return bl.face ? &facecoeffs0[ce.facecoeffsindex[FaceIndex::FromNr0(bl.nr)]]
                     : &edgecoeffs0[ce.edgecoeffsindex[EdgeIndex::FromNr0(bl.nr)]];
    }

    auto Elements (const Block & bl) const
    {
      return bl.face ? face2el[FaceIndex::FromNr0(bl.nr)] : edge2el[EdgeIndex::FromNr0(bl.nr)];
    }

    double EntitySize (const Block & bl) const
    {
      if (bl.face)
        {
          auto v = top.GetFaceVertices (FaceIndex::FromNr0(bl.nr));
          Vec<3> e1 = mesh[v[1]] - mesh[v[0]], e2 = mesh[v[2]] - mesh[v[0]];
          return sqrt (Cross(e1,e2).Length());
        }
      auto [p1, p2] = top.GetEdgeVertices (EdgeIndex::FromNr0(bl.nr));
      return Dist (mesh[p1], mesh[p2]);
    }

    double DepthWeight (ElementIndex ei) const
    {
      return pow (params.beta, max(0, eldepth[ei]));
    }

    void AddBezierData (ElemEval & ev, int pos, int nc, double w, double dets, BlockWork & work) const
    {
      int p = order;
      const auto & tp = TetBernstein::Get(p);
      int n1 = TetBernstein::Get(p-1).n;
      const auto & Vinv = BernsteinVandermondeInverse(p);
      int ndof = ev.info.ndof;
      // map without the block and the block shapes at the lattice
      Array<double> shapes(ndof), vals((3+nc)*tp.n), bern(tp.n);
      for (int l = 0; l < tp.n; l++)
        {
          ce.CalcElementShapes (ev.info, tp.Lattice(l), TFlatVector<double>(ndof, &shapes[0]));
          Vec<3> x(0,0,0);
          for (int j = 0; j < ndof; j++)
            if (j < pos || j >= pos+nc) x += shapes[j] * ev.coefs[j];
          for (int r = 0; r < 3; r++) vals[r*tp.n+l] = x(r);
          for (int j = 0; j < nc; j++) vals[(3+j)*tp.n+l] = shapes[pos+j];
        }
      size_t oJ = work.bzJ.Size(), odJ = work.bzdJ.Size();
      work.bzJ.SetSize (oJ + 9*n1);
      work.bzdJ.SetSize (odJ + 3*nc*n1);
      for (int m = 0; m < 3+nc; m++)
        {
          for (int i = 0; i < tp.n; i++)
            {
              double v = 0;
              for (int l = 0; l < tp.n; l++) v += Vinv(i,l) * vals[m*tp.n+l];
              bern[i] = v;
            }
          if (m < 3) BernsteinGradient (p, bern, work.bzJ.Range(oJ+3*m*n1, oJ+3*(m+1)*n1));
          else BernsteinGradient (p, bern, work.bzdJ.Range(odJ+3*(m-3)*n1, odJ+3*(m-2)*n1));
        }
      work.bzw.Append (w / TetBernstein::Get(3*(p-1)).n);
      work.bzinv.Append (1/dets);
    }

    /*
      penalty sum w (eps - b)^2 over the Bernstein coefficients b < eps of
      det J_c / det J_s, positive coefficients certify det J_c > 0 on the
      whole element. Gradient and Gauss-Newton Hessian added to g, H.
    */
    double BezierPenalty (const BlockWork & work, int nc, const Vec<3> * c,
                          FlatArray<double> g = {}, FlatArray<double> H = {}) const
    {
      if (work.bzw.Size() == 0) return 0;
      int p = order, n = 3*nc;
      int n1 = TetBernstein::Get(p-1).n, n2 = TetBernstein::Get(2*(p-1)).n, nd = TetBernstein::Get(3*(p-1)).n;
      const auto & p12 = BernsteinProduct::Get(p-1, 2*(p-1));
      bool deriv = g.Size() > 0;
      Array<double> J(9*n1), cof(9*n2), det(nd), dd(n*nd), gb(n);
      double P = 0;
      for (int e = 0; e < work.bzw.Size(); e++)
        {
          auto dJ = work.bzdJ.Range(3*nc*n1*e, 3*nc*n1*(e+1));
          J = work.bzJ.Range(9*n1*e, 9*n1*(e+1));
          for (int j = 0; j < nc; j++)
            for (int r = 0; r < 3; r++)
              for (int k = 0; k < 3; k++)
                for (int i = 0; i < n1; i++)
                  J[(3*r+k)*n1+i] += c[j](r) * dJ[(3*j+k)*n1+i];
          BernsteinCofactors (p, J, deriv ? 3 : 1, cof);
          det = 0.0;
          for (int k = 0; k < 3; k++)
            p12.MultAdd (J.Range(k*n1, (k+1)*n1), cof.Range(k*n2, (k+1)*n2), 1.0, det);

          double w = work.bzw[e], s = work.bzinv[e];
          bool active = false;
          for (int a = 0; a < nd; a++)
            {
              double v = params.eps - s*det[a];
              if (v > 0) { P += w*v*v; active = true; }
            }
          if (!deriv || !active) continue;

          // d det / d c_j(r) = sum_k dJ_jk C(r,k)
          dd = 0.0;
          for (int j = 0; j < nc; j++)
            for (int r = 0; r < 3; r++)
              for (int k = 0; k < 3; k++)
                p12.MultAdd (dJ.Range((3*j+k)*n1, (3*j+k+1)*n1), cof.Range((3*r+k)*n2, (3*r+k+1)*n2),
                             1.0, dd.Range((3*j+r)*nd, (3*j+r+1)*nd));
          for (int a = 0; a < nd; a++)
            {
              double v = params.eps - s*det[a];
              if (v <= 0) continue;
              for (int i = 0; i < n; i++) gb[i] = s * dd[i*nd+a];
              for (int i = 0; i < n; i++)
                {
                  g[i] -= 2*w*v * gb[i];
                  for (int l = 0; l < n; l++) H[i*n+l] += 2*w * gb[i]*gb[l];
                }
            }
        }
      return P;
    }

    // energy of the block's elements for coefficients c, F = Ffix + sum c_j g_j^T
    double BlockEnergy (const BlockWork & work, int nc, const Vec<3> * c, const Distortion & dist) const
    {
      double E = BezierPenalty (work, nc, c);
      for (auto & ip : work.ips)
        {
          Mat<3> F = ip.Ffix;
          for (int j = 0; j < nc; j++)
            for (int a = 0; a < 3; a++)
              for (int b = 0; b < 3; b++)
                F(a,b) += c[j](a) * ip.g[j](b);
          double f = dist.Eval (F, nullptr, nullptr);
          if (f >= 1e300) return 1e300;
          E += ip.w * f;
        }
      return E;
    }

    double Fidelity (const Block & bl, int nc, const Vec<3> * c, double & scale) const
    {
      if (!(bl.face ? freleased[FaceIndex::FromNr0(bl.nr)] : ereleased[EdgeIndex::FromNr0(bl.nr)]))
        { scale = 0; return 0; }
      double h = EntitySize (bl), vol = 0;
      for (auto ei : Elements(bl)) vol += elvol[ei];
      double lam = bl.face ? flambda[FaceIndex::FromNr0(bl.nr)] : elambda[EdgeIndex::FromNr0(bl.nr)];
      scale = lam * vol / (h*h);
      const Vec<3> * c0 = Coefs0 (bl);
      double s = 0;
      for (int j = 0; j < nc; j++) s += (c[j]-c0[j]).Length2();
      return scale * s;
    }

    // one damped Newton step with line search, returns max relative change
    double BlockStep (const Block & bl, BlockWork & work)
    {
      int nc = NCoefs (bl);
      if (nc == 0 || nc > 8) return 0;
      Vec<3> * c = Coefs (bl);
      double h = EntitySize (bl);

      // collect quadrature data
      double t0 = WallTime();
      work.ips.SetSize(0);
      work.bzJ.SetSize(0); work.bzdJ.SetSize(0); work.bzw.SetSize(0); work.bzinv.SetSize(0);
      double smin = 1e300;
      ElemEval ev;
      for (auto ei : Elements(bl))
        {
          SetupEval (ei, ev);
          int pos = CoefficientOffset (ev.info, bl.face, bl.nr);
          const auto & ir = GetIntRule (ev.type, quad_n);
          double wdepth = DepthWeight (ei);
          Mat<3> jc, js, jsinv;
          for (int i = 0; i < ir.pts.Size(); i++)
            {
              ev.Jacobians (ce, ir.pts[i], jc, js);
              CalcInverse (js, jsinv);
              IPData ip;
              Mat<3> F = jc * jsinv;
              for (int j = 0; j < nc; j++)
                {
                  Vec<3> dphi (ev.DShape(pos+j,0), ev.DShape(pos+j,1), ev.DShape(pos+j,2));
                  ip.g[j] = Trans(jsinv) * dphi;
                }
              smin = min (smin, Det(F));
              for (int j = 0; j < nc; j++)
                for (int a = 0; a < 3; a++)
                  for (int b = 0; b < 3; b++)
                    F(a,b) -= c[j](a) * ip.g[j](b);
              ip.Ffix = F;
              ip.w = ir.w[i] * Det(js) * wdepth;
              work.ips.Append (ip);
            }
          if (ev.type == TET && order >= 2)
            AddBezierData (ev, pos, nc, bezier_weight * wdepth * fabs (elvol[ei]), Det(js), work);
        }

      // regularisation (Escobar et al.): h(smin) = epsreg, vanishes for valid blocks
      Distortion dist;
      const double epsreg = 1e-2;
      dist.delta2 = (smin < epsreg) ? epsreg * (epsreg - smin) : 0.0;
      double t1 = WallTime();
      AtomicAdd (t_collect, t1-t0);

      // gradient and Hessian
      int n = 3*nc;
      Array<double> g(n), H(n*n), d(n), Hd(n*n);
      g = 0.0; H = 0.0;
      double D2[81];
      Mat<3> G;
      for (auto & ip : work.ips)
        {
          Mat<3> F = ip.Ffix;
          for (int j = 0; j < nc; j++)
            for (int a = 0; a < 3; a++)
              for (int b = 0; b < 3; b++)
                F(a,b) += c[j](a) * ip.g[j](b);
          dist.Eval (F, &G, D2);
          for (int j = 0; j < nc; j++)
            {
              Vec<3> gj = G * ip.g[j];
              for (int k = 0; k < 3; k++) g[3*j+k] += ip.w * gj(k);
            }
          for (int j = 0; j < nc; j++)
            for (int l = 0; l < nc; l++)
              for (int k = 0; k < 3; k++)
                for (int m2 = 0; m2 < 3; m2++)
                  {
                    double v = 0;
                    for (int m = 0; m < 3; m++)
                      for (int nn = 0; nn < 3; nn++)
                        v += ip.g[j](m) * D2[(3*k+m)*9 + 3*m2+nn] * ip.g[l](nn);
                    H[(3*j+k)*n + 3*l+m2] += ip.w * v;
                  }
        }
      BezierPenalty (work, nc, c, g, H);
      double fidscale;
      ArrayMem<Vec<3>,8> c0(nc), ctrial(nc);
      for (int j = 0; j < nc; j++) c0[j] = c[j];
      double E0 = BlockEnergy (work, nc, c0.Data(), dist) + Fidelity (bl, nc, c0.Data(), fidscale);
      if (fidscale > 0)
        {
          const Vec<3> * cg = Coefs0 (bl);
          for (int j = 0; j < nc; j++)
            for (int k = 0; k < 3; k++)
              {
                g[3*j+k] += 2*fidscale * (c0[j](k) - cg[j](k));
                H[(3*j+k)*n + 3*j+k] += 2*fidscale;
              }
        }

      double t2 = WallTime();
      AtomicAdd (t_assemble, t2-t1);

      // damped Newton direction
      double trace = 0;
      for (int i = 0; i < n; i++) trace += H[i*n+i];
      double mu = 0;
      Array<double> mg(n);
      for (int i = 0; i < n; i++) mg[i] = -g[i];
      bool ok = false;
      for (int it = 0; it < 10 && !ok; it++)
        {
          Hd = H;
          for (int i = 0; i < n; i++) Hd[i*n+i] += mu;
          ok = CholeskySolve (n, Hd, mg, d);
          if (!ok) mu = (mu == 0) ? 1e-6 * trace / n + 1e-300 : 10*mu;
        }
      cnt_steps++;
      if (!ok) { cnt_nochol++; return 0; }
      double gd = 0;
      for (int i = 0; i < n; i++) gd += g[i]*d[i];
      if (gd >= 0) { cnt_nodesc++; return 0; }

      // line search
      double alpha = 1, change = 0;
      for (int it = 0; it < 15; it++, alpha *= 0.5)
        {
          for (int j = 0; j < nc; j++)
            for (int k = 0; k < 3; k++)
              ctrial[j](k) = c0[j](k) + alpha * d[3*j+k];
          double E = BlockEnergy (work, nc, ctrial.Data(), dist) + Fidelity (bl, nc, ctrial.Data(), fidscale);
          if (E <= E0 + 1e-4 * alpha * gd)
            {
              for (int j = 0; j < nc; j++)
                {
                  change = max (change, (ctrial[j]-c0[j]).Length() / h);
                  c[j] = ctrial[j];
                }
              if (change > 1e-12)
                for (auto ei : Elements(bl)) touched[ei] = true;
              break;
            }
          if (it == 14) cnt_nols++;
        }
      AtomicAdd (t_solve, WallTime()-t2);
      return change;
    }

    // ------------------------------------------------------------
    //   sweeps
    // ------------------------------------------------------------

    // all active entities, or (release) the released ones and the free
    // interior entities of their elements
    void BuildBlocks (bool near_released = false)
    {
      blocks.SetSize(0);
      if (near_released)
        {
          Array<bool, EdgeIndex> ue(efree.Size()); ue = false;
          Array<bool, FaceIndex> uf(ffree.Size()); uf = false;
          auto mark = [&](auto els)
            {
              for (auto ei : els)
                {
                  for (auto e : top.GetEdges(ei)) if (efree[e]) ue[e] = true;
                  for (auto f : top.GetFaces(ei)) if (ffree[f]) uf[f] = true;
                }
            };
          for (auto e : efree.Range()) if (ereleased[e]) { ue[e] = true; mark (edge2el[e]); }
          for (auto f : ffree.Range()) if (freleased[f]) { uf[f] = true; mark (face2el[f]); }
          for (auto e : ue.Range()) if (ue[e]) blocks.Append ({false, e.Nr0()});
          for (auto f : uf.Range()) if (uf[f]) blocks.Append ({true, f.Nr0()});
        }
      else
        {
          for (auto e : efree.Range())
            if (efree[e] || ereleased[e]) blocks.Append ({false, e.Nr0()});
          for (auto f : ffree.Range())
            if (ffree[f] || freleased[f]) blocks.Append ({true, f.Nr0()});
        }

      Array<int> col(blocks.Size());
      int ncol = ngcore::ComputeColoring (col, mesh.GetNE(), [&](int i) { return Elements(blocks[i]); });
      colors = CreateSortedTable<int,int> (Range(blocks.Size()),
                                     [&](auto & table, int i) { table.Add (col[i], i); }, ncol);

      eblock.SetSize (efree.Size()); eblock = -1;
      fblock.SetSize (ffree.Size()); fblock = -1;
      for (int i = 0; i < blocks.Size(); i++)
        if (blocks[i].face) fblock[FaceIndex::FromNr0(blocks[i].nr)] = i;
        else eblock[EdgeIndex::FromNr0(blocks[i].nr)] = i;
      dirty.SetSize (blocks.Size());
      dirty = true;
    }

    void MarkDirty ()
    { dirty = true; }

    void MarkElementDirty (ElementIndex ei)
    {
      for (auto e : top.GetEdges(ei)) if (eblock[e] >= 0) dirty[eblock[e]] = true;
      for (auto f : top.GetFaces(ei)) if (fblock[f] >= 0) dirty[fblock[f]] = true;
    }

    void MarkNeighboursDirty (const Block & bl)
    {
      for (auto ei : Elements(bl))
        MarkElementDirty (ei);
    }

    double Sweep ()
    {
      double maxchange = 0;
      std::mutex m;
      static Timer tblocks("ExtendCurving - blocks of one color");
      for (auto icol : Range(colors.Size()))
        ParallelForRange (Range(colors[icol].Size()), [&](auto r)
          {
            RegionTracer reg(ngcore::TaskManager::GetThreadId(), tblocks, r.Size());
            BlockWork work;
            double mc = 0;
            for (auto i : r)
              {
                int bi = colors[icol][i];
                if (!dirty[bi]) { cnt_skipped++; continue; }
                dirty[bi] = false;
                double ch = BlockStep (blocks[bi], work);
                if (ch > 1e-7) MarkNeighboursDirty (blocks[bi]);
                mc = max (mc, ch);
              }
            std::lock_guard<std::mutex> lock(m);
            maxchange = max (maxchange, mc);
          });
      return maxchange;
    }

    // elements adjacent to active blocks
    void RegionElements (Array<ElementIndex> & region) const
    {
      Array<bool, ElementIndex> mark(mesh.GetNE());
      mark = false;
      for (auto & bl : blocks)
        for (auto ei : Elements(bl)) mark[ei] = true;
      region.SetSize(0);
      for (auto ei : mark.Range()) if (mark[ei]) region.Append (ei);
    }

    void UpdateBounds (FlatArray<ElementIndex> els, double t, int depth)
    {
      ParallelForRange (Range(els.Size()), [&](auto r)
        {
          for (auto i : r) ElementBounds (els[i], t, lower[els[i]], upper[els[i]], depth);
        });
    }

    // sweeps until all targets reach t or the coefficients stop moving,
    // returns the largest relative change of the last sweep
    double Optimize (double t, FlatArray<ElementIndex> targets, int maxsweeps)
    {
      static Timer tsweep("ExtendCurving - sweep"), tbounds("ExtendCurving - bounds");
      if (blocks.Size() == 0) return 0;
      touched.SetSize (mesh.GetNE());
      double change = 0;
      for (int sweep = 0; sweep < maxsweeps; sweep++)
        {
          touched = false;
          cnt_steps = cnt_nochol = cnt_nodesc = cnt_nols = cnt_skipped = 0;
          tsweep.Start();
          change = Sweep ();
          tsweep.Stop();
          Array<ElementIndex> upd;
          for (auto ei : touched.Range()) if (touched[ei]) upd.Append (ei);
          tbounds.Start();
          UpdateBounds (upd, t, 1);
          tbounds.Stop();
          size_t nbelow = 0;
          double minlow = 1e300;
          for (auto ei : targets)
            { if (lower[ei] < t) nbelow++; minlow = min (minlow, lower[ei]); }
          PrintMessage (6, "  sweep ", sweep, ": change ", change, ", blocks ", blocks.Size(),
                        " skipped ", size_t(cnt_skipped), " steps ", size_t(cnt_steps), " nochol ", size_t(cnt_nochol), " nodesc ", size_t(cnt_nodesc),
                        " nols ", size_t(cnt_nols), ", targets below ", nbelow, " min lower ", minlow);
          if (nbelow == 0 || change < 1e-4) break;
        }
      // sharpen the undecided ones
      Array<ElementIndex> region, amb;
      RegionElements (region);
      for (auto ei : region) if (lower[ei] < t && upper[ei] >= t) amb.Append (ei);
      tbounds.Start();
      UpdateBounds (amb, t, params.bezier_depth);
      tbounds.Stop();
      return change;
    }

    // ------------------------------------------------------------
    //   unlocking by need
    // ------------------------------------------------------------

    // a boundary face whose own parametrisation folds makes the adjacent
    // element invalid at that face for any interior bending
    bool SurfaceFolded (ElementIndex ei, double t)
    {
      int n = order+2;
      for (auto f : top.GetFaces(ei))
        {
          if (!ffixed[f] || freleased[f]) continue;
          auto sei = top.GetFace2SurfaceElement (f);
          if (!sei.IsValid()) continue;
          auto sel = mesh[sei];
          if (sel.GetType() != TRIG) continue;
          Vec<3> ns = Cross (mesh[sel[0]]-mesh[sel[2]], mesh[sel[1]]-mesh[sel[2]]);
          for (int i = 0; i <= n; i++)
            for (int j = 0; i+j <= n; j++)
              {
                Mat<3,2> dxdxi;
                ce.CalcSurfaceTransformation (Point<2>(double(i)/n, double(j)/n), sei, dxdxi);
                Vec<3> nc = Cross (Vec<3>(dxdxi(0,0), dxdxi(1,0), dxdxi(2,0)),
                                   Vec<3>(dxdxi(0,1), dxdxi(1,1), dxdxi(2,1)));
                if ((nc*ns) < t * (ns*ns)) return true;
              }
        }
      return false;
    }

    bool UnfixableByBending (ElementIndex ei, double t)
    {
      return SurfaceFolded (ei, t) || FixedCorner (ei, t);
    }

    // a corner with det below t whose edges are all fixed: face and interior
    // bubbles vanish at the corner, bending cannot repair it
    bool FixedCorner (ElementIndex ei, double t)
    {
      ElemEval ev;
      SetupEval (ei, ev);
      if (ev.type != TET) return false;
      auto el = mesh[ei];
      auto edges = MeshTopology::GetEdges (TET);
      auto elegdes = top.GetEdges (ei);
      Point<3> vxi[4] = { Point<3>(1,0,0), Point<3>(0,1,0), Point<3>(0,0,1), Point<3>(0,0,0) };
      for (int v = 0; v < 4; v++)
        {
          Mat<3> jc, js;
          ev.Jacobians (ce, vxi[v], jc, js);
          if (Det(jc)/Det(js) >= t) continue;
          bool allfixed = true;
          for (int i = 0; i < 6; i++)
            if (edges[i][0] == v || edges[i][1] == v)
              if (!efixed[elegdes[i]] || ereleased[elegdes[i]]) allfixed = false;
          if (allfixed) return true;
        }
      return false;
    }

    // an inner vertex at or behind a curved boundary face of the element:
    // the face bulges at least up to the height of the vertex over the
    // straight face, no bending repairs that, the vertex has to move
    bool Pierced (ElementIndex ei)
    {
      auto el = mesh[ei];
      if (el.GetType() != TET) return false;
      auto lfaces = MeshTopology::GetFaces (TET);
      auto faces = top.GetFaces (ei);
      int n = order+2;
      for (int i = 0; i < 4; i++)
        {
          if (!ffixed[faces[i]]) continue;
          auto sei = top.GetFace2SurfaceElement (faces[i]);
          if (!sei.IsValid()) continue;
          int vl = 6 - lfaces[i][0] - lfaces[i][1] - lfaces[i][2];
          if (vfixed[el[vl]]) continue;
          Point<3> a = mesh[el[lfaces[i][0]]], b = mesh[el[lfaces[i][1]]], c = mesh[el[lfaces[i][2]]];
          Vec<3> nrm = Cross (b-a, c-a);
          nrm.Normalize();
          if ((mesh[el[vl]] - a) * nrm < 0) nrm *= -1;
          double hv = (mesh[el[vl]] - a) * nrm;
          for (int ix = 0; ix <= n; ix++)
            for (int iy = 0; ix+iy <= n; iy++)
              {
                Point<3> x;
                ce.CalcSurfaceTransformation (Point<2>(double(ix)/n, double(iy)/n), sei, x);
                if ((x - a) * nrm >= hv) return true;
              }
        }
      return false;
    }

    void Classify ()
    {
      info.causes.SetSize (5);
      info.causes = 0;
      for (auto ei : lower.Range())
        {
          if (lower[ei] >= params.eps) continue;
          bool allbnd = true;
          for (auto pi : mesh[ei].PNums()) if (mesh[pi].Type() == INNERPOINT) allbnd = false;
          int cause = 4;
          if (SurfaceFolded (ei, params.eps)) cause = 0;
          else if (allbnd) cause = 1;
          else if (Pierced (ei)) cause = 2;
          else if (FixedCorner (ei, params.eps)) cause = 3;
          info.causes[cause]++;
        }
      PrintMessage (5, "extend curving, still invalid after bending: ", info.causes[0], " folded boundary face, ",
                    info.causes[1], " all vertices on boundary, ", info.causes[2], " inner vertex behind curved face, ",
                    info.causes[3], " corner with fixed edges, ", info.causes[4], " other");
    }

    // ------------------------------------------------------------
    //   vertex relaxation before the bending
    // ------------------------------------------------------------

    /*
      Distortion energy of the curved map (boundary curved, interior as it
      is) w.r.t. the straight element of the original vertices, as function
      of the position x of one vertex: J = J_rest + x dlambda^T for the
      curved and the straight map. Newton with a difference Hessian.
    */
    struct VertexIP { Mat<3> jc, js, jrefinv; Vec<3> dl; double w, detref; };

    double VertexStep (PointIndex v, FlatArray<Point<3>, PointIndex> p0,
                       const Table<ElementIndex, PointIndex> & point2el, Array<VertexIP> & ips)
    {
      Vec<3> x0 (mesh[v]);
      double h = 0;
      for (auto ei : point2el[v]) h += cbrt (fabs (elvol[ei]));
      h /= max (size_t(1), point2el[v].Size());
      ips.SetSize(0);
      ElemEval ev;
      double smin = 1e300;
      for (auto ei : point2el[v])
        {
          if (mesh[ei].GetType() != TET) return 0;
          SetupEval (ei, ev);
          auto el = mesh[ei];
          int lv = -1;
          for (int j = 0; j < ev.info.nv; j++) if (el[j] == v) lv = j;
          const auto & ir = GetIntRule (TET, quad_n);
          for (int i = 0; i < ir.pts.Size(); i++)
            {
              VertexIP ip;
              Mat<3> jc, js, jref;
              ev.Jacobians (ce, ir.pts[i], jc, js);
              jref = 0.0;
              for (int j = 0; j < ev.info.nv; j++)
                for (int a = 0; a < 3; a++)
                  for (int k = 0; k < 3; k++)
                    jref(a,k) += p0[el[j]](a) * ev.DShape(j,k);
              CalcInverse (jref, ip.jrefinv);
              ip.detref = Det (jref);
              ip.dl = Vec<3> (ev.DShape(lv,0), ev.DShape(lv,1), ev.DShape(lv,2));
              smin = min (smin, Det (jc * ip.jrefinv));
              for (int a = 0; a < 3; a++)
                for (int k = 0; k < 3; k++)
                  { jc(a,k) -= x0(a) * ip.dl(k); js(a,k) -= x0(a) * ip.dl(k); }
              ip.jc = jc; ip.js = js;
              ip.w = ir.w[i] * fabs (ip.detref);
              ips.Append (ip);
            }
        }
      Distortion dist;
      const double epsreg = 1e-2;
      dist.delta2 = (smin < epsreg) ? epsreg * (epsreg - smin) : 0.0;

      auto energy = [&] (const Vec<3> & x)
        {
          double E = 0;
          for (auto & ip : ips)
            {
              Mat<3> jc = ip.jc, js = ip.js;
              for (int a = 0; a < 3; a++)
                for (int k = 0; k < 3; k++)
                  { jc(a,k) += x(a) * ip.dl(k); js(a,k) += x(a) * ip.dl(k); }
              // the straight element must keep its orientation
              if (Det (js) / ip.detref < 0.1) return 1e300;
              double f = dist.Eval (jc * ip.jrefinv, nullptr, nullptr);
              if (f >= 1e300) return 1e300;
              E += ip.w * f;
            }
          return E;
        };

      double E0 = energy (x0);
      if (E0 >= 1e300) return 0;
      double d = 1e-4 * h;
      Vec<3> g;
      Mat<3> Hm;
      for (int i = 0; i < 3; i++)
        {
          Vec<3> xp = x0, xm = x0;
          xp(i) += d; xm(i) -= d;
          double Ep = energy (xp), Em = energy (xm);
          if (Ep >= 1e300 || Em >= 1e300) return 0;
          g(i) = (Ep - Em) / (2*d);
          Hm(i,i) = (Ep - 2*E0 + Em) / (d*d);
        }
      for (int i = 0; i < 3; i++)
        for (int j = i+1; j < 3; j++)
          {
            Vec<3> xpp = x0, xpm = x0, xmp = x0, xmm = x0;
            xpp(i) += d; xpp(j) += d;  xpm(i) += d; xpm(j) -= d;
            xmp(i) -= d; xmp(j) += d;  xmm(i) -= d; xmm(j) -= d;
            Hm(i,j) = Hm(j,i) = (energy(xpp) - energy(xpm) - energy(xmp) + energy(xmm)) / (4*d*d);
          }

      ArrayMem<double,9> Hd(9);
      ArrayMem<double,3> mg(3), dir(3);
      double trace = Hm(0,0) + Hm(1,1) + Hm(2,2), mu = 0;
      bool ok = false;
      for (int it = 0; it < 10 && !ok; it++)
        {
          for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
              Hd[3*i+j] = Hm(i,j) + (i == j ? mu : 0);
          for (int i = 0; i < 3; i++) mg[i] = -g(i);
          ok = CholeskySolve (3, Hd, mg, dir);
          if (!ok) mu = (mu == 0) ? 1e-6 * fabs(trace) / 3 + 1e-300 : 10*mu;
        }
      if (!ok) return 0;
      double gd = g(0)*dir[0] + g(1)*dir[1] + g(2)*dir[2];
      if (gd >= 0) return 0;

      double alpha = 1;
      for (int it = 0; it < 15; it++, alpha *= 0.5)
        {
          Vec<3> x = x0 + alpha * Vec<3>(dir[0], dir[1], dir[2]);
          if (energy (x) <= E0 + 1e-4 * alpha * gd)
            {
              mesh[v] = Point<3> (x);
              return (x - x0).Length() / h;
            }
        }
      return 0;
    }

    void MakeRoom ()
    {
      static Timer t("ExtendCurving - move vertices"); RegionTimer reg(t);
      Array<bool, PointIndex> cand(mesh.GetNP());
      cand = false;
      bool any = false;
      for (auto ei : lower.Range())
        if (lower[ei] < params.quality)
          for (auto pi : mesh[ei].PNums()) { cand[pi] = true; any = true; }
      if (!any) return;

      auto point2el = CreateSortedTable<ElementIndex, PointIndex>
        (mesh.VolumeElements().Range(),
         [&](auto & table, ElementIndex ei)
         { for (auto pi : mesh[ei].PNums()) table.Add (pi, ei); }, mesh.GetNP());
      for (int l = 0; l < params.vertex_layers; l++)
        {
          Array<bool, PointIndex> next(cand);
          for (auto pi : cand.Range())
            if (cand[pi])
              for (auto ei : point2el[pi])
                for (auto pj : mesh[ei].PNums()) next[pj] = true;
          cand = std::move(next);
        }
      Array<PointIndex> verts;
      for (auto pi : cand.Range())
        if (cand[pi] && !vfixed[pi]) verts.Append (pi);
      if (verts.Size() == 0) return;

      Array<Point<3>, PointIndex> p0(mesh.GetNP());
      for (auto pi : mesh.Points().Range()) p0[pi] = mesh[pi];

      Array<int> col(verts.Size());
      int ncol = ngcore::ComputeColoring (col, mesh.GetNE(), [&](int i) { return point2el[verts[i]]; });
      auto vcolors = CreateSortedTable<int,int> (Range(verts.Size()),
                                                 [&](auto & table, int i) { table.Add (col[i], i); }, ncol);
      for (int sweep = 0; sweep < params.vertex_sweeps; sweep++)
        {
          double maxchange = 0;
          std::mutex m;
          for (auto icol : Range(vcolors.Size()))
            ParallelForRange (Range(vcolors[icol].Size()), [&](auto r)
              {
                Array<VertexIP> ips;
                double mc = 0;
                for (auto i : r)
                  mc = max (mc, VertexStep (verts[vcolors[icol][i]], p0, point2el, ips));
                std::lock_guard<std::mutex> lock(m);
                maxchange = max (maxchange, mc);
              });
          PrintMessage (6, "  vertex sweep ", sweep, ": max move ", maxchange);
          if (maxchange < 1e-4) break;
        }

      // the straight elements changed: volumes and bounds
      Array<bool, ElementIndex> changed(mesh.GetNE());
      changed = false;
      for (auto pi : verts)
        if (Dist (mesh[pi], p0[pi]) > 0)
          {
            info.nmoved++;
            double h = 0;
            for (auto ei : point2el[pi]) { changed[ei] = true; h += cbrt (fabs (elvol[ei])); }
            h /= max (size_t(1), point2el[pi].Size());
            info.max_vertex_move = max (info.max_vertex_move, Dist (mesh[pi], p0[pi]) / h);
          }
      Array<ElementIndex> upd;
      for (auto ei : changed.Range()) if (changed[ei]) upd.Append (ei);
      ParallelForRange (Range(upd.Size()), [&](auto r)
        {
          for (auto i : r)
            {
              elvol[upd[i]] = mesh[upd[i]].Volume (mesh.Points());
              ElementBounds (upd[i], params.quality, lower[upd[i]], upper[upd[i]]);
            }
        });
      PrintMessage (5, "extend curving: ", info.nmoved, " vertices moved, max ", info.max_vertex_move, " h");
    }

    // free the interior entities of the elements, raise their order
    void Unlock (FlatArray<ElementIndex> els, int depth)
    {
      Array<EdgeIndex> newe;
      Array<FaceIndex> newf;
      for (auto ei : els)
        {
          if (eldepth[ei] == -1) eldepth[ei] = depth;
          for (auto e : top.GetEdges(ei))
            if (!efixed[e] && !efree[e]) { efree[e] = true; edepth[e] = depth; newe.Append(e); }
          for (auto f : top.GetFaces(ei))
            if (!ffixed[f] && !ffree[f]) { ffree[f] = true; fdepth[f] = depth; newf.Append(f); }
        }
      for (auto e : newe)
        for (auto ei : edge2el[e])
          if (eldepth[ei] == -1) eldepth[ei] = depth+1;
      for (auto f : newf)
        for (auto ei : face2el[f])
          if (eldepth[ei] == -1) eldepth[ei] = depth+1;
      ce.RaiseOrder (newe, newf);
    }

    void BendInterior ()
    {
      Array<ElementIndex> bad;
      for (auto ei : lower.Range())
        if (lower[ei] < params.quality) bad.Append (ei);
      if (bad.Size() == 0) return;

      Array<bool, ElementIndex> stalled(mesh.GetNE());
      stalled = false;
      size_t nunfix = 0;
      for (auto ei : bad)
        if (UnfixableByBending (ei, params.quality)) { stalled[ei] = true; nunfix++; }
      PrintMessage (5, "extend curving: ", bad.Size(), " bad elements, ", nunfix,
                    " not repairable by bending");
      info.depth_histogram.SetSize (params.max_depth+1);
      info.depth_histogram = 0;

      for (int depth = 0; depth <= params.max_depth; depth++)
        {
          Array<ElementIndex> active;
          for (auto ei : bad) if (!stalled[ei]) active.Append (ei);
          if (active.Size() == 0) break;
          info.depth_histogram[depth] = active.Size();

          Array<double, ElementIndex> before(lower);

          if (depth == 0)
            Unlock (active, 0);
          else
            {
              // next layer: elements sharing an edge with the still-bad ones
              Array<ElementIndex> layer;
              Array<bool, ElementIndex> mark(mesh.GetNE());
              mark = false;
              for (auto ei : active)
                for (auto e : top.GetEdges(ei))
                  for (auto nb : edge2el[e])
                    if (!mark[nb]) { mark[nb] = true; layer.Append (nb); }
              Unlock (layer, depth);
            }
          BuildBlocks ();
          Optimize (params.quality, active, params.max_sweeps);

          Array<ElementIndex> stillbad;
          for (auto ei : active)
            {
              if (lower[ei] >= params.quality) continue;
              stillbad.Append (ei);
              // a layer that closed less than 1% of the gap will not be helped by the next one
              if (lower[ei] - before[ei] < 0.01 * (params.quality - before[ei]))
                stalled[ei] = true;
            }
          PrintMessage (5, "extend curving depth ", depth, ": ", active.Size(),
                        " bad elements, ", stillbad.Size(), " remaining");
          bad = std::move(stillbad);
        }
    }

    // ------------------------------------------------------------
    //   release of boundary curving
    // ------------------------------------------------------------

    void StoreSurfaceSamples (EdgeIndex e)
    {
      for (auto sei : edge2sel[e]) StoreSurfaceSamples (sei);
    }
    void StoreSurfaceSamples (FaceIndex f)
    {
      auto sei = top.GetFace2SurfaceElement (f);
      if (sei.IsValid()) StoreSurfaceSamples (sei);
    }
    void StoreSurfaceSamples (SurfaceElementIndex sei)
    {
      if (surf_samples.count(sei.Nr0())) return;
      Array<Point<3>> pts;
      SurfaceSamples (sei, pts);
      surf_samples[sei.Nr0()] = std::move(pts);
    }
    void SurfaceSamples (SurfaceElementIndex sei, Array<Point<3>> & pts)
    {
      int n = 2*order+1;
      bool quad = mesh[sei].GetType() == QUAD || mesh[sei].GetType() == QUAD8;
      for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
          {
            if (!quad && i+j > n) continue;
            Point<3> x;
            ce.CalcSurfaceTransformation (Point<2>(double(i)/n, double(j)/n), sei, x);
            pts.Append (x);
          }
    }

    void ReleaseBoundary ()
    {
      edgecoeffs0 = ce.edgecoeffs;
      facecoeffs0 = ce.facecoeffs;
      elambda.SetSize (efixed.Size()); elambda = 0.0;
      flambda.SetSize (ffixed.Size()); flambda = 0.0;

      for (int pass = 0; pass < params.max_release_passes; pass++)
        {
          Array<ElementIndex> bad;
          for (auto ei : lower.Range())
            if (lower[ei] < params.eps) bad.Append (ei);
          if (bad.Size() == 0) break;

          // boundary entities of the bad elements and `pass` layers around them
          Array<ElementIndex> layer(bad);
          Array<bool, ElementIndex> mark(mesh.GetNE());
          mark = false;
          for (auto ei : bad) mark[ei] = true;
          for (int l = 0; l < pass; l++)
            {
              Array<ElementIndex> next;
              for (auto ei : layer)
                for (auto pi : mesh[ei].PNums())
                  for (auto nb : top.GetVertexElements(pi))
                    if (!mark[nb]) { mark[nb] = true; next.Append (nb); }
              layer = std::move(next);
            }
          Array<ElementIndex> marked;
          for (auto ei : mark.Range()) if (mark[ei]) marked.Append (ei);

          // free the interior of these elements too; RaiseOrder reallocates,
          // so the geometry reference of already released entities is remapped
          Array<int, EdgeIndex> oldeidx(ce.edgecoeffsindex);
          Array<int, FaceIndex> oldfidx(ce.facecoeffsindex);
          Array<Vec<3>> olde0(edgecoeffs0), oldf0(facecoeffs0);
          Unlock (marked, 0);
          edgecoeffs0 = ce.edgecoeffs;
          facecoeffs0 = ce.facecoeffs;
          for (auto e : ereleased.Range())
            if (ereleased[e])
              for (int j = 0; j < min (oldeidx[e+1]-oldeidx[e], ce.edgecoeffsindex[e+1]-ce.edgecoeffsindex[e]); j++)
                edgecoeffs0[ce.edgecoeffsindex[e]+j] = olde0[oldeidx[e]+j];
          for (auto f : freleased.Range())
            if (freleased[f])
              for (int j = 0; j < min (oldfidx[f+1]-oldfidx[f], ce.facecoeffsindex[f+1]-ce.facecoeffsindex[f]); j++)
                facecoeffs0[ce.facecoeffsindex[f]+j] = oldf0[oldfidx[f]+j];

          size_t nnew = 0;
          for (auto ei : marked)
              {
                for (auto e : top.GetEdges(ei))
                  if (efixed[e] && !ereleased[e] && ce.edgecoeffsindex[e+1] > ce.edgecoeffsindex[e])
                    { StoreSurfaceSamples (e); ereleased[e] = true; nnew++; }
                for (auto f : top.GetFaces(ei))
                  if (ffixed[f] && !freleased[f] && ce.facecoeffsindex[f+1] > ce.facecoeffsindex[f])
                    { StoreSurfaceSamples (f); freleased[f] = true; nnew++; }
              }
          if (nnew == 0 && pass > 0) break;
          info.nreleased += nnew;
          BuildBlocks (true);

          Array<ElementIndex> region;
          RegionElements (region);
          // continuation in the fidelity weight: large steps while nothing
          // moves, halving once the boundary starts to give way, the full
          // sweep budget at the floor
          bool done = false;
          const double lam_floor = 1e-5;
          const double lam_start = 1e8;
          for (auto e : ereleased.Range()) if (ereleased[e] && elambda[e] == 0) elambda[e] = lam_start;
          for (auto f : freleased.Range()) if (freleased[f] && flambda[f] == 0) flambda[f] = lam_start;

          Array<bool, ElementIndex> isbad(mesh.GetNE());
          bool last = false;
          int radius = 1;
          double lam_min = lam_start;
          while (true)
            {
              MarkDirty ();
              double change = Optimize (params.eps, bad, last ? params.max_sweeps : min (params.max_sweeps, 3));
              // released entities may have hurt neighbours, they join the targets
              bad.SetSize(0);
              for (auto ei : region) if (lower[ei] < params.eps) bad.Append (ei);
              done = (bad.Size() == 0);
              if (done || last) break;

              // lower the weight only on entities within `radius` vertex layers
              // of a still invalid element (it may be invalid through a
              // neighbour's entity)
              isbad = false;
              Array<ElementIndex> front(bad);
              for (auto ei : bad) isbad[ei] = true;
              for (int l = 0; l < radius; l++)
                {
                  Array<ElementIndex> next;
                  for (auto ei : front)
                    for (auto pi : mesh[ei].PNums())
                      for (auto nb : top.GetVertexElements(pi))
                        if (!isbad[nb]) { isbad[nb] = true; next.Append (nb); }
                  front = std::move(next);
                }
              double factor = (change < 1e-3) ? 10 : 2;
              size_t nlowered = 0;
              lam_min = lam_start;
              auto lower_lam = [&](double & lam, auto els)
                {
                  if (lam > lam_floor)
                    for (auto ei : els)
                      if (isbad[ei]) { lam = max (lam_floor, lam/factor); nlowered++; break; }
                  lam_min = min (lam_min, lam);
                };
              for (auto e : ereleased.Range()) if (ereleased[e]) lower_lam (elambda[e], edge2el[e]);
              for (auto f : freleased.Range()) if (freleased[f]) lower_lam (flambda[f], face2el[f]);
              {
                double bmin = 1e300;
                for (auto ei : bad) bmin = min (bmin, lower[ei]);
                PrintMessage (6, "  release step: ", bad.Size(), " invalid, min ", bmin, ", lowered ", nlowered,
                              ", radius ", radius, ", change ", change);
              }
              // everything near the invalid elements is at the floor: widen,
              // then one last full optimization
              if (nlowered == 0)
                {
                  if (radius < 3) radius++;
                  else last = true;
                }
            }
          if (!done)
            for (auto ei : bad)
              {
                int ef = 0, er = 0, efr = 0, ff = 0, fr = 0, ffr = 0, inblk = 0;
                for (auto e : top.GetEdges(ei))
                  { ef += efixed[e]; er += ereleased[e]; efr += efree[e]; inblk += (eblock[e] >= 0); }
                for (auto f : top.GetFaces(ei))
                  { ff += ffixed[f]; fr += freleased[f]; ffr += ffree[f]; inblk += (fblock[f] >= 0); }
                PrintMessage (6, "   stuck el ", ei, " q ", lower[ei], " edges fixed/released/free ", ef, "/", er, "/", efr,
                              " faces ", ff, "/", fr, "/", ffr, " in blocks ", inblk, " unfixable ", UnfixableByBending (ei, params.eps));
              }
          PrintMessage (5, "extend curving release pass ", pass, ": ", nnew,
                        " entities released, min lambda ", lam_min, done ? ", valid" : ", still invalid");
          if (done) break;
        }

      // geometric error on the touched surface elements
      for (auto & [seinr, old] : surf_samples)
        {
          Array<Point<3>> now;
          SurfaceSamples (SurfaceElementIndex::FromNr0(seinr), now);
          for (int i = 0; i < now.Size(); i++)
            info.geometry_error = max (info.geometry_error, Dist (now[i], old[i]));
        }
    }

    // ------------------------------------------------------------
    //   straightening
    // ------------------------------------------------------------

    void Straighten ()
    {
      struct Cand { Block bl; int depth; double rel; };
      Array<Cand> cands;
      for (auto e : efree.Range())
        if (efree[e])
          {
            Block bl{false, e.Nr0()};
            double mx = 0;
            for (int j = 0; j < NCoefs(bl); j++) mx = max (mx, Coefs(bl)[j].Length());
            if (mx > 0) cands.Append ({bl, edepth[e], mx / EntitySize(bl)});
          }
      for (auto f : ffree.Range())
        if (ffree[f])
          {
            Block bl{true, f.Nr0()};
            double mx = 0;
            for (int j = 0; j < NCoefs(bl); j++) mx = max (mx, Coefs(bl)[j].Length());
            if (mx > 0) cands.Append ({bl, fdepth[f], mx / EntitySize(bl)});
          }
      if (params.straighten_tol > 0)
        {
          for (auto e : efixed.Range())
            if (efixed[e] && !ereleased[e])
              {
                Block bl{false, e.Nr0()};
                double mx = 0;
                for (int j = 0; j < NCoefs(bl); j++) mx = max (mx, Coefs(bl)[j].Length());
                if (mx > 0 && mx < params.straighten_tol * EntitySize(bl))
                  cands.Append ({bl, -1, mx / EntitySize(bl)});
              }
          for (auto f : ffixed.Range())
            if (ffixed[f] && !freleased[f])
              {
                Block bl{true, f.Nr0()};
                double mx = 0;
                for (int j = 0; j < NCoefs(bl); j++) mx = max (mx, Coefs(bl)[j].Length());
                if (mx > 0 && mx < params.straighten_tol * EntitySize(bl))
                  cands.Append ({bl, -1, mx / EntitySize(bl)});
              }
        }

      std::sort (cands.begin(), cands.end(), [](const Cand & a, const Cand & b)
                 { return a.depth > b.depth || (a.depth == b.depth && a.rel < b.rel); });

      ArrayMem<Vec<3>,8> save;
      for (auto & cand : cands)
        {
          int nc = NCoefs (cand.bl);
          Vec<3> * c = Coefs (cand.bl);
          save.SetSize (nc);
          for (int j = 0; j < nc; j++) { save[j] = c[j]; c[j] = Vec<3>(0,0,0); }

          bool ok = true;
          ArrayMem<double,16> newlo, newup;
          for (auto ei : Elements(cand.bl))
            {
              double t = min (params.quality, lower[ei] - 1e-12);
              double lo, up;
              ElementBounds (ei, t, lo, up);
              newlo.Append (lo); newup.Append (up);
              if (lo < t) { ok = false; break; }
            }
          if (ok)
            {
              int i = 0;
              for (auto ei : Elements(cand.bl))
                { lower[ei] = newlo[i]; upper[ei] = newup[i]; i++; }
              if (cand.depth == -1)
                for (int j = 0; j < nc; j++)
                  info.geometry_error = max (info.geometry_error, save[j].Length());
            }
          else
            for (int j = 0; j < nc; j++) c[j] = save[j];
        }
    }

    // ------------------------------------------------------------
    //   driver
    // ------------------------------------------------------------

    void Run ()
    {
      static Timer t("ExtendCurving"); RegionTimer reg(t);
      static Timer tsetup("ExtendCurving - setup");
      static Timer tbend("ExtendCurving - bend");
      static Timer trel("ExtendCurving - release");
      static Timer tstr("ExtendCurving - straighten");

      tsetup.Start();
      Setup ();
      tsetup.Stop();
      for (auto ei : mesh.VolumeElements().Range())
        {
          if (ce.IsElementCurved (ei)) info.ncurved_before++;
          if (lower[ei] < params.quality) info.nbad++;
        }

      if (params.move_vertices && info.nbad > 0)
        MakeRoom ();

      tbend.Start();
      BendInterior ();
      tbend.Stop();

      Classify ();

      if (params.reduce && info.nbad > 0)
        {
          trel.Start();
          ReleaseBoundary ();
          trel.Stop();
        }

      if (info.nbad > 0)
        {
          tstr.Start();
          Straighten ();
          tstr.Stop();
        }

      PrintMessage (5, "extend curving block time (thread sum): collect ", double(t_collect),
                    " s, assemble ", double(t_assemble), " s, solve+linesearch ", double(t_solve), " s");
      ce.CompressCoefficients ();
      for (auto ei : mesh.VolumeElements().Range())
        {
          if (lower[ei] < params.eps)
            {
              if (upper[ei] >= params.eps)
                {
                  ElementBounds (ei, params.eps, lower[ei], upper[ei], params.bezier_depth + 3);
                  if (lower[ei] >= params.eps) { if (ce.IsElementCurved (ei)) info.ncurved_after++; continue; }
                  if (upper[ei] >= params.eps) info.nundecided++;
                }
              info.nremaining++;
            }
          if (ce.IsElementCurved (ei)) info.ncurved_after++;
        }
    }
  };



  void CurvingOptimizer :: ElementBounds (ElementIndex ei, double t, double & lo, double & up, int depth) const
  { ElementBoundsImpl (ei, t, lo, up, depth); }


  ExtendCurvingInfo ExtendCurving (Mesh & mesh, const ExtendCurvingParameters & params)
  {
    if (mesh.GetDimension() != 3)
      throw NgException ("ExtendCurving: only 3D volume meshes are supported");
    auto & ce = mesh.GetCurvedElements();
    if (ce.GetOrder() <= 1 || ce.IsRational())
      return ExtendCurvingInfo();

    const_cast<Mesh&>(mesh).UpdateTopology();
    CurvingOptimizer opt (mesh, params);
    opt.Run ();

    PrintMessage (3, "ExtendCurving: ", opt.info.nbad, " bad elements, ",
                  opt.info.nreleased, " boundary entities released, geometry error ",
                  opt.info.geometry_error, ", ", opt.info.nremaining, " elements remaining");
    return opt.info;
  }


  void CurvedElementQuality (const Mesh & mesh,
                             Array<double, ElementIndex> & lower,
                             Array<double, ElementIndex> & upper,
                             int bezier_depth)
  {
    ExtendCurvingParameters params;
    params.bezier_depth = bezier_depth;
    const_cast<Mesh&>(mesh).UpdateTopology();
    CurvingOptimizer opt (const_cast<Mesh&>(mesh), params);
    lower.SetSize (mesh.GetNE());
    upper.SetSize (mesh.GetNE());
    ParallelForRange (mesh.VolumeElements().Range(), [&](auto r)
      {
        for (auto ei : r)
          opt.ElementBounds (ei, std::numeric_limits<double>::quiet_NaN(), lower[ei], upper[ei]);
      });
  }


  ExtendCurvingInfo Mesh :: ExtendCurving (const ExtendCurvingParameters & params)
  {
    auto info = netgen::ExtendCurving (*this, params);
    for (SurfaceElementIndex sei : SurfaceElements().Range())
      (*this)[sei].SetCurved (GetCurvedElements().IsCurved (sei));
    for (ElementIndex ei : VolumeElements().Range())
      (*this)[ei].SetCurved (GetCurvedElements().IsCurved (ei));
    SetNextMajorTimeStamp();
    return info;
  }

}

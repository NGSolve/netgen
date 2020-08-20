#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <set>

#include "csg2d.hpp"

// Polygon clipping algorithm based on:
// Foster, Erich & Hormann, Kai & Popa, Romeo. (2019). Clipping Simple Polygons with Degenerate Intersections. Computers & Graphics: X. 2. 100007. 10.1016/j.cagx.2019.100007.
// extended to handle quadratic spline segments

namespace netgen
{

constexpr static double EPSILON=0.000000001;

void ToggleLabel(EntryExitLabel& status)
{
  if (status == ENTRY)
  {
    status = EXIT;
    return;
  }
  if (status == EXIT)
  {
    status = ENTRY;
    return;
  }
}

Spline Split( const Spline & s, double t0, double t1 )
{
  if(t0==0.0 && t1==1.0) return s;

  Point<2> a = s.StartPI();
  if(t0!=0.0)
    a = s.GetPoint(t0);

  Point<2> c = s.EndPI();
  if(t1!=1.0)
    c = s.GetPoint(t1);

  // Find new midpoints by cutting the tangents at the new end points
  auto tang0 = s.GetTangent(t0);
  auto tang1 = s.GetTangent(t1);

  netgen::Mat<2,2> m, minv;
  m(0,0) = tang0[0];
  m(1,0) = tang0[1];
  m(0,1) = -tang1[0];
  m(1,1) = -tang1[1];

  CalcInverse(m, minv);

  Vec<2> lam = minv*(c-a);

  Point<2> b = a+lam[0]*tang0;

  auto res = Spline{a, b, c};

  // compute weight of new spline such that p lies on it
  Point<2> p = s.GetPoint(0.5*(t0+t1));
  double A = (p[1]-a[1])*(b[0]-p[0]) - (p[0]-a[0])*(b[1]-p[1]);
  double B = (p[1]-c[1])*(b[0]-p[0]) - (p[0]-c[0])*(b[1]-p[1]);
  double det = sqrt(-A*B);
  double tt = (B-det)/(A+det);
  auto v = b-p;
  int dim = fabs(v[0]) > fabs(v[1]) ? 0 : 1;
  double weight = fabs(tt*(p[dim]-a[dim])/v[dim] + 1.0/tt*(p[dim]-c[dim])/v[dim]);
  res.SetWeight(weight);
  return res;
}

Vertex * Vertex :: Insert(Point<2> p, double lam)
{
  auto vnew = make_unique<Vertex>(p);
  vnew->lam = lam;

  Vertex * current = this;

  if(lam > -1.0)
  {
    do {
      current = current->next;
    } while (!current->is_source && current->lam < lam);
  }
  else
    current = current->next;

  auto pre = current->prev;
  vnew->bc = pre->bc;

  pre->next = vnew.get();
  vnew->prev = pre;
  vnew->next = current;

  vnew->pnext = std::move(current->prev->pnext);

  current->prev = vnew.get();

  pre->pnext = std::move(vnew);

  return pre->next;
}

IntersectionType ClassifyNonOverlappingIntersection( double alpha, double beta )
{
  // classify alpha
  bool alpha_is_0 = false;
  bool alpha_in_0_1 = false;

  if ( (alpha > EPSILON) && (alpha < 1.0-EPSILON) )
    alpha_in_0_1 = true;
  else
    if (fabs(alpha) <= EPSILON)
      alpha_is_0 = true;

  // classify beta
  bool beta_is_0 = false;
  bool beta_in_0_1 = false;

  if ( (beta > EPSILON) && (beta < 1.0-EPSILON) )
    beta_in_0_1 = true;
  else
    if (fabs(beta) <= EPSILON)
      beta_is_0 = true;

  // distinguish intersection types
  if (alpha_in_0_1 && beta_in_0_1)
    return (X_INTERSECTION);

  if (alpha_is_0 && beta_in_0_1)
    return (T_INTERSECTION_Q);

  if (beta_is_0 && alpha_in_0_1)
    return (T_INTERSECTION_P);

  if (alpha_is_0 && beta_is_0)
    return (V_INTERSECTION);

  return NO_INTERSECTION;
}

IntersectionType ClassifyOverlappingIntersection( double alpha, double beta )
{
  // classify alpha
  bool alpha_is_0 = false;
  bool alpha_in_0_1 = false;
  bool alpha_not_in_0_1 = false;

  if ( (alpha > EPSILON) && (alpha < 1.0-EPSILON) )
    alpha_in_0_1 = true;
  else
    if (fabs(alpha) <= EPSILON)
      alpha_is_0 = true;
    else
      alpha_not_in_0_1 = true;

  // classify beta
  bool beta_is_0 = false;
  bool beta_in_0_1 = false;
  bool beta_not_in_0_1 = false;

  if ( (beta > EPSILON) && (beta < 1.0-EPSILON) )
    beta_in_0_1 = true;
  else
    if (fabs(alpha) <= EPSILON)
      beta_is_0 = true;
    else
      beta_not_in_0_1 = true;

  // distinguish intersection types
  if (alpha_in_0_1 && beta_in_0_1)
    return (X_OVERLAP);

  if (alpha_not_in_0_1 && beta_in_0_1)
    return (T_OVERLAP_Q);

  if (beta_not_in_0_1 && alpha_in_0_1)
    return (T_OVERLAP_P);

  if (alpha_is_0 && beta_is_0)
    return (V_OVERLAP);

  return NO_INTERSECTION;
}

IntersectionType intersect(const Point<2> P1, const Point<2> P2, const Point<2> Q1, const Point<2> Q2, double& alpha, double& beta)
{
  double AP1 = Area(P1,Q1,Q2);
  double AP2 = Area(P2,Q1,Q2);

  if (fabs(AP1-AP2) > EPSILON)
  {
    // (P1,P2) and (Q1,Q2) are not parallel

    double AQ1 = Area(Q1,P1,P2);
    double AQ2 = Area(Q2,P1,P2);

    alpha = AP1 / (AP1-AP2);
    beta  = AQ1 / (AQ1-AQ2);

    return ClassifyNonOverlappingIntersection(alpha, beta);
  }
  else
    if (fabs(AP1) < EPSILON)
    {
      // (P1,P2) and (Q1,Q2) are collinear

      auto dP = P2-P1;
      auto dQ = Q2-Q1;
      auto PQ = Q1-P1;

      alpha = (PQ*dP) / (dP*dP);
      beta = -(PQ*dQ) / (dQ*dQ);

      return ClassifyOverlappingIntersection(alpha, beta);
    }
  return NO_INTERSECTION;
}

IntersectionType IntersectSplineSegment( const Spline & s, const Point<2> & r0, const Point<2> & r1, double& alpha, double& beta )
{
  Point<2> p0 = s.StartPI();
  Point<2> p1 = s.TangentPoint();
  Point<2> p2 = s.EndPI();

  auto vr = r1-r0;
  double a0 = vr[1]*(p0[0] - r0[0]) - vr[0]*(p0[1] - r0[1]);
  double a1 = vr[1]*(p1[0] - r0[0]) - vr[0]*(p1[1] - r0[1]);
  double a2 = vr[1]*(p2[0] - r0[0]) - vr[0]*(p2[1] - r0[1]);
  a1 *= s.GetWeight();

  double a_ = a0-a1+a2;
  double b_ = a1-2*a0;
  double c_ = a0;

  double det =  b_*b_ - 4*a_*c_;
  if(det<0.0)
    return NO_INTERSECTION;
  double sqrt_det =  sqrt(det);
  double t1 = 1.0/(2*a_) * (-b_ + sqrt_det);
  double t2 = 1.0/(2*a_) * (-b_ - sqrt_det);

  double t = min(t1,t2);
  if(t<alpha)
    t = max(t1,t2);

  if(t+EPSILON<alpha)
    return NO_INTERSECTION;

  alpha = t;

  int dim = fabs(vr[0]) > fabs(vr[1]) ? 0 : 1;
  beta = 1.0/vr[dim] * (s.GetPoint(t)[dim] - r0[dim]);

  return ClassifyNonOverlappingIntersection(alpha, beta);
}

IntersectionType IntersectSplineSegment1( const Spline & s, const Point<2> & r0, const Point<2> & r1, double& alpha, double& beta )
{
  Point<2> p0 = s.StartPI();
  Point<2> p1 = s.TangentPoint();
  Point<2> p2 = s.EndPI();

  auto vr = r1-r0;
  double a0 = vr[1]*(p0[0] - r0[0]) - vr[0]*(p0[1] - r0[1]);
  double a1 = vr[1]*(p1[0] - r0[0]) - vr[0]*(p1[1] - r0[1]);
  double a2 = vr[1]*(p2[0] - r0[0]) - vr[0]*(p2[1] - r0[1]);
  a1 *= s.GetWeight();

  double a_ = a0-a1+a2;
  double b_ = a1-2*a0;
  double c_ = a0;

  double det =  b_*b_ - 4*a_*c_;
  if(det<0.0)
    return NO_INTERSECTION;
  double sqrt_det =  sqrt(det);
  double vbeta[2];
  vbeta[0] = 1.0/(2*a_) * (-b_ + sqrt_det);
  vbeta[1] = 1.0/(2*a_) * (-b_ - sqrt_det);

  int dim = fabs(vr[0]) > fabs(vr[1]) ? 0 : 1;
  double valpha[2];
  valpha[0] = 1.0/vr[dim] * (s.GetPoint(vbeta[0])[dim] - r0[dim]);
  valpha[1] = 1.0/vr[dim] * (s.GetPoint(vbeta[1])[dim] - r0[dim]);


  IntersectionType vtype[2];
  vtype[0] = ClassifyNonOverlappingIntersection(valpha[0], vbeta[0]);
  vtype[1] = ClassifyNonOverlappingIntersection(valpha[1], vbeta[1]);

  if(valpha[0]>valpha[1])
  {
    swap(valpha[0], valpha[1]);
    swap(vbeta[0], vbeta[1]);
    swap(vtype[0], vtype[1]);
  }

  int choice = 0;
  if(vtype[0]==NO_INTERSECTION && vtype[1]!=NO_INTERSECTION)
    choice = 1;

  if(valpha[0] < alpha+EPSILON)
    choice = 1;

  if(valpha[choice] < alpha+EPSILON)
    return NO_INTERSECTION;

  alpha = valpha[choice];
  beta = vbeta[choice];
  return vtype[choice];
}

bool IsOverlapping( Spline p, Spline s, double & alpha, double & beta, IntersectionType & type )
{

  auto p_mid = Center(p.StartPI(), p.EndPI());
  auto s_mid = Center(s.StartPI(), s.EndPI());

  double lam0 = -1e3*EPSILON;
  double lam1 = -1e3*EPSILON;
  alpha=-1e8;
  beta=-1e8;

  // Check if s.p0 lies on p and vice versa, also check if tangents are in same direction (TODO: TEST)
  // If so, assume overlapping splines
  // TODO: Better checks! False positives could happen here!
  IntersectSplineSegment1( p, s.StartPI(), p_mid, lam0, alpha );
  IntersectSplineSegment1( s, p.StartPI(), s_mid, lam1, beta );
  auto tang0 = s.GetTangent(0.);
  auto tang1 = p.GetTangent(alpha);
  double err = tang0*tang1;
  err*=err;
  err *= 1.0/(tang0.Length2()*tang1.Length2());

  if(fabs(lam0) < 1e3*EPSILON && fabs(lam1) < 1e3*EPSILON /*&& err < EPSILON*/)
  {
    type = ClassifyOverlappingIntersection( alpha, beta );
    return true;
  }
  return false;
}

bool IsInsideTrig( const array<Point<2>,3> & t, Point<2> r )
{
  int w = 0;
  Point<2> trig[4] = {t[0],t[1],t[2],t[0]};
  for(auto i : Range(3))
    w += CalcSide(trig[i], trig[i+1], r);
  return ( (w % 2) != 0 );
}


IntersectionType IntersectTrig( Point<2> p0, Point<2> p1, const array<Point<2>,3> & trig)
{
  Point<2> lt[4] = { trig[0], trig[1], trig[2], trig[0] };

  double alpha, beta;
  for(auto i : IntRange(3))
  {
    auto type = intersect(p0, p1, lt[i], lt[i+1], alpha, beta);
    if(type != NO_INTERSECTION)
      return type;
  }

  return NO_INTERSECTION;
}

bool IntersectTrigs( const array<Point<2>,3> & trig0, const array<Point<2>,3> & trig1)
{
  Point<2> lt0[4] = { trig0[0], trig0[1], trig0[2], trig0[0] };

  for(auto i : IntRange(3))
  {
    if(IntersectTrig(lt0[i], lt0[i+1], trig1))
      return true;
    if(IsInsideTrig(trig0, trig1[i]))
      return true;
    if(IsInsideTrig(trig1, trig0[i]))
      return true;
  }
  return false;
}

bool BisectIntersect( Spline p, Spline s, double &t0, double &t1, double &s0, double &s1, int depth=-50)
{
  if(depth==0)
  {
    s0 = s1;
    t0 = t1;
    return true;
  }

  bool side = depth%2==0;

  double & lam0 = side ? t0 : s0;
  double & lam1 = side ? t1 : s1;
  Spline & spline = side ? p : s;
  Spline & spline_other = side ? s : p;

  double lam_mid = 0.5*(lam0+lam1);
  auto left = Split(spline, lam0, lam_mid);
  auto right = Split(spline, lam_mid, lam1);

  double & lam0_other = side ? s0 : t0;
  double & lam1_other = side ? s1 : t1;
  auto curr = Split(spline_other, lam0_other, lam1_other);

  bool left_hull_intersecting = IntersectTrigs( {left.StartPI(), left.TangentPoint(), left.EndPI()}, {curr.StartPI(), curr.TangentPoint(), curr.EndPI()});
  bool right_hull_intersecting = IntersectTrigs( {right.StartPI(), right.TangentPoint(), right.EndPI()}, {curr.StartPI(), curr.TangentPoint(), curr.EndPI()});

  // TODO: Additionaly check if one spline intersects with convex hull of other?
  //   // Check if one spline intersects with convex hull of spline
  //   if(left_hull_intersecting)
  //   {
  //     double a,b;
  //     left_hull_intersecting  = left.Intersect( curr.p0, curr.p1, a, b );
  //     left_hull_intersecting |= left.Intersect( curr.p1, curr.p2, a, b );
  //     left_hull_intersecting |= left.Intersect( curr.p2, curr.p0, a, b );
  //   }
  //
  //   if(right_hull_intersecting)
  //   {
  //     double a,b;
  //     right_hull_intersecting  = right.Intersect( curr.p0, curr.p1, a, b );
  //     right_hull_intersecting |= right.Intersect( curr.p1, curr.p2, a, b );
  //     right_hull_intersecting |= right.Intersect( curr.p2, curr.p0, a, b );
  //   }


  if(!left_hull_intersecting && !right_hull_intersecting)
    return false;

  if(left_hull_intersecting && right_hull_intersecting)
  {
    // cout << "intersect both sides " << endl;
    double temp_lam;
    temp_lam = lam1;
    lam1 = lam_mid;

    double t0_ = t0;
    double t1_ = t1;
    double s0_ = s0;
    double s1_ = s1;

    // cout << "recursive bisect " << t0 << ',' << t1 << ',' << s0 << ',' << s1 << endl;
    bool first_intersecting = BisectIntersect(p, s, t0_, t1_, s0_, s1_, depth+1);
    if(first_intersecting)
    {
      t0 = t0_;
      t1 = t1_;
      s0 = s0_;
      s1 = s1_;
      return true;
    }
    else
    {
      // cout << "search other side " << endl;
      // no first intersection -> search other side
      lam1 = temp_lam;
      left_hull_intersecting = false;
    }
  }

  if(left_hull_intersecting)
    lam1 = lam_mid;
  else
    lam0 = lam_mid;

  return BisectIntersect(p, s, t0, t1, s0, s1, depth+1);
}

bool NewtonIntersect( Spline p, Spline s, double & alpha, double & beta )
{

  Point<2> p0, s0;
  Vec<2> dp, ds, ddp, dds;

  p.GetDerivatives(alpha, p0, dp, ddp);
  s.GetDerivatives(beta,  s0, ds, dds);

  netgen::Mat<2,2> m, minv;

  m(0,0) = dp[0];
  m(1,0) = dp[1];
  m(0,1) = -ds[0];
  m(1,1) = -ds[1];

  CalcInverse(m, minv);

  Vec<2> res = s0-p0;
  Vec<2> h = minv*res;
  alpha +=h[0];
  beta  +=h[1];
  return true;
}


IntersectionType Intersect( Spline p, Spline s, double &alpha, double &beta)
{
  bool is_convex_hull_intersecting = IntersectTrigs( {p.StartPI(), p.TangentPoint(), p.EndPI()}, {s.StartPI(), s.TangentPoint(), s.EndPI()});
  if(!is_convex_hull_intersecting)
    return NO_INTERSECTION;

  {
    // Check if splines overlap
    double alpha_ = alpha;
    double beta_ = beta;
    IntersectionType overlap_type;
    bool have_overlap = IsOverlapping( p, s, alpha_, beta_, overlap_type );
    if(have_overlap)
    {
      alpha = alpha_;
      beta = beta_;
      return overlap_type;
    }
  }

  // Bisection
  double t1 = 1.0;
  double s1 = 1.0;

  bool have_intersection = false;
  if(alpha>0.0) // alpha > 0 means, we have found one intersection already
  {
    // reverse parametrization of first spline to make sure, we find the second intersection first
    auto p_ = Spline{p.EndPI(), p.TangentPoint(), p.StartPI(), p.GetWeight()};
    t1 = 1.0-alpha;
    alpha = 0.0;
    beta = 0.0;

    have_intersection = BisectIntersect(p_,s,alpha,t1,beta,s1);
    alpha = 1.0-alpha;
  }
  else
    have_intersection = BisectIntersect(p,s,alpha,t1,beta,s1);

  if(have_intersection)
  {
    for(auto i : IntRange(10))
      NewtonIntersect(p, s, alpha, beta);
    return ClassifyNonOverlappingIntersection( alpha, beta );
  }

  return NO_INTERSECTION;
}


IntersectionType intersect(const Edge& edgeP, const Edge& edgeQ, double& alpha, double& beta)
{
  const Point<2>& P1 = *edgeP.v0;
  const Point<2>& P2 = *edgeP.v1;
  const Point<2>& Q1 = *edgeQ.v0;
  const Point<2>& Q2 = *edgeQ.v1;

  if(edgeP.v0->spline)
  {
    if(edgeQ.v0->spline)
      return Intersect(*edgeP.v0->spline, *edgeQ.v0->spline, alpha, beta);
    else
      return IntersectSplineSegment(*edgeP.v0->spline, Q1, Q2, alpha, beta);
  }
  else
  {
    if(edgeQ.v0->spline)
      return IntersectSplineSegment1(*edgeQ.v0->spline, P1, P2, alpha, beta);
    else
      return intersect(P1, P2, Q1, Q2, alpha, beta);
  }
}

void AddIntersectionPoint(Edge edgeP, Edge edgeQ, IntersectionType i, double alpha, double beta)
{
  Point<2> I;
  Vertex* I_P;
  Vertex* I_Q;

  Vertex* P1 = edgeP.v0;
  Vertex* Q1 = edgeQ.v0;

  switch(i)
  {
    case X_INTERSECTION:
      if(edgeP.v0->spline)
        I = edgeP.v0->spline->GetPoint(alpha);
      else
        I = *edgeP.v0 + alpha*(*edgeP.v1 - *edgeP.v0);
      I_P = edgeP.v0->Insert(I, alpha);
      I_Q = edgeQ.v0->Insert(I, beta);
      I_P->Link(I_Q);
      break;

    case X_OVERLAP:
      I_Q = edgeQ.v0->Insert(*P1, beta);
      P1->Link( I_Q);

      I_P = edgeP.v0->Insert(*Q1, alpha);
      I_P->Link( Q1);
      break;

    case T_INTERSECTION_Q:
    case T_OVERLAP_Q:
      I_Q = edgeQ.v0->Insert(*P1, beta);
      P1->Link( I_Q);
      break;

    case T_INTERSECTION_P:
    case T_OVERLAP_P:
      I_P = edgeP.v0->Insert(*Q1, alpha);
      I_P->Link( Q1);
      break;

    case V_INTERSECTION:
    case V_OVERLAP:
      P1->Link(Q1);
      break;
    default:
      break;
  }
}


void ComputeIntersections(Solid2d & sp, Solid2d & sq)
{
  auto & PP = sp.polys;
  auto & QQ = sq.polys;

  for (Polygon2d& P : PP)
    for (Edge edgeP : P.Edges(SOURCE))
      for (Polygon2d& Q : QQ)
        for (Edge edgeQ : Q.Edges(SOURCE))
        {
          double alpha = 0.0;
          double beta = 0.0;
          IntersectionType i = intersect(edgeP, edgeQ, alpha, beta);
          AddIntersectionPoint(edgeP, edgeQ, i, alpha, beta);
          if(i==X_INTERSECTION && (edgeP.v0->spline || edgeQ.v0->spline))
          {
            double alpha1 = alpha+1e2*EPSILON;
            double beta1 = 0.0; //beta+1e2*EPSILON;

            // search for possible second intersection
            i = intersect(edgeP, edgeQ, alpha1, beta1);
            // cout << "second intersection " << i << ',' << alpha1 << ',' << beta1 << ',' << alpha1-alpha << ',' << beta1-beta << endl;
            if(i!=NO_INTERSECTION && alpha+EPSILON<alpha1)
            {
              // Add midpoint of two intersection points to avoid false overlap detection of splines
              // TODO: Check if this is really necessary
              auto alpha_mid = 0.5*(alpha+alpha1);
              auto beta_mid = 0.5*(beta+beta1);
              Point<2> MP;
              if(edgeP.v0->spline)
              {
                MP = edgeP.v0->spline->GetPoint(alpha_mid);
                edgeP.v0->Insert(MP, alpha_mid);
              }
              else
                MP = edgeQ.v0->spline->GetPoint(beta_mid);

              if(edgeQ.v0->spline)
                edgeQ.v0->Insert(MP, beta_mid);

              AddIntersectionPoint(edgeP, edgeQ, i, alpha1, beta1);
            }
          }
        }

  // Split splines at new vertices
  auto split_spline_at_vertex = [](Vertex *v)
  {
    if(!v->spline)
      return;
    Spline ori{*v->spline};
    Vertex * curr = v;
    do
    {
      auto next = curr->next;
      if(!curr->is_source || !next->is_source)
      {
        double t0 = curr->is_source ? 0.0 : curr->lam;
        double t1 = next->is_source ? 1.0 : next->lam;
        curr->spline = Split(ori, t0, t1);
      }
      curr = next;
    } while(!curr->is_source);
  };

  for (Polygon2d& P : PP)
    for (Vertex* v : P.Vertices(SOURCE))
      split_spline_at_vertex(v);
  for (Polygon2d& Q : QQ)
    for (Vertex* v : Q.Vertices(SOURCE))
      split_spline_at_vertex(v);
}

enum RelativePositionType
{
  LEFT,
  RIGHT,
  IS_P_m,
  IS_P_p
};

RelativePositionType oracle(bool prev, Vertex* P1, Vertex* P2, Vertex* P3)
{
  Vertex* Q;
  Point<2> q;
  if(prev)
  {
    Q = P2->neighbour->prev;
    q = *Q;
    if(Q->spline)
      q = Q->spline->TangentPoint();
  }
  else
  {
    Q = P2->neighbour->next;
    q = *Q;
    if(P2->neighbour->spline)
      q = P2->neighbour->spline->TangentPoint();
  }

  // is Q linked to P1 ?
  if ( P1->is_intersection && (P1->neighbour == Q) )
    return(IS_P_m);

  // is Q linked to P2 ?
  if ( P3->is_intersection && (P3->neighbour == Q) )
    return(IS_P_p);

  Point<2> p1 = *P1;
  Point<2> p2 = *P2;
  Point<2> p3 = *P3;

  if(P1->spline)
    p1 = P1->spline->TangentPoint();
  if(P2->spline)
    p3 = P2->spline->TangentPoint();

  // check relative position of Q with respect to chain (P1,P2,P3)
  double s1 = Area(  q, p1, p2);
  double s2 = Area(  q, p2, p3);
  double s3 = Area( p1, p2, p3);

  if (s3 > 0)
  {
    // chain makes a left turn
    if (s1 > 0 && s2 > 0)
      return(LEFT);
    else
      return(RIGHT);
  }
  else
  {
    // chain makes a right turn (or is straight)
    if (s1 < 0 && s2 < 0)
      return(RIGHT);
    else
      return(LEFT);
  }
}

void LabelIntersections(Solid2d & sp, Solid2d & sq, Solid2d & sr, bool UNION)
{
  auto & PP = sp.polys;
  auto & QQ = sq.polys;
  auto & RR = sr.polys;

  // 1) initial classification
  for (Polygon2d& P : PP)
    for (Vertex* I : P.Vertices(INTERSECTION))
    {

      // determine local configuration at this intersection vertex
      Vertex* P_m = I->prev;
      Vertex* P_p = I->next;

      // check positions of Q- and Q+ relative to (P-, I, P+)
      RelativePositionType Q_m_type = oracle(true,  P_m, I, P_p);
      RelativePositionType Q_p_type = oracle(false, P_m, I, P_p);

      // check non-overlapping cases
      if ((Q_m_type == LEFT  && Q_p_type == RIGHT) ||
          (Q_m_type == RIGHT && Q_p_type == LEFT ))
      {
        I->label = CROSSING;
      }

      if ((Q_m_type == LEFT  && Q_p_type == LEFT ) ||
          (Q_m_type == RIGHT && Q_p_type == RIGHT))
      {
        I->label = BOUNCING;
      }

      // check overlapping cases
      if ( ( (Q_p_type == IS_P_p) && (Q_m_type == RIGHT) ) ||
          ( (Q_m_type == IS_P_p) && (Q_p_type == RIGHT) ) )
        I->label = LEFT_ON;

      if ( ( (Q_p_type == IS_P_p) && (Q_m_type == LEFT) ) ||
          ( (Q_m_type == IS_P_p) && (Q_p_type == LEFT) ) )
        I->label = RIGHT_ON;

      if ( ( (Q_p_type == IS_P_p) && (Q_m_type == IS_P_m) ) ||
          ( (Q_m_type == IS_P_p) && (Q_p_type == IS_P_m) ) )
        I->label = ON_ON;

      if ( ( (Q_m_type == IS_P_m) && (Q_p_type == RIGHT) ) ||
          ( (Q_p_type == IS_P_m) && (Q_m_type == RIGHT) ) )
        I->label = ON_LEFT;

      if ( ( (Q_m_type == IS_P_m) && (Q_p_type == LEFT) ) ||
          ( (Q_p_type == IS_P_m) && (Q_m_type == LEFT) ) )
        I->label = ON_RIGHT;
    }

  // 2) classify intersection chains
  for (Polygon2d& P : PP)
    for (Vertex* I : P.Vertices(INTERSECTION))
    {

      // start of an intersection chain ?
      if (I->label == LEFT_ON ||
          I->label == RIGHT_ON)
      {

        // remember status of the first chain vertex and vertex itself
        RelativePositionType x;
        if (I->label == LEFT_ON)
          x = LEFT;
        else
          x = RIGHT;
        Vertex* X = I;

        // proceed to end of intersection chain and mark all visited vertices as NONE
        do {
          I->label = NONE;
          I = I->next;
        } while (I->label == ON_ON);

        RelativePositionType y;
        if (I->label == ON_LEFT)
          y = LEFT;
        else
          y = RIGHT;

        // determine type of intersection chain
        IntersectionLabel chainType;
        if (x != y)
          chainType = DELAYED_CROSSING;
        else
          chainType = DELAYED_BOUNCING;

        // mark both ends of an intersection chain with chainType (i.e., as DELAYED_*)
        X->label = chainType;
        I->label = chainType;
      }
    }

  // 3) copy labels from P to Q
  // loop over intersection vertices of P
  for (Polygon2d& P : PP)
    for (Vertex* I : P.Vertices(INTERSECTION))
      I->neighbour->label = I->label;

  // 3.5) check for special cases

  set<Polygon2d*> noIntersection[2];
  set<Polygon2d*> identical[2];

  for (int i=0; i<2; ++i)
  {
    Array<Polygon2d>* P_or_Q = &PP;      // if i=0, then do it for P w.r.t. Q
    Array<Polygon2d>* Q_or_P = &QQ;

    if (i==1) {                         // if i=1, then do it for Q w.r.t. P
      P_or_Q = &QQ;
      Q_or_P = &PP;
    }

    // loop over all components of P (or Q)
    for (Polygon2d& P : *P_or_Q)
      if (P.noCrossingVertex(UNION))
      {
        // P_ has no crossing vertex (but may have bounces or delayed bounces, except for UNION),
        // hence it does not intersect with Q_or_P
        noIntersection[i].insert(&P);   // remember component, and ignore it later in step 4

        // is P identical to some component of and Q_or_P?
        if (P.allOnOn())
        {
          identical[i].insert(&P);      // -> remember for further processing below
        }
        else
        {
          // is P inside Q_or_P?
          bool isInside = false;
          auto p = P.getNonIntersectionPoint();
          for (Polygon2d& Q : *Q_or_P)
            if ( Q.IsInside(p) )
              isInside = !isInside;
          if (isInside ^ UNION)
            RR.Append(P);             // -> add P to the result
        }
      }
  }

  // handle components of P that are identical to some component of Q
  for (Polygon2d* P : identical[0])
  {
    // is P a hole?
    bool P_isHole = false;
    for (Polygon2d& P_ : PP)
      if ( ( P_.first.get() != P->first.get() ) && (P_.IsInside(*P->first)) )
        P_isHole = !P_isHole;

    for (Polygon2d* Q : identical[1])
      for (Vertex* V : Q->Vertices(ALL))
        if (V == P->first->neighbour) {  // found Q that matches P
          // is Q a hole?
          bool Q_isHole = false;
          for (Polygon2d& Q_ : QQ)
            if ( ( Q_.first.get() != Q->first.get() ) && (Q_.IsInside(*Q->first)) )
              Q_isHole = !Q_isHole;

          // if P and Q are both holes or both are not holes
          if (P_isHole == Q_isHole)
            RR.Append(*P);           // -> add P to the result
          goto next_P;
        }
next_P: ;
  }

  // 4) set entry/exit flags
  set<Vertex*> split[2];                // split vertex candidates for P and Q
  set<Vertex*> crossing[2];             // CROSSING vertex candidates for P and Q

  for (int i=0; i<2; ++i)
  {
    Array<Polygon2d>* P_or_Q = &PP;      // if i=0, then do it for P w.r.t. Q
    Array<Polygon2d>* Q_or_P = &QQ;

    if (i==1) {                         // if i=1, then do it for Q w.r.t. P
      P_or_Q = &QQ;
      Q_or_P = &PP;
    }

    // loop over all components of P (or Q)
    for (Polygon2d& P : *P_or_Q)
    {

      // ignore P if it does not intersect with Q_or_P (detected in step 3.5 above)
      if(noIntersection[i].find(&P) != noIntersection[i].end())
        continue;

      // start at a non-intersection vertex of P
      Vertex* V = P.getNonIntersectionVertex();

      // check if it is inside or outside Q (or P)
      // and set ENTRY/EXIT status accordingly
      EntryExitLabel status = ENTRY;
      for (Polygon2d& Q : *Q_or_P)
        if (Q.IsInside(*V))
          ToggleLabel(status);

      // starting at V, loop over those vertices of P, that are either
      // a crossing intersection or marked as ends of an intersection chain
      bool first_chain_vertex = true;     // needed for dealing with crossing chains

      for (Vertex* I : P.Vertices(INTERSECTION, V))
      {
        // in the case of normal crossings, we...
        if (I->label == CROSSING)
        {
          // mark vertex with current ENTRY/EXIT status
          I->enex = status;
          // toggle status from ENTRY to EXIT or vice versa
          ToggleLabel(status);
        }

        // identify split vertex candidates (INTERIOR bouncing vertices)
        if ( (I->label == BOUNCING) && ((status == EXIT) ^ UNION) )
          split[i].insert(I);

        //
        // in the case of a delayed crossing chain, we
        // mark both end points of the chain with the current ENTRY/EXIT status,
        // toggling the status only at the end last chain vertex,
        // and, in case of a delayed EXIT  crossing, the first vertex
        //  or, in case of a delayed ENTRY crossing, the last  vertex,
        // of the chain as CROSSING
        //
        if (I->label == DELAYED_CROSSING)
        {
          // mark vertex with current ENTRY/EXIT status
          I->enex = status;

          if (first_chain_vertex) {       // are we at the first vertex of a delayed crossing chain?
            if ((status == EXIT) ^ UNION)
              I->label = CROSSING;        // mark first vertex as CROSSING
            first_chain_vertex = false;
          }
          else {                          // here we are at the last vertex of a delayed crossing chain
            if ((status == ENTRY) ^ UNION)
              I->label = CROSSING;        // mark last vertex as CROSSING
            first_chain_vertex = true;

            // toggle status from ENTRY to EXIT or vice versa (only for last chain vertex)
            ToggleLabel(status);
          }
        }

        //
        // in the case of a delayed bouncing chain, we
        // mark both end points of the chain with the current ENTRY/EXIT status
        // toggling the status at both end points of the chain,
        // and, in case of a delayed INTERIOR bouncing, both end points
        // of the chain as CROSSING candidates
        //
        if (I->label == DELAYED_BOUNCING)
        {
          // mark vertex with current ENTRY/EXIT status
          I->enex = status;

          if (first_chain_vertex) {       // are we at the first vertex of a delayed crossing chain?
            if ((status == EXIT) ^ UNION)
              crossing[i].insert(I);      // mark first EXIT vertex as CROSSING candidate
            first_chain_vertex = false;
          }
          else {                          // here we are at the last vertex of a delayed crossing chain
            if ((status == ENTRY) ^ UNION)
              crossing[i].insert(I);      // mark last ENTRY vertex as CROSSING candidate
            first_chain_vertex = true;

          }
          // toggle status from ENTRY to EXIT or vice versa (for first AND last chain vertex)
          ToggleLabel(status);
        }
      }
    }
  }

  // 5) handle split vertex pairs
  // loop over P's split candidates
  for (Vertex* I_P : split[0])
  {
    Vertex* I_Q = I_P->neighbour;

    // check if the neighbour on Q is also a split candidate
    if (split[1].find(I_Q) != split[1].end())
    {
      // compute areas to compare local orientation
      double sP = Area( *I_P->prev, *I_P, *I_P->next);
      double sQ = Area( *I_Q->prev, *I_Q, *I_Q->next);

      // add duplicate vertices to P and Q
      auto V_P = I_P->Insert(*I_P);
      V_P->spline = I_P->spline;
      auto V_Q = I_Q->Insert(*I_Q);
      V_Q->spline = I_Q->spline;

      // link vertices correctly
      if (sP*sQ > 0) {                  // same local orientation
        I_P->Link( V_Q);
        I_Q->Link( V_P);
      }
      else {                            // different local orientation
        V_P->Link( V_Q);
      }

      // mark all four vertices correctly
      if (!UNION)
      {
        I_P->enex = EXIT;
        V_P->enex = ENTRY;
        I_Q->enex = EXIT;
        V_Q->enex = ENTRY;
      }
      else
      {
        I_P->enex = ENTRY;
        V_P->enex = EXIT;
        I_Q->enex = ENTRY;
        V_Q->enex = EXIT;
      }

      I_P->label = CROSSING;
      V_P->label = CROSSING;
      I_Q->label = CROSSING;
      V_Q->label = CROSSING;
    }
  }

  // 6) handle CROSSING vertex candidates
  // loop over P's CROSSING candidates
  for (Vertex* I_P : crossing[0])
  {
    Vertex* I_Q = I_P->neighbour;

    // check if the neighbour on Q is also a CROSSING candidate
    if (crossing[1].find(I_Q) != crossing[1].end())
    {
      // mark CROSSING candidate pair as such
      I_P->label = CROSSING;
      I_Q->label = CROSSING;
    }
  }
}

void CreateResult(Solid2d & sp, Solid2d & sr, bool UNION)
{
  auto & PP = sp.polys;
  auto & RR = sr.polys;
  //
  // for all crossing vertices
  //
  // NOTE: all crossing vertices that are visited while contructing a
  //       component of the result polygon are marked as "not intersection",
  //       so that they cannot serve as start vertex of another component
  //

  for (Polygon2d& P : PP)
  {
    for (Vertex* I : P.Vertices(CROSSING_INTERSECTION))
    {
      Polygon2d R;                         // result polygon component

      Vertex* V = I;                      // start traversal at I
      V->is_intersection = false;            // mark visited vertices

      do {
        EntryExitLabel status = V->enex;
        ToggleLabel(status);
        while ( !(V->enex == status))    // ... we arrive at a vertex with opposite entry/exit flag, or
        {
          auto & vnew = R.AppendVertex(*V);
          if ((status == EXIT) ^ UNION)
          {
            vnew.bc = V->bc;
            if(V->spline)
              vnew.spline = *V->spline;
            else
              vnew.spline = nullopt;
            V = V->next;                  // move forward  from an ENTRY vertex to the next EXIT  vertex
            V->is_intersection = false;        // mark visited vertices
          }
          else
          {
            V = V->prev;                  // move backward from an EXIT  vertex to the next ENTRY vertex
            if(V->spline)
            {
              auto & s = *V->spline;
              vnew.spline = Spline{s.EndPI(), s.TangentPoint(), s.StartPI(), s.GetWeight()};
            }
            else
              vnew.spline = nullopt;
            vnew.bc = V->bc;
            V->is_intersection = false;        // mark visited vertices
          }
          if(V == I)
            break;
        }

        if (V != I)
        {
          V = V->neighbour;               // switch from P to Q or vice versa
          V->is_intersection = false;        // mark visited vertices
        }
      } while (V != I);                   // the result polygon component is complete,
      // if we are back to the initial vertex I
      RR.Append(R);
    }
  }
}

void CleanUpResult(Solid2d & sr)
{
  auto & RR = sr.polys;
  for (Polygon2d& R : RR)
  {
    while ( (R.first.get() != NULL) && (fabs(Area(*R.first->prev,*R.first,*R.first->next)) < EPSILON) )
      R.Remove(R.first.get());

    if (R.first.get() != NULL)
      for (Vertex* V : R.Vertices(ALL))
        if (!V->spline && !V->prev->spline && fabs(Area(*V->prev,*V,*V->next)) < EPSILON)
        {
          R.Remove(V);
        }
  }
  for (int i = RR.Size()-1; i>=0; i--)
    if(RR[i].Size()==0)
      RR.RemoveElement(i);
}

void RemoveDuplicates(Solid2d & sr)
{
  for(auto & poly : sr.polys)
  {
    if(poly.first==nullptr) continue;
    Vertex * last = poly.first->prev;
    for(auto v : poly.Vertices(ALL))
    {
      if(Dist2(*v, *last)<EPSILON*EPSILON)
        poly.Remove(last);
      last = v;
    }
  }
}

Polygon2d RectanglePoly(double x0, double x1, double y0, double y1, string bc)
{
  Polygon2d r;
  r.Append( {x0, y0} );
  r.Append( {x1, y0} );
  r.Append( {x1, y1} );
  r.Append( {x0, y1} );
  r.SetBC(bc);
  return r;
}

Solid2d Rectangle(double x0, double x1, double y0, double y1, string name, string bc)
{
  Solid2d s;
  s.name = name;
  s.polys.Append(RectanglePoly(x0,x1,y0,y1, bc));
  s.SetBC(bc);
  return s;
}

Solid2d Circle(double x, double y, double r, string name, string bc)
{
  Solid2d s;
  s.name = name;
  Polygon2d poly;

  Point<2> ps[] =
  {
    {x+r, y+0},
    {x+r, y+r},
    {x+0, y+r},
    {x-r, y+r},
    {x-r, y+0},
    {x-r, y-r},
    {x+0, y-r},
    {x+r, y-r}
  };

  for (auto i : IntRange(4))
  {
    int i0 = 2*i;
    int i1 = (i0+1)%8;
    int i2 = (i0+2)%8;
    auto & v0 = poly.Append( ps[i0] );
    v0.spline = { ps[i0], ps[i1], ps[i2] };
  }

  s.polys.Append(poly);
  s.SetBC(bc);
  return s;
}

Solid2d AddIntersectionPoints ( Solid2d s1, Solid2d s2 )
{
  ComputeIntersections(s1, s2);
  RemoveDuplicates(s1);
  return s1;
}

Solid2d ClipSolids ( Solid2d s1, Solid2d s2, bool intersect)
{
  static Timer t1("intersection");
  static Timer t2("label");
  static Timer t3("cut");
  static Timer t4("cleanup");

  for(auto & poly : s1.polys)
    for(auto v : poly.Vertices(ALL))
    {
      v->is_source = true;
      v->neighbour = nullptr;
      v->lam = -1.0;
      v->is_intersection = false;
      v->label = NONE;
      v->enex = NEITHER;
    }

  for(auto & poly : s2.polys)
    for(auto v : poly.Vertices(ALL))
    {
      v->is_source = true;
      v->neighbour = nullptr;
      v->lam = -1.0;
      v->is_intersection = false;
      v->label = NONE;
      v->enex = NEITHER;
    }

  Solid2d res;
  res.name = s1.name;

  t1.Start();
  ComputeIntersections(s1, s2);
  t1.Stop();

  t2.Start();
  LabelIntersections(s1, s2, res, !intersect);
  t2.Stop();

  t3.Start();
  CreateResult(s1, res, !intersect);
  t3.Stop();

  t4.Start();
  CleanUpResult(res);
  RemoveDuplicates(res);
  t4.Stop();

  return res;
}

Solid2d Solid2d :: operator+(Solid2d & other)
{
  if(polys.Size()==0)
    return other;

  auto res = ClipSolids(*this, other, false);
  res.name = name;
  return res;
}

Solid2d Solid2d :: operator*(Solid2d & other)
{
  auto res = ClipSolids(*this, other, true);
  res.name = name;
  return res;
}

Solid2d Solid2d :: operator-(Solid2d other)
{
  // TODO: Check dimensions of solids with bounding box
  other.Append(RectanglePoly(-1e8, 1e8, -1e8, 1e8, "JUST_FOR_CLIPPING"));
  auto res = ClipSolids(*this, other);

  for (auto i : Range(other.polys))
  {
    auto & first = *other.polys[i].first;
    if(first[0] == -1e8)
      other.polys.DeleteElement(i);
  }
  res.name = name;
  return res;
}

bool Solid2d :: IsInside( Point<2> r ) const
{
  int w = 0;
  for(auto & poly : polys)
    for(auto v : poly.Vertices(ALL))
      w += CalcSide(*v, *v->next, r);
  return ( (w % 2) != 0 );
}

bool Solid2d :: IsLeftInside( const Vertex & p0 )
{
  auto & p1 = *p0.next;
  auto v = p1-p0;
  auto n = Vec<2>{v[1], -v[0]};
  auto q = p0 + 0.5*v + 1e-6*n;
  return IsInside(q);
}

bool Solid2d :: IsRightInside( const Vertex & p0 )
{
  auto & p1 = *p0.next;
  auto v = p1-p0;
  auto n = Vec<2>{-v[1], v[0]};
  auto q = p0 + 0.5*v + 1e-6*n;
  return IsInside(q);
}


shared_ptr<netgen::SplineGeometry2d> CSG2d :: GenerateSplineGeometry()
{
  static Timer t_intersections("CSG2d - AddIntersections()");
  static Timer tall("CSG2d - GenerateSplineGeometry()");
  RegionTimer rt(tall);

  struct Seg
  {
    int p0;
    int p1;
    int left;
    int right;
    int bc;
    int p2;
    double weight;
  };

  auto geo = std::make_shared<netgen::SplineGeometry2d>();
  std::map<std::tuple<int,int,int>, Seg> seg_map;
  std::map<string, int> bcmap;
  Array<int> points;

  // Cut each solid with each other one to add all possible intersection points and have conforming edges from both domains
  // TODO: OPTIMIZE!!!
  // Idea: Find edges with just one neighbor (either leftdomain or rightdomain unset after the marking below) -> just cut those edges with each other
  t_intersections.Start();
  for(auto & s1 : solids)
    for(auto & s2 : solids)
      if(&s1!=&s2)
        s1 = AddIntersectionPoints(s1,s2);
  t_intersections.Stop();

  // Add geometry points to SplineGeometry

  netgen::Box<2> box(netgen::Box<2>::EMPTY_BOX);
  for(auto & s : solids)
    for(auto & poly : s.polys)
      for(auto v : poly.Vertices(ALL))
        box.Add(*v);

  netgen::BoxTree <2, int> ptree(box);

  auto getPoint = [&](Point<2> p )
  {
    int res = -1;
    ptree.GetFirstIntersecting(p, p, [&] (int pi)
        {
        res = pi;
        return true;
        });
    return res;
  };

  auto insertPoint = [&](Point<2> p )
  {
    int pi = getPoint(p);
    if(pi==-1)
    {
      // not found -> insert to tree
      netgen::GeomPoint<2> gp(p);
      gp.name = "";
      geo->geompoints.Append(gp);
      ptree.Insert(p,p,geo->geompoints.Size()-1);
    }
  };

  for(auto & s : solids)
    for(auto & poly : s.polys)
      for(auto v : poly.Vertices(ALL))
      {
        box.Add(*v);
        insertPoint(*v);
        if(v->spline)
          insertPoint(v->spline->TangentPoint());
      }


  // Generate segments from polygon edges and find left/right domain of each segment
  int dom = 0;
  for(auto & s : solids)
  {
    dom++;
    geo->SetMaterial(dom, s.name);
    for(auto & poly : s.polys)
    {
      for(auto v : poly.Vertices(ALL))
      {
        auto & p0 = *v;
        auto & p1 = *v->next;

        auto pi0 = getPoint(p0);
        auto pi1 = getPoint(p1);
        int pi2 = -1;
        double weight = 0.0;

        if(v->spline)
        {
          auto p2 = v->spline->TangentPoint();
          pi2 = getPoint(p2);
          weight = v->spline->GetWeight();
        }

        bool flip = false;
        if(pi1<pi0)
        {
          flip = true;
          Swap(pi1,pi0);
        }

        auto li = s.IsLeftInside(p0);
        auto ri = s.IsRightInside(p0);

        if(li!=ri)
        {
          auto & ls = seg_map[{pi0,pi1,pi2}];
          ls.p0 = pi0;
          ls.p1 = pi1;
          ls.p2 = pi2;
          ls.weight = weight;
          if(s.IsLeftInside(p0) == flip)
            ls.left = dom;
          else
            ls.right = dom;
          if(bcmap.count(p0.bc)==0)
            bcmap[p0.bc] = bcmap.size()+1;
          ls.bc = bcmap[p0.bc];
        }
      }
    }
  }

  for(auto & [name, bc] : bcmap)
  {
    geo->SetBCName(bc, name);
  }

  for(auto const &m : seg_map)
  {
    auto ls = m.second;
    netgen::SplineSegExt * seg;
    if(ls.p2!=-1)
    {
      // spline segment
      auto * seg3 = new netgen::SplineSeg3<2>( geo->GetPoint(ls.p0), geo->GetPoint(ls.p2), geo->GetPoint(ls.p1), ls.weight );
      seg = new netgen::SplineSegExt(*seg3);
    }
    else
    {
      // line segment
      auto * l = new netgen::LineSeg<2>(geo->GetPoint(ls.p0), geo->GetPoint(ls.p1));
      seg = new netgen::SplineSegExt(*l);
    }

    seg->leftdom = ls.left;
    seg->rightdom = ls.right;
    seg->bc = ls.bc;
    seg->reffak = 1;
    seg->copyfrom = -1;
    seg->hmax = 1e99;
    geo->AppendSegment(seg);
  }
  return geo;
}
}

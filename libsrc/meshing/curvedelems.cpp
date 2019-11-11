#include <mystdlib.h>

#include "meshing.hpp"

#include "../general/autodiff.hpp"


namespace netgen
{
  
  //   bool rational = true;

  
  static void ComputeGaussRule (int n, NgArray<double> & xi, NgArray<double> & wi)
  {
    xi.SetSize (n);
    wi.SetSize (n);
    
    int m = (n+1)/2;
    double p1, p2, p3;
    double pp, z, z1;
    for (int i = 1; i <= m; i++)
      {
	z = cos ( M_PI * (i - 0.25) / (n + 0.5));
	while(1)
	  {
	    p1 = 1; p2 = 0;
	    for (int j = 1; j <= n; j++)
	      {
		p3 = p2; p2 = p1;
		p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
	      }
	    // p1 is legendre polynomial
	    
	    pp = n * (z*p1-p2) / (z*z - 1);
	    z1 = z;
	    z = z1-p1/pp;
	    
	    if (fabs (z - z1) < 1e-14) break;
	  }
	
	xi[i-1] = 0.5 * (1 - z);
	xi[n-i] = 0.5 * (1 + z);
	wi[i-1] = wi[n-i] = 1.0 / ( (1  - z * z) * pp * pp);
      }
  }
  
  

  // compute edge bubbles up to order n, x \in (-1, 1)
  template <typename T>
  static void CalcEdgeShape (int n, T x, T * shape)
  {
    T p1 = x, p2 = -1, p3 = 0;
    for (int j=2; j<=n; j++)
      {
	p3=p2; p2=p1;
	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	shape[j-2] = p1;
      } 
  }
  template <typename T, typename FUNC>
  static void CalcEdgeShapeLambda (int n, T x, FUNC func)
  {
    T p1(x), p2(-1.0), p3(0.0);
    for (int j=2; j<=n; j++)
      {
	p3=p2; p2=p1;
	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	func(j-2, p1);
      } 
  }


  template <typename T>
  static void CalcEdgeDx (int n, T x, T * dshape)
  {
    T p1 = x, p2 = -1, p3 = 0;
    T p1dx = 1, p2dx = 0, p3dx = 0;

    for (int j=2; j<=n; j++)
      {
	p3=p2; p2=p1;
	p3dx = p2dx; p2dx = p1dx;

	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	p1dx = ( (2*j-3) * (x * p2dx + p2) - (j-3) * p3dx) / j;

	dshape[j-2] = p1dx;
      }    
  }

  template <typename T>
  static void CalcEdgeShapeDx (int n, T x, T * shape, T * dshape)
  {
    T p1 = x, p2 = -1, p3 = 0;
    T p1dx = 1, p2dx = 0, p3dx = 0;

    for (int j=2; j<=n; j++)
      {
	p3=p2; p2=p1;
	p3dx = p2dx; p2dx = p1dx;

	p1=( (2*j-3) * x * p2 - (j-3) * p3) / j;
	p1dx = ( (2*j-3) * (x * p2dx + p2) - (j-3) * p3dx) / j;

	shape[j-2] = p1;
	dshape[j-2] = p1dx;
      }    
  }

  // compute L_i(x/t) * t^i
  template <typename T>
  static void CalcScaledEdgeShape (int n, T x, T t, T * shape)
  {
    static bool init = false;
    static double coefs[100][2];
    if (!init)
      {
        for (int j = 0; j < 100; j++)
          {
            coefs[j][0] = double(2*j+1)/(j+2);
            coefs[j][1] = -double(j-1)/(j+2);
          }
        init = true;
      }

    T p1 = x, p2 = -1, p3 = 0;
    T tt = t*t;
    for (int j=0; j<=n-2; j++)
      {
	p3=p2; p2=p1;
        p1= coefs[j][0] * x * p2 + coefs[j][1] * tt*p3;
	// p1=( (2*j+1) * x * p2 - t*t*(j-1) * p3) / (j+2);
	shape[j] = p1;
      }    
  }
  
  template <typename T, typename FUNC>
  static void CalcScaledEdgeShapeLambda (int n, T x, T t, FUNC func)
  {
    static bool init = false;
    static double coefs[100][2];
    if (!init)
      {
        for (int j = 0; j < 100; j++)
          {
            coefs[j][0] = double(2*j+1)/(j+2);
            coefs[j][1] = -double(j-1)/(j+2);
          }
        init = true;
      }

    T p1(x), p2(-1.0), p3(0.0);
    T tt = t*t;
    for (int j=0; j<=n-2; j++)
      {
	p3=p2; p2=p1;
        p1= coefs[j][0] * x * p2 + coefs[j][1] * tt*p3;
	// p1=( (2*j+1) * x * p2 - t*t*(j-1) * p3) / (j+2);
	func(j, p1);
      }    
  }


  
  template <int DIST, typename T>
  static void CalcScaledEdgeShapeDxDt (int n, T x, T t, T * dshape)
  {
    T p1 = x, p2 = -1, p3 = 0;
    T p1dx = 1, p1dt = 0;
    T p2dx = 0, p2dt = 0;
    T p3dx = 0, p3dt = 0;
     
    for (int j=0; j<=n-2; j++)
      {
	p3=p2; p3dx=p2dx; p3dt = p2dt;
	p2=p1; p2dx=p1dx; p2dt = p1dt;

	p1   = ( (2*j+1) * x * p2 - t*t*(j-1) * p3) / (j+2);
	p1dx = ( (2*j+1) * (x * p2dx + p2) - t*t*(j-1) * p3dx) / (j+2);
	p1dt = ( (2*j+1) * x * p2dt - (j-1)* (t*t*p3dt+2*t*p3)) / (j+2);

	// shape[j] = p1;
	dshape[DIST*j  ] = p1dx;
	dshape[DIST*j+1] = p1dt;
      }    
  }


  template <class Tx, class Tres>
  static void LegendrePolynomial (int n, Tx x, Tres * values)
  {
    switch (n)
      {
      case 0:
	values[0] = 1;
	break;
      case 1:
	values[0] = 1;
	values[1] = x;
	break;

      default:

	if (n < 0) return;

	Tx p1 = 1.0, p2 = 0.0, p3;
	
	values[0] = 1.0;
	for (int j=1; j<=n; j++)
	  {
	    p3 = p2; p2 = p1;
	    p1 = ((2.0*j-1.0)*x*p2 - (j-1.0)*p3) / j;
	    values[j] = p1;
	  }
      }
  }

  template <class Tx, class Tt, class Tres>
  static void ScaledLegendrePolynomial (int n, Tx x, Tt t, Tres * values)
  {
    switch (n)
      {
      case 0:
	values[0] = 1.0;
	break;

      case 1:
	values[0] = 1.0;
	values[1] = x;
	break;

      default:

	if (n < 0) return;

	Tx p1 = 1.0, p2 = 0.0, p3;
	values[0] = 1.0;
	for (int j=1; j<=n; j++)
	  {
	    p3 = p2; p2 = p1;
	    p1 = ((2.0*j-1.0)*x*p2 - t*t*(j-1.0)*p3) / j;
	    values[j] = p1;
	  }
      }
  }

  class RecPol
  {
  protected:
    int maxorder;
    double *a, *b, *c;
  public:
    RecPol (int amaxorder)
    {
      maxorder = amaxorder;
      a = new double[maxorder+1];
      b = new double[maxorder+1];
      c = new double[maxorder+1];
    }
    ~RecPol ()
    {
      delete [] a;
      delete [] b;
      delete [] c;
    }
    
    template <class S, class T>
    void Evaluate (int n, S x, T * values)
    {
      S p1(1.0), p2(0.0), p3;
      
      if (n >= 0) 
	p2 = values[0] = 1.0;
      if (n >= 1) 
	p1 = values[1] = a[0]+b[0]*x;
      
      for (int i  = 1; i < n; i++)
	{
	  p3 = p2; p2=p1;
	  p1 = (a[i]+b[i]*x)*p2-c[i]*p3;
	  values[i+1] = p1;
	}
    }

    template <class S, class T>
    void EvaluateScaled (int n, S x, S y, T * values)
    {
      S p1(1.0), p2(0.0), p3;
      
      if (n >= 0) 
	p2 = values[0] = 1.0;
      if (n >= 1) 
	p1 = values[1] = a[0]*y+b[0]*x;
      
      for (int i  = 1; i < n; i++)
	{
	  p3 = p2; p2=p1;
	  p1 = (a[i]*y+b[i]*x)*p2-c[i]*y*y*p3;
	  values[i+1] = p1;
	}
    }

    template <class S, class FUNC>
    void EvaluateScaledLambda (int n, S x, S y, FUNC func)
    {
      S p1(1.0), p2(0.0), p3;
      
      if (n >= 0)
        {
          p2 = 1.0;
          func(0, p2);
        }
      if (n >= 1)
        {
          p1 = a[0]*y+b[0]*x;
          func(1, p1);
        }
      
      for (int i = 1; i < n; i++)
	{
	  p3 = p2; p2=p1;
	  p1 = (a[i]*y+b[i]*x)*p2-c[i]*y*y*p3;
	  func(i+1, p1);
	}
    }


  };

  class JacobiRecPol : public RecPol
  {
  public:
    JacobiRecPol (int amo, double al, double be)
      : RecPol (amo)
    {
      for (int i = 0; i <= maxorder; i++)
	{
	  double den = 2*(i+1)*(i+al+be+1)*(2*i+al+be);
	  a[i] = (2*i+al+be+1)*(al*al-be*be) / den;
	  b[i] = (2*i+al+be)*(2*i+al+be+1)*(2*i+al+be+2) / den;
	  c[i] = 2*(i+al)*(i+be)*(2*i+al+be+2) / den;
	}
    }
  };



  template <class S, class T>
  inline void JacobiPolynomial (int n, S x, double alpha, double beta, T * values)
  {
    S p1 = 1.0, p2 = 0.0, p3;

    if (n >= 0) 
      p2 = values[0] = 1.0;
    if (n >= 1) 
      p1 = values[1] = 0.5 * (2*(alpha+1)+(alpha+beta+2)*(x-1));

    for (int i  = 1; i < n; i++)
      {
        p3 = p2; p2=p1;
        p1 =
          1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
          ( 
           ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) + 
             (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
           * p2
           - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * p3
           );
        values[i+1] = p1;
      }
  }

  
  

  template <class S, class St, class T>
  inline void ScaledJacobiPolynomial (int n, S x, St t, double alpha, double beta, T * values)
  {
    /*
      S p1 = 1.0, p2 = 0.0, p3;

      if (n >= 0) values[0] = 1.0;
    */

    S p1(1.0), p2(0.0), p3;

    if (n >= 0) 
      p2 = values[0] = 1.0;
    if (n >= 1) 
      p1 = values[1] = 0.5 * (2*(alpha+1)*t+(alpha+beta+2)*(x-t));

    for (int i=1; i < n; i++)
      {
        p3 = p2; p2=p1;
        p1 =
          1.0 / ( 2 * (i+1) * (i+alpha+beta+1) * (2*i+alpha+beta) ) *
          ( 
           ( (2*i+alpha+beta+1)*(alpha*alpha-beta*beta) * t + 
             (2*i+alpha+beta)*(2*i+alpha+beta+1)*(2*i+alpha+beta+2) * x) 
           * p2
           - 2*(i+alpha)*(i+beta) * (2*i+alpha+beta+2) * t * t * p3
           );
        values[i+1] = p1;
      }
  }

  
  static NgArray<shared_ptr<RecPol>> jacpols2;

  void CurvedElements::buildJacPols()
  {
    if (!jacpols2.Size())
      {
	jacpols2.SetSize (100);
	for (int i = 0; i < 100; i++)
	  jacpols2[i] = make_shared<JacobiRecPol> (100, i, 2);
      }
  }

  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
  template <class Tx, class Ty, class Ts>
  static void CalcTrigShape (int n, Tx x, Ty y, Ts * shape)
  {
    // cout << "calc trig shape" << endl;
    if (n < 3) return;
    Tx hx[50], hy[50*50];

    jacpols2[2] -> EvaluateScaled (n-3, x, 1-y, hx);

    for (int ix = 0; ix <= n-3; ix++)
      jacpols2[2*ix+5] -> Evaluate (n-3, 2*y-1, hy+50*ix);

    int ii = 0;

    Tx bub = (1+x-y)*y*(1-x-y);
    for (int ix = 0; ix <= n-3; ix++)
      hx[ix] *= bub;

    /*
    for (int iy = 0; iy <= n-3; iy++)
      for (int ix = 0; ix <= n-3-iy; ix++)
	shape[ii++] = hx[ix]*hy[iy+50*ix];
    */
    // change loops:
    for (int ix = 0; ix <= n-3; ix++)
      for (int iy = 0; iy <= n-3-ix; iy++)
	shape[ii++] = hx[ix]*hy[iy+50*ix];
  }

  template <typename T>
  static void CalcTrigShapeDxDy (int n, T x, T y, T * dshape)
  { 
    if (n < 3) return;

    AutoDiff<2,T> adx(x, 0);
    AutoDiff<2,T> ady(y, 1);
    AutoDiff<2,T> res[2000];
    CalcTrigShape (n, adx, ady, &res[0]);
    int ndof = (n-1)*(n-2)/2;
    for (int i = 0; i < ndof; i++)
      {
	dshape[2*i] = res[i].DValue(0);
	dshape[2*i+1] = res[i].DValue(1);
      }
  }


  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
  template <class Tx, class Ty, class Tt, class Tr>
  static void CalcScaledTrigShape (int n, Tx x, Ty y, Tt t, Tr * shape)
  {
    if (n < 3) return;

    Tx hx[50], hy[50];
    ScaledJacobiPolynomial (n-3, x, t-y, 2, 2, hx);

    int ii = 0;
    Tx bub = (t+x-y)*y*(t-x-y);
    for (int ix = 0; ix <= n-3; ix++)
      {
        jacpols2[2*ix+5] -> EvaluateScaled (n-3, 2*y-1, t, hy);
        for (int iy = 0; iy <= n-3-ix; iy++)
          shape[ii++] = bub * hx[ix]*hy[iy];
      }
  }

  template <class Tx, class Ty, class Tt, typename FUNC>
  static void CalcScaledTrigShapeLambda (int n, Tx x, Ty y, Tt t, FUNC func)
  {
    if (n < 3) return;
    int ii = 0;
    Tx bub = (t+x-y)*y*(t-x-y);
    jacpols2[2]->EvaluateScaledLambda
      (n-3, x, t-y,
       [&](int ix, Tx valx)
       {
         jacpols2[2*ix+5] -> EvaluateScaledLambda (n-3-ix, 2*y-1, t, [&](int iy, Ty valy)
                                                   {
                                                     func(ii++, bub*valx*valy);
                                                   });
       });
  }

  
  // compute face bubbles up to order n, 0 < y, y-x < 1, x+y < 1
  template <typename T>
  static void CalcScaledTrigShapeDxDyDt (int n, T x, T y, T t, T * dshape)
  {
    /*
    if (n < 3) return;
    AutoDiff<3,T> adx(x, 0);
    AutoDiff<3,T> ady(y, 1);
    AutoDiff<3,T> adt(t, 2);
    AutoDiff<3,T> res[2000];
    CalcScaledTrigShape (n, adx, ady, adt, &res[0]);
    int ndof = (n-1)*(n-2)/2;
    for (int i = 0; i < ndof; i++)
      {
	dshape[3*i] = res[i].DValue(0);
	dshape[3*i+1] = res[i].DValue(1);
	dshape[3*i+2] = res[i].DValue(2);
      }
    */
    if (n < 3) return;
    AutoDiff<3,T> adx(x, 0);
    AutoDiff<3,T> ady(y, 1);
    AutoDiff<3,T> adt(t, 2);
    CalcScaledTrigShapeLambda (n, adx, ady, adt,
                               [&] (int i, AutoDiff<3,T> shape)
                               {
                                 dshape[3*i] = shape.DValue(0);
                                 dshape[3*i+1] = shape.DValue(1);
                                 dshape[3*i+2] = shape.DValue(2);
                               });
  }

      

  CurvedElements :: CurvedElements (const Mesh & amesh)
    : mesh(amesh)
  {
    order = 1;
    rational = 0;
    ishighorder = 0;
  }


  CurvedElements :: ~CurvedElements()
  {
  }


  void CurvedElements :: BuildCurvedElements(const Refinement * ref, int aorder,
                                             bool arational)
  {
    auto & geo = *mesh.GetGeometry();

    ishighorder = 0;
    order = 1;

    // MPI_Comm curve_comm;
    const auto & curve_comm = mesh.GetCommunicator();
#ifdef PARALLEL
    enum { MPI_TAG_CURVE = MPI_TAG_MESH+20 };

    const ParallelMeshTopology & partop = mesh.GetParallelTopology ();
    // MPI_Comm_dup (mesh.GetCommunicator(), &curve_comm);      
    NgArray<int> procs;
#else
    // curve_comm = mesh.GetCommunicator();
#endif
    int id = curve_comm.Rank();
    int ntasks = curve_comm.Size();

    bool working = (ntasks == 1) || (id > 0);

    if (working)
      order = aorder;

    if (mesh.coarsemesh)
      {
	mesh.coarsemesh->GetCurvedElements().BuildCurvedElements (ref, aorder, arational);
        order = aorder;
        rational = arational;
        ishighorder = (order > 1);
	return;
      }


    PrintMessage (1, "Curve elements, order = ", aorder);
    if (rational) PrintMessage (1, "curved elements with rational splines");

    // if (working)
    const_cast<Mesh&> (mesh).UpdateTopology();
    const MeshTopology & top = mesh.GetTopology();

    rational = arational;

    NgArray<int> edgenrs;
    int nedges = top.GetNEdges();
    int nfaces = top.GetNFaces();

    edgeorder.SetSize (nedges);
    faceorder.SetSize (nfaces);

    edgeorder = 1;
    faceorder = 1;

    if (rational)
      {
        edgeweight.SetSize (nedges);
        edgeweight = 1.0;
      }

    
    if (aorder <= 1) 
      {
	for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	  if (mesh[ei].GetType() == TET10)
	    ishighorder = 1;
	return; 
      }


    if (rational) aorder = 2;

    if (working)
      {
	if (mesh.GetDimension() == 3)
	  for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
	    {
	      top.GetEdges (i, edgenrs);
	      for (int j = 0; j < edgenrs.Size(); j++)
		edgeorder[edgenrs[j]] = aorder;
	      faceorder[top.GetFace (i)] = aorder;
	    }
	for (SegmentIndex i = 0; i < mesh.GetNSeg(); i++)
	  edgeorder[top.GetEdge (i)] = aorder;
      }

    if (rational)
      {
        edgeorder = 2;
        faceorder = 1;
      }


#ifdef PARALLEL
    TABLE<int> send_orders(ntasks), recv_orders(ntasks);

    if (ntasks > 1 && working)
      {
	for (int e = 0; e < edgeorder.Size(); e++)
	  {
	    partop.GetDistantEdgeNums (e+1, procs);
	    for (int j = 0; j < procs.Size(); j++)
	      send_orders.Add (procs[j], edgeorder[e]);
	  }
	for (int f = 0; f < faceorder.Size(); f++)
	  {
	    partop.GetDistantFaceNums (f+1, procs);
	    for (int j = 0; j < procs.Size(); j++)
	      send_orders.Add (procs[j], faceorder[f]);
	  }
      }

    if (ntasks > 1)
      MyMPI_ExchangeTable (send_orders, recv_orders, MPI_TAG_CURVE, curve_comm);

    if (ntasks > 1 && working)
      {
	NgArray<int> cnt(ntasks);
	cnt = 0;
	for (int e = 0; e < edgeorder.Size(); e++)
	  {
	    partop.GetDistantEdgeNums (e+1, procs);
	    for (int j = 0; j < procs.Size(); j++)
	      edgeorder[e] = max(edgeorder[e], recv_orders[procs[j]][cnt[procs[j]]++]);
	  }
	for (int f = 0; f < faceorder.Size(); f++)
	  {
	    partop.GetDistantFaceNums (f+1, procs);
	    for (int j = 0; j < procs.Size(); j++)
	      faceorder[f] = max(faceorder[f], recv_orders[procs[j]][cnt[procs[j]]++]);
	  }
      }
#endif


    edgecoeffsindex.SetSize (nedges+1);
    int nd = 0;
    for (int i = 0; i < nedges; i++)
      {
	edgecoeffsindex[i] = nd;
	nd += max (0, edgeorder[i]-1);
      }
    edgecoeffsindex[nedges] = nd;

    edgecoeffs.SetSize (nd);
    edgecoeffs = Vec<3> (0,0,0);
    

    facecoeffsindex.SetSize (nfaces+1);
    nd = 0;
    for (int i = 0; i < nfaces; i++)
      {
	facecoeffsindex[i] = nd;
	if (top.GetFaceType(i+1) == TRIG)
	  nd += max2 (0, (faceorder[i]-1)*(faceorder[i]-2)/2);
	else
	  nd += max2 (0, sqr(faceorder[i]-1));
      }
    facecoeffsindex[nfaces] = nd;

    facecoeffs.SetSize (nd);
    facecoeffs = Vec<3> (0,0,0);


    if (!ref || aorder <= 1) 
      {
        order = aorder;
	return; 
      }
    
    NgArray<double> xi, weight;

    ComputeGaussRule (aorder+4, xi, weight);  // on (0,1)

    buildJacPols();
    PrintMessage (3, "Curving edges");

    if (mesh.GetDimension() == 3 || rational)
      {
        static Timer tce("curve edges"); RegionTimer reg(tce);
	NgArray<int> surfnr(nedges);
	NgArray<PointGeomInfo> gi0(nedges);
	NgArray<PointGeomInfo> gi1(nedges);
	surfnr = -1;

	if (working)
	  for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
	    {
	      top.GetEdges (i, edgenrs);
	      const Element2d & el = mesh[i];
	      const ELEMENT_EDGE * edges = MeshTopology::GetEdges0 (el.GetType());

	      for (int i2 = 0; i2 < edgenrs.Size(); i2++)
		{
		  // PointIndex pi1 = el[edges[i2][0]];
		  // PointIndex pi2 = el[edges[i2][1]];

		  // bool swap = pi1 > pi2;
		
		  // Point<3> p1 = mesh[pi1];
		  // Point<3> p2 = mesh[pi2];
		
		  // int order1 = edgeorder[edgenrs[i2]];
		  // int ndof = max (0, order1-1);

		  surfnr[edgenrs[i2]] = mesh.GetFaceDescriptor(el.GetIndex()).SurfNr();
		  gi0[edgenrs[i2]] = el.GeomInfoPi(edges[i2][0]+1);
		  gi1[edgenrs[i2]] = el.GeomInfoPi(edges[i2][1]+1);
		}
	    }


#ifdef PARALLEL
	if (ntasks > 1)
	  {
	    // distribute it ...
	    TABLE<double> senddata(ntasks), recvdata(ntasks);
	    if (working)
	      for (int e = 0; e < nedges; e++)
		{
		  partop.GetDistantEdgeNums (e+1, procs);
		  for (int j = 0; j < procs.Size(); j++)
		    {
		      senddata.Add (procs[j], surfnr[e]);
		      if (surfnr[e] != -1)
			{
			  senddata.Add (procs[j], gi0[e].trignum);
			  senddata.Add (procs[j], gi0[e].u);
			  senddata.Add (procs[j], gi0[e].v);
			  senddata.Add (procs[j], gi1[e].trignum);
			  senddata.Add (procs[j], gi1[e].u);
			  senddata.Add (procs[j], gi1[e].v);
			}
		    }
		}
	    
	    MyMPI_ExchangeTable (senddata, recvdata, MPI_TAG_CURVE, curve_comm);
	    

	    NgArray<int> cnt(ntasks);
	    cnt = 0;
	    if (working)
	      for (int e = 0; e < nedges; e++)
		{
		  partop.GetDistantEdgeNums (e+1, procs);
		  for (int j = 0; j < procs.Size(); j++)
		    {
		      int surfnr1 = recvdata[procs[j]][cnt[procs[j]]++];
		      if (surfnr1 != -1)
			{
			  surfnr[e] = surfnr1; 
			  gi0[e].trignum = int (recvdata[procs[j]][cnt[procs[j]]++]);
			  gi0[e].u = recvdata[procs[j]][cnt[procs[j]]++];
			  gi0[e].v = recvdata[procs[j]][cnt[procs[j]]++];
			  gi1[e].trignum = int (recvdata[procs[j]][cnt[procs[j]]++]);
			  gi1[e].u = recvdata[procs[j]][cnt[procs[j]]++];
			  gi1[e].v = recvdata[procs[j]][cnt[procs[j]]++];
			}
		    }
		}
	    
	  }
#endif    


	if (working)
	  for (int e = 0; e < surfnr.Size(); e++)
	    {
	      if (surfnr[e] == -1) continue;
	      SetThreadPercent(double(e)/surfnr.Size()*100.);

	      PointIndex pi1, pi2;
	      top.GetEdgeVertices (e+1, pi1, pi2);
	      bool swap = (pi1 > pi2);

	      Point<3> p1 = mesh[pi1];
	      Point<3> p2 = mesh[pi2];

	      int order1 = edgeorder[e];
	      int ndof = max (0, order1-1);

	      if (rational && order1 >= 2)
		{
		  Point<3> pm = Center (p1, p2);

		  Vec<3> n1 = geo.GetNormal (surfnr[e], p1, &gi0[e]);
		  Vec<3> n2 = geo.GetNormal (surfnr[e], p2, &gi1[e]);

		  // p3 = pm + alpha1 n1 + alpha2 n2
		
		  Mat<2> mat, inv;
		  Vec<2> rhs, sol;
		
		  mat(0,0) = n1*n1;
		  mat(0,1) = mat(1,0) = n1*n2;
		  mat(1,1) = n2*n2;
                
		  rhs(0) = n1 * (p1-pm);
		  rhs(1) = n2 * (p2-pm);
                  

		  Point<3> p3;
		
		  if (fabs (Det (mat)) > 1e-10)
		    {
		      CalcInverse (mat, inv);
		      sol = inv * rhs;
		    
		      p3 = pm + sol(0) * n1 + sol(1) * n2;
		    }
		  else
		    p3 = pm;
		
		  edgecoeffs[edgecoeffsindex[e]] = Vec<3> (p3);
		

		  double wold = 1, w = 1, dw = 0.1;
		  double dold = 1e99;
		  while (fabs (dw) > 1e-12)
		    {
		      Vec<3> v05 = 0.25 * Vec<3> (p1) + 0.5*w* Vec<3>(p3) + 0.25 * Vec<3> (p2);
		      v05 /= 1 + (w-1) * 0.5;
		      Point<3> p05 (v05), pp05(v05);
		      geo.ProjectPointGI(surfnr[e], pp05, gi0[e]);
		      double d = Dist (pp05, p05);
                    
		      if (d < dold)
			{
			  dold = d;
			  wold = w;
			  w += dw;
			}
		      else
			{
			  dw *= -0.7;
			  w = wold + dw;
			}
		    }
		
		  edgeweight[e] = w;
		  continue;
		}
	    
	      Vector shape(ndof);
	      DenseMatrix mat(ndof, ndof), inv(ndof, ndof),
		rhs(ndof, 3), sol(ndof, 3);
	    
	      rhs = 0.0;
	      mat = 0.0;
	      for (int j = 0; j < xi.Size(); j++)
		{
		  Point<3> p;
		  Point<3> pp;
		  PointGeomInfo ppgi;
		
		  if (swap)
		    {
		      p = p1 + xi[j] * (p2-p1);
		      geo.PointBetween (p1, p2, xi[j],
                                        surfnr[e], gi0[e], gi1[e],
                                        pp, ppgi);
		    }
		  else
		    {
		      p = p2 + xi[j] * (p1-p2);
		      geo.PointBetween (p2, p1, xi[j],
                                        surfnr[e], gi1[e], gi0[e],
                                        pp, ppgi);
		    }
		
		  Vec<3> dist = pp - p;
		
		  CalcEdgeShape (order1, 2*xi[j]-1, &shape(0));
		
		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < ndof; l++)
		      mat(k,l) += weight[j] * shape(k) * shape(l);
		
		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < 3; l++)
		      rhs(k,l) += weight[j] * shape(k) * dist(l);
		}
	    
	      CalcInverse (mat, inv);
	      Mult (inv, rhs, sol);
	    
	      int first = edgecoeffsindex[e];
	      for (int j = 0; j < ndof; j++)
		for (int k = 0; k < 3; k++)
		  edgecoeffs[first+j](k) = sol(j,k);
	    }
      }


    NgArray<int> use_edge(nedges);
    NgArray<int> edge_surfnr1(nedges);
    NgArray<int> edge_surfnr2(nedges);
    NgArray<int> swap_edge(nedges);
    NgArray<EdgePointGeomInfo> edge_gi0(nedges);
    NgArray<EdgePointGeomInfo> edge_gi1(nedges);
    use_edge = 0;

    if (working)
      for (SegmentIndex i = 0; i < mesh.GetNSeg(); i++)
	{
	  const Segment & seg = mesh[i];
	  int edgenr = top.GetEdge (i);
	  use_edge[edgenr] = 1;
	  edge_surfnr1[edgenr] = seg.surfnr1;
	  edge_surfnr2[edgenr] = seg.surfnr2;
	  edge_gi0[edgenr] = seg.epgeominfo[0];
	  edge_gi1[edgenr] = seg.epgeominfo[1];
	  swap_edge[edgenr] = int (seg[0] > seg[1]);
	}

#ifdef PARALLEL
    if (ntasks > 1)
      {
	// distribute it ...
	TABLE<double> senddata(ntasks), recvdata(ntasks);
	if (working)
	  for (int e = 0; e < nedges; e++)
	    {
	      partop.GetDistantEdgeNums (e+1, procs);
	      for (int j = 0; j < procs.Size(); j++)
		{
		  senddata.Add (procs[j], use_edge[e]);
		  if (use_edge[e])
		    {
		      senddata.Add (procs[j], edge_surfnr1[e]);
		      senddata.Add (procs[j], edge_surfnr2[e]);
		      senddata.Add (procs[j], edge_gi0[e].edgenr);
		      senddata.Add (procs[j], edge_gi0[e].body);
		      senddata.Add (procs[j], edge_gi0[e].dist);
		      senddata.Add (procs[j], edge_gi0[e].u);
		      senddata.Add (procs[j], edge_gi0[e].v);
		      senddata.Add (procs[j], edge_gi1[e].edgenr);
		      senddata.Add (procs[j], edge_gi1[e].body);
		      senddata.Add (procs[j], edge_gi1[e].dist);
		      senddata.Add (procs[j], edge_gi1[e].u);
		      senddata.Add (procs[j], edge_gi1[e].v);
		      senddata.Add (procs[j], swap_edge[e]);
		    }
		}
	    }
	MyMPI_ExchangeTable (senddata, recvdata, MPI_TAG_CURVE, curve_comm);
	NgArray<int> cnt(ntasks);
	cnt = 0;
	if (working)
	  for (int e = 0; e < edge_surfnr1.Size(); e++)
	    {
	      partop.GetDistantEdgeNums (e+1, procs);
	      for (int j = 0; j < procs.Size(); j++)
		{
		  int get_edge = int(recvdata[procs[j]][cnt[procs[j]]++]);
		  if (get_edge)
		    {
		      use_edge[e] = 1;
		      edge_surfnr1[e] = int (recvdata[procs[j]][cnt[procs[j]]++]);
		      edge_surfnr2[e] = int (recvdata[procs[j]][cnt[procs[j]]++]);
		      edge_gi0[e].edgenr = int (recvdata[procs[j]][cnt[procs[j]]++]);
		      edge_gi0[e].body = int (recvdata[procs[j]][cnt[procs[j]]++]);
		      edge_gi0[e].dist = recvdata[procs[j]][cnt[procs[j]]++];
		      edge_gi0[e].u = recvdata[procs[j]][cnt[procs[j]]++];
		      edge_gi0[e].v = recvdata[procs[j]][cnt[procs[j]]++];
		      edge_gi1[e].edgenr = int (recvdata[procs[j]][cnt[procs[j]]++]);
		      edge_gi1[e].body = int (recvdata[procs[j]][cnt[procs[j]]++]);
		      edge_gi1[e].dist = recvdata[procs[j]][cnt[procs[j]]++];
		      edge_gi1[e].u = recvdata[procs[j]][cnt[procs[j]]++];
		      edge_gi1[e].v = recvdata[procs[j]][cnt[procs[j]]++];
		      swap_edge[e] = recvdata[procs[j]][cnt[procs[j]]++];
		    }
		}
	    }

      }
#endif    

    if (working)
      for (int edgenr = 0; edgenr < use_edge.Size(); edgenr++)
	{
	  int segnr = edgenr;
	  if (!use_edge[edgenr]) continue;

	  SetThreadPercent(double(edgenr)/edge_surfnr1.Size()*100.);

	  PointIndex pi1, pi2;
	  top.GetEdgeVertices (edgenr+1, pi1, pi2);

	  bool swap = swap_edge[edgenr]; // (pi1 > pi2);
	  if (swap) Swap (pi1, pi2);

	  Point<3> p1 = mesh[pi1];
	  Point<3> p2 = mesh[pi2];

	  int order1 = edgeorder[segnr];
	  int ndof = max (0, order1-1);

	  if (rational)
	    {
	      Vec<3> tau1 = geo.GetTangent(p1, edge_surfnr2[edgenr], edge_surfnr1[edgenr],
                                           edge_gi0[edgenr]);
	      Vec<3> tau2 = geo.GetTangent(p2, edge_surfnr2[edgenr], edge_surfnr1[edgenr],
                                           edge_gi1[edgenr]);
	      // p1 + alpha1 tau1 = p2 + alpha2 tau2;

	      Mat<3,2> mat;
	      Mat<2,3> inv;
	      Vec<3> rhs;
	      Vec<2> sol;
	      for (int j = 0; j < 3; j++)
		{
		  mat(j,0) = tau1(j); 
		  mat(j,1) = -tau2(j); 
		  rhs(j) = p2(j)-p1(j); 
		}
	      CalcInverse (mat, inv);
	      sol = inv * rhs;

	      Point<3> p3 = p1+sol(0) * tau1;
	      edgecoeffs[edgecoeffsindex[segnr]] = Vec<3> (p3);

	      double wold = 1, w = 1, dw = 0.1;
	      double dold = 1e99;
	      while (fabs (dw) > 1e-12)
		{
		  Vec<3> v05 = 0.25 * Vec<3> (p1) + 0.5*w* Vec<3>(p3) + 0.25 * Vec<3> (p2);
		  v05 /= 1 + (w-1) * 0.5;
		  Point<3> p05 (v05), pp05(v05);
		  geo.ProjectPointEdge(edge_surfnr1[edgenr], edge_surfnr2[edgenr], pp05,
                                       &edge_gi0[edgenr]);
		  double d = Dist (pp05, p05);

		  if (d < dold)
		    {
		      dold = d;
		      wold = w;
		      w += dw;
		    }
		  else
		    {
		      dw *= -0.7;
		      w = wold + dw;
		    }
		  // *testout << "w = " << w << ", dw = " << dw << endl;
		}

	      // cout << "wopt = " << w << ", dopt = " << dold << endl;
	      edgeweight[segnr] = w;
            
	      //             cout << "p1 = " << p1 << ", tau1 = " << tau1 << ", alpha1 = " << sol(0) << endl;
	      //             cout << "p2 = " << p2 << ", tau2 = " << tau2 << ", alpha2 = " << -sol(1) << endl;
	      //             cout << "p+alpha tau = " << p1 + sol(0) * tau1 
	      //                  << " =?= " << p2 +sol(1) * tau2 << endl;
            
	    }

	  else
          
	    {
	      Vector shape(ndof);
	      DenseMatrix mat(ndof, ndof), inv(ndof, ndof),
		rhs(ndof, 3), sol(ndof, 3);

	      rhs = 0.0;
	      mat = 0.0;
	      for (int j = 0; j < xi.Size(); j++)
		{
		  Point<3> p, pp;
		  EdgePointGeomInfo ppgi;
	    
		  if (swap)
		    {
		      p = p1 + xi[j] * (p2-p1);
		      geo.PointBetweenEdge(p1, p2, xi[j],
                                           edge_surfnr2[edgenr], edge_surfnr1[edgenr],
                                           edge_gi0[edgenr], edge_gi1[edgenr],
                                           pp, ppgi);
		    }
		  else
		    {
		      p = p2 + xi[j] * (p1-p2);
		      geo.PointBetweenEdge(p2, p1, xi[j],
					   edge_surfnr2[edgenr], edge_surfnr1[edgenr],
					   edge_gi1[edgenr], edge_gi0[edgenr],
					   pp, ppgi);
		    }
	    
		  Vec<3> dist = pp - p;

		  CalcEdgeShape (order1, 2*xi[j]-1, &shape(0));

		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < ndof; l++)
		      mat(k,l) += weight[j] * shape(k) * shape(l);

		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < 3; l++)
		      rhs(k,l) += weight[j] * shape(k) * dist(l);
		}


	      CalcInverse (mat, inv);
	      Mult (inv, rhs, sol);

	      int first = edgecoeffsindex[segnr];
	      for (int j = 0; j < ndof; j++)
		for (int k = 0; k < 3; k++)
		  edgecoeffs[first+j](k) = sol(j,k);
	    }
	}

   
    
    PrintMessage (3, "Curving faces");

    NgArray<int> surfnr(nfaces);
    surfnr = -1;

    if (working)
      for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
	surfnr[top.GetFace(i)] = 
	  mesh.GetFaceDescriptor(mesh[i].GetIndex()).SurfNr();

#ifdef PARALLEL
    TABLE<int> send_surfnr(ntasks), recv_surfnr(ntasks);

    if (ntasks > 1 && working)
      {
	for (int f = 0; f < nfaces; f++)
	  {
	    partop.GetDistantFaceNums (f+1, procs);
	    for (int j = 0; j < procs.Size(); j++)
	      send_surfnr.Add (procs[j], surfnr[f]);
	  }
      }

    if (ntasks > 1)
      MyMPI_ExchangeTable (send_surfnr, recv_surfnr, MPI_TAG_CURVE, curve_comm);

    if (ntasks > 1 && working)
      {
	NgArray<int> cnt(ntasks);
	cnt = 0;
	for (int f = 0; f < nfaces; f++)
	  {
	    partop.GetDistantFaceNums (f+1, procs);
	    for (int j = 0; j < procs.Size(); j++)
	      surfnr[f] = max(surfnr[f], recv_surfnr[procs[j]][cnt[procs[j]]++]);
	  }
      }
#endif

    if (mesh.GetDimension() == 3 && working)
      {
        static Timer tcf("curve faces"); RegionTimer reg(tcf);
	for (int f = 0; f < nfaces; f++)
	  {
	    int facenr = f;
	    if (surfnr[f] == -1) continue;
	    // if (el.GetType() == TRIG && order >= 3)
	    if (top.GetFaceType(facenr+1) == TRIG && order >= 3)
	      {
		NgArrayMem<int, 3> verts(3);
		top.GetFaceVertices (facenr+1, verts);

		int fnums[] = { 0, 1, 2 };
		/*
		if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
		if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		*/
		if (verts[fnums[0]] > verts[fnums[1]]) swap (fnums[0], fnums[1]);
		if (verts[fnums[1]] > verts[fnums[2]]) swap (fnums[1], fnums[2]);
		if (verts[fnums[0]] > verts[fnums[1]]) swap (fnums[0], fnums[1]);

		int order1 = faceorder[facenr];
		int ndof = max (0, (order1-1)*(order1-2)/2);
	    
		Vector shape(ndof), dmat(ndof);
		MatrixFixWidth<3> rhs(ndof), sol(ndof);
	    
		rhs = 0.0;
		dmat = 0.0;

		int np = sqr(xi.Size());
		NgArray<Point<2> > xia(np);
		NgArray<Point<3> > xa(np);

		for (int jx = 0, jj = 0; jx < xi.Size(); jx++)
		  for (int jy = 0; jy < xi.Size(); jy++, jj++)
		    xia[jj] = Point<2> ((1-xi[jy])*xi[jx], xi[jy]);

		// CalcMultiPointSurfaceTransformation (&xia, i, &xa, NULL);

		NgArray<int> edgenrs;
		top.GetFaceEdges (facenr+1, edgenrs);
		for (int k = 0; k < edgenrs.Size(); k++) edgenrs[k]--;

		for (int jj = 0; jj < np; jj++)
		  {
		    Point<3> pp(0,0,0);
		    double lami[] = { xia[jj](0), xia[jj](1), 1-xia[jj](0)-xia[jj](1)};

		    for (int k = 0; k < verts.Size(); k++)
		      pp += lami[k] * Vec<3> (mesh.Point(verts[k]));

		    // const ELEMENT_EDGE * edges = MeshTopology::GetEdges0 (TRIG);
		    for (int k = 0; k < edgenrs.Size(); k++)
		      {
			int eorder = edgeorder[edgenrs[k]];
			if (eorder < 2) continue;

			int first = edgecoeffsindex[edgenrs[k]];
			Vector eshape(eorder-1);
			int vi1, vi2;
			top.GetEdgeVertices (edgenrs[k]+1, vi1, vi2);
			if (vi1 > vi2) swap (vi1, vi2);
			int v1 = -1, v2 = -1;
			for (int j = 0; j < 3; j++)
			  {
			    if (verts[j] == vi1) v1 = j;
			    if (verts[j] == vi2) v2 = j;
			  }

			CalcScaledEdgeShape (eorder, lami[v1]-lami[v2], lami[v1]+lami[v2], &eshape(0));
			for (int n = 0; n < eshape.Size(); n++)
			  pp += eshape(n) * edgecoeffs[first+n];
		      }
		    xa[jj] = pp;
		  }

		for (int jx = 0, jj = 0; jx < xi.Size(); jx++)
		  for (int jy = 0; jy < xi.Size(); jy++, jj++)
		    {
		      double y = xi[jy];
		      double x = (1-y) * xi[jx];
		      double lami[] = { x, y, 1-x-y };
		      double wi = weight[jx]*weight[jy]*(1-y);
 
		      Point<3> pp = xa[jj];
		      // ref -> ProjectToSurface (pp, mesh.GetFaceDescriptor(el.GetIndex()).SurfNr());
		      /**
			 with MPI and an interior surface element between volume elements assigned to different
			 procs, only one of them has the surf-el
		      **/
                      SurfaceElementIndex sei = top.GetFace2SurfaceElement (f+1)-1;
		      if (sei != SurfaceElementIndex(-1)) {
			PointGeomInfo gi = mesh[sei].GeomInfoPi(1);
			geo.ProjectPointGI(surfnr[facenr], pp, gi);
		      }
		      else
			{ geo.ProjectPoint(surfnr[facenr], pp); }
		      Vec<3> dist = pp-xa[jj];
		
		      CalcTrigShape (order1, lami[fnums[1]]-lami[fnums[0]],
				     1-lami[fnums[1]]-lami[fnums[0]], &shape(0));

		      for (int k = 0; k < ndof; k++)
			dmat(k) += wi * shape(k) * shape(k);

		      dist *= wi;
		      for (int k = 0; k < ndof; k++)
			for (int l = 0; l < 3; l++)
			  rhs(k,l) += shape(k) * dist(l);
		    }

		for (int i = 0; i < ndof; i++)
		  for (int j = 0; j < 3; j++)
		    sol(i,j) = rhs(i,j) / dmat(i);   // Orthogonal basis !

		int first = facecoeffsindex[facenr];
		for (int j = 0; j < ndof; j++)
		  for (int k = 0; k < 3; k++)
		    facecoeffs[first+j](k) = sol(j,k);
	      }
	  }
      }


    // compress edge and face tables
    int newbase = 0;
    for (int i = 0; i < edgeorder.Size(); i++)
      {
	bool curved = 0;
	int oldbase = edgecoeffsindex[i];
	int nd = edgecoeffsindex[i+1] - edgecoeffsindex[i];

	for (int j = 0; j < nd; j++)
	  if (edgecoeffs[oldbase+j].Length() > 1e-12)
	    curved = 1;
	if (rational) curved = 1;

	if (curved && newbase != oldbase)
	  for (int j = 0; j < nd; j++)
	    edgecoeffs[newbase+j] = edgecoeffs[oldbase+j];

	edgecoeffsindex[i] = newbase;
	if (!curved) edgeorder[i] = 1;
	if (curved) newbase += nd;
      }
    edgecoeffsindex.Last() = newbase;


    newbase = 0;
    for (int i = 0; i < faceorder.Size(); i++)
      {
	bool curved = 0;
	int oldbase = facecoeffsindex[i];
	int nd = facecoeffsindex[i+1] - facecoeffsindex[i];

	for (int j = 0; j < nd; j++)
	  if (facecoeffs[oldbase+j].Length() > 1e-12)
	    curved = 1;

	if (curved && newbase != oldbase)
	  for (int j = 0; j < nd; j++)
	    facecoeffs[newbase+j] = facecoeffs[oldbase+j];

	facecoeffsindex[i] = newbase;
	if (!curved) faceorder[i] = 1;
	if (curved) newbase += nd;
      }
    facecoeffsindex.Last() = newbase;
    
    if (working)
      ishighorder = (order > 1);
    // (*testout) << "edgecoeffs = " << endl << edgecoeffs << endl;
    // (*testout) << "facecoeffs = " << endl << facecoeffs << endl;


#ifdef PARALLEL
    curve_comm.Barrier();
    // MPI_Comm_free (&curve_comm);      
#endif
  }










  // ***********************  Transform edges *****************************

  
  bool CurvedElements ::  IsSegmentCurved (SegmentIndex elnr) const
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	return mesh.coarsemesh->GetCurvedElements().IsSegmentCurved (hpref_el.coarse_elnr);
      }

    SegmentInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = 2;
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	info.edgenr = top.GetSegmentEdge (elnr+1)-1;	
	info.ndof += edgeorder[info.edgenr]-1;
      }

    return (info.ndof > info.nv);
  }


 
  
  template <typename T>
  void CurvedElements :: 
  CalcSegmentTransformation (T xi, SegmentIndex elnr,
			     Point<3,T> * x, Vec<3,T> * dxdxi, bool * curved)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	// xi umrechnen
	T lami[2] = { xi, 1-xi };
	T dlami[2] = { 1, -1 };

	T coarse_xi = 0;
	T trans = 0;
	for (int i = 0; i < 2; i++)
	  {
	    coarse_xi += hpref_el.param[i][0] * lami[i];
	    trans += hpref_el.param[i][0] * dlami[i];
	  }

	mesh.coarsemesh->GetCurvedElements().CalcSegmentTransformation (coarse_xi, hpref_el.coarse_elnr, x, dxdxi, curved);
	if (dxdxi) *dxdxi *= trans;
	
	return;
      }
    


    // TVector<T> shapes, dshapes;
    //     NgArray<Vec<3> > coefs;

    SegmentInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = 2;

    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	info.edgenr = top.GetSegmentEdge (elnr+1)-1;	
	info.ndof += edgeorder[info.edgenr]-1;
      }

    NgArrayMem<Vec<3>,100> coefs(info.ndof);
    NgArrayMem<T, 100> shapes_mem(info.ndof);
    TFlatVector<T> shapes(info.ndof, &shapes_mem[0]);
    NgArrayMem<T, 200> dshapes_mem(info.ndof);
    TFlatVector<T> dshapes(info.ndof, &dshapes_mem[0]);

    
    CalcElementShapes (info, xi, shapes);
    GetCoefficients (info, coefs);

    *x = 0;
    for (int i = 0; i < shapes.Size(); i++)
      // *x += shapes(i) * coefs[i];
      for (int j = 0; j < 3; j++)
        (*x)(j) += shapes(i) * coefs[i](j);


    if (dxdxi)
      {
	CalcElementDShapes (info, xi, dshapes);
	
	*dxdxi = 0;
	for (int i = 0; i < shapes.Size(); i++)
	  for (int j = 0; j < 3; j++)
	    (*dxdxi)(j) += dshapes(i) * coefs[i](j);
      }

    if (curved)
      *curved = (info.order > 1);

    // cout << "Segment, |x| = " << Abs2(Vec<3> (*x) ) << endl;
  }


  template <typename T>
  void CurvedElements :: 
  CalcElementShapes (SegmentInfo & info, T xi, TFlatVector<T> shapes) const
  {
    /*
    if (rational && info.order == 2)
      {
	shapes.SetSize(3);
	double w = edgeweight[info.edgenr];
	shapes(0) = xi*xi;
	shapes(1) = (1-xi)*(1-xi);
	shapes(2) = 2*w*xi*(1-xi);
	shapes *= 1.0 / (1 + (w-1) *2*xi*(1-xi));
	return;
      }
    */

    // shapes.SetSize(info.ndof);
    shapes(0) = xi;
    shapes(1) = 1-xi;

    if (info.order >= 2)
      {
	if (mesh[info.elnr][0] > mesh[info.elnr][1])
	  xi = 1-xi;
	CalcEdgeShape (edgeorder[info.edgenr], 2*xi-1, &shapes(2));
      }
  }

  template <typename T>
  void CurvedElements :: 
  CalcElementDShapes (SegmentInfo & info, T xi, TFlatVector<T> dshapes) const
  {
    /*
    if (rational && info.order == 2)
      {
	dshapes.SetSize(3);
	double wi = edgeweight[info.edgenr];
	double shapes[3];
	shapes[0] = xi*xi;
	shapes[1] = (1-xi)*(1-xi);
	shapes[2] = 2*wi*xi*(1-xi);
	double w = 1 + (wi-1) *2*xi*(1-xi);
	double dw = (wi-1) * (2 - 4*xi);
        
	dshapes(0) = 2*xi;
	dshapes(1) = 2*(xi-1);
	dshapes(2) = 2*wi*(1-2*xi);

	for (int j = 0;j < 3; j++)
	  dshapes(j) = dshapes(j) / w - shapes[j] * dw / (w*w);
	return;
      }
    */



    // dshapes.SetSize(info.ndof);
    dshapes = 0;
    dshapes(0) = 1;
    dshapes(1) = -1;

    // int order = edgeorder[info.edgenr];

    if (info.order >= 2)
      {
	T fac = 2;
	if (mesh[info.elnr][0] > mesh[info.elnr][1])
	  {
	    xi = 1-xi; 
	    fac *= -1;
	  }
	CalcEdgeDx (edgeorder[info.edgenr], 2*xi-1, &dshapes(2));
	for (int i = 2; i < dshapes.Size(); i++)
	  dshapes(i) *= fac;
      }

    // ??? not implemented ????
  }

  void CurvedElements :: 
  GetCoefficients (SegmentInfo & info, NgArray<Vec<3> > & coefs) const
  {
    const Segment & el = mesh[info.elnr];

    coefs.SetSize(info.ndof);

    coefs[0] = Vec<3> (mesh[el[0]]);
    coefs[1] = Vec<3> (mesh[el[1]]);

    if (info.order >= 2)
      {
	int first = edgecoeffsindex[info.edgenr]; 
	int next = edgecoeffsindex[info.edgenr+1]; 
	for (int i = 0; i < next-first; i++)
	  coefs[i+2] = edgecoeffs[first+i];
      }
  }













  // ********************** Transform surface elements *******************


  bool CurvedElements :: IsSurfaceElementCurved (SurfaceElementIndex elnr) const
  {
    if (mesh[elnr].GetType() != TRIG) return true;
    if (!IsHighOrder()) return false;

    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	return mesh.coarsemesh->GetCurvedElements().IsSurfaceElementCurved (hpref_el.coarse_elnr);
      }

    const Element2d & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();
    
    SurfaceElementInfo info;
    info.elnr = elnr;
    info.order = order;

    switch (type)
      {
      case TRIG : info.nv = 3; break;
      case QUAD : info.nv = 4; break;
      case TRIG6: return true;
      default:
	cerr << "undef element in CalcSurfaceTrafo" << endl;
      }
    info.ndof = info.nv;

    // info.ndof = info.nv = ( (type == TRIG) || (type == TRIG6) ) ? 3 : 4;
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.edgenrs[i]--;
	info.facenr = top.GetSurfaceElementFace (elnr+1)-1;

	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
      }

    return (info.ndof > info.nv);
  }
  
  void CurvedElements :: 
  CalcSurfaceTransformation (Point<2> xi, SurfaceElementIndex elnr,
			     Point<3> * x, Mat<3,2> * dxdxi, bool * curved)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	// xi umrechnen
	double lami[4];
	FlatVector vlami(4, lami);
	vlami = 0;
	mesh[elnr].GetShapeNew (xi, vlami);
	
	Mat<2,2> trans;
	Mat<3,2> dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<2> dlami(4);
	    dlami = 0;
	    mesh[elnr].GetDShapeNew (xi, dlami);	  
	    
	    trans = 0;
	    for (int k = 0; k < 2; k++)
	      for (int l = 0; l < 2; l++)
		for (int i = 0; i < hpref_el.np; i++)
		  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
	  }
	
	Point<2> coarse_xi(0,0);
	for (int i = 0; i < hpref_el.np; i++)
	  for (int j = 0; j < 2; j++)
	    coarse_xi(j) += hpref_el.param[i][j] * lami[i];
	
	mesh.coarsemesh->GetCurvedElements().CalcSurfaceTransformation (coarse_xi, hpref_el.coarse_elnr, x, &dxdxic, curved);
	
	if (dxdxi)
	  *dxdxi = dxdxic * trans;
	
	return;
      }
    



    const Element2d & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    SurfaceElementInfo info;
    info.elnr = elnr;
    info.order = order;

    switch (type)
      {
      case TRIG : info.nv = 3; break;
      case QUAD : info.nv = 4; break;
      case TRIG6: info.nv = 6; break;
      case QUAD8 : info.nv = 8; break;
      default:
	cerr << "undef element in CalcSurfaceTrafo" << endl;
      }
    info.ndof = info.nv;

    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.edgenrs[i]--;
	info.facenr = top.GetSurfaceElementFace (elnr+1)-1;


	bool firsttry = true;
	bool problem = false;

	while(firsttry || problem)
	  {
	    problem = false;

	    for (int i = 0; !problem && i < info.edgenrs.Size(); i++)
	      {
		if(info.edgenrs[i]+1 >= edgecoeffsindex.Size())
		  problem = true;
		else
		  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	      }
	    if(info.facenr+1 >= facecoeffsindex.Size())
	      problem = true;
	    else
	      info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];

	    if(problem && !firsttry)
	      throw NgException("something wrong with curved elements");
	    
	    if(problem)
	      BuildCurvedElements(NULL,order,rational);

	    firsttry = false;
	  }
      }

    
    Point<2> _xi(xi);
    Point<3> _x;
    Mat<3,2> _dxdxi;
    if (EvaluateMapping (info, _xi, _x, _dxdxi))
      {
        if (x) *x = _x;
        if (dxdxi) *dxdxi = _dxdxi;
        return;
      }

    
    NgArrayMem<Vec<3>,100> coefs(info.ndof);
    NgArrayMem<double, 100> shapes_mem(info.ndof);
    TFlatVector<double> shapes(info.ndof, &shapes_mem[0]);
    NgArrayMem<double, 200> dshapes_mem(2*info.ndof);
    MatrixFixWidth<2> dshapes(info.ndof, &dshapes_mem[0]);


    CalcElementShapes (info, xi, shapes);
    GetCoefficients (info, coefs);

    *x = 0;
    for (int i = 0; i < coefs.Size(); i++)
      *x += shapes(i) * coefs[i];

    if (dxdxi)
      {
	CalcElementDShapes (info, xi, dshapes);
	
	*dxdxi = 0;
	for (int i = 0; i < coefs.Size(); i++)
	  for (int j = 0; j < 3; j++)
	    for (int k = 0; k < 2; k++)
	      (*dxdxi)(j,k) += dshapes(i,k) * coefs[i](j);
      }

    if (curved)
      *curved = (info.ndof > info.nv);
  }



  template <typename T>
  void CurvedElements :: 
  CalcElementShapes (SurfaceElementInfo & info, const Point<2,T> xi, TFlatVector<T> shapes) const
  {
    const Element2d & el = mesh[info.elnr];
    // shapes.SetSize(info.ndof);
    
    if (rational && info.order >= 2)
      {
	// shapes.SetSize(6);
	T w(1);
	T lami[3] = { xi(0), xi(1), 1-xi(0)-xi(1) };
	for (int j = 0; j < 3; j++)
	  shapes(j) = lami[j] * lami[j];

	const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TRIG);
	for (int j = 0; j < 3; j++)
	  {
	    T wi = edgeweight[info.edgenrs[j]];
	    shapes(j+3) = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	    w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	  }

	shapes *= 1.0 / w;
	return;
      }

    switch (el.GetType())
      {
      case TRIG:
	{
	  shapes(0) = xi(0);
	  shapes(1) = xi(1);
	  shapes(2) = 1-xi(0)-xi(1);

	  if (info.order == 1) return;

	  int ii = 3;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges0 (TRIG);
	  
	  for (int i = 0; i < 3; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0], vi2 = edges[i][1];
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShape (eorder, shapes(vi1)-shapes(vi2), shapes(vi1)+shapes(vi2), &shapes(ii));
		  ii += eorder-1;
		}
	    }

	  int forder = faceorder[info.facenr];
	  if (forder >= 3)
	    {
	      int fnums[] = { 0, 1, 2 };
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      
	      CalcTrigShape (forder, 
			     shapes(fnums[1])-shapes(fnums[0]),
			     1-shapes(fnums[1])-shapes(fnums[0]), &shapes(ii));
	    }
	  break;
	}

      case TRIG6:
	{
	  if (shapes.Size() == 3)
	    {
	      shapes(0) = xi(0);
	      shapes(1) = xi(1);
	      shapes(2) = 1-xi(0)-xi(1);
	    }
	  else
	    {
	      T x = xi(0);
	      T y = xi(1);
	      T lam3 = 1-x-y;
	      
	      shapes(0) = x * (2*x-1);
	      shapes(1) = y * (2*y-1);
	      shapes(2) = lam3 * (2*lam3-1);
	      shapes(3) = 4 * y * lam3;
	      shapes(4) = 4 * x * lam3;
	      shapes(5) = 4 * x * y;
	    }
	  break;
	}

      case QUAD:
	{
	  shapes(0) = (1-xi(0))*(1-xi(1));
	  shapes(1) =    xi(0) *(1-xi(1));
	  shapes(2) =    xi(0) *   xi(1) ;
	  shapes(3) = (1-xi(0))*   xi(1) ;

	  if (info.order == 1) return;
	  
	  T mu[4] = { 
	    1 - xi(0) + 1 - xi(1), 
	    xi(0) + 1 - xi(1), 
	    xi(0) +     xi(1), 
	    1 - xi(0) +     xi(1), 
	  };
	    
	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (QUAD);
	  
	  for (int i = 0; i < 4; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcEdgeShape (eorder, mu[vi1]-mu[vi2], &shapes(ii));
		  T lame = shapes(vi1)+shapes(vi2);
		  for (int j = 0; j < order-1; j++)
		    shapes(ii+j) *= lame;
		  ii += eorder-1;
		}
	    }
	  
	  for (int i = ii; i < info.ndof; i++)
	    shapes(i) = 0;

	  break;
	}

      case QUAD8:
	{
          auto x = xi(0), y = xi(1);
	  shapes(0) = (1-x)*(1-y);
	  shapes(1) = x*(1-y);
	  shapes(2) = x*y;
	  shapes(3) = (1-x)*y;
          shapes(4) = 4*(1-x)*x*(1-y);
          shapes(5) = 4*(1-x)*x*y;
          shapes(6) = 4*(1-y)*y*(1-x);
          shapes(7) = 4*(1-y)*y*x;
          shapes(0) -= 0.5*(shapes(4)+shapes(6));
          shapes(1) -= 0.5*(shapes(4)+shapes(7));
          shapes(2) -= 0.5*(shapes(5)+shapes(7));
          shapes(3) -= 0.5*(shapes(5)+shapes(6));
          break;
        }
        
      default:
	throw NgException("CurvedElements::CalcShape 2d, element type not handled");
      };
  }

  template <typename T>
  void CurvedElements :: 
  CalcElementDShapes (SurfaceElementInfo & info, const Point<2,T> xi, MatrixFixWidth<2,T> dshapes) const
  {
    const Element2d & el = mesh[info.elnr];
    ELEMENT_TYPE type = el.GetType();

    T lami[4];

    dshapes.SetSize(info.ndof);
    // dshapes = 0;	  

    // *testout << "calcelementdshapes, info.ndof = " << info.ndof << endl;

    if (rational && info.order >= 2)
      {
	T w = 1;
	T dw[2] = { 0, 0 };


	lami[0] = xi(0); lami[1] = xi(1); lami[2] = 1-xi(0)-xi(1);
	T dlami[3][2] = { { 1, 0 }, { 0, 1 }, { -1, -1 }};
	T shapes[6];

	for (int j = 0; j < 3; j++)
	  {
	    shapes[j] = lami[j] * lami[j];
	    dshapes(j,0) = 2 * lami[j] * dlami[j][0];
	    dshapes(j,1) = 2 * lami[j] * dlami[j][1];
	  }

	const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TRIG);
	for (int j = 0; j < 3; j++)
	  {
	    T wi = edgeweight[info.edgenrs[j]];

	    shapes[j+3] = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	    for (int k = 0; k < 2; k++)
	      dshapes(j+3,k) = 2*wi* (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
				      lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);

	    w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	    for (int k = 0; k < 2; k++)
	      dw[k] += 2*(wi-1) * (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
				   lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);
	  }
	// shapes *= 1.0 / w;
	dshapes *= 1.0 / w;
	for (int i = 0; i < 6; i++)
	  for (int j = 0; j < 2; j++)
	    dshapes(i,j) -= shapes[i] * dw[j] / (w*w);
	return;
      }





    switch (type)
      {
      case TRIG:
	{
	  dshapes(0,0) = 1;
	  dshapes(0,1) = 0.0;
	  dshapes(1,0) = 0.0;
	  dshapes(1,1) = 1;
	  dshapes(2,0) = -1;
	  dshapes(2,1) = -1;
	  
	  if (info.order == 1) return;

	  // *testout << "info.order = " << info.order << endl;


	  lami[0] = xi(0);
	  lami[1] = xi(1);
	  lami[2] = 1-xi(0)-xi(1);

	  int ii = 3;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TRIG);
	  
	  for (int i = 0; i < 3; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShapeDxDt<2> (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &dshapes(ii,0));

		  Mat<2,2,T> trans;
		  for (int j = 0; j < 2; j++)
		    {
		      trans(0,j) = dshapes(vi1,j)-dshapes(vi2,j);
		      trans(1,j) = dshapes(vi1,j)+dshapes(vi2,j);
		    }
		  
		  for (int j = 0; j < eorder-1; j++)
		    {
		      T ddx = dshapes(ii+j,0);
		      T ddt = dshapes(ii+j,1);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		    }

		  ii += eorder-1;
		}
	    }

	  int forder = faceorder[info.facenr];
	  // *testout << "forder = " << forder << endl;
	  if (forder >= 3)
	    {
	      int fnums[] = { 0, 1, 2 };
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
	      if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
	      
	      CalcTrigShapeDxDy (forder, 
				 lami[fnums[1]]-lami[fnums[0]],
				 1-lami[fnums[1]]-lami[fnums[0]], &dshapes(ii,0));

	      int nd = (forder-1)*(forder-2)/2;
	      Mat<2,2,T> trans;
	      for (int j = 0; j < 2; j++)
		{
		  trans(0,j) = dshapes(fnums[1],j)-dshapes(fnums[0],j);
		  trans(1,j) = -dshapes(fnums[1],j)-dshapes(fnums[0],j);
		}

	      for (int j = 0; j < nd; j++)
		{
		  T ddx = dshapes(ii+j,0);
		  T ddt = dshapes(ii+j,1);
		  dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		  dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		}
	    }

	  break;
	}

      case TRIG6:
	{
	  if (dshapes.Height() == 3)
	    {
	      dshapes = T(0.0);
	      dshapes(0,0) = 1;
	      dshapes(1,1) = 1;
	      dshapes(2,0) = -1;
	      dshapes(2,1) = -1;	    
	    }
	  else
	    {
	      AutoDiff<2,T> x(xi(0), 0);
	      AutoDiff<2,T> y(xi(1), 1);
	      AutoDiff<2,T> lam3 = 1-x-y;
	      AutoDiff<2,T> shapes[6];
	      shapes[0] = x * (2*x-1);
	      shapes[1] = y * (2*y-1);
	      shapes[2] = lam3 * (2*lam3-1);
	      shapes[3] = 4 * y * lam3;
	      shapes[4] = 4 * x * lam3;
	      shapes[5] = 4 * x * y;

	      for (int i = 0; i < 6; i++)
		{
		  dshapes(i,0) = shapes[i].DValue(0);
		  dshapes(i,1) = shapes[i].DValue(1);
		}
	      
	    }
	  break;
	}

      case QUAD:
	{
	  dshapes(0,0) = -(1-xi(1));
	  dshapes(0,1) = -(1-xi(0));
	  dshapes(1,0) =  (1-xi(1));
	  dshapes(1,1) =    -xi(0);
	  dshapes(2,0) =     xi(1);
	  dshapes(2,1) =     xi(0);
	  dshapes(3,0) =    -xi(1);
	  dshapes(3,1) =  (1-xi(0));

	  if (info.order == 1) return;

	  T shapes[4] = {
	    (1-xi(0))*(1-xi(1)),
	    xi(0) *(1-xi(1)),
	    xi(0) *   xi(1) ,
	    (1-xi(0))*   xi(1) 
	  };

	  T mu[4] = { 
	    1 - xi(0) + 1 - xi(1), 
	    xi(0) + 1 - xi(1), 
	    xi(0) +     xi(1), 
	    1 - xi(0) +     xi(1), 
	  };

	  T dmu[4][2] = {
	    { -1, -1 },
	    { 1, -1 },
	    { 1, 1 },
	    { -1, 1 } };
	    
	  // double hshapes[20], hdshapes[20];
	  NgArrayMem<T, 20> hshapes(order+1), hdshapes(order+1);

	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (QUAD);
	  
	  for (int i = 0; i < 4; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcEdgeShapeDx (eorder, mu[vi1]-mu[vi2], &hshapes[0], &hdshapes[0]);

		  T lame = shapes[vi1]+shapes[vi2];
		  T dlame[2] = {
		    dshapes(vi1, 0) + dshapes(vi2, 0),
		    dshapes(vi1, 1) + dshapes(vi2, 1) };
		    
		  for (int j = 0; j < eorder-1; j++)
		    for (int k = 0; k < 2; k++)
		      dshapes(ii+j, k) = 
			lame * hdshapes[j] * (dmu[vi1][k]-dmu[vi2][k])
			+ dlame[k] * hshapes[j];

		  ii += eorder-1;
		}
	    }

	  /*	  
	   *testout << "quad, dshape = " << endl << dshapes << endl;
	   for (int i = 0; i < 2; i++)
	   {
	   Point<2> xil = xi, xir = xi;
	   Vector shapesl(dshapes.Height()), shapesr(dshapes.Height());
	   xil(i) -= 1e-6;
	   xir(i) += 1e-6;
	   CalcElementShapes (info, xil, shapesl);
	   CalcElementShapes (info, xir, shapesr);
	      
	   for (int j = 0; j < dshapes.Height(); j++)
	   dshapes(j,i) = 1.0 / 2e-6 * (shapesr(j)-shapesl(j));
	   }
	  
	   *testout << "quad, num dshape = " << endl << dshapes << endl;
	   */
	  break;
	}
      default:
	throw NgException("CurvedElements::CalcDShape 2d, element type not handled");

      };
  }

  template <int DIM_SPACE, typename T>
  bool CurvedElements ::
  EvaluateMapping (SurfaceElementInfo & info, const Point<2,T> xi, Point<DIM_SPACE,T> & mx, Mat<DIM_SPACE,2,T> & jac) const
  {
    const Element2d & el = mesh[info.elnr];
    if (rational && info.order >= 2) return false; // not supported     

    AutoDiff<2,T> x(xi(0), 0);
    AutoDiff<2,T> y(xi(1), 1);

    AutoDiff<2,T> mapped_x[DIM_SPACE];
    for (int i = 0; i < DIM_SPACE; i++)
      mapped_x[i] = AutoDiff<2,T>(0.0);
    
    switch (el.GetType())
      {
      case TRIG6:
        {
          AutoDiff<2,T> lam3 = 1-x-y;
          AutoDiff<2,T> lami[6] = { x * (2*x-1), y * (2*y-1), lam3 * (2*lam3-1),
                                    4 * y * lam3, 4 * x * lam3, 4 * x * y };
          for (int j = 0; j < 6; j++)
            {
              Point<3> p = mesh[el[j]];
              for (int k = 0; k < DIM_SPACE; k++)
                mapped_x[k] += p(k) * lami[j];
            }
          break;
        }
        
      case TRIG:
        {
          // if (info.order >= 2) return false; // not yet supported
          AutoDiff<2,T> lami[4] = { x, y, 1-x-y };
          for (int j = 0; j < 3; j++)
            {
              Point<3> p = mesh[el[j]];
              for (int k = 0; k < DIM_SPACE; k++)
                mapped_x[k] += p(k) * lami[j];
            }
          if (info.order == 1) break;
          
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TRIG);
	  for (int i = 0; i < 3; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
                  int first = edgecoeffsindex[info.edgenrs[i]];
                  
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShapeLambda (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2],
                                             [&](int i, AutoDiff<2,T> shape)
                                             {
                                               for (int k = 0; k < DIM_SPACE; k++)
                                                 mapped_x[k] += edgecoeffs[first+i](k) * shape;
                                             });
		}              
	    }
          
          int forder = faceorder[info.facenr];
          if (forder >= 3)
            {
              int first = facecoeffsindex[info.facenr];
              
              int fnums[] = { 0, 1, 2 };
              if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
              if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
              if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
              
              CalcScaledTrigShapeLambda (forder, 
                                         lami[fnums[1]]-lami[fnums[0]], lami[fnums[2]], AutoDiff<2,T>(1.0),
                                         [&](int i, AutoDiff<2,T> shape)
                                         {
                                           for (int k = 0; k < DIM_SPACE; k++)
                                             mapped_x[k] += facecoeffs[first+i](k) * shape;
                                         });
            }
          break;
        }
      case QUAD: 
        {
          if (info.order >= 2) return false; // not yet supported
          AutoDiff<2,T> lami[4] = { (1-x)*(1-y), x*(1-y), x*y, (1-x)*y };
          for (int j = 0; j < 4; j++)
            {
              Point<3> p = mesh[el[j]];
              for (int k = 0; k < DIM_SPACE; k++)
                mapped_x[k] += p(k) * lami[j];
            }
          break;
        }
      case QUAD8:
        {
          // AutoDiff<2,T> lami[4] = { (1-x)*(1-y), x*(1-y), x*y, (1-x)*y };
          AutoDiff<2,T> lami[8] =
            { (1-x)*(1-y),
              x*(1-y),
              x*y,
              (1-x)*y,
              4*(1-x)*x*(1-y),
              4*(1-x)*x*y,
              4*(1-y)*y*(1-x), 
              4*(1-y)*y*x };
          
          lami[0] -= 0.5*(lami[4]+lami[6]);
          lami[1] -= 0.5*(lami[4]+lami[7]);
          lami[2] -= 0.5*(lami[5]+lami[7]);
          lami[3] -= 0.5*(lami[5]+lami[6]);
          
          for (int j = 0; j < 8; j++)
            {
              Point<3> p = mesh[el[j]];
              for (int k = 0; k < DIM_SPACE; k++)
                mapped_x[k] += p(k) * lami[j];
            }
          break;
        }
        
      default:
        return false;
      }
        
    for (int i = 0; i < DIM_SPACE; i++)
      {
        mx(i) = mapped_x[i].Value();
        for (int j = 0; j < 2; j++)
          jac(i,j) = mapped_x[i].DValue(j);
      }
    return true;
  }

  template <int DIM_SPACE>
  void CurvedElements :: 
  GetCoefficients (SurfaceElementInfo & info, NgArray<Vec<DIM_SPACE> > & coefs) const
  {
    const Element2d & el = mesh[info.elnr];
    coefs.SetSize (info.ndof);
    
    for (int i = 0; i < info.nv; i++)
      {
	Point<3> hv = mesh[el[i]];
	for (int j = 0; j < DIM_SPACE; j++)
	  coefs[i](j) = hv(j);
      }
    
    if (info.order == 1) return;

    int ii = info.nv;
	  
    for (int i = 0; i < info.edgenrs.Size(); i++)
      {
	int first = edgecoeffsindex[info.edgenrs[i]];
	int next = edgecoeffsindex[info.edgenrs[i]+1];
	for (int j = first; j < next; j++, ii++)
	  for (int k = 0; k < DIM_SPACE; k++)
	    coefs[ii](k) = edgecoeffs[j](k);
      }
    
    int first = facecoeffsindex[info.facenr];
    int next = facecoeffsindex[info.facenr+1];
    for (int j = first; j < next; j++, ii++)
      for (int k = 0; k < DIM_SPACE; k++)
	coefs[ii](k) = facecoeffs[j](k);
  }


  template void CurvedElements :: 
  GetCoefficients<2> (SurfaceElementInfo & info, NgArray<Vec<2> > & coefs) const;

  template void CurvedElements :: 
  GetCoefficients<3> (SurfaceElementInfo & info, NgArray<Vec<3> > & coefs) const;





  // ********************** Transform volume elements *******************


  bool CurvedElements :: IsElementCurved (ElementIndex elnr) const
  {
    if (mesh[elnr].GetType() != TET) return true;
    
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	return mesh.coarsemesh->GetCurvedElements().IsElementCurved (hpref_el.coarse_elnr);
      }

    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    int nfaces = MeshTopology::GetNFaces (type);
    if (nfaces > 4)
      { // not a tet
	const ELEMENT_FACE * faces = MeshTopology::GetFaces0 (type);
	for (int j = 0; j < nfaces; j++)
	  {
	    if (faces[j][3] != -1)
	      {  // a quad face
		Point<3> pts[4];
		for (int k = 0; k < 4; k++)
		  pts[k] = mesh.Point(el[faces[j][k]]);
		Vec<3> twist = (pts[1] - pts[0]) - (pts[2]-pts[3]);
		if (twist.Length() > 1e-8 * (pts[1]-pts[0]).Length())
		  return true;
	      }
	  }
      }
      
    

    ElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = MeshTopology::GetNPoints (type);
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	  info.edgenrs[i]--;

	info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	  info.facenrs[i]--;

	for (int i = 0; i < info.nedges; i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	for (int i = 0; i < info.nfaces; i++)
	  info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
      }

    return (info.ndof > info.nv);
  }


  bool CurvedElements :: IsElementHighOrder (ElementIndex elnr) const
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	return mesh.coarsemesh->GetCurvedElements().IsElementHighOrder (hpref_el.coarse_elnr);
      }

    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    ElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = MeshTopology::GetNPoints (type);
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++) info.edgenrs[i]--;

	info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++) info.facenrs[i]--;

	for (int i = 0; i < info.nedges; i++)
          if (edgecoeffsindex[info.edgenrs[i]+1] > edgecoeffsindex[info.edgenrs[i]]) return true;
	for (int i = 0; i < info.nfaces; i++)
          if (facecoeffsindex[info.facenrs[i]+1] > facecoeffsindex[info.facenrs[i]]) return true;
      }
    return false;
  }






  void CurvedElements :: 
  CalcElementTransformation (Point<3> xi, ElementIndex elnr,
			     Point<3> * x, Mat<3,3> * dxdxi, //  bool * curved,
			     void * buffer, bool valid)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	  
	// xi umrechnen
	double lami[8];
	FlatVector vlami(8, lami);
	vlami = 0;
	mesh[elnr].GetShapeNew<double> (xi, vlami);

	Mat<3,3> trans, dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<3> dlami(8);
	    dlami = 0;
	    mesh[elnr].GetDShapeNew (xi, dlami);	  
	      
	    trans = 0;
	    for (int k = 0; k < 3; k++)
	      for (int l = 0; l < 3; l++)
		for (int i = 0; i < hpref_el.np; i++)
		  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
	  }

	Point<3> coarse_xi(0,0,0);
	for (int i = 0; i < hpref_el.np; i++)
	  for (int j = 0; j < 3; j++)
	    coarse_xi(j) += hpref_el.param[i][j] * lami[i];

	mesh.coarsemesh->GetCurvedElements().CalcElementTransformation (coarse_xi, hpref_el.coarse_elnr, x, &dxdxic /* , curved */);

	if (dxdxi)
	  *dxdxi = dxdxic * trans;

	return;
      }


    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    ElementInfo hinfo;
    ElementInfo & info = (buffer) ? *static_cast<ElementInfo*> (buffer) : hinfo;
    

    if (!valid)
      {
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = MeshTopology::GetNPoints (type);
	if (info.order > 1)
	  {
	    const MeshTopology & top = mesh.GetTopology();
            
	    info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	    for (int i = 0; i < info.nedges; i++)
	      info.edgenrs[i]--;
            
	    info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	    for (int i = 0; i < info.nfaces; i++)
	      info.facenrs[i]--;
            
	    for (int i = 0; i < info.nedges; i++)
	      info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	    for (int i = 0; i < info.nfaces; i++)
	      info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
	  }
      }

    NgArrayMem<double,100> mem(info.ndof);
    TFlatVector<double> shapes(info.ndof, &mem[0]);
    NgArrayMem<double,100> dshapes_mem(info.ndof*3);
    MatrixFixWidth<3> dshapes(info.ndof, &dshapes_mem[0]);
    
    CalcElementShapes (info, xi, shapes);

    Vec<3> * coefs =  (info.ndof <= 10) ? 
      &info.hcoefs[0] : new Vec<3> [info.ndof];

    if (info.ndof > 10 || !valid)
      GetCoefficients (info, coefs);

    if (x)
      {
	*x = 0;
	for (int i = 0; i < shapes.Size(); i++)
	  *x += shapes(i) * coefs[i];
      }

    if (dxdxi)
      {
	if (valid && info.order == 1 && info.nv == 4)   // a linear tet
	  {
	    *dxdxi = info.hdxdxi;
	  }
	else
	  {
	    CalcElementDShapes (info, xi, dshapes);
            
	    *dxdxi = 0;
	    for (int i = 0; i < shapes.Size(); i++)
	      for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
		  (*dxdxi)(j,k) += dshapes(i,k) * coefs[i](j);
            
	    info.hdxdxi = *dxdxi;
	  }
      }

    // *testout << "curved_elements, dshapes = " << endl << dshapes << endl;

    //    if (curved) *curved = (info.ndof > info.nv);

    if (info.ndof > 10) delete [] coefs;
  }



  template <typename T>
  void CurvedElements :: CalcElementShapes (ElementInfo & info, Point<3,T> xi, TFlatVector<T> shapes) const
  {
    const Element & el = mesh[info.elnr];

    if (rational && info.order >= 2)
      {
	// shapes.SetSize(10);
	T w = 1;
	T lami[4] = { xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2) };
	for (int j = 0; j < 4; j++)
	  shapes(j) = lami[j] * lami[j];

	const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TET);
	for (int j = 0; j < 6; j++)
	  {
	    double wi = edgeweight[info.edgenrs[j]];
	    shapes(j+4) = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	    w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	  }

	shapes *= 1.0 / w;
	return;
      }

    // shapes.SetSize(info.ndof);
    
    switch (el.GetType())
      {
      case TET:
	{
	  shapes(0) = xi(0);
	  shapes(1) = xi(1);
	  shapes(2) = xi(2);
	  shapes(3) = 1-xi(0)-xi(1)-xi(2);

	  if (info.order == 1) return;

	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TET);
	  for (int i = 0; i < 6; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShape (eorder, shapes(vi1)-shapes(vi2), shapes(vi1)+shapes(vi2), &shapes(ii));
		  ii += eorder-1;
		}
	    }
	  const ELEMENT_FACE * faces = MeshTopology::GetFaces1 (TET);
	  for (int i = 0; i < 4; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if (forder >= 3)
		{
		  int fnums[] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 }; 
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		  if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);

		  CalcScaledTrigShape (forder, 
				       shapes(fnums[1])-shapes(fnums[0]), shapes(fnums[2]), 
				       shapes(fnums[0])+shapes(fnums[1])+shapes(fnums[2]), &shapes(ii));
		  ii += (forder-1)*(forder-2)/2;
		}
	    }

	  break;
	}
        
      case TET10:
	{
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
	  T lam4 = 1 - x - y - z;
	  /*
	    shapes(0) = xi(0);
	    shapes(1) = xi(1);
	    shapes(2) = xi(2);
	    shapes(3) = 1-xi(0)-xi(1)-xi(2);
	  */
          
	  shapes(0) = 2 * x * x - x;  
	  shapes(1) = 2 * y * y - y;
	  shapes(2) = 2 * z * z - z;
	  shapes(3) = 2 * lam4 * lam4 - lam4;
          
	  shapes(4) = 4 * x * y;
	  shapes(5) = 4 * x * z;
	  shapes(6) = 4 * x * lam4;
	  shapes(7) = 4 * y * z;
	  shapes(8) = 4 * y * lam4;
	  shapes(9) = 4 * z * lam4;

	  break;
	}

      case PRISM:
	{
	  T lami[6] = { xi(0), xi(1), 1-xi(0)-xi(1), xi(0), xi(1), 1-xi(0)-xi(1) };
	  T lamiz[6] = { 1-xi(2), 1-xi(2), 1-xi(2), xi(2), xi(2), xi(2) };
	  for (int i = 0; i < 6; i++)
	    shapes(i) = lami[i] * lamiz[i]; 
	  for (int i = 6; i < info.ndof; i++)
	    shapes(i) = 0;

	  if (info.order == 1) return;


	  int ii = 6;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (PRISM);
	  for (int i = 0; i < 6; i++)    // horizontal edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShape (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &shapes(ii));
		  T facz = (i < 3) ? (1-xi(2)) : xi(2);
		  for (int j = 0; j < eorder-1; j++)
		    shapes(ii+j) *= facz;

		  ii += eorder-1;
		}
	    }

	  for (int i = 6; i < 9; i++)    // vertical edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  T bubz = lamiz[vi1]*lamiz[vi2];
		  T polyz = lamiz[vi1] - lamiz[vi2];
		  T bubxy = lami[vi1];

		  for (int j = 0; j < eorder-1; j++)
		    {
		      shapes(ii+j) = bubxy * bubz;
		      bubz *= polyz;
		    }
		  ii += eorder-1;
		}
	    }

	  // FACE SHAPES
	  const ELEMENT_FACE * faces = MeshTopology::GetFaces1 (PRISM);
	  for (int i = 0; i < 2; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if ( forder < 3 ) continue;
	      int fav[3] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 };
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 
	      if(el[fav[1]] > el[fav[2]]) swap(fav[1],fav[2]);
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 	

	      CalcTrigShape (forder, 
			     lami[fav[2]]-lami[fav[1]], lami[fav[0]],
			     &shapes(ii));
	      
	      int ndf = (forder+1)*(forder+2)/2 - 3 - 3*(forder-1);
	      for ( int j = 0; j < ndf; j++ )
		shapes(ii+j) *= lamiz[fav[1]];
	      ii += ndf;
	    }
	  break;
	}

      case PRISM15:
        {
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
          T lam = 1-x-y;
          T lamz = 1-z;
          shapes[0] = (2*x*x-x) * (2*lamz*lamz-lamz);
          shapes[1] = (2*y*y-y) * (2*lamz*lamz-lamz);
          shapes[2] = (2*lam*lam-lam) * (2*lamz*lamz-lamz);
          shapes[3] = (2*x*x-x) * (2*z*z-z);
          shapes[4] = (2*y*y-y) * (2*z*z-z);
          shapes[5] = (2*lam*lam-lam) * (2*z*z-z);
          shapes[6] = 4 * x * y * (2*lamz*lamz-lamz);
          shapes[7] = 4 * x * lam * (2*lamz*lamz-lamz);
          shapes[8] = 4 * y * lam * (2*lamz*lamz-lamz);
          shapes[9] = x * 4 * z * (1-z);
          shapes[10] = y * 4 * z * (1-z);
          shapes[11] = lam * 4 * z * (1-z);
          shapes[12] = 4 * x * y * (2*z*z-z);
          shapes[13] = 4 * x * lam * (2*z*z-z);
          shapes[14] = 4 * y * lam * (2*z*z-z);
          break;
        }

      case PYRAMID:
	{
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
	  
	  // if (z == 1.) z = 1-1e-10;
          z *= (1-1e-12);
	  shapes[0] = (1-z-x)*(1-z-y) / (1-z);
	  shapes[1] = x*(1-z-y) / (1-z);
	  shapes[2] = x*y / (1-z);
	  shapes[3] = (1-z-x)*y / (1-z);
	  shapes[4] = z;
          
	  if (info.order == 1) return;

          T sigma[4] =
            {
              sigma[0] = ( (1-z-x) + (1-z-y) ),
              sigma[1] = (       x + (1-z-y) ),
              sigma[2] = (       x +       y ),
              sigma[3] = ( (1-z-x) +       y ),
            };

	  int ii = 5;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (PYRAMID);
	  for (int i = 0; i < 4; i++)    // horizontal edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = (edges[i][0]-1), vi2 = (edges[i][1]-1);
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

                  CalcScaledEdgeShape (eorder, sigma[vi1]-sigma[vi2], 1-z, &shapes(ii));
		  T fac = (shapes[vi1]+shapes[vi2]) / (1-z);
		  for (int j = 0; j < eorder-1; j++)
		    shapes(ii+j) *= fac;

		  ii += eorder-1;
		}
	    }



	  break;
	}

      case PYRAMID13:
        {
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
          z *= 1-1e-12;
          shapes[0] = (-z + z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (-2*x - z + 2)*(-2*y - z + 2))*(-0.5*x - 0.5*y - 0.5*z + 0.25);
          shapes[1] = (0.5*x - 0.5*y - 0.25)*(-z - z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (2*x + z)*(-2*y - z + 2));
          shapes[2] = (-z + z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (2*x + z)*(2*y + z))*(0.5*x + 0.5*y + 0.5*z - 0.75);
          shapes[3] = (-0.5*x + 0.5*y - 0.25)*(-z - z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + (2*y + z)*(-2*x - z + 2));
          shapes[4] = z*(2*z - 1);
          shapes[5] = 2*x*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-2*z + 2);
          shapes[6] = 4*x*y*(-2*x - 2*z + 2)/(-2*z + 2);
          shapes[7] = 2*y*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-2*z + 2);
          shapes[8] = 4*x*y*(-2*y - 2*z + 2)/(-2*z + 2);
          shapes[9] = z*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-z + 1);
          shapes[10] = 2*x*z*(-2*y - 2*z + 2)/(-z + 1);
          shapes[11] = 4*x*y*z/(-z + 1);
          shapes[12] = 2*y*z*(-2*x - 2*z + 2)/(-z + 1);
          break;
        }

      case HEX:
	{
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
	  
	  shapes[0] = (1-x)*(1-y)*(1-z);
	  shapes[1] =    x *(1-y)*(1-z);
	  shapes[2] =    x *   y *(1-z);
	  shapes[3] = (1-x)*   y *(1-z);
	  shapes[4] = (1-x)*(1-y)*(z);
	  shapes[5] =    x *(1-y)*(z);
	  shapes[6] =    x *   y *(z);
	  shapes[7] = (1-x)*   y *(z);

	  if (info.order == 1) return;
	  
	  T mu[8] = {
            (1-x)+(1-y)+(1-z),
            x    +(1-y)+(1-z),
            x    +   y +(1-z),
            (1-x)+   y +(1-z),
            (1-x)+(1-y)+(z),
            x    +(1-y)+(z),
            x    +   y +(z),
            (1-x)+   y +(z),
          };
	    
	  int ii = 8;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (HEX);
	  
	  for (int i = 0; i < 8; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcEdgeShape (eorder, mu[vi1]-mu[vi2], &shapes(ii));
		  T lame = shapes(vi1)+shapes(vi2);
		  for (int j = 0; j < order-1; j++)
		    shapes(ii+j) *= lame;
		  ii += eorder-1;
		}
	    }

          
	  break;
        }
        
      case HEX20:
	{
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
	  
	  shapes[0] = (1-x)*(1-y)*(1-z);
	  shapes[1] =    x *(1-y)*(1-z);
	  shapes[2] =    x *   y *(1-z);
	  shapes[3] = (1-x)*   y *(1-z);
	  shapes[4] = (1-x)*(1-y)*(z);
	  shapes[5] =    x *(1-y)*(z);
	  shapes[6] =    x *   y *(z);
	  shapes[7] = (1-x)*   y *(z);

          T sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                      (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 

          static const int e[12][2] =
            {
              { 0, 1 }, { 2, 3 }, { 3, 0 }, { 1, 2 },
              { 4, 5 }, { 6, 7 }, { 7, 4 }, { 5, 6 },
              { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
            };
          
          for (int i = 0; i < 12; i++)
            {
              T lame = shapes[e[i][0]]+shapes[e[i][1]];
              T xi = sigma[e[i][1]]-sigma[e[i][0]];
              shapes[8+i] = (1-xi*xi)*lame;
            }
          for (int i = 0; i < 12; i++)
            {
              shapes[e[i][0]] -= 0.5 * shapes[8+i];
              shapes[e[i][1]] -= 0.5 * shapes[8+i];
            }
          break;
	}

      default:
	throw NgException("CurvedElements::CalcShape 3d, element type not handled");

      };
  }


  template <typename T>
  void CurvedElements :: 
  CalcElementDShapes (ElementInfo & info, const Point<3,T> xi, MatrixFixWidth<3,T> dshapes) const
  {
    // static int timer = NgProfiler::CreateTimer ("calcelementdshapes");
    
    const Element & el = mesh[info.elnr];

    // dshapes.SetSize(info.ndof);
    if ( (long int)(&dshapes(0,0)) % alignof(T) != 0)
      throw NgException ("alignment problem");
    if (dshapes.Height() != info.ndof)
      throw NgException ("wrong height");
    if (rational && info.order >= 2)
      {
	T w = 1;
	T dw[3] = { 0, 0, 0 };

	T lami[4] = { xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2) };
	T dlami[4][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { -1, -1, -1 }};
	T shapes[10];

	for (int j = 0; j < 4; j++)
	  {
	    shapes[j] = lami[j] * lami[j];
	    dshapes(j,0) = 2 * lami[j] * dlami[j][0];
	    dshapes(j,1) = 2 * lami[j] * dlami[j][1];
	    dshapes(j,2) = 2 * lami[j] * dlami[j][2];
	  }

	const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TET);
	for (int j = 0; j < 6; j++)
	  {
	    T wi = edgeweight[info.edgenrs[j]];

	    shapes[j+4] = 2 * wi * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	    for (int k = 0; k < 3; k++)
	      dshapes(j+4,k) = 2*wi* (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
				      lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);

	    w += (wi-1) * 2 * lami[edges[j][0]-1] * lami[edges[j][1]-1];
	    for (int k = 0; k < 3; k++)
	      dw[k] += 2*(wi-1) * (lami[edges[j][0]-1] * dlami[edges[j][1]-1][k] +
				   lami[edges[j][1]-1] * dlami[edges[j][0]-1][k]);
	  }
	// shapes *= 1.0 / w;
	dshapes *= 1.0 / w;
	for (int i = 0; i < 10; i++)
	  for (int j = 0; j < 3; j++)
	    dshapes(i,j) -= shapes[i] * dw[j] / (w*w);
	return;
      }

    /*
    if (typeid(T) == typeid(SIMD<double>))
      {
        if (el.GetType() == HEX)
          dshapes = T(0.0);
        return;
      }
    */
    switch (el.GetType())
      {
      case TET:
	{
          // if (typeid(T) == typeid(SIMD<double>)) return;
          
          dshapes = T(0.0);
          
	  dshapes(0,0) = 1;
	  dshapes(1,1) = 1;
	  dshapes(2,2) = 1;
	  dshapes(3,0) = -1;
	  dshapes(3,1) = -1;
	  dshapes(3,2) = -1;

	  if (info.order == 1) return;

	  T lami[] = { xi(0), xi(1), xi(2), 1-xi(0)-xi(1)-xi(2) };
	  int ii = 4;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TET);
	  for (int i = 0; i < 6; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShapeDxDt<3> (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &dshapes(ii,0));

		  Mat<2,3,T> trans;
		  for (int j = 0; j < 3; j++)
		    {
		      trans(0,j) = dshapes(vi1,j)-dshapes(vi2,j);
		      trans(1,j) = dshapes(vi1,j)+dshapes(vi2,j);
		    }
		  
		  for (int j = 0; j < order-1; j++)
		    {
		      T ddx = dshapes(ii+j,0);
		      T ddt = dshapes(ii+j,1);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		      dshapes(ii+j,2) = ddx * trans(0,2) + ddt * trans(1,2);
		    }

		  ii += eorder-1;
		}
	    }

	  const ELEMENT_FACE * faces = MeshTopology::GetFaces1 (TET);
	  for (int i = 0; i < 4; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if (forder >= 3)
		{
		  int fnums[] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 }; 
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		  if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);

		  CalcScaledTrigShapeDxDyDt (forder, 
					     lami[fnums[1]]-lami[fnums[0]], 
					     lami[fnums[2]], lami[fnums[0]]+lami[fnums[1]]+lami[fnums[2]],
					     &dshapes(ii,0));

		  Mat<3,3,T> trans;
		  for (int j = 0; j < 3; j++)
		    {
		      trans(0,j) = dshapes(fnums[1],j)-dshapes(fnums[0],j);
		      trans(1,j) = dshapes(fnums[2],j);
		      trans(2,j) = dshapes(fnums[0],j)+dshapes(fnums[1],j)+dshapes(fnums[2],j);
		    }
		  
		  int nfd = (forder-1)*(forder-2)/2;
		  for (int j = 0; j < nfd; j++)
		    {
		      T ddx = dshapes(ii+j,0);
		      T ddy = dshapes(ii+j,1);
		      T ddt = dshapes(ii+j,2);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddy * trans(1,0) + ddt * trans(2,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddy * trans(1,1) + ddt * trans(2,1);
		      dshapes(ii+j,2) = ddx * trans(0,2) + ddy * trans(1,2) + ddt * trans(2,2);
		    }

		  ii += nfd;
		}
	    }

	  break;
	}

      case TET10:
	{
          // if (typeid(T) == typeid(SIMD<double>)) return;
          
	  if (dshapes.Height() == 4)
	    {
	      dshapes = T(0.0);

	      dshapes(0,0) = 1;
	      dshapes(1,1) = 1;
	      dshapes(2,2) = 1;
	      dshapes(3,0) = -1;
	      dshapes(3,1) = -1;
	      dshapes(3,2) = -1;
	    }
	  else
	    {
	      AutoDiff<3,T> x(xi(0), 0);
	      AutoDiff<3,T> y(xi(1), 1);
	      AutoDiff<3,T> z(xi(2), 2);
	      AutoDiff<3,T> lam4 = 1-x-y-z;
	      AutoDiff<3,T> shapes[10];
              
	      shapes[0] = 2 * x * x - x;  
	      shapes[1] = 2 * y * y - y;
	      shapes[2] = 2 * z * z - z;
	      shapes[3] = 2 * lam4 * lam4 - lam4;
              
	      shapes[4] = 4 * x * y;
	      shapes[5] = 4 * x * z;
	      shapes[6] = 4 * x * lam4;
	      shapes[7] = 4 * y * z;
	      shapes[8] = 4 * y * lam4;
	      shapes[9] = 4 * z * lam4;

	      for (int i = 0; i < 10; i++)
		{
		  dshapes(i,0) = shapes[i].DValue(0);
		  dshapes(i,1) = shapes[i].DValue(1);
		  dshapes(i,2) = shapes[i].DValue(2);
		}
	      
	    }
	  break;

	  break;
	}


      case PRISM:
	{
	  T lami[6] = { xi(0), xi(1), 1-xi(0)-xi(1), xi(0), xi(1), 1-xi(0)-xi(1)  };
	  T lamiz[6] = { 1-xi(2), 1-xi(2), 1-xi(2), xi(2), xi(2), xi(2) };
	  T dlamiz[6] = { -1, -1, -1, 1, 1, 1 };
	  T dlami[6][2] = 
	    { { 1, 0, },
	      { 0, 1, },
	      { -1, -1 },
	      { 1, 0, },
	      { 0, 1, },
	      { -1, -1 } };
	  for (int i = 0; i < 6; i++)
	    {
	      // shapes(i) = lami[i%3] * ( (i < 3) ? (1-xi(2)) : xi(2) );
	      dshapes(i,0) = dlami[i%3][0] * ( (i < 3) ? (1-xi(2)) : xi(2) );
	      dshapes(i,1) = dlami[i%3][1] * ( (i < 3) ? (1-xi(2)) : xi(2) );
	      dshapes(i,2) = lami[i%3] * ( (i < 3) ? -1 : 1 );
	    }

	  int ii = 6;

	  if (info.order == 1) return;

          
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (PRISM);
	  for (int i = 0; i < 6; i++)    // horizontal edges
	    {
	      int order = edgeorder[info.edgenrs[i]];
	      if (order >= 2)
		{
		  int vi1 = (edges[i][0]-1), vi2 = (edges[i][1]-1);
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);
		  vi1 = vi1 % 3;
		  vi2 = vi2 % 3;

                  NgArrayMem<T,20> shapei_mem(order+1);
		  TFlatVector<T> shapei(order+1, &shapei_mem[0]);
		  CalcScaledEdgeShapeDxDt<3> (order, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &dshapes(ii,0) );
		  CalcScaledEdgeShape(order, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2], &shapei(0) );

		  Mat<2,2,T> trans;
		  for (int j = 0; j < 2; j++)
		    {
		      trans(0,j) = dlami[vi1][j]-dlami[vi2][j];
		      trans(1,j) = dlami[vi1][j]+dlami[vi2][j];
		    }
		  
		  for (int j = 0; j < order-1; j++)
		    {
		      T ddx = dshapes(ii+j,0);
		      T ddt = dshapes(ii+j,1);
		      dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		      dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		    }



		  T facz = (i < 3) ? (1-xi(2)) : xi(2);
		  T dfacz = (i < 3) ? (-1) : 1;
		  for (int j = 0; j < order-1; j++)
		    {
		      dshapes(ii+j,0) *= facz;
		      dshapes(ii+j,1) *= facz;
		      dshapes(ii+j,2) = shapei(j) * dfacz;
		    }

		  ii += order-1;
		}
	    }

          // if (typeid(T) == typeid(SIMD<double>)) return;


	  for (int i = 6; i < 9; i++)    // vertical edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = (edges[i][0]-1), vi2 = (edges[i][1]-1);
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  T bubz = lamiz[vi1] * lamiz[vi2];
		  T dbubz = dlamiz[vi1]*lamiz[vi2] + lamiz[vi1]*dlamiz[vi2];
		  T polyz = lamiz[vi1] - lamiz[vi2];
		  T dpolyz = dlamiz[vi1] - dlamiz[vi2];
		  T bubxy = lami[(vi1)%3];
		  T dbubxydx = dlami[(vi1)%3][0];
		  T dbubxydy = dlami[(vi1)%3][1];

		  for (int j = 0; j < eorder-1; j++)
		    {
		      dshapes(ii+j,0) = dbubxydx * bubz;
		      dshapes(ii+j,1) = dbubxydy * bubz;
		      dshapes(ii+j,2) = bubxy * dbubz;

		      dbubz = bubz * dpolyz + dbubz * polyz;
		      bubz *= polyz;
		    }
		  ii += eorder-1;
		}
	    }


	  if (info.order == 2) return;
	  // FACE SHAPES
	  const ELEMENT_FACE * faces = MeshTopology::GetFaces1 (PRISM);
	  for (int i = 0; i < 2; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];

	      if ( forder < 3 ) continue;
	      int ndf = (forder+1)*(forder+2)/2 - 3 - 3*(forder-1);

	      int fav[3] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 };
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 
	      if(el[fav[1]] > el[fav[2]]) swap(fav[1],fav[2]);
	      if(el[fav[0]] > el[fav[1]]) swap(fav[0],fav[1]); 	

              NgArrayMem<T,2*20> dshapei_mem(ndf);
              NgArrayMem<T,20> shapei_mem(ndf);
	      MatrixFixWidth<2,T> dshapei(ndf, &dshapei_mem[0]);
	      TFlatVector<T> shapei(ndf, &shapei_mem[0]);

	      CalcTrigShapeDxDy (forder, 
				 lami[fav[2]]-lami[fav[1]], lami[fav[0]],
				 &dshapei(0,0));
	      CalcTrigShape (forder, lami[fav[2]]-lami[fav[1]], lami[fav[0]],
			     &shapei(0));
	      
	      Mat<2,2,T> trans;
	      for (int j = 0; j < 2; j++)
		{
		  trans(0,j) = dlami[fav[2]][j]-dlami[fav[1]][j];
		  trans(1,j) = dlami[fav[0]][j];
		}
		  
	      for (int j = 0; j < ndf; j++)
		{
		  // double ddx = dshapes(ii+j,0);
		  // double ddt = dshapes(ii+j,1);
		  T ddx = dshapei(j,0);
		  T ddt = dshapei(j,1);
		  dshapes(ii+j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		  dshapes(ii+j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		}

	      for ( int j = 0; j < ndf; j++ )
		{
		  dshapes(ii+j,0) *= lamiz[fav[1]];
		  dshapes(ii+j,1) *= lamiz[fav[1]];
		  dshapes(ii+j,2) = shapei(j) * dlamiz[fav[1]];
		}
	      ii += ndf;
	    }

	  break;

	}

      case PRISM15:
        {
          AutoDiff<3,T> x(xi(0), 0);
          AutoDiff<3,T> y(xi(1), 1);
          AutoDiff<3,T> z(xi(2), 2);
          AutoDiff<3,T> ad[15];
          AutoDiff<3,T> lam = 1-x-y;
          AutoDiff<3,T> lamz = 1-z;

          ad[0] = (2*x*x-x) * (2*lamz*lamz-lamz);
          ad[1] = (2*y*y-y) * (2*lamz*lamz-lamz);
          ad[2] = (2*lam*lam-lam) * (2*lamz*lamz-lamz);
          ad[3] = (2*x*x-x) * (2*z*z-z);
          ad[4] = (2*y*y-y) * (2*z*z-z);
          ad[5] = (2*lam*lam-lam) * (2*z*z-z);
          ad[6] = 4 * x * y * (2*lamz*lamz-lamz);
          ad[7] = 4 * x * lam * (2*lamz*lamz-lamz);
          ad[8] = 4 * y * lam * (2*lamz*lamz-lamz);
          ad[9] = x * 4 * z * (1-z);
          ad[10] = y * 4 * z * (1-z);
          ad[11] = lam * 4 * z * (1-z);
          ad[12] = 4 * x * y * (2*z*z-z);
          ad[13] = 4 * x * lam * (2*z*z-z);
          ad[14] = 4 * y * lam * (2*z*z-z);

          for(int i=0; i<15; i++)
            for(int j=0; j<3; j++)
              dshapes(i,j) = ad[i].DValue(j);
          break;
        }
      case PYRAMID:
	{
          // if (typeid(T) == typeid(SIMD<double>)) return;
          
	  dshapes = T(0.0);
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
	  
	  // if (z == 1.) z = 1-1e-10;
          z *= 1-1e-12;
	  T z1 = 1-z;
	  T z2 = z1*z1;
	  
	  dshapes(0,0) = -(z1-y)/z1;
	  dshapes(0,1) = -(z1-x)/z1;
	  dshapes(0,2) = ((x+y+2*z-2)*z1+(z1-y)*(z1-x))/z2;

	  dshapes(1,0) = (z1-y)/z1;
	  dshapes(1,1) = -x/z1;
	  dshapes(1,2) = (-x*z1+x*(z1-y))/z2;

	  dshapes(2,0) = y/z1;
	  dshapes(2,1) = x/z1;
	  dshapes(2,2) = x*y/z2;

	  dshapes(3,0) = -y/z1;
	  dshapes(3,1) = (z1-x)/z1;
	  dshapes(3,2) = (-y*z1+y*(z1-x))/z2;

	  dshapes(4,0) = 0;
	  dshapes(4,1) = 0;
	  dshapes(4,2) = 1;

	  if (info.order == 1) return;

	  int ii = 5;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (PYRAMID);
	  // if (z == 1.) z = 1-1e-10;
          z *= 1-1e-12;
          T shapes[5];
	  shapes[0] = (1-z-x)*(1-z-y) / (1-z);
	  shapes[1] = x*(1-z-y) / (1-z);
	  shapes[2] = x*y / (1-z);
	  shapes[3] = (1-z-x)*y / (1-z);
	  shapes[4] = z;

          T sigma[4] =
            {
              ( (1-z-x) + (1-z-y) ),
              (       x + (1-z-y) ),
              (       x +       y ),
              ( (1-z-x) +       y ),
            };
          T dsigma[4][3] =
            {
              { -1, -1, -2 },
              { 1, -1, -1 },
              { 1, 1, 0 },
              { -1, 1, -1 }
            };
          T dz[3] = { 0, 0, 1 };
	  for (int i = 0; i < 4; i++)    // horizontal edges
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = (edges[i][0]-1), vi2 = (edges[i][1]-1);
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

                  NgArrayMem<T,20> shapei_mem(eorder+1);
		  TFlatVector<T> shapei(eorder+1,&shapei_mem[0]);
		  CalcScaledEdgeShapeDxDt<3> (eorder, sigma[vi1]-sigma[vi2], 1-z, &dshapes(ii,0) );
		  CalcScaledEdgeShape(eorder, sigma[vi1]-sigma[vi2], 1-z, &shapei(0) );
		  T fac = (shapes[vi1]+shapes[vi2]) / (1-z);
                  T dfac[3];
                  for (int k = 0; k < 3; k++)
                    dfac[k] = ( (dshapes(vi1,k)+dshapes(vi2,k)) * (1-z) -
                                (shapes[vi1]+shapes[vi2]) *(-dshapes(4,k)) )
                      / sqr(1-z);
                      
		  for (int j = 0; j < eorder-1; j++)
                    {
                      T ddx = dshapes(ii+j,0);
                      T ddt = dshapes(ii+j,1);
                      for (int k = 0; k < 3; k++)
                        dshapes(ii+j,k) = fac * (ddx * (dsigma[vi1][k]-dsigma[vi2][k]) - ddt*dz[k])
                          + dfac[k] * shapei(j);
                    }

		  ii += eorder-1;
		}
	    }
          
	  break;
	}

      case PYRAMID13:
        {
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
          z *= 1-1e-12;
          dshapes(0,0) = 0.5*z - 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) - 0.5*(-2*x - z + 2)*(-2*y - z + 2) + (-0.5*x - 0.5*y - 0.5*z + 0.25)*(4*y + 2*z + 2*z*(2*y + z - 1)/(-z + 1) - 4);
          dshapes(0,1) = 0.5*z - 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) - 0.5*(-2*x - z + 2)*(-2*y - z + 2) + (-0.5*x - 0.5*y - 0.5*z + 0.25)*(4*x + 2*z + 2*z*(2*x + z - 1)/(-z + 1) - 4);
          dshapes(0,2) = 0.5*z - 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) - 0.5*(-2*x - z + 2)*(-2*y - z + 2) + (-0.5*x - 0.5*y - 0.5*z + 0.25)*(2*x + 2*y + 2*z + z*(2*x + z - 1)/(-z + 1) + z*(2*y + z - 1)/(-z + 1) + z*(2*x + z - 1)*(2*y + z - 1)/((-z + 1)*(-z + 1)) - 5 + (2*x + z - 1)*(2*y + z - 1)/(-z + 1));
          dshapes(1,0) = -0.5*z - 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + 0.5*(2*x + z)*(-2*y - z + 2) + (0.5*x - 0.5*y - 0.25)*(-4*y - 2*z - 2*z*(2*y + z - 1)/(-z + 1) + 4);
          dshapes(1,1) = 0.5*z + 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) - 0.5*(2*x + z)*(-2*y - z + 2) + (-4*x - 2*z - 2*z*(2*x + z - 1)/(-z + 1))*(0.5*x - 0.5*y - 0.25);
          dshapes(1,2) = (0.5*x - 0.5*y - 0.25)*(-2*x - 2*y - 2*z - z*(2*x + z - 1)/(-z + 1) - z*(2*y + z - 1)/(-z + 1) - z*(2*x + z - 1)*(2*y + z - 1)/((-z + 1)*(-z + 1)) + 1 - (2*x + z - 1)*(2*y + z - 1)/(-z + 1));
          dshapes(2,0) = -0.5*z + 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + 0.5*(2*x + z)*(2*y + z) + (4*y + 2*z + 2*z*(2*y + z - 1)/(-z + 1))*(0.5*x + 0.5*y + 0.5*z - 0.75);
          dshapes(2,1) = -0.5*z + 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + 0.5*(2*x + z)*(2*y + z) + (4*x + 2*z + 2*z*(2*x + z - 1)/(-z + 1))*(0.5*x + 0.5*y + 0.5*z - 0.75);
          dshapes(2,2) = -0.5*z + 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + 0.5*(2*x + z)*(2*y + z) + (0.5*x + 0.5*y + 0.5*z - 0.75)*(2*x + 2*y + 2*z + z*(2*x + z - 1)/(-z + 1) + z*(2*y + z - 1)/(-z + 1) + z*(2*x + z - 1)*(2*y + z - 1)/((-z + 1)*(-z + 1)) - 1 + (2*x + z - 1)*(2*y + z - 1)/(-z + 1));
          dshapes(3,0) = 0.5*z + 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) - 0.5*(2*y + z)*(-2*x - z + 2) + (-0.5*x + 0.5*y - 0.25)*(-4*y - 2*z - 2*z*(2*y + z - 1)/(-z + 1));
          dshapes(3,1) = -0.5*z - 0.5*z*(2*x + z - 1)*(2*y + z - 1)/(-z + 1) + 0.5*(2*y + z)*(-2*x - z + 2) + (-0.5*x + 0.5*y - 0.25)*(-4*x - 2*z - 2*z*(2*x + z - 1)/(-z + 1) + 4);
          dshapes(3,2) = (-0.5*x + 0.5*y - 0.25)*(-2*x - 2*y - 2*z - z*(2*x + z - 1)/(-z + 1) - z*(2*y + z - 1)/(-z + 1) - z*(2*x + z - 1)*(2*y + z - 1)/((-z + 1)*(-z + 1)) + 1 - (2*x + z - 1)*(2*y + z - 1)/(-z + 1));
          dshapes(4,0) = 0;
          dshapes(4,1) = 0;
          dshapes(4,2) = 4*z - 1;
          dshapes(5,0) = -4*x*(-2*y - 2*z + 2)/(-2*z + 2) + 2*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-2*z + 2);
          dshapes(5,1) = -4*x*(-2*x - 2*z + 2)/(-2*z + 2);
          dshapes(5,2) = -4*x*(-2*x - 2*z + 2)/(-2*z + 2) - 4*x*(-2*y - 2*z + 2)/(-2*z + 2) + 4*x*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/((-2*z + 2)*(-2*z + 2));
          dshapes(6,0) = -8*x*y/(-2*z + 2) + 4*y*(-2*x - 2*z + 2)/(-2*z + 2);
          dshapes(6,1) = 4*x*(-2*x - 2*z + 2)/(-2*z + 2);
          dshapes(6,2) = -8*x*y/(-2*z + 2) + 8*x*y*(-2*x - 2*z + 2)/((-2*z + 2)*(-2*z + 2));
          dshapes(7,0) = -4*y*(-2*y - 2*z + 2)/(-2*z + 2);
          dshapes(7,1) = -4*y*(-2*x - 2*z + 2)/(-2*z + 2) + 2*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-2*z + 2);
          dshapes(7,2) = -4*y*(-2*x - 2*z + 2)/(-2*z + 2) - 4*y*(-2*y - 2*z + 2)/(-2*z + 2) + 4*y*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/((-2*z + 2)*(-2*z + 2));
          dshapes(8,0) = 4*y*(-2*y - 2*z + 2)/(-2*z + 2);
          dshapes(8,1) = -8*x*y/(-2*z + 2) + 4*x*(-2*y - 2*z + 2)/(-2*z + 2);
          dshapes(8,2) = -8*x*y/(-2*z + 2) + 8*x*y*(-2*y - 2*z + 2)/((-2*z + 2)*(-2*z + 2));
          dshapes(9,0) = -2*z*(-2*y - 2*z + 2)/(-z + 1);
          dshapes(9,1) = -2*z*(-2*x - 2*z + 2)/(-z + 1);
          dshapes(9,2) = -2*z*(-2*x - 2*z + 2)/(-z + 1) - 2*z*(-2*y - 2*z + 2)/(-z + 1) + z*(-2*x - 2*z + 2)*(-2*y - 2*z + 2)/((-z + 1)*(-z + 1)) + (-2*x - 2*z + 2)*(-2*y - 2*z + 2)/(-z + 1);
          dshapes(10,0) = 2*z*(-2*y - 2*z + 2)/(-z + 1);
          dshapes(10,1) = -4*x*z/(-z + 1);
          dshapes(10,2) = -4*x*z/(-z + 1) + 2*x*z*(-2*y - 2*z + 2)/((-z + 1)*(-z + 1)) + 2*x*(-2*y - 2*z + 2)/(-z + 1);
          dshapes(11,0) = 4*y*z/(-z + 1);
          dshapes(11,1) = 4*x*z/(-z + 1);
          dshapes(11,2) = 4*x*y*z/((-z + 1)*(-z + 1)) + 4*x*y/(-z + 1);
          dshapes(12,0) = -4*y*z/(-z + 1);
          dshapes(12,1) = 2*z*(-2*x - 2*z + 2)/(-z + 1);
          dshapes(12,2) = -4*y*z/(-z + 1) + 2*y*z*(-2*x - 2*z + 2)/((-z + 1)*(-z + 1)) + 2*y*(-2*x - 2*z + 2)/(-z + 1);
          break;
        }

      case HEX:
	{
          // if (typeid(T) == typeid(SIMD<double>)) return;
          
          // NgProfiler::StartTimer(timer);
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);

	  // shapes[0] = (1-x)*(1-y)*(1-z);
	  dshapes(0,0) = - (1-y)*(1-z);
	  dshapes(0,1) = (1-x) * (-1) * (1-z);
	  dshapes(0,2) = (1-x) * (1-y) * (-1);

	  // shapes[1] =    x *(1-y)*(1-z);
	  dshapes(1,0) = (1-y)*(1-z);
	  dshapes(1,1) = -x * (1-z);
	  dshapes(1,2) = -x * (1-y);

	  // shapes[2] =    x *   y *(1-z);
	  dshapes(2,0) = y * (1-z);
	  dshapes(2,1) = x * (1-z);
	  dshapes(2,2) = -x * y;

	  // shapes[3] = (1-x)*   y *(1-z);
	  dshapes(3,0) = -y * (1-z);
	  dshapes(3,1) = (1-x) * (1-z);
	  dshapes(3,2) = -(1-x) * y;

	  // shapes[4] = (1-x)*(1-y)*z;
	  dshapes(4,0) = - (1-y)*z;
	  dshapes(4,1) = (1-x) * (-1) * z;
	  dshapes(4,2) = (1-x) * (1-y) * 1;

	  // shapes[5] =    x *(1-y)*z;
	  dshapes(5,0) = (1-y)*z;
	  dshapes(5,1) = -x * z;
	  dshapes(5,2) = x * (1-y);

	  // shapes[6] =    x *   y *z;
	  dshapes(6,0) = y * z;
	  dshapes(6,1) = x * z;
	  dshapes(6,2) = x * y;

	  // shapes[7] = (1-x)*   y *z;
	  dshapes(7,0) = -y * z;
	  dshapes(7,1) = (1-x) * z;
	  dshapes(7,2) = (1-x) * y;

          // NgProfiler::StopTimer(timer);

	  if (info.order == 1) return;          

	  T shapes[8] = {
            (1-x)*(1-y)*(1-z),
               x *(1-y)*(1-z),
               x *   y *(1-z),
            (1-x)*   y *(1-z),
            (1-x)*(1-y)*(z),
               x *(1-y)*(z),
               x *   y *(z),
            (1-x)*   y *(z),
	  };

	  T mu[8] = {
            (1-x)+(1-y)+(1-z),
            x    +(1-y)+(1-z),
            x    +   y +(1-z),
            (1-x)+   y +(1-z),
            (1-x)+(1-y)+(z),
            x    +(1-y)+(z),
            x    +   y +(z),
            (1-x)+   y +(z)
	  };

	  T dmu[8][3] = {
	    { -1, -1, -1 },
	    { 1, -1, -1 },
	    { 1, 1, -1 },
	    { -1, 1, -1 },
	    { -1, -1, 1 },
	    { 1, -1, 1 },
	    { 1, 1, 1 },
	    { -1, 1, 1 }
          };
	    
	  NgArrayMem<T, 20> hshapes(order+1), hdshapes(order+1);

	  int ii = 8;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (HEX);
	  
	  for (int i = 0; i < 8; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcEdgeShapeDx (eorder, mu[vi1]-mu[vi2], &hshapes[0], &hdshapes[0]);

		  T lame = shapes[vi1]+shapes[vi2];
		  T dlame[3] = {
		    dshapes(vi1, 0) + dshapes(vi2, 0),
		    dshapes(vi1, 1) + dshapes(vi2, 1),
                    dshapes(vi1, 2) + dshapes(vi2, 2)
                  };
		    
		  for (int j = 0; j < eorder-1; j++)
		    for (int k = 0; k < 3; k++)
		      dshapes(ii+j, k) = 
			lame * hdshapes[j] * (dmu[vi1][k]-dmu[vi2][k])
			+ dlame[k] * hshapes[j];

		  ii += eorder-1;
		}
	    }

	  /*	  
	   *testout << "quad, dshape = " << endl << dshapes << endl;
	   for (int i = 0; i < 2; i++)
	   {
	   Point<2> xil = xi, xir = xi;
	   Vector shapesl(dshapes.Height()), shapesr(dshapes.Height());
	   xil(i) -= 1e-6;
	   xir(i) += 1e-6;
	   CalcElementShapes (info, xil, shapesl);
	   CalcElementShapes (info, xir, shapesr);
	      
	   for (int j = 0; j < dshapes.Height(); j++)
	   dshapes(j,i) = 1.0 / 2e-6 * (shapesr(j)-shapesl(j));
	   }
	  
	   *testout << "quad, num dshape = " << endl << dshapes << endl;
	   */
	  break;
	}
      case HEX20:
        {
          AutoDiff<3,T> x(xi(0), 0);
          AutoDiff<3,T> y(xi(1), 1);
          AutoDiff<3,T> z(xi(2), 2);
          AutoDiff<3,T> ad[20];
          
          ad[0] = (1-x)*(1-y)*(1-z);
	  ad[1] =    x *(1-y)*(1-z);
	  ad[2] =    x *   y *(1-z);
	  ad[3] = (1-x)*   y *(1-z);
	  ad[4] = (1-x)*(1-y)*(z);
	  ad[5] =    x *(1-y)*(z);
	  ad[6] =    x *   y *(z);
	  ad[7] = (1-x)*   y *(z);
          
          AutoDiff<3,T> sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                                  (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
          
          static const int e[12][2] =
            {
              { 0, 1 }, { 2, 3 }, { 3, 0 }, { 1, 2 },
              { 4, 5 }, { 6, 7 }, { 7, 4 }, { 5, 6 },
              { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
            };
          
          for (int i = 0; i < 12; i++)
            {
              auto lame = ad[e[i][0]]+ad[e[i][1]];
              auto xi = sigma[e[i][1]]-sigma[e[i][0]];
              ad[8+i] = (1-xi*xi)*lame;
            }
          for (int i = 0; i < 12; i++)
            {
              ad[e[i][0]] -= 0.5 * ad[8+i];
              ad[e[i][1]] -= 0.5 * ad[8+i];
            }
          for (int i = 0; i < 20; i++)
            for (int j = 0; j < 3; j++)
              dshapes(i,j) = ad[i].DValue(j);
          break;
        }
      default:
	throw NgException("CurvedElements::CalcDShape 3d, element type not handled");
      }
    
    /*
      DenseMatrix dshapes2 (info.ndof, 3);
      Vector shapesl(info.ndof); 
      Vector shapesr(info.ndof); 
    
      double eps = 1e-6;
      for (int i = 0; i < 3; i++)
      {
      Point<3> xl = xi;
      Point<3> xr = xi;
	
      xl(i) -= eps;
      xr(i) += eps;
      CalcElementShapes (info, xl, shapesl);
      CalcElementShapes (info, xr, shapesr);
	
      for (int j = 0; j < info.ndof; j++)
      dshapes2(j,i) = (shapesr(j)-shapesl(j)) / (2*eps);
      }
      (*testout) << "dshapes = " << endl << dshapes << endl;
      (*testout) << "dshapes2 = " << endl << dshapes2 << endl;
      dshapes2 -= dshapes;
      (*testout) << "diff = " << endl << dshapes2 << endl;
    */
  }

  // extern int mappingvar;
  template <typename T>
  bool CurvedElements ::
  EvaluateMapping (ElementInfo & info, Point<3,T> xi, Point<3,T> & mx, Mat<3,3,T> & jac) const
  {
    const Element & el = mesh[info.elnr];
    if (rational && info.order >= 2) return false; // not supported     

    AutoDiff<3,T> x(xi(0), 0);
    AutoDiff<3,T> y(xi(1), 1);
    AutoDiff<3,T> z(xi(2), 2);

    AutoDiff<3,T> mapped_x[3] = { T(0.0), T(0.0), T(0.0) } ;
    
    switch (el.GetType())
      {
      case TET:
        {
          // if (info.order >= 2) return false; // not yet supported
          AutoDiff<3,T> lami[4] = { x, y, z, 1-x-y-z };
          for (int j = 0; j < 4; j++)
            {
              Point<3> p = mesh[el[j]];
              for (int k = 0; k < 3; k++)
                mapped_x[k] += p(k) * lami[j];
            }
          if (info.order == 1) break;

	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (TET);
	  for (int i = 0; i < 6; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
                  int first = edgecoeffsindex[info.edgenrs[i]];
                  
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

		  CalcScaledEdgeShapeLambda (eorder, lami[vi1]-lami[vi2], lami[vi1]+lami[vi2],
                                             [&](int i, AutoDiff<3,T> shape)
                                             {
                                               Vec<3> coef = edgecoeffs[first+i];
                                               for (int k = 0; k < 3; k++)
                                                 mapped_x[k] += coef(k) * shape;
                                             });
		}              
	    }
          
	  const ELEMENT_FACE * faces = MeshTopology::GetFaces1 (TET);
	  for (int i = 0; i < 4; i++)
	    {
	      int forder = faceorder[info.facenrs[i]];
	      if (forder >= 3)
		{
                  int first = facecoeffsindex[info.facenrs[i]];
                  
		  int fnums[] = { faces[i][0]-1, faces[i][1]-1, faces[i][2]-1 }; 
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		  if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
		  if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);

		  CalcScaledTrigShapeLambda (forder, 
                                             lami[fnums[1]]-lami[fnums[0]], lami[fnums[2]], 
                                             lami[fnums[0]]+lami[fnums[1]]+lami[fnums[2]],
                                             [&](int i, AutoDiff<3,T> shape)
                                             {
                                               Vec<3> coef = facecoeffs[first+i];
                                               for (int k = 0; k < 3; k++)
                                                 mapped_x[k] += coef(k) * shape;
                                             });
                }
	    }
          
          break;
        }
      case HEX:
        {
          if (info.order >= 2) return false; // not yet supported          
          AutoDiff<3,T> lami[8] =
            { (1-x)*(1-y)*(1-z),
              (  x)*(1-y)*(1-z),
              (  x)*   y *(1-z),
              (1-x)*   y *(1-z),
              (1-x)*(1-y)*(z),
              (  x)*(1-y)*(z),
              (  x)*   y *(z),
              (1-x)*   y *(z) };

          for (int j = 0; j < 8; j++)
            {
              Point<3> p = mesh[el[j]];
              for (int k = 0; k < 3; k++)
                mapped_x[k] += p(k) * lami[j];
            }

	  if (info.order == 1) break;

	  AutoDiff<3,T> mu[8] = {
            (1-x)+(1-y)+(1-z),
            x    +(1-y)+(1-z),
            x    +   y +(1-z),
            (1-x)+   y +(1-z),
            (1-x)+(1-y)+(z),
            x    +(1-y)+(z),
            x    +   y +(z),
            (1-x)+   y +(z),
          };
	  // int ii = 8;
	  const ELEMENT_EDGE * edges = MeshTopology::GetEdges1 (HEX);
	  
	  for (int i = 0; i < 8; i++)
	    {
	      int eorder = edgeorder[info.edgenrs[i]];
	      if (eorder >= 2)
		{
                  int first = edgecoeffsindex[info.edgenrs[i]];                  
		  int vi1 = edges[i][0]-1, vi2 = edges[i][1]-1;
		  if (el[vi1] > el[vi2]) swap (vi1, vi2);

                  AutoDiff<3,T> lame = lami[vi1]+lami[vi2];
		  CalcEdgeShapeLambda (eorder, mu[vi1]-mu[vi2], 
                                       [&](int i, AutoDiff<3,T> shape)
                                       {
                                         Vec<3> coef = edgecoeffs[first+i];
                                         for (int k = 0; k < 3; k++)
                                           mapped_x[k] += coef(k) * (lame*shape);
                                       });
                  
		}
	    }
          
          break;
        }
      default:
        return false;
      }

    for (int i = 0; i < 3; i++)
      {
        mx(i) = mapped_x[i].Value();
        for (int j = 0; j < 3; j++)
          jac(i,j) = mapped_x[i].DValue(j);
      }
    return true;
  }
  



  
  void CurvedElements :: 
  GetCoefficients (ElementInfo & info, Vec<3> * coefs) const
  {
    const Element & el = mesh[info.elnr];

    for (int i = 0; i < info.nv; i++)
      coefs[i] = Vec<3> (mesh[el[i]]);

    if (info.order == 1) return;

    int ii = info.nv;
	  
    for (int i = 0; i < info.nedges; i++)
      {
	int first = edgecoeffsindex[info.edgenrs[i]];
	int next = edgecoeffsindex[info.edgenrs[i]+1];
	for (int j = first; j < next; j++, ii++)
	  coefs[ii] = edgecoeffs[j];
      }
    for (int i = 0; i < info.nfaces; i++)
      {
	int first = facecoeffsindex[info.facenrs[i]];
	int next = facecoeffsindex[info.facenrs[i]+1];
	for (int j = first; j < next; j++, ii++)
	  coefs[ii] = facecoeffs[j];
      }
  }






























  /*
  void CurvedElements :: 
  CalcMultiPointSegmentTransformation (NgArray<double> * xi, SegmentIndex segnr,
				       NgArray<Point<3> > * x,
				       NgArray<Vec<3> > * dxdxi)
  {
    ;
  }
  */

  template <int DIM_SPACE, typename T>
  void CurvedElements :: 
  CalcMultiPointSegmentTransformation (SegmentIndex elnr, int n,
				       const T * xi, size_t sxi,
				       T * x, size_t sx,
				       T * dxdxi, size_t sdxdxi)
  {
    for (int ip = 0; ip < n; ip++)
      {
	Point<3,T> xg;
	Vec<3,T> dx;

	// mesh->GetCurvedElements().
	CalcSegmentTransformation<T> (xi[ip*sxi], elnr, &xg, &dx);
      
	if (x)
	  for (int i = 0; i < DIM_SPACE; i++)
	    x[ip*sx+i] = xg(i);
	  
	if (dxdxi)
	  for (int i=0; i<DIM_SPACE; i++)
	    dxdxi[ip*sdxdxi+i] = dx(i);
      }
  }


  template void CurvedElements :: 
  CalcMultiPointSegmentTransformation<2> (SegmentIndex elnr, int npts,
					  const double * xi, size_t sxi,
					  double * x, size_t sx,
					  double * dxdxi, size_t sdxdxi);

  template void CurvedElements :: 
  CalcMultiPointSegmentTransformation<3> (SegmentIndex elnr, int npts,
					  const double * xi, size_t sxi,
					  double * x, size_t sx,
					  double * dxdxi, size_t sdxdxi);

  template void CurvedElements :: 
  CalcMultiPointSegmentTransformation<2> (SegmentIndex elnr, int npts,
					  const SIMD<double> * xi, size_t sxi,
					  SIMD<double> * x, size_t sx,
					  SIMD<double> * dxdxi, size_t sdxdxi);

  template void CurvedElements :: 
  CalcMultiPointSegmentTransformation<3> (SegmentIndex elnr, int npts,
					  const SIMD<double> * xi, size_t sxi,
					  SIMD<double> * x, size_t sx,
					  SIMD<double> * dxdxi, size_t sdxdxi);

  template void CurvedElements :: 
  CalcSegmentTransformation<double> (double xi, SegmentIndex elnr,
                                     Point<3,double> * x, Vec<3,double> * dxdxi, bool * curved);


  void CurvedElements :: 
  CalcMultiPointSurfaceTransformation (NgArray< Point<2> > * xi, SurfaceElementIndex elnr,
				       NgArray< Point<3> > * x,
				       NgArray< Mat<3,2> > * dxdxi)
  {
    double * px = (x) ? &(*x)[0](0) : NULL;
    double * pdxdxi = (dxdxi) ? &(*dxdxi)[0](0) : NULL;

    CalcMultiPointSurfaceTransformation <3> (elnr, xi->Size(),
					     &(*xi)[0](0), 2, 
					     px, 3,
					     pdxdxi, 6);
  }




  template <int DIM_SPACE, typename T>
  void CurvedElements :: 
  CalcMultiPointSurfaceTransformation (SurfaceElementIndex elnr, int npts,
				       const T * xi, size_t sxi,
				       T * x, size_t sx,
				       T * dxdxi, size_t sdxdxi)
  {
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	// xi umrechnen
	T lami[4];
	TFlatVector<T> vlami(4, lami);

	NgArrayMem<Point<2,T>, 50> coarse_xi (npts);
	
	for (int pi = 0; pi < npts; pi++)
	  {
	    vlami = 0;
	    Point<2,T> hxi(xi[pi*sxi], xi[pi*sxi+1]);
	    mesh[elnr].GetShapeNew ( hxi, vlami);
	    
	    Point<2,T> cxi(0,0);
	    for (int i = 0; i < hpref_el.np; i++)
	      for (int j = 0; j < 2; j++)
		cxi(j) += hpref_el.param[i][j] * lami[i];

	    coarse_xi[pi] = cxi;
	  }

	mesh.coarsemesh->GetCurvedElements().
	  CalcMultiPointSurfaceTransformation<DIM_SPACE,T> (hpref_el.coarse_elnr, npts,
                                                            &coarse_xi[0](0), &coarse_xi[1](0)-&coarse_xi[0](0), 
                                                            x, sx, dxdxi, sdxdxi);

	// Mat<3,2> dxdxic;
	if (dxdxi)
	  {
            T mem_dlami[8]; // avoid alignment problems if T is SIMD
	    MatrixFixWidth<2,T> dlami(4, mem_dlami);
	    dlami = T(0.0);

	    for (int pi = 0; pi < npts; pi++)
	      {
		Point<2,T> hxi(xi[pi*sxi], xi[pi*sxi+1]);
		mesh[elnr].GetDShapeNew ( hxi, dlami);	  
		
		Mat<2,2,T> trans;
		trans = 0;
		for (int k = 0; k < 2; k++)
		  for (int l = 0; l < 2; l++)
		    for (int i = 0; i < hpref_el.np; i++)
		      trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
		
		Mat<DIM_SPACE,2,T> hdxdxic, hdxdxi;
		for (int k = 0; k < 2*DIM_SPACE; k++)
		  hdxdxic(k) = dxdxi[pi*sdxdxi+k];

		hdxdxi = hdxdxic * trans;

		for (int k = 0; k < 2*DIM_SPACE; k++)
		  dxdxi[pi*sdxdxi+k] = hdxdxi(k);
                    
		// dxdxic = (*dxdxi)[pi];
		// (*dxdxi)[pi] = dxdxic * trans;
	      }
	  }	

	return;
      }


    const Element2d & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    SurfaceElementInfo info;
    info.elnr = elnr;
    info.order = order;
    switch (type)
      {
      case TRIG : info.nv = 3; break;
      case QUAD : info.nv = 4; break;
      case TRIG6: info.nv = 6; break;
      case QUAD8 : info.nv = 8; break;
      default:
	cerr << "undef element in CalcMultPointSurfaceTrafo" << endl;
      }
    info.ndof = info.nv;

    // if (info.order > 1)
    //   {
    //     const MeshTopology & top = mesh.GetTopology();
	
    //     top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
    //     for (int i = 0; i < info.edgenrs.Size(); i++)
    //       info.edgenrs[i]--;
    //     info.facenr = top.GetSurfaceElementFace (elnr+1)-1;

    //     for (int i = 0; i < info.edgenrs.Size(); i++)
    //       info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
    //     info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
    //   }

// Michael Woopen: THESE FOLLOWING LINES ARE COPIED FROM CurvedElements::CalcSurfaceTransformation

    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	  info.edgenrs[i]--;
	info.facenr = top.GetSurfaceElementFace (elnr+1)-1;


	bool firsttry = true;
	bool problem = false;

	while(firsttry || problem)
	  {
	    problem = false;

	    for (int i = 0; !problem && i < info.edgenrs.Size(); i++)
	      {
		if(info.edgenrs[i]+1 >= edgecoeffsindex.Size())
		  problem = true;
		else
		  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	      }
	    if(info.facenr+1 >= facecoeffsindex.Size())
	      problem = true;
	    else
	      info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];

	    if(problem && !firsttry)
	      throw NgException("something wrong with curved elements");
	    
	    if(problem)
	      BuildCurvedElements(NULL,order,rational);

	    firsttry = false;
	  }
      }



    bool ok = true;
    for (int i = 0; i < npts; i++)
      {
        Point<2,T> _xi(xi[i*sxi], xi[i*sxi+1]);
        Point<DIM_SPACE,T> _x;
        Mat<DIM_SPACE,2,T> _dxdxi;
        if (!EvaluateMapping (info, _xi, _x, _dxdxi))
          { ok = false; break; }
        // *testout << "x = " << _x << ", dxdxi = " << _dxdxi << endl;
        if (x)
          for (int j = 0; j < DIM_SPACE; j++)
            x[i*sx+j] = _x[j];
        if (dxdxi)
          for (int j = 0; j < DIM_SPACE; j++)
            for (int k = 0; k < 2; k++)
              dxdxi[i*sdxdxi+2*j+k] = _dxdxi(j,k);
      }
    if (ok) return;

    
// THESE LAST LINES ARE COPIED FROM CurvedElements::CalcSurfaceTransformation

    NgArrayMem<Vec<DIM_SPACE>,100> coefs(info.ndof);
    GetCoefficients (info, coefs);
    
    NgArrayMem<T, 100> shapes_mem(info.ndof);
    TFlatVector<T> shapes(info.ndof, &shapes_mem[0]);

    NgArrayMem<T, 100> dshapes_mem(info.ndof*2);
    MatrixFixWidth<2,T> dshapes(info.ndof,&shapes_mem[0]);



    if (x)
      {
	if (info.order == 1 && type == TRIG)
	  {
	    for (int j = 0; j < npts; j++)
	      {
		Point<2,T> vxi(xi[j*sxi], xi[j*sxi+1]);

		Point<DIM_SPACE,T> val;
                for (int k = 0; k < DIM_SPACE; k++)
                  val(k) = coefs[2](k) + (coefs[0](k)-coefs[2](k)) * vxi(0) + (coefs[1](k)-coefs[2](k)) * vxi(1);
                /*
                (coefs[2]);
		val += (coefs[0]-coefs[2]) * vxi(0);
		val += (coefs[1]-coefs[2]) * vxi(1);
                */
		for (int k = 0; k < DIM_SPACE; k++)
		  x[j*sx+k] = val(k);
	      }
	  }
	else
	  for (int j = 0; j < npts; j++)
	    {
	      Point<2,T> vxi(xi[j*sxi], xi[j*sxi+1]);
	      CalcElementShapes (info, vxi, shapes);
	      
	      Point<DIM_SPACE,T> val = T(0.0);
	      for (int i = 0; i < coefs.Size(); i++)
                for (int k = 0; k < DIM_SPACE; k++)
                  val(k) += shapes(i) * coefs[i](k);
	      
	      for (int k = 0; k < DIM_SPACE; k++)
		x[j*sx+k] = val(k);
	    }
      }

    if (dxdxi)
      {
	if (info.order == 1 && type == TRIG)
	  {
	    Point<2,T> xij(xi[0], xi[1]);
	    CalcElementDShapes (info, xij, dshapes);
	    
	    Mat<3,2,T> dxdxij;
	    dxdxij = 0.0;
	    for (int i = 0; i < coefs.Size(); i++)
	      for (int j = 0; j < DIM_SPACE; j++)
		for (int k = 0; k < 2; k++)
		  dxdxij(j,k) += dshapes(i,k) * coefs[i](j);
	    
		
	    for (int ip = 0; ip < npts; ip++)
	      for (int j = 0; j < DIM_SPACE; j++)
		for (int k = 0; k < 2; k++)
		  dxdxi[ip*sdxdxi+2*j+k] = dxdxij(j,k);
	  }
	else
	  {
	    for (int j = 0; j < npts; j++)
	      {
		Point<2,T> vxi(xi[j*sxi], xi[j*sxi+1]);
		CalcElementDShapes (info, vxi, dshapes);
		
		Mat<DIM_SPACE,2,T> ds;
		ds = 0.0;
		for (int i = 0; i < coefs.Size(); i++)
		  for (int j = 0; j < DIM_SPACE; j++)
		    for (int k = 0; k < 2; k++)
		      ds(j,k) += dshapes(i,k) * coefs[i](j);
		// (*dxdxi)[ip] = ds;
		
		for (int k = 0; k < 2*DIM_SPACE; k++)
		  dxdxi[j*sdxdxi+k] = ds(k);
	      }
	  }
      }
  }



  template void CurvedElements :: 
  CalcMultiPointSurfaceTransformation<2> (SurfaceElementIndex elnr, int npts,
					  const double * xi, size_t sxi,
					  double * x, size_t sx,
					  double * dxdxi, size_t sdxdxi);

  template void CurvedElements :: 
  CalcMultiPointSurfaceTransformation<3> (SurfaceElementIndex elnr, int npts,
					  const double * xi, size_t sxi,
					  double * x, size_t sx,
					  double * dxdxi, size_t sdxdxi);


  template void CurvedElements :: 
  CalcMultiPointSurfaceTransformation<2> (SurfaceElementIndex elnr, int npts,
					  const SIMD<double> * xi, size_t sxi,
					  SIMD<double> * x, size_t sx,
					  SIMD<double> * dxdxi, size_t sdxdxi);

  template void CurvedElements :: 
  CalcMultiPointSurfaceTransformation<3> (SurfaceElementIndex elnr, int npts,
					  const SIMD<double> * xi, size_t sxi,
					  SIMD<double> * x, size_t sx,
					  SIMD<double> * dxdxi, size_t sdxdxi);











  void CurvedElements :: 
  CalcMultiPointElementTransformation (NgArray< Point<3> > * xi, ElementIndex elnr,
				       NgArray< Point<3> > * x,
				       NgArray< Mat<3,3> > * dxdxi)
  {
    double * px = (x) ? &(*x)[0](0) : NULL;
    double * pdxdxi = (dxdxi) ? &(*dxdxi)[0](0) : NULL;

    CalcMultiPointElementTransformation (elnr, xi->Size(),
					 &(*xi)[0](0), 3, 
					 px, 3,
					 pdxdxi, 9);
    
    return;
#ifdef OLD

    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	// xi umrechnen
	double lami[8];
	FlatVector vlami(8, lami);


	NgArrayMem<Point<3>, 50> coarse_xi (xi->Size());
	
	for (int pi = 0; pi < xi->Size(); pi++)
	  {
	    vlami = 0;
	    mesh[elnr].GetShapeNew ( (*xi)[pi], vlami);
	    
	    Point<3> cxi(0,0,0);
	    for (int i = 0; i < hpref_el.np; i++)
	      for (int j = 0; j < 3; j++)
		cxi(j) += hpref_el.param[i][j] * lami[i];

	    coarse_xi[pi] = cxi;
	  }

	mesh.coarsemesh->GetCurvedElements().
	  CalcMultiPointElementTransformation (&coarse_xi, hpref_el.coarse_elnr, x, dxdxi);


	Mat<3,3> trans, dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<3> dlami(8);
	    dlami = 0;

	    for (int pi = 0; pi < xi->Size(); pi++)
	      {
		mesh[elnr].GetDShapeNew ( (*xi)[pi], dlami);	  
		
		trans = 0;
		for (int k = 0; k < 3; k++)
		  for (int l = 0; l < 3; l++)
		    for (int i = 0; i < hpref_el.np; i++)
		      trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
		
		dxdxic = (*dxdxi)[pi];
		(*dxdxi)[pi] = dxdxic * trans;
	      }
	  }	

	return;
      }








    Vector shapes;
    MatrixFixWidth<3> dshapes;


    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();

    ElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = MeshTopology::GetNPoints (type);
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	  info.edgenrs[i]--;

	info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	  info.facenrs[i]--;

	for (int i = 0; i < info.nedges; i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	for (int i = 0; i < info.nfaces; i++)
	  info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
	// info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
      }

    NgArray<Vec<3> > coefs(info.ndof);
    GetCoefficients (info, &coefs[0]);
    if (x)
      {
	for (int j = 0; j < xi->Size(); j++)
	  {
	    CalcElementShapes (info, (*xi)[j], shapes);
	    (*x)[j] = 0;
	    for (int i = 0; i < coefs.Size(); i++)
	      (*x)[j] += shapes(i) * coefs[i];
	  }
      }

    if (dxdxi)
      {
	if (info.order == 1 && type == TET)
	  {
	    if (xi->Size() > 0)
	      {
		CalcElementDShapes (info, (*xi)[0], dshapes);
		Mat<3,3> ds;
		ds = 0;
		for (int i = 0; i < coefs.Size(); i++)
		  for (int j = 0; j < 3; j++)
		    for (int k = 0; k < 3; k++)
		      ds(j,k) += dshapes(i,k) * coefs[i](j);
	    
		for (int ip = 0; ip < xi->Size(); ip++)
		  (*dxdxi)[ip] = ds;
	      }
	  }
	else
	  for (int ip = 0; ip < xi->Size(); ip++)
	    {
	      CalcElementDShapes (info, (*xi)[ip], dshapes);
	      
	      Mat<3,3> ds;
	      ds = 0;
	      for (int i = 0; i < coefs.Size(); i++)
		for (int j = 0; j < 3; j++)
		  for (int k = 0; k < 3; k++)
		    ds(j,k) += dshapes(i,k) * coefs[i](j);
	      (*dxdxi)[ip] = ds;
	    }
      }
#endif
  }


  // extern int multipointtrafovar;
  template <typename T>
  void  CurvedElements :: 
  CalcMultiPointElementTransformation (ElementIndex elnr, int n,
				       const T * xi, size_t sxi,
				       T * x, size_t sx,
				       T * dxdxi, size_t sdxdxi)
  {
    // multipointtrafovar++;
    /*
    static int timer = NgProfiler::CreateTimer ("calcmultipointelementtrafo");
    static int timer1 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 1");
    static int timer2 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 2");
    static int timer3 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 3");
    static int timer4 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 4");
    static int timer5 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 5");
    NgProfiler::RegionTimer reg(timer);
    */
    // NgProfiler::StartTimer (timer);
    // NgProfiler::StartTimer (timer1);
    if (mesh.coarsemesh)
      {
	const HPRefElement & hpref_el =
	  (*mesh.hpelements) [mesh[elnr].hp_elnr];
	
	// xi umrechnen
	T lami[8];
	TFlatVector<T> vlami(8, &lami[0]);


	NgArrayMem<T, 100> coarse_xi (3*n);
	
	for (int pi = 0; pi < n; pi++)
	  {
	    vlami = 0;
	    Point<3,T> pxi;
	    for (int j = 0; j < 3; j++)
	      pxi(j) = xi[pi*sxi+j];

	    mesh[elnr].GetShapeNew (pxi, vlami);
	    
	    Point<3,T> cxi(0,0,0);
	    for (int i = 0; i < hpref_el.np; i++)
	      for (int j = 0; j < 3; j++)
		cxi(j) += hpref_el.param[i][j] * lami[i];

	    for (int j = 0; j < 3; j++)
	      coarse_xi[3*pi+j] = cxi(j);
	  }

	mesh.coarsemesh->GetCurvedElements().
	  CalcMultiPointElementTransformation (hpref_el.coarse_elnr, n, 
					       &coarse_xi[0], 3, 
					       x, sx, 
					       dxdxi, sdxdxi);

	Mat<3,3,T> trans, dxdxic;
	if (dxdxi)
	  {
	    MatrixFixWidth<3,T> dlami(8);
	    dlami = T(0);

	    for (int pi = 0; pi < n; pi++)
	      {
		Point<3,T> pxi;
		for (int j = 0; j < 3; j++)
		  pxi(j) = xi[pi*sxi+j];

		mesh[elnr].GetDShapeNew (pxi, dlami);	  
		
		trans = 0;
		for (int k = 0; k < 3; k++)
		  for (int l = 0; l < 3; l++)
		    for (int i = 0; i < hpref_el.np; i++)
		      trans(l,k) += hpref_el.param[i][l] * dlami(i, k);

		Mat<3,3,T> mat_dxdxic, mat_dxdxi;
		for (int j = 0; j < 3; j++)
		  for (int k = 0; k < 3; k++)
		    mat_dxdxic(j,k) = dxdxi[pi*sdxdxi+3*j+k];
		
		mat_dxdxi = mat_dxdxic * trans;

		for (int j = 0; j < 3; j++)
		  for (int k = 0; k < 3; k++)
		    dxdxi[pi*sdxdxi+3*j+k] = mat_dxdxi(j,k);

		// dxdxic = (*dxdxi)[pi];
		// (*dxdxi)[pi] = dxdxic * trans;
	      }
	  }	
	return;
      }

    // NgProfiler::StopTimer (timer1);
    // NgProfiler::StartTimer (timer2);


    const Element & el = mesh[elnr];
    ELEMENT_TYPE type = el.GetType();


    ElementInfo info;
    info.elnr = elnr;
    info.order = order;
    info.ndof = info.nv = MeshTopology::GetNPoints (type);
    if (info.order > 1)
      {
	const MeshTopology & top = mesh.GetTopology();
	
	info.nedges = top.GetElementEdges (elnr+1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	  info.edgenrs[i]--;

	info.nfaces = top.GetElementFaces (elnr+1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	  info.facenrs[i]--;

	for (int i = 0; i < info.nedges; i++)
	  info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	for (int i = 0; i < info.nfaces; i++)
	  info.ndof += facecoeffsindex[info.facenrs[i]+1] - facecoeffsindex[info.facenrs[i]];
	// info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
      }

    // NgProfiler::StopTimer (timer2);
    // NgProfiler::StartTimer (timer3);


    bool ok = true;
    for (int i = 0; i < n; i++)
      {
        Point<3,T> _xi(xi[i*sxi], xi[i*sxi+1], xi[i*sxi+2]);
        Point<3,T> _x;
        Mat<3,3,T> _dxdxi;
        if (!EvaluateMapping (info, _xi, _x, _dxdxi))
          { ok = false; break; }
        // cout << "x = " << _x << ", dxdxi = " << _dxdxi << endl;
        if (x)
          for (int j = 0; j < 3; j++)
            x[i*sx+j] = _x[j];
        if (dxdxi)
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
              dxdxi[i*sdxdxi+3*j+k] = _dxdxi(j,k);
      }
    if (ok) return;

    NgArrayMem<Vec<3>,100> coefs(info.ndof);
    NgArrayMem<T,500> shapes_mem(info.ndof);
    
    TFlatVector<T> shapes(info.ndof, &shapes_mem[0]);

    NgArrayMem<T,1500> dshapes_mem(3*info.ndof);
    MatrixFixWidth<3,T> dshapes(info.ndof, &dshapes_mem[0]);

    // NgProfiler::StopTimer (timer3);
    // NgProfiler::StartTimer (timer4);

    GetCoefficients (info, &coefs[0]);
    if (x)
      {
	for (int j = 0; j < n; j++)
	  {
	    Point<3,T> xij, xj;
	    for (int k = 0; k < 3; k++)
	      xij(k) = xi[j*sxi+k];
	    CalcElementShapes (info, xij, shapes);
	    xj = T(0.0);
	    for (int i = 0; i < coefs.Size(); i++)
              for (int k = 0; k < 3; k++)
                xj(k) += shapes(i) * coefs[i](k);

            // cout << "old, xj = " << xj << endl;

	    for (int k = 0; k < 3; k++)
	      x[j*sx+k] = xj(k);
	  }
      }


    // NgProfiler::StopTimer (timer4);
    // NgProfiler::StartTimer (timer5);
                
    if (dxdxi)
      {
	if (info.order == 1 && type == TET)
	  {
	    if (n > 0)
	      {

		Point<3,T> xij;
		for (int k = 0; k < 3; k++)
		  xij(k) = xi[k];
		
		CalcElementDShapes (info, xij, dshapes);
		
		Mat<3,3,T> dxdxij;
		dxdxij = 0.0;
		for (int i = 0; i < coefs.Size(); i++)
		  for (int j = 0; j < 3; j++)
		    for (int k = 0; k < 3; k++)
		      dxdxij(j,k) += dshapes(i,k) * coefs[i](j);
		
		
		for (int ip = 0; ip < n; ip++)
		  for (int j = 0; j < 3; j++)
		    for (int k = 0; k < 3; k++)
		      dxdxi[ip*sdxdxi+3*j+k] = dxdxij(j,k);
	      }
	  }
	else
	  {
	    for (int ip = 0; ip < n; ip++)
	      {
		Point<3,T> xij;
		for (int k = 0; k < 3; k++)
		  xij(k) = xi[ip*sxi+k];

                CalcElementDShapes (info, xij, dshapes);


		Mat<3,3,T> dxdxij;
		dxdxij = 0.0;
		for (int i = 0; i < coefs.Size(); i++)
		  for (int j = 0; j < 3; j++)
		    for (int k = 0; k < 3; k++)
		      dxdxij(j,k) += dshapes(i,k) * coefs[i](j);
                
                // cout << "old, jac = " << dxdxij << endl;

		for (int j = 0; j < 3; j++)
		  for (int k = 0; k < 3; k++)
		    dxdxi[ip*sdxdxi+3*j+k] = dxdxij(j,k);

                /*
                T dxdxi00 = T(0.0);
                T dxdxi01 = T(0.0);
                T dxdxi02 = T(0.0);
                T dxdxi10 = T(0.0);
                T dxdxi11 = T(0.0);
                T dxdxi12 = T(0.0);
                T dxdxi20 = T(0.0);
                T dxdxi21 = T(0.0);
                T dxdxi22 = T(0.0);
                
		for (int i = 0; i < coefs.Size(); i++)
                  {
                    T ds0 = dshapes(i,0);
                    T ds1 = dshapes(i,1);
                    T ds2 = dshapes(i,2);
                    T cf0 = coefs[i](0);
                    T cf1 = coefs[i](1);
                    T cf2 = coefs[i](2);

                    dxdxi00 += ds0*cf0;
                    dxdxi01 += ds1*cf0;
                    dxdxi02 += ds2*cf0;
                    dxdxi10 += ds0*cf1;
                    dxdxi11 += ds1*cf1;
                    dxdxi12 += ds2*cf1;
                    dxdxi20 += ds0*cf2;
                    dxdxi21 += ds1*cf2;
                    dxdxi22 += ds2*cf2;
                  }

                dxdxi[ip*sdxdxi+3*0+0] = dxdxi00;
                dxdxi[ip*sdxdxi+3*0+1] = dxdxi01;
                dxdxi[ip*sdxdxi+3*0+2] = dxdxi02;

                dxdxi[ip*sdxdxi+3*1+0] = dxdxi10;
                dxdxi[ip*sdxdxi+3*1+1] = dxdxi11;
                dxdxi[ip*sdxdxi+3*1+2] = dxdxi12;
                
                dxdxi[ip*sdxdxi+3*2+0] = dxdxi20;
                dxdxi[ip*sdxdxi+3*2+1] = dxdxi21;
                dxdxi[ip*sdxdxi+3*2+2] = dxdxi22;
                */
	      }
	  }
      }
    // NgProfiler::StopTimer (timer5);
    // NgProfiler::StopTimer (timer);    
  }


  template
  void  CurvedElements :: 
  CalcMultiPointElementTransformation
  (ElementIndex elnr, int n,
   const double * xi, size_t sxi,
   double * x, size_t sx,
   double * dxdxi, size_t sdxdxi);

  template
  void  CurvedElements :: 
  CalcMultiPointElementTransformation
  (ElementIndex elnr, int n,
   const SIMD<double> * xi, size_t sxi,
   SIMD<double> * x, size_t sx,
   SIMD<double> * dxdxi, size_t sdxdxi);

};

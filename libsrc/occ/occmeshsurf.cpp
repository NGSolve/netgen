#ifdef OCCGEOMETRY

#include <mystdlib.h>

#include <occgeom.hpp>
#include <meshing.hpp>
#include <GeomLProp_SLProps.hxx>
#include <ShapeAnalysis_Surface.hxx>


namespace netgen
{
#include "occmeshsurf.hpp"


  bool glob_testout(false);

  void OCCSurface :: GetNormalVector (const Point<3> & p, 
				      const PointGeomInfo & geominfo,
				      Vec<3> & n) const
  {
    gp_Pnt pnt;
    gp_Vec du, dv;

    /*
      double gu = geominfo.u;
      double gv = geominfo.v;

      if (fabs (gu) < 1e-3) gu = 0;
      if (fabs (gv) < 1e-3) gv = 0;

      occface->D1(gu,gv,pnt,du,dv);
    */

    /*
      occface->D1(geominfo.u,geominfo.v,pnt,du,dv);

      n = Cross (Vec<3>(du.X(), du.Y(), du.Z()),
      Vec<3>(dv.X(), dv.Y(), dv.Z()));
      n.Normalize();
    */



    GeomLProp_SLProps lprop(occface,geominfo.u,geominfo.v,1,1e-5);
    double setu=geominfo.u,setv=geominfo.v;

    if(lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5)
      {
	double ustep = 0.01*(umax-umin);
	double vstep = 0.01*(vmax-vmin);

	n=0;

	while(setu < umax && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setu += ustep;
	if(setu < umax)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setu = geominfo.u;
	while(setu > umin && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setu -= ustep;
	if(setu > umin)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setu = geominfo.u;

	while(setv < vmax && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setv += ustep;
	if(setv < vmax)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setv = geominfo.v;
	while(setv > vmin && (lprop.D1U().Magnitude() < 1e-5 || lprop.D1V().Magnitude() < 1e-5))
	  setv -= ustep;
	if(setv > vmin)
	  {
	    lprop.SetParameters(setu,setv);
	    n(0)+=lprop.Normal().X();
	    n(1)+=lprop.Normal().Y();
	    n(2)+=lprop.Normal().Z();
	  }
	setv = geominfo.v;

	n.Normalize();
      }
    else
      {
	n(0)=lprop.Normal().X();
	n(1)=lprop.Normal().Y();
	n(2)=lprop.Normal().Z();
      }

    if(glob_testout)
      {
	(*testout) << "u " << geominfo.u << " v " << geominfo.v 
		   << " du " << lprop.D1U().X() << " "<< lprop.D1U().Y() << " "<< lprop.D1U().Z()
		   << " dv " << lprop.D1V().X() << " "<< lprop.D1V().Y() << " "<< lprop.D1V().Z() << endl;
      }



    if (orient == TopAbs_REVERSED) n = -1*n;
    //  (*testout) << "GetNormalVector" << endl;
  }


  void OCCSurface :: DefineTangentialPlane (const Point<3> & ap1,
					    const PointGeomInfo & geominfo1,
					    const Point<3> & ap2,
					    const PointGeomInfo & geominfo2)
  {
    if (projecttype == PLANESPACE)
      {
	p1 = ap1; p2 = ap2;

	//cout << "p1 = " << p1 << endl;
	//cout << "p2 = " << p2 << endl;
      
	GetNormalVector (p1, geominfo1, ez);
      
	ex = p2 - p1;
	ex -= (ex * ez) * ez;
	ex.Normalize();
	ey = Cross (ez, ex); 

	GetNormalVector (p2, geominfo2, n2);
  
	nmid = 0.5*(n2+ez);
      
	ez = nmid;
	ez.Normalize(); 
      
	ex = (p2 - p1).Normalize();
	ez -= (ez * ex) * ex;
	ez.Normalize();
	ey = Cross (ez, ex);
	nmid = ez;
	//cout << "ex " << ex << " ey " << ey << " ez " << ez << endl;
      }
    else
      {
	if ( (geominfo1.u < umin) ||
	     (geominfo1.u > umax) ||
	     (geominfo2.u < umin) ||
	     (geominfo2.u > umax) ||
	     (geominfo1.v < vmin) ||
	     (geominfo1.v > vmax) ||
	     (geominfo2.v < vmin) ||
	     (geominfo2.v > vmax) ) throw UVBoundsException();
	  

	p1 = ap1; p2 = ap2;
	psp1 = Point<2>(geominfo1.u, geominfo1.v);
	psp2 = Point<2>(geominfo2.u, geominfo2.v);
      
	Vec<3> n;
	GetNormalVector (p1, geominfo1, n);

	gp_Pnt pnt;
	gp_Vec du, dv;
	occface->D1 (geominfo1.u, geominfo1.v, pnt, du, dv);

        // static Timer t("occ-defintangplane calculations");
        // RegionTimer reg(t);

        Mat<3,2> D1_;
	D1_(0,0) = du.X(); D1_(1,0) = du.Y(); D1_(2,0) = du.Z();
	D1_(0,1) = dv.X(); D1_(1,1) = dv.Y(); D1_(2,1) = dv.Z();
        auto D1T_ = Trans(D1_);
	auto D1TD1_ = D1T_*D1_;
	if (Det (D1TD1_) == 0) throw SingularMatrixException();
        Mat<2,2> DDTinv_;
        CalcInverse (D1TD1_, DDTinv_);

        Mat<3,2> Y_;
	Vec<3> y1_ = (ap2-ap1).Normalize();
	Vec<3> y2_ = Cross(n, y1_).Normalize();
	for (int i = 0; i < 3; i++)
	  {
	    Y_(i,0) = y1_(i);
	    Y_(i,1) = y2_(i);
	  }

        auto A_ = DDTinv_ * D1T_ * Y_;
	Mat<2,2> Ainv_;
	if (Det(A_) == 0) throw SingularMatrixException();
	CalcInverse (A_, Ainv_);

	Vec<2> temp_ = Ainv_ * (psp2-psp1);
	double r_ = temp_.Length();
        Mat<2,2> R_;
        R_(0,0) = temp_(0)/r_;
        R_(1,0) = temp_(1)/r_;
        R_(0,1) = -R_(1,0);
        R_(1,1) = R_(0,0);

        A_ = A_ * R_;
        Ainv_ = Trans(R_) * Ainv_;

        Amat = A_;
        Amatinv = Ainv_;
        
	// temp = Amatinv * (psp2-psp1);
        

#ifdef OLD
	DenseMatrix D1(3,2), D1T(2,3), DDTinv(2,2);
	D1(0,0) = du.X(); D1(1,0) = du.Y(); D1(2,0) = du.Z();
	D1(0,1) = dv.X(); D1(1,1) = dv.Y(); D1(2,1) = dv.Z();

	/*
	  (*testout) << "DefineTangentialPlane" << endl
	  << "---------------------" << endl;
	  (*testout) << "D1 = " << endl << D1 << endl;
	*/

	Transpose (D1, D1T);
	DenseMatrix D1TD1(3,3);

	D1TD1 = D1T*D1;
	if (D1TD1.Det() == 0) throw SingularMatrixException();
      
	CalcInverse (D1TD1, DDTinv);
        // cout << " =?= inv = " << DDTinv << endl;
        
	DenseMatrix Y(3,2);
	Vec<3> y1 = (ap2-ap1).Normalize();
	Vec<3> y2 = Cross(n, y1).Normalize();
	for (int i = 0; i < 3; i++)
	  {
	    Y(i,0) = y1(i);
	    Y(i,1) = y2(i);
	  }

	DenseMatrix A(2,2);
	A = DDTinv * D1T * Y;
	DenseMatrix Ainv(2,2);

	if (A.Det() == 0) throw SingularMatrixException();

	CalcInverse (A, Ainv);

	for (int i = 0; i < 2; i++)
	  for (int j = 0; j < 2; j++)
	    {
	      Amat(i,j) = A(i,j);
	      Amatinv(i,j) = Ainv(i,j);
	    }

	Vec<2> temp = Amatinv * (psp2-psp1);
      

	double r = temp.Length();
	//      double alpha = -acos (temp(0)/r);
	double alpha = -atan2 (temp(1),temp(0));
	DenseMatrix R(2,2);
	R(0,0) = cos (alpha);
	R(1,0) = -sin (alpha);
	R(0,1) = sin (alpha);
	R(1,1) = cos (alpha);

        // cout << "=?= R = " << R << endl;

	A = A*R;

	if (A.Det() == 0) throw SingularMatrixException();

	CalcInverse (A, Ainv);
    

	for (int i = 0; i < 2; i++)
	  for (int j = 0; j < 2; j++)
	    {
	      Amat(i,j) = A(i,j);
	      Amatinv(i,j) = Ainv(i,j);
	    }
        // cout << "=?= Ainv = " << endl << Ainv << endl;
	temp = Amatinv * (psp2-psp1);
        cout << " =?= Amatinv = " << Amatinv << endl;        
#endif
      };
 
  }


  void OCCSurface :: ToPlane (const Point<3> & p3d,
			      const PointGeomInfo & geominfo,
			      Point<2> & pplane, 
			      double h, int & zone) const
  {
    if (projecttype == PLANESPACE)
      {
	Vec<3> p1p, n;
	GetNormalVector (p3d, geominfo, n);
      
	p1p = p3d - p1;
	pplane(0) = (p1p * ex) / h;
	pplane(1) = (p1p * ey) / h;
      
	if (n * nmid < 0)
	  zone = -1;
	else
	  zone = 0;

	/*
	  if(zone == -1)
	  {
	  (*testout) << "zone = -1 for " << p3d << " 2D: " << pplane << " n " << n << " nmid " << nmid << endl;
	  glob_testout = true;
	  GetNormalVector (p3d, geominfo, n);
	  glob_testout = false;
	  }
	*/
      }
    else
      {
	pplane = Point<2>(geominfo.u, geominfo.v);
	//      (*testout) << "(u,v) = " << geominfo.u << ", " << geominfo.v << endl;
	pplane = Point<2> (1/h * (Amatinv * (pplane-psp1)));
	//      pplane = Point<2> (h * (Amatinv * (pplane-psp1)));
	//      pplane = Point<2> (1/h * ((pplane-psp1)));

	zone = 0;
      };
  }	


  void OCCSurface :: FromPlane (const Point<2> & pplane, 
				Point<3> & p3d,
				PointGeomInfo & gi,
				double h) 
  {
    static Timer t("FromPlane"); RegionTimer reg(t);
    
    if (projecttype == PLANESPACE)
      {
	//      cout << "2d   : " << pplane << endl;
	p3d = p1 + (h * pplane(0)) * ex + (h * pplane(1)) * ey;
	//      cout << "3d   : " << p3d << endl;
	Project (p3d, gi);  
	//      cout << "proj : " << p3d << endl;
      }
    else
      {
	//      Point<2> pspnew = Point<2>(1/h * (Amat * Vec<2>(pplane)) + Vec<2>(psp1));
	Point<2> pspnew = Point<2>(h * (Amat * Vec<2>(pplane)) + Vec<2>(psp1));
	//      Point<2> pspnew = Point<2>(h * (Vec<2>(pplane)) + Vec<2>(psp1));
	gi.u = pspnew(0);
	gi.v = pspnew(1);
	gi.trignum = 1;
	gp_Pnt val = occface->Value (gi.u, gi.v);
	p3d = Point<3> (val.X(), val.Y(), val.Z());
      };
  }



  void OCCSurface :: Project (Point<3> & ap, PointGeomInfo & gi)
  {
    static Timer t("OccSurface::Project"); RegionTimer reg(t);
    static Timer t2("OccSurface::Project actural"); 


    // try Newton's method ...
    
    gp_Pnt p(ap(0), ap(1), ap(2));

    double u = gi.u;
    double v = gi.v;
    gp_Pnt x = occface->Value (u,v);

    if (p.SquareDistance(x) <= sqr(PROJECTION_TOLERANCE)) return;

    gp_Vec du, dv;
    occface->D1(u,v,x,du,dv);

    int count = 0;
    
    gp_Pnt xold;
    gp_Vec n;
    double det, lambda, mu;
    
    do
      {
        count++;
        
        n = du^dv;
        
        det = Det3 (n.X(), du.X(), dv.X(),
                    n.Y(), du.Y(), dv.Y(),
                    n.Z(), du.Z(), dv.Z());
        
        if (det < 1e-15)
          break;
        
        lambda = Det3 (n.X(), p.X()-x.X(), dv.X(),
                       n.Y(), p.Y()-x.Y(), dv.Y(),
                       n.Z(), p.Z()-x.Z(), dv.Z())/det;
        
        mu     = Det3 (n.X(), du.X(), p.X()-x.X(),
                       n.Y(), du.Y(), p.Y()-x.Y(),
                       n.Z(), du.Z(), p.Z()-x.Z())/det;
        
        u += lambda;
        v += mu;
        
        xold = x;
        occface->D1(u,v,x,du,dv);
        
        if (xold.SquareDistance(x) < sqr(PROJECTION_TOLERANCE))
          {
            ap = Point<3> (x.X(), x.Y(), x.Z());
            gi.u = u;
            gi.v = v;
            return;
          }
      }
    while (count < 20);
    

    // Newton did not converge, use OCC projection

    
    //   static int cnt = 0;
    //  if (cnt++ % 1000 == 0) cout << "********************************************** OCCSurfce :: Project, cnt = " << cnt << endl;
  
    gp_Pnt pnt = p;  // (p(0), p(1), p(2));

    //(*testout) << "pnt = " << pnt.X() << ", " << pnt.Y() << ", " << pnt.Z() << endl;


    /*
    GeomAPI_ProjectPointOnSurf proj(pnt, occface, umin, umax, vmin, vmax);

    if (!proj.NbPoints())
      {
	cout << "Project Point on Surface FAIL" << endl;
	throw UVBoundsException();
      }
    */

    



    /*
      cout << "NP = " << proj.NbPoints() << endl;

      for (int i = 1; i <= proj.NbPoints(); i++)
      {
      gp_Pnt pnt2 = proj.Point(i);
      Point<3> p2 = Point<3> (pnt2.X(), pnt2.Y(), pnt2.Z());
      cout << i << ". p = " << p2 << ", dist = " << (p2-p).Length() << endl;
      }
    */

    /*
    pnt = proj.NearestPoint();
    proj.LowerDistanceParameters (gi.u, gi.v);
    */

    // double u,v;
    Handle( ShapeAnalysis_Surface ) su = new ShapeAnalysis_Surface( occface );
    auto toltool =  BRep_Tool::Tolerance( topods_face );

    // gp_Pnt2d suval = su->ValueOfUV ( pnt, toltool);
    t2.Start();
    gp_Pnt2d suval = su->NextValueOfUV (gp_Pnt2d(u,v), pnt, toltool);
    t2.Stop();
    suval.Coord( u, v);
    pnt = occface->Value( u, v );
    
    //(*testout) << "pnt(proj) = " << pnt.X() << ", " << pnt.Y() << ", " << pnt.Z() << endl;
    gi.u = u;
    gi.v = v;
    gi.trignum = 1;

    ap = Point<3> (pnt.X(), pnt.Y(), pnt.Z());
  } 


  Meshing2OCCSurfaces :: Meshing2OCCSurfaces (const NetgenGeometry& geo,
                                              const TopoDS_Shape & asurf,
					      const Box<3> & abb, int aprojecttype,
                                              const MeshingParameters & mparam)
    : Meshing2(geo, mparam, Box<3>(abb.PMin(), abb.PMax())),
      surface(TopoDS::Face(asurf), aprojecttype)
  {
    ;
  }


  void Meshing2OCCSurfaces :: DefineTransformation (const Point<3> & p1, const Point<3> & p2,
						    const PointGeomInfo * geominfo1,
						    const PointGeomInfo * geominfo2)
  {
    ((OCCSurface&)surface).DefineTangentialPlane (p1, *geominfo1, p2, *geominfo2);
  }
 
  void Meshing2OCCSurfaces :: TransformToPlain (const Point<3>& locpoint, 
						const MultiPointGeomInfo & geominfo,
						Point<2> & planepoint, 
						double h, int & zone)
  {
    surface.ToPlane (locpoint, geominfo.GetPGI(1), planepoint, h, zone);
  }

  int Meshing2OCCSurfaces :: TransformFromPlain (const Point<2> & planepoint,
						 Point<3> & locpoint,
						 PointGeomInfo & gi,
						 double h)
  {
    surface.FromPlane (planepoint, locpoint, gi, h);
    return 0;
  }



  double Meshing2OCCSurfaces :: CalcLocalH (const Point<3> & p, double gh) const
  {
    return gh;
  }

  /*
    inline double Det3 (double a00, double a01, double a02,
    double a10, double a11, double a12,
    double a20, double a21, double a22)
    {
    return a00*a11*a22 + a01*a12*a20 + a10*a21*a02 - a20*a11*a02 - a10*a01*a22 - a21*a12*a00;
    }

    bool ProjectToSurface (gp_Pnt & p, Handle(Geom_Surface) surface, double& u, double& v)
    {
    gp_Pnt x = surface->Value (u,v);

    if (p.SquareDistance(x) <= sqr(PROJECTION_TOLERANCE)) return true;

    gp_Vec du, dv;

    surface->D1(u,v,x,du,dv);

    int count = 0;

    gp_Pnt xold;
    gp_Vec n;
    double det, lambda, mu;

    do {
    count++;

    n = du^dv;

    det = Det3 (n.X(), du.X(), dv.X(),
    n.Y(), du.Y(), dv.Y(),
    n.Z(), du.Z(), dv.Z());

    if (det < 1e-15) return false; 

    lambda = Det3 (n.X(), p.X()-x.X(), dv.X(),
    n.Y(), p.Y()-x.Y(), dv.Y(),
    n.Z(), p.Z()-x.Z(), dv.Z())/det;

    mu     = Det3 (n.X(), du.X(), p.X()-x.X(),
    n.Y(), du.Y(), p.Y()-x.Y(),
    n.Z(), du.Z(), p.Z()-x.Z())/det;
  
    u += lambda;
    v += mu;

    xold = x;
    surface->D1(u,v,x,du,dv);

    } while (xold.SquareDistance(x) > sqr(PROJECTION_TOLERANCE) || count > 50);

    if (count > 50) return false;

    p = x;

    return true;
    }
  */
}


#endif

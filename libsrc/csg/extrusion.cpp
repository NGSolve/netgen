#include <mystdlib.h>
#include <core/register_archive.hpp>

#include <linalg.hpp>
#include <csg.hpp>


namespace netgen
{

  NgArray<Point<3> > project1, project2;



  void ExtrusionFace :: Init(void)
  {
    p0.SetSize(path->GetNSplines());
    x_dir.SetSize(path->GetNSplines());
    y_dir.SetSize(path->GetNSplines());
    z_dir.SetSize(path->GetNSplines());
    loc_z_dir.SetSize(path->GetNSplines());
    spline3_path.SetSize(path->GetNSplines());
    line_path.SetSize(path->GetNSplines());

    for(int i=0; i<path->GetNSplines(); i++)
      {
	spline3_path[i] = dynamic_cast < const SplineSeg3<3>* >(&path->GetSpline(i));
	line_path[i] = dynamic_cast < const LineSeg<3>* >(&path->GetSpline(i));
	
	if(line_path[i])
	  {
	    y_dir[i] = line_path[i]->EndPI() - line_path[i]->StartPI();
	    y_dir[i].Normalize();
	    z_dir[i] = glob_z_direction;
	    Orthogonalize(y_dir[i],z_dir[i]);
	    x_dir[i] = Cross(y_dir[i],z_dir[i]);
	    loc_z_dir[i] = z_dir[i];
	  }
	else
	  {
	    z_dir[i] = glob_z_direction;
	    loc_z_dir[i] = glob_z_direction;
	  }
      }

    double cum_angle = 0.;
    for(auto i : Range(path->GetSplines()))
      {
        const auto& sp = path->GetSpline(i);
        auto t1 = sp.GetTangent(0.);
        t1.Normalize();
        auto t2 = sp.GetTangent(1.);
        t2.Normalize();
        cum_angle += acos(t1 * t2);
        angles.Append(cum_angle);
      }
    
    profile->GetCoeff(profile_spline_coeff);
    latest_point3d = -1.111e30;
  }

  
  ExtrusionFace :: ExtrusionFace(const SplineSeg<2> * profile_in,
				 const SplineGeometry<3> * path_in,
				 const Vec<3> & z_direction) :
    profile(profile_in), path(path_in), glob_z_direction(z_direction)
  {
    deletable = false;

    Init();
  }

  ExtrusionFace :: ExtrusionFace(const NgArray<double> & raw_data)
  {
    deletable = true;

    int pos=0;

    NgArray< Point<2> > p(3);

    int ptype = int(raw_data[pos]); pos++;

    for(int i=0; i<ptype; i++)
      {
	p[i](0) = raw_data[pos]; pos++;
	p[i](1) = raw_data[pos]; pos++;
      }
    if(ptype == 2)
      {
	profile = new LineSeg<2>(GeomPoint<2>(p[0],1),
				 GeomPoint<2>(p[1],1));
      }
    else if(ptype == 3)
      {
	profile = new SplineSeg3<2>(GeomPoint<2>(p[0],1),
				    GeomPoint<2>(p[1],1),
				    GeomPoint<2>(p[2],1));
      }

    path = new SplineGeometry<3>;
    pos = const_cast< SplineGeometry<3> *>(path)->Load(raw_data,pos);

    for(int i = 0; i < 3; i++)
      {
	glob_z_direction(i) = raw_data[pos]; 
	pos++;
      }
    
    Init();
  }

  ExtrusionFace :: ~ExtrusionFace()
  {
    if(deletable)
      {
	delete profile;
	delete path;
      }
  }

  
  int ExtrusionFace :: IsIdentic (const Surface & s2, int & inv, double eps) const
  {
    const ExtrusionFace * ext2 = dynamic_cast<const ExtrusionFace*>(&s2);

    if(!ext2) return 0;

    if(ext2 == this)
      return 1;

    return 0;
  } 
  
  void ExtrusionFace :: Orthogonalize(const Vec<3> & v1, Vec<3> & v2) const
  {
    v2 -= (v1*v2)*v1;
    v2.Normalize();
  }


  void ExtrusionFace :: CalcProj(const Point<3> & point3d, Point<2> & point2d,
				 int & seg, double & t) const
  {
    static mutex set_latest_point;

    auto eps = 1e-25 * Dist2(path->GetSpline(0).StartPI(), path->GetSpline(0).EndPI());
    if (Dist2 (point3d, latest_point3d) < eps)
      {
        std::lock_guard<std::mutex> guard(set_latest_point);
        if (Dist2 (point3d, latest_point3d) < eps)
        {
            point2d = latest_point2d;
            seg = latest_seg;
            t = latest_t;
            return;
        }
      }
    

    double cutdist = -1;
    

    NgArray<double> mindist(path->GetNSplines());

    for(int i = 0; i < path->GetNSplines(); i++)
      {
	double auxcut = -1;
	double auxmin = -1;

	if(spline3_path[i])
	  {
	    Point<3> startp(path->GetSpline(i).StartPI());
	    Point<3> endp(path->GetSpline(i).EndPI());
	    Point<3> tanp(spline3_path[i]->TangentPoint());
            
            // lower bound for dist
            auxmin = sqrt (MinDistTP2 (startp, endp, tanp, point3d)); 
            
            // upper bound for dist
            auxcut = min2 (Dist (startp, point3d), Dist (endp, point3d));
	  }
	else if(line_path[i])
	  {
            auxmin = auxcut = sqrt (MinDistLP2 (path->GetSpline(i).StartPI(),
                                                path->GetSpline(i).EndPI(),
                                                point3d));
	  }
	
	mindist[i] = auxmin;
	
	if(i==0 || auxcut < cutdist)
	  cutdist = auxcut;
      }
	


    Point<2> testpoint2d;
    Point<3> testpoint3d;
    
    double minproj(-1);
    bool minproj_set(false);



    for(int i=0; i<path->GetNSplines(); i++)
      {
	if(mindist[i] > cutdist*(1+1e-10)) continue;

	double thist = CalcProj(point3d,testpoint2d,i);

	testpoint3d = p0[i] + testpoint2d(0)*x_dir[i] + testpoint2d(1)*loc_z_dir[i];
	double d = Dist2(point3d,testpoint3d);


	if(!minproj_set || d < minproj)
	  {
	    minproj_set = true;
	    minproj = d;
	    point2d = testpoint2d;
	    t = thist;
	    seg = i;
	  }
      }
    std::lock_guard<std::mutex> guard(set_latest_point);
    latest_seg = seg;
    latest_t = t;
    latest_point2d = point2d;
    latest_point3d = point3d;
  }

  double ExtrusionFace :: CalcProj(const Point<3> & point3d, Point<2> & point2d,
				   int seg) const
  {
    double t = -1;

    if(line_path[seg])
      {
	point2d(0) = (point3d-line_path[seg]->StartPI())*x_dir[seg];
	point2d(1) = (point3d-line_path[seg]->StartPI())*z_dir[seg];
	double l = Dist(line_path[seg]->StartPI(),
			line_path[seg]->EndPI());
	t = min2(max2((point3d - line_path[seg]->StartPI()) * y_dir[seg],0.),
		 l);	
	p0[seg] = line_path[seg]->StartPI() + t*y_dir[seg];
	t *= 1./l;
      }
    else if(spline3_path[seg])
      {
	spline3_path[seg]->Project(point3d,p0[seg],t);
	
	y_dir[seg] = spline3_path[seg]->GetTangent(t); 
        y_dir[seg].Normalize();
	loc_z_dir[seg] = z_dir[seg];
	Orthogonalize(y_dir[seg],loc_z_dir[seg]);
	x_dir[seg] = Cross(y_dir[seg],loc_z_dir[seg]);
	Vec<3> dir = point3d-p0[seg];
	point2d(0) = x_dir[seg]*dir;
	point2d(1) = loc_z_dir[seg]*dir;	
      }
    return t;
  }



  double ExtrusionFace :: CalcFunctionValue (const Point<3> & point) const
  {
    Point<2> p;

    double dummyd;
    int dummyi;

    CalcProj(point, p, dummyi, dummyd);

    return 
      profile_spline_coeff(0)*p(0)*p(0) + 
      profile_spline_coeff(1)*p(1)*p(1) + 
      profile_spline_coeff(2)*p(0)*p(1) + 
      profile_spline_coeff(3)*p(0) + 
      profile_spline_coeff(4)*p(1) + 
      profile_spline_coeff(5);    
  }




  void ExtrusionFace :: CalcGradient (const Point<3> & point, Vec<3> & grad) const
  {
    Point<2> p2d;

    double t_path;
    int seg;
    CalcProj (point, p2d, seg, t_path);


    Point<3> phi;
    Vec<3> phip, phipp, phi_minus_point;

    path->GetSpline(seg).GetDerivatives(t_path, phi, phip, phipp);
    phi_minus_point = phi-point;

    Vec<3> grad_t = (1.0/(phipp*phi_minus_point + phip*phip)) * phip;
    Vec<3> grad_xbar, grad_ybar;

    Vec<3> hex, hey, hez, dex, dey, dez;
    CalcLocalCoordinatesDeriv (seg, t_path, hex, hey, hez, dex, dey, dez);

    grad_xbar = hex - (phi_minus_point*dex + hex*phip) * grad_t;
    grad_ybar = hez - (phi_minus_point*dez + hez*phip) * grad_t;

    double dFdxbar = 2.*profile_spline_coeff(0)*p2d(0) +
      profile_spline_coeff(2)*p2d(1) + profile_spline_coeff(3);

    double dFdybar = 2.*profile_spline_coeff(1)*p2d(1) +
      profile_spline_coeff(2)*p2d(0) + profile_spline_coeff(4);
    

    grad = dFdxbar * grad_xbar + dFdybar * grad_ybar;    
  }

  void ExtrusionFace :: CalcHesse (const Point<3> & point, Mat<3> & hesse) const
  {
    const double eps = 1e-7*Dist(path->GetSpline(0).StartPI(),path->GetSpline(0).EndPI());
    

    Point<3> auxpoint1(point),auxpoint2(point);
    Vec<3> auxvec,auxgrad1,auxgrad2;

    for(int i=0; i<3; i++)
      {
	auxpoint1(i) -= eps;
	auxpoint2(i) += eps;
	CalcGradient(auxpoint1,auxgrad1);
	CalcGradient(auxpoint2,auxgrad2);
	auxvec = (1./(2.*eps)) * (auxgrad2-auxgrad1);
	for(int j=0; j<3; j++)
	  hesse(i,j) = auxvec(j);
	auxpoint1(i) = point(i);
	auxpoint2(i) = point(i);
      }

    /*
    Vec<3> grad;
    CalcGradient(point,grad);

    Point<3> auxpoint(point);
    Vec<3> auxvec,auxgrad;

    for(int i=0; i<3; i++)
      {
	auxpoint(i) -= eps;
	CalcGradient(auxpoint,auxgrad);
	auxvec = (1./eps) * (grad-auxgrad);
	for(int j=0; j<3; j++)
	  hesse(i,j) = auxvec(j);
	auxpoint(i) = point(i);
      }
    */

    
    for(int i=0; i<3; i++)
      for(int j=i+1; j<3; j++)
	hesse(i,j) = hesse(j,i) = 0.5*(hesse(i,j)+hesse(j,i));
  }
  


  double ExtrusionFace :: HesseNorm () const
  {
    return fabs(profile_spline_coeff(0) + profile_spline_coeff(1)) +
      sqrt(pow(profile_spline_coeff(0)+profile_spline_coeff(1),2)+4.*pow(profile_spline_coeff(2),2));
  }

  double ExtrusionFace :: MaxCurvature () const
  {
    double retval,actmax;
    
    retval = profile->MaxCurvature();
    for(int i=0; i<path->GetNSplines(); i++)
      {
	actmax = path->GetSpline(i).MaxCurvature();
	if(actmax > retval)
	  retval = actmax;
      }

    return 2.*retval;
  }


  void ExtrusionFace :: Project (Point<3> & p) const
  {
    double dummyt;
    int seg;
    Point<2> p2d;

    CalcProj(p,p2d,seg,dummyt);

    profile->Project(p2d,p2d,profile_par);
    
    p = p0[seg] + p2d(0)*x_dir[seg] + p2d(1)*loc_z_dir[seg];

    Vec<2> tangent2d = profile->GetTangent(profile_par);
    profile_tangent = tangent2d(0)*x_dir[seg] + tangent2d(1)*y_dir[seg];
  }


  
  Point<3> ExtrusionFace :: GetSurfacePoint () const
  {
    p0[0] = path->GetSpline(0).GetPoint(0.5);
    if(!line_path[0])
      {
	y_dir[0] = path->GetSpline(0).GetTangent(0.5);
	y_dir[0].Normalize();
	loc_z_dir[0] = z_dir[0];
	Orthogonalize(y_dir[0],loc_z_dir[0]);
	x_dir[0] = Cross(y_dir[0],loc_z_dir[0]);
      }

    Point<2> locpoint = profile->GetPoint(0.5);

    return p0[0] + locpoint(0)*x_dir[0] + locpoint(1)*loc_z_dir[0];
  }
  

  bool ExtrusionFace :: BoxIntersectsFace(const Box<3> & box) const
  {
    Point<3> center = box.Center();

    Project(center);

    //(*testout) << "box.Center() " << box.Center() << " projected " << center << " diam " << box.Diam() 
    //       << " dist " << Dist(box.Center(),center) << endl;

    return (Dist(box.Center(),center) < 0.5*box.Diam());
  }


  bool ExtrusionFace :: PointInFace (const Point<3> & p, const double eps) const
  {
    Point<3> hp = p;
    Project(hp);
    return Dist2(p,hp) < sqr(eps);
  }
  

  void ExtrusionFace :: LineIntersections ( const Point<3> & p,
					    const Vec<3> & v,
					    const double eps,
					    int & before,
					    int & after,
					    bool & intersecting ) const
  {
    Point<2> p2d;
    Vec<2> v2d;

    intersecting = false;

    double segt;
    int seg;

    CalcProj(p,p2d,seg,segt);

    if(seg == 0 && segt < 1e-20)
      {
	Vec<3> v1,v2;
	v1 = path->GetSpline(0).GetTangent(0);
	v2 = p-p0[seg];
	if(v1*v2 < -eps)
	  return;
      }
    if(seg == path->GetNSplines()-1 && 1.-segt < 1e-20)
      {
	Vec<3> v1,v2;
	v1 = path->GetSpline(seg).GetTangent(1);
	v2 = p-p0[seg];
	if(v1*v2 > eps)
	  return;
      }

    v2d(0) = v * x_dir[seg];
    v2d(1) = v * loc_z_dir[seg];
    
    Vec<2> n(v2d(1),-v2d(0));
    NgArray < Point<2> > ips;


    profile->LineIntersections(v2d(1),
			      -v2d(0),
			      -v2d(1)*p2d(0) + v2d(0)*p2d(1),
			      ips,eps);
    int comp;

    if(fabs(v2d(0)) >= fabs(v2d(1)))
      comp = 0;
    else
      comp = 1;

    //(*testout) << "p2d " << p2d;

    for(int i=0; i<ips.Size(); i++)
      {
	//(*testout) << " ip " << ips[i];

	double t = (ips[i](comp)-p2d(comp))/v2d(comp);

	if(t < -eps)
	  before++;
	else if(t > eps)
	  after++;
	else
	  intersecting = true;
      }
    //(*testout) << endl;
  }

  void ExtrusionFace :: Print (ostream & str) const{}

  INSOLID_TYPE ExtrusionFace :: VecInFace ( const Point<3> & p,
					    const Vec<3> & v,
					    const double eps ) const
  {
    
    Vec<3> normal1;
    CalcGradient(p,normal1); normal1.Normalize();

    double d1 = normal1*v;


    if(d1 > eps)
      return IS_OUTSIDE;
    if(d1 < -eps)
      return IS_INSIDE;
    

    return DOES_INTERSECT;

    /*
    Point<2> p2d;

    double t_path;
    int seg;
    CalcProj(p,p2d,seg,t_path);

    double t;
    profile.Project(p2d,p2d,t);


    
    Vec<2> profile_tangent = profile.GetTangent(t);

    double d;

    Vec<3> normal1;
    CalcGradient(p,normal1); normal1.Normalize();

    double d1 = normal1*v;

    Vec<2> v2d;

    v2d(0) = v*x_dir[seg];
    v2d(1) = v*loc_z_dir[seg];

			    	    
    Vec<2> normal(-profile_tangent(1),profile_tangent(0));
    
    //d = normal*v2d;
    

    d = d1;


    if(d > eps)
      return IS_OUTSIDE;
    if(d < -eps)
      return IS_INSIDE;
    

    return DOES_INTERSECT;
    */
  }


  void ExtrusionFace :: GetTriangleApproximation (TriangleApproximation & tas, 
						  const Box<3> & boundingbox, 
						  double facets) const
  {
    int n = int(facets) + 1;

    for(int k = 0; k < path -> GetNSplines(); k++)
      {
	for(int i = 0; i <= n; i++)
	  {
	    Point<3> origin = path -> GetSpline(k).GetPoint(double(i)/double(n));
	    if(!line_path[k])
	      {
		y_dir[k] = path->GetSpline(k).GetTangent(double(i)/double(n));
		y_dir[k].Normalize();
	      }
	    loc_z_dir[k] = z_dir[k];
	    Orthogonalize(y_dir[k],loc_z_dir[k]);
	    if(!line_path[k])
	      x_dir[k] = Cross(y_dir[k],loc_z_dir[k]);
	    
	    for(int j = 0; j <= n; j++)
	      {
		Point<2> locp = profile->GetPoint(double(j)/double(n));
		tas.AddPoint(origin + locp(0)*x_dir[k] + locp(1)*loc_z_dir[k]);
	      }
	  }
      }
    
    for(int k = 0; k < path->GetNSplines(); k++)
      for(int i = 0; i < n; i++)
	for(int j = 0; j < n; j++)
	  {
	    int pi = k*(n+1)*(n+1) + (n+1)*i +j;
	  
	    tas.AddTriangle( TATriangle (0, pi,pi+1,pi+n+1) );
	    tas.AddTriangle( TATriangle (0, pi+1,pi+n+1,pi+n+2) );
	  }
  }
  

  void ExtrusionFace :: GetRawData(NgArray<double> & data) const
  {
    data.DeleteAll();
    profile->GetRawData(data);
    path->GetRawData(data);
    for(int i=0; i<3; i++)
      data.Append(glob_z_direction[i]);
  }


  void ExtrusionFace :: 
  CalcLocalCoordinates (int seg, double t, 
                        Vec<3> & ex, Vec<3> & ey, Vec<3> & ez) const
  {
    ey = path->GetSpline(seg).GetTangent(t); 
    ey /= ey.Length();
    ex = Cross (ey, glob_z_direction);
    ex /= ex.Length();
    ez = Cross (ex, ey);
  }

  void ExtrusionFace :: 
  CalcLocalCoordinatesDeriv (int seg, double t, 
                             Vec<3> & ex, Vec<3> & ey, Vec<3> & ez,
                             Vec<3> & dex, Vec<3> & dey, Vec<3> & dez) const
  {
    Point<3> point;
    Vec<3> first, second;
    path->GetSpline(seg).GetDerivatives (t, point, first, second);

    ey = first;
    ex = Cross (ey, glob_z_direction);
    ez = Cross (ex, ey);
    
    dey = second;
    dex = Cross (dey, glob_z_direction);
    dez = Cross (dex, ey) + Cross (ex, dey);
    
    double lenx = ex.Length();
    double leny = ey.Length();
    double lenz = ez.Length();

    ex /= lenx;
    ey /= leny;
    ez /= lenz;
    
    dex /= lenx;
    dex -= (dex * ex) * ex;

    dey /= leny;
    dey -= (dey * ey) * ey;

    dez /= lenz;
    dez -= (dez * ez) * ez;
  }

  void ExtrusionFace :: DefineTangentialPlane(const Point<3>& ap1,
                                              const Point<3>& ap2)
  {
    Surface::DefineTangentialPlane(ap1, ap2);
    tangential_plane_seg = latest_seg;
  }

  void ExtrusionFace :: ToPlane(const Point<3>& p3d, Point<2>& p2d,
                                double h, int& zone) const
  {
    Surface::ToPlane(p3d, p2d, h, zone);
    double angle = angles[tangential_plane_seg] - angles[latest_seg];
    if(fabs(angle) > 3.14/2.)
      zone = -1;
  }

  Extrusion :: Extrusion(shared_ptr<SplineGeometry<3>> path_in,
			 shared_ptr<SplineGeometry<2>> profile_in,
			 const Vec<3> & z_dir) :
    path(path_in), profile(profile_in), z_direction(z_dir)
  {
    surfaceactive.SetSize(0);
    surfaceids.SetSize(0);

    for(int j=0; j<profile->GetNSplines(); j++)
      {
	ExtrusionFace * face = new ExtrusionFace(&(profile->GetSpline(j)),
						 path.get(),
						 z_direction);
	faces.Append(face);
	surfaceactive.Append(true);
	surfaceids.Append(0);
      }

  }


  Extrusion :: ~Extrusion()
  {
    for(int i=0; i<faces.Size(); i++)
      delete faces[i];
  }





  INSOLID_TYPE Extrusion :: BoxInSolid (const BoxSphere<3> & box) const
  {
    for(int i=0; i<faces.Size(); i++)
      {
	if(faces[i]->BoxIntersectsFace(box))
	  return DOES_INTERSECT;
      }

    return PointInSolid(box.Center(),0);
  }


  INSOLID_TYPE Extrusion :: PointInSolid (const Point<3> & p,
					  const double eps,
					  NgArray<int> * const facenums) const
  {
    Vec<3> random_vec(-0.4561,0.7382,0.4970247);

    int before(0), after(0);
    bool intersects(false);
    bool does_intersect(false);

    for(int i=0; i<faces.Size(); i++)
      {
	faces[i]->LineIntersections(p,random_vec,eps,before,after,intersects);

	//(*testout) << "intersects " << intersects << " before " << before << " after " << after << endl;
	if(intersects)
	  {
	    if(facenums)
	      {
		facenums->Append(i);
		does_intersect = true;
	      }
	    else
	      return DOES_INTERSECT;
	  }
      }

    if(does_intersect)
      return DOES_INTERSECT;


    if(before % 2 == 0)
      return IS_OUTSIDE;

    return IS_INSIDE;
  }


  INSOLID_TYPE Extrusion :: PointInSolid (const Point<3> & p,
					  double eps) const
  {
    return PointInSolid(p,eps,NULL);    
  }

  void Extrusion :: GetTangentialSurfaceIndices (const Point<3> & p, 
                                                 NgArray<int> & surfind, double eps) const
  {
    for (int j = 0; j < faces.Size(); j++)
      if (faces[j] -> PointInFace(p, eps))
        if (!surfind.Contains (GetSurfaceId(j)))
          surfind.Append (GetSurfaceId(j));
  }

  
  INSOLID_TYPE Extrusion :: VecInSolid (const Point<3> & p,
					const Vec<3> & v,
					double eps) const
  {
    NgArray<int> facenums;
    INSOLID_TYPE pInSolid = PointInSolid(p,eps,&facenums);

    if(pInSolid != DOES_INTERSECT)
      return pInSolid;


    double d(0);

    if(facenums.Size() == 1)
      {
	Vec<3> normal;
	faces[facenums[0]]->CalcGradient(p,normal);
	normal.Normalize();
	d = normal*v;
	
	latestfacenum = facenums[0];
      }
    else if (facenums.Size() == 2)
      {
	Vec<3> checkvec;

	Point<3> dummy(p);
	faces[facenums[0]]->Project(dummy);
	if(fabs(faces[facenums[0]]->GetProfilePar()) < 0.1)
	  {
	    int aux = facenums[0];
	    facenums[0] = facenums[1]; facenums[1] = aux;
	  }
	
	checkvec = faces[facenums[0]]->GetYDir();
     
	Vec<3> n0, n1;
	faces[facenums[0]]->CalcGradient(p,n0);
	faces[facenums[1]]->CalcGradient(p,n1);
	n0.Normalize();
	n1.Normalize();
	

	Vec<3> t = Cross(n0,n1);
	if(checkvec*t < 0) t*= (-1.);
	
	Vec<3> t0 = Cross(n0,t);
	Vec<3> t1 = Cross(t,n1);
	
	t0.Normalize();
	t1.Normalize();
	

	const double t0v = t0*v;
	const double t1v = t1*v;

	if(t0v > t1v)
	  {
	    latestfacenum = facenums[0];
	    d = n0*v;
	  }
	else
	  {
	    latestfacenum = facenums[1];
	    d = n1*v;
	  }

	if(fabs(t0v) < eps && fabs(t1v) < eps)
	  latestfacenum = -1;
      }

    else
      {
	cerr << "WHY ARE THERE " << facenums.Size() << " FACES?" << endl;
      }

    if(d > eps)
      return IS_OUTSIDE;
    if(d < -eps)
      return IS_INSIDE;
      
    return DOES_INTERSECT;
  }



  // checks if lim s->0 lim t->0  p + t(v1 + s v2) in solid
  INSOLID_TYPE Extrusion :: VecInSolid2 (const Point<3> & p,
					 const Vec<3> & v1,
					 const Vec<3> & v2,
					 double eps) const
  {
    INSOLID_TYPE retval;
    retval = VecInSolid(p,v1,eps);

    // *testout << "extr, vecinsolid=" << int(retval) << endl;

    if(retval != DOES_INTERSECT)
      return retval;

    if(latestfacenum >= 0)
      return faces[latestfacenum]->VecInFace(p,v2,eps);
    else
      return VecInSolid(p,v2,eps);
  }

  
  int Extrusion :: GetNSurfaces() const
  {
    return faces.Size();
  }

  Surface & Extrusion :: GetSurface (int i)
  {
    return *faces[i];
  }

  const Surface & Extrusion :: GetSurface (int i) const
  {
    return *faces[i];
  }


  void Extrusion :: Reduce (const BoxSphere<3> & box)
  {
    for(int i = 0; i < faces.Size(); i++)
      surfaceactive[i] = faces[i]->BoxIntersectsFace(box);
  }

  void Extrusion :: UnReduce ()
  {
    for(int i = 0; i < faces.Size(); i++)
      surfaceactive[i] = true;
  }

  RegisterClassForArchive<ExtrusionFace, Surface> regexf;
  RegisterClassForArchive<Extrusion, Primitive> regextr;
}

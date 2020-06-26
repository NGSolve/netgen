/* *************************************************************************/
/* File:   surfacegeom.cpp                                                 */
/* Author: Michael Neunteufel                                              */
/* Date:   Jun. 2020                                                       */
/* *************************************************************************/

#include <meshing.hpp>

namespace netgen
{
  SurfaceGeometry :: SurfaceGeometry()
  {
    //identity
    func = [](Point<2> p) { return Vec<3>(p[0],p[1],0.0); };
  }

  SurfaceGeometry :: SurfaceGeometry(function<Vec<3>(Point<2>)> _func) : func(_func)
  {
    ;
  }

  SurfaceGeometry :: SurfaceGeometry(const SurfaceGeometry& geom) : func(geom.func), eps(geom.eps)
  {
    ;
  }
  
  void SurfaceGeometry :: CalcHesse(double u, double v, Vec<3>& f_uu, Vec<3>& f_vv, Vec<3>& f_uv) const
  {
    Point<2> p = Point<2>(u,v);
    double pr  = p[0]+eps;
    double pl  = p[0]-eps;
    double prr = p[0]+2*eps;
    double pll = p[0]-2*eps;
    
    auto dr = GetTangentVectors( pr, v );
    auto dl = GetTangentVectors( pl, v );
    auto drr = GetTangentVectors( prr, v );
    auto dll = GetTangentVectors( pll, v );
    
    f_uu = (1.0/(12.0*eps)) * (8.0*dr[0]-8.0*dl[0]-drr[0]+dll[0]);
    f_uv = (1.0/(12.0*eps)) * (8.0*dr[1]-8.0*dl[1]-drr[1]+dll[1]);

    pr  = p[1]+eps;
    pl  = p[1]-eps;
    prr = p[1]+2*eps;
    pll = p[1]-2*eps;

    dr = GetTangentVectors(u, pr);
    dl = GetTangentVectors(u, pl);
    drr = GetTangentVectors(u, prr);
    dll = GetTangentVectors(u, pll);
    
    f_vv = (1.0/(12.0*eps)) * (8.0*dr[1]-8.0*dl[1]-drr[1]+dll[1]);
  }
  
  Array<Vec<3>> SurfaceGeometry :: GetTangentVectors(double u, double v) const
  {
    Array<Vec<3>> tang(2);
    
    Point<2> pru  = Point<2>(u+eps,v);
    Point<2> plu  = Point<2>(u-eps,v);
    Point<2> prru = Point<2>(u+2*eps,v);
    Point<2> pllu = Point<2>(u-2*eps,v);
    
    Point<2> prv  = Point<2>(u,v+eps);
    Point<2> plv  = Point<2>(u,v-eps);
    Point<2> prrv = Point<2>(u,v+2*eps);
    Point<2> pllv = Point<2>(u,v-2*eps);

    
    tang[0] = 1/(12.0*eps)*( 8.0*func(pru) - 8.0*func(plu) - func(prru) + func(pllu) );
    tang[1] = 1/(12.0*eps)*( 8.0*func(prv) - 8.0*func(plv) - func(prrv) + func(pllv) );
    
    return tang;
  }

  Vec<3> SurfaceGeometry :: GetNormal(int surfind, const Point<3> & p, const PointGeomInfo* gi) const
  {
    Array<Vec<3>> tang = GetTangentVectors(gi->u, gi->v);
    auto normal = Cross(tang[0], tang[1]);
    return Cross(tang[0], tang[1]);
  }


  PointGeomInfo SurfaceGeometry :: ProjectPoint(int surfind, Point<3> & p) const
  {
    throw Exception("In SurfaceGeometry::ProjectPoint");
  }
  
  void SurfaceGeometry :: ProjectPointEdge (int surfind, int surfind2, Point<3> & p,
                                            EdgePointGeomInfo* gi) const
  {
    if (gi == nullptr)
      throw Exception("In SurfaceGeometry::ProjectPointEdge: gi is nullptr");
    throw Exception("In SurfaceGeometry::ProjectPointEdge: not implemented");
  }
    
  bool  SurfaceGeometry :: ProjectPointGI (int surfind, Point<3> & p, PointGeomInfo & gi) const
  {
    Array<Vec<3>> tangs;
    Vec<3> diff, f_uu, f_vv, f_uv;
    Vec<2> r, dx;
    double norm_r, det, energy=0.0, new_energy=0.0, alpha=2.0,u=0.0,v=0.0;
    Mat<2,2> mat, inv;
    int num=0, maxit=20;
    double damping=0.2;


    //Solve minimization problem
    //   argmin_(u,v) 0.5*\| f(u,v)-p\|^2
    //via Neton's method:
    //  F(u,v) = ( (f(u,v)-p)*f_u(u,v), (f(u,v)-p)*f_v(u,v))^T = (0,0)^T
    //Stiffness matrix
    //  F'(u,v) = ( f_u*f_u + (f-p)*f_uu, f_v*f_u + (f-p)*f_uv, f_v*f_u + (f-p)*f_uv, f_v*f_v + (f-p)*f_vv )
    do
      {
        num++;
        tangs = GetTangentVectors(gi.u, gi.v);
        diff = func(Point<2>(gi.u, gi.v)) - Vec<3>(p);
        energy = diff.Length2();
        r = Vec<2>( diff*tangs[0], diff*tangs[1] );
        norm_r = r.Length2();

        CalcHesse(gi.u, gi.v, f_uu, f_vv, f_uv);
        
       
        mat(0,0) = tangs[0]*tangs[0] + diff*f_uu;
        mat(1,0) = mat(0,1) = tangs[0]*tangs[1]+diff*f_uv;
        mat(1,1) = tangs[1]*tangs[1]+diff*f_vv;
    
        CalcInverse(mat,inv);

        dx = inv*r;

        //Linesearch 
        alpha = 2.0;
        do
          {
            alpha /= 2.0;
            u = gi.u - min(1.0,alpha*damping*num)*dx[0];
            v = gi.v - min(1.0,alpha*damping*num)*dx[1];

            diff = func(Point<2>(u, v)) - Vec<3>(p);
            new_energy = diff.Length2();
          }
        while (alpha > 1e-10 && new_energy > energy+1e-14);
        if (alpha <= 1e-10)
          throw Exception("In SurfaceGeometry::ProjectPointGI: Linesearch min alpha reached!");
        gi.u = u;
        gi.v = v;
        

      }
    while ( norm_r > 1e-12 && num < maxit);

    //Stay in reference domain [0,1]^2
    if (gi.u < 0 || gi.u > 1 || gi.v < 0 || gi.v > 1)
      {
        cout << "Warning: Projected point outside [0,1]^2: u=" << gi.u << ",v=" << gi.v <<". Setting back." << endl;
        gi.u = min(max(gi.u,0.0),1.0);
        gi.v = min(max(gi.v,0.0),1.0);
      }

    p = Point<3>(func(Point<2>(gi.u,gi.v)));

    if (num == maxit)
      {
        //cout << "In SurfaceGeometry::ProjectPointGI: Newton did not converge" << endl;
        throw Exception("In SurfaceGeometry::ProjectPointGI: Newton did not converge");
      }
    return true;
  }

  bool  SurfaceGeometry :: CalcPointGeomInfo(int surfind, PointGeomInfo& gi, const Point<3> & p3) const
  {
    throw Exception("In SurfaceGeometry::CalcPointGeomInfo: not implemented");
    return false;
  }

  void  SurfaceGeometry :: PointBetweenEdge(const Point<3> & p1, const Point<3> & p2, double secpoint, int surfi1, int surfi2, const EdgePointGeomInfo & ap1, const EdgePointGeomInfo & ap2, Point<3> & newp, EdgePointGeomInfo & newgi) const
  {
    newp = p1+secpoint*(p2-p1);

    PointGeomInfo pgi;
    pgi.u = ap1.u+secpoint*(ap2.u-ap1.u);
    pgi.v = ap1.v+secpoint*(ap2.v-ap1.v);

    ProjectPointGI(surfi1, newp, pgi);

    newgi.u = pgi.u;
    newgi.v = pgi.v;
    newgi.edgenr = ap1.edgenr;
    newgi.body = -1;
    newgi.dist = -1.0;
  }
  
  void SurfaceGeometry :: PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint,
                                       int surfi, 
                                       const PointGeomInfo & gi1, 
                                       const PointGeomInfo & gi2,
                                       Point<3> & newp, PointGeomInfo & newgi) const
  {
    newp = p1+secpoint*(p2-p1);
    
    newgi.u = gi1.u+secpoint*(gi2.u-gi1.u);
    newgi.v = gi1.v+secpoint*(gi2.v-gi1.v);  
    newgi.trignum = -1;

    ProjectPointGI(surfi, newp, newgi);
  }

  int SurfaceGeometry :: GenerateMesh(shared_ptr<Mesh> & mesh, bool quads, int nx, int ny, bool flip_triangles, const Array<Point<3>>& bbbpts, const Array<string>& bbbnames)
  {
    mesh->SetDimension(3);

    Array<bool> found(bbbpts.Size());
    found = false;
    Array<PointIndex> indbbbpts(bbbpts.Size());

    
    Array<PointIndex> pids;
    Array<PointGeomInfo> pgis;
    for(int i=0; i <= ny; i++)
      for(int j=0; j <= nx; j++)
        {
          PointGeomInfo pgi;
          pgi.trignum = -1;
          pgi.u = double(j)/nx;
          pgi.v = double(i)/ny;

          Point<3> pnt = Point<3>(func(Point<2>(pgi.u,pgi.v)));
          pids.Append(mesh->AddPoint(pnt));
          pgis.Append(pgi);
          
          for (int k = 0; k < bbbpts.Size(); k++)
            {
              auto diff = pnt - bbbpts[k];
              if(diff.Length2() < 1e-14)
                {
                  found[k] = true;
                  indbbbpts[k] = pids[pids.Size()-1];
                }
            }
        }

    for (bool f : found)
      if (!f)
        throw Exception("In SurfaceGeometry :: GenerateMesh: bbbpts not resolved in mesh.");

    FaceDescriptor fd;
    fd.SetSurfNr(1);
    fd.SetDomainIn(1);
    fd.SetDomainOut(0);
    fd.SetBCProperty(1);
    mesh->AddFaceDescriptor(fd);


    for(int i=0; i < ny; i++)
      {
        for(int j=0; j < nx; j++)
          {
            int base = i * (nx+1) + j;
            if (quads)
              {
                int pnum[4] = {base,base+1,base+nx+2,base+nx+1};
                Element2d el = Element2d(QUAD);
                for (int i = 0; i < 4; i++)
                  {
                    el[i] = pids[pnum[i]];
                    el.GeomInfoPi(i+1) = pgis[pnum[i]];
                  }
                el.SetIndex(1);
            
                mesh->AddSurfaceElement(el);
              }
            else
              {
                Array<int> pnum1(3);
                Array<int> pnum2(3);
                if (flip_triangles)
                  {
                    pnum1[0] = base;
                    pnum1[1] = base+1;
                    pnum1[2] = base+nx+2;
                    pnum2[0] = base;
                    pnum2[1] = base+nx+2;
                    pnum2[2] = base+nx+1;
                  }
                else
                  {
                    pnum1[0] = base;
                    pnum1[1] = base+1;
                    pnum1[2] = base+nx+1;
                    pnum2[0] = base+1;
                    pnum2[1] = base+nx+2;
                    pnum2[2] = base+nx+1;
                  }

                Element2d el = Element2d(TRIG);
                for (int i = 0; i < 3; i++)
                  {
                    el[i] = pids[pnum1[i]];
                    el.GeomInfoPi(i+1) = pgis[pnum1[i]];
                  }
                el.SetIndex(1);
            
                mesh->AddSurfaceElement(el);
                for (int i = 0; i < 3; i++)
                  {
                    el[i] = pids[pnum2[i]];
                    el.GeomInfoPi(i+1) = pgis[pnum2[i]];
                  }
                mesh->AddSurfaceElement(el);
              }
          }
      }

    Segment seg;
    seg.si = 1;
    seg.edgenr = 1;
    seg.epgeominfo[0].edgenr = 1;
    seg.epgeominfo[1].edgenr = 1;
    // needed for codim2 in 3d
    seg.edgenr = 1;
    for(int i=0; i < nx; i++)
      {
        seg[0] = pids[i];
        seg[1] = pids[i+1];
        
        seg.geominfo[0] = pgis[i];
        seg.geominfo[1] = pgis[i+1];
        seg.epgeominfo[0].u = pgis[i].u;
        seg.epgeominfo[0].v = pgis[i].v;
        seg.epgeominfo[0].edgenr = seg.edgenr;
        seg.epgeominfo[1].u = pgis[i+1].u;
        seg.epgeominfo[1].v = pgis[i+1].v;
        seg.epgeominfo[1].edgenr = seg.edgenr;
        
        mesh->AddSegment(seg);
      }

    seg.si = 2;
    seg.edgenr = 2;
    for(int i=0; i<ny; i++)
      {
        seg[0] = pids[i*(nx+1)+nx];
        seg[1] = pids[(i+1)*(nx+1)+nx];

        seg.geominfo[0] = pgis[i*(nx+1)+nx];
        seg.geominfo[1] = pgis[(i+1)*(nx+1)+nx];
        seg.epgeominfo[0].u = pgis[i*(nx+1)+nx].u;
        seg.epgeominfo[0].v = pgis[i*(nx+1)+nx].v;
        seg.epgeominfo[0].edgenr = seg.edgenr;
        seg.epgeominfo[1].u = pgis[(i+1)*(nx+1)+nx].u;
        seg.epgeominfo[1].v = pgis[(i+1)*(nx+1)+nx].v;
        seg.epgeominfo[1].edgenr = seg.edgenr;

        mesh->AddSegment(seg);
      }

    seg.si = 3;
    seg.edgenr = 3;
    for(int i=0; i<nx; i++)
      {
        seg[0] = pids[ny*(nx+1)+i+1];
        seg[1] = pids[ny*(nx+1)+i];

        seg.geominfo[0] = pgis[ny*(nx+1)+i+1];
        seg.geominfo[1] = pgis[ny*(nx+1)+i];
        seg.epgeominfo[0].u = pgis[ny*(nx+1)+i+1].u;
        seg.epgeominfo[0].v = pgis[ny*(nx+1)+i+1].v;
        seg.epgeominfo[0].edgenr = seg.edgenr;
        seg.epgeominfo[1].u = pgis[ny*(nx+1)+i].u;
        seg.epgeominfo[1].v = pgis[ny*(nx+1)+i].v;
        seg.epgeominfo[1].edgenr = seg.edgenr;
        
        mesh->AddSegment(seg);
      }

    seg.si = 4;
    seg.edgenr = 4;
    for(int i=0; i<ny; i++)
      {
        seg[0] = pids[(i+1)*(nx+1)];
        seg[1] = pids[i*(nx+1)];

        seg.geominfo[0] = pgis[(i+1)*(nx+1)];
        seg.geominfo[1] = pgis[i*(nx+1)];
        seg.epgeominfo[0].u = pgis[(i+1)*(nx+1)].u;
        seg.epgeominfo[0].v = pgis[(i+1)*(nx+1)].v;
        seg.epgeominfo[0].edgenr = seg.edgenr;
        seg.epgeominfo[1].u = pgis[i*(nx+1)].u;
        seg.epgeominfo[1].v = pgis[i*(nx+1)].v;
        seg.epgeominfo[1].edgenr = seg.edgenr;

        mesh->AddSegment(seg);
      }

    mesh->SetCD2Name(1, "bottom");
    mesh->SetCD2Name(2, "right");
    mesh->SetCD2Name(3, "top");
    mesh->SetCD2Name(4, "left");

    for (int i = 0; i < bbbpts.Size(); i++)
      {
        Element0d el;
        el.pnum = indbbbpts[i];
        el.index = i+1;
        mesh->pointelements.Append(el);
        mesh->SetCD3Name(i+1, bbbnames[i]);
      }

    mesh->Compress();
    mesh->UpdateTopology();

    return 0;
  }
  
};

#include <mystdlib.h>

#include "meshing.hpp"
#include <csg.hpp>
#include <geometry2d.hpp>

namespace netgen
{

  int CheckSurfaceMesh (const Mesh & mesh)
  {
    PrintMessage (3, "Check Surface mesh");

    int nf = mesh.GetNSE();
    ClosedHashTable<PointIndices<2>, int> edges(nf+2);
    int cnt1 = 0, cnt2 = 0;

    for (SurfaceElementIndex i : T_Range<SurfaceElementIndex>(nf))
      for (int j = 1; j <= 3; j++)
        {
          PointIndex pi1 = mesh[i].PNumMod(j);
          PointIndex pi2 = mesh[i].PNumMod(j+1);
          if (edges.Used ( { pi1, pi2 } ))
            {
              int hi = edges.Get ( { pi1, pi2 } );
              if (hi != 1) 
                PrintSysError ("CheckSurfaceMesh, hi = ", hi);
              edges.Set ( { pi1, pi2 }, 2);
              cnt2++;
            }
          else
            {
              edges.Set ( { pi2, pi1 }, 1);
              cnt1++;
            }
        }
  

    if (cnt1 != cnt2)
      {
        PrintUserError ("Surface mesh not consistent");
        //      MyBeep(2);
        //      (*mycout) << "cnt1 = " << cnt1 << " cnt2 = " << cnt2 << endl;
        return 0;
      }
    return 1;
  }



  int CheckSurfaceMesh2 (const Mesh & mesh)
  {
    const Point<3> *tri1[3], *tri2[3];

    for (int i = 1; i <= mesh.GetNOpenElements(); i++)
      {
        PrintDot ();
        for (int j = 1; j < i; j++)
          {
            for (int k = 0; k < 3; k++)
              {
                tri1[k] = &mesh.Point (mesh.OpenElement(i)[k]);
                tri2[k] = &mesh.Point (mesh.OpenElement(j)[k]);
              }
            if (IntersectTriangleTriangle (&tri1[0], &tri2[0]))
              {
                PrintSysError ("Surface elements are intersecting");
                (*testout) << "Intersecting: " << endl;
                for (int k = 0; k <= 2; k++)
                  (*testout) << *tri1[k] << "   ";
                (*testout) << endl;
                for (int k = 0; k <= 2; k++)
                  (*testout) << *tri2[k] << "   ";
                (*testout) << endl;
              }

          }
      }
    return 0;
  }





  static double TriangleQualityInst (const Point<3> & p1, const Point<3> & p2,
                                     const Point<3> & p3)
  {
    // quality 0 (worst) .. 1 (optimal)

    Vec<3> v1, v2, v3;
    double s1, s2, s3;
    double an1, an2, an3;

    v1 = p2 - p1;
    v2 = p3 - p1;
    v3 = p3 - p2;

    an1 = Angle (v1, v2);
    v1 *= -1;
    an2 = Angle (v1, v3);
    an3 = Angle (v2, v3);

    s1 = sin (an1/2);
    s2 = sin (an2/2);
    s3 = sin (an3/2);

    return 8 * s1 * s2 * s3;
  }














  void MeshQuality2d (const Mesh & mesh)
  {
    int ncl = 20;
    Array<int> incl(ncl);

    incl = 0;

    for (auto & el : mesh.SurfaceElements())
      {
        double qual = TriangleQualityInst (mesh[el[0]],
                                           mesh[el[1]],
                                           mesh[el[2]]);

        int cl = int ( (ncl-1e-3) * qual ) + 1;
        incl[cl-1]++;
      }

    (*testout) << endl << endl;

    (*testout) << "Points:           " << mesh.GetNP() << endl;
    (*testout) << "Surface Elements: " << mesh.GetNSE() << endl;

    (*testout) << endl;
    (*testout) << "Elements in qualityclasses:" << endl;
    // (*testout).precision(2);
    (*testout) << setprecision(2);
    for (int i = 1; i <= ncl; i++)
      {
        (*testout) << setw(4) << double (i-1)/ncl << " - "
                   << setw(4) << double (i) / ncl << ": "
                   << incl[i-1] << endl;
      }
  }


  static double TetElementQuality (const Point<3> & p1, const Point<3> & p2,
                                   const Point<3> & p3, const Point<3> & p4)
  {
    double vol, l, l4, l5, l6;


    Vec<3> v1 = p2 - p1;
    Vec<3> v2 = p3 - p1;
    Vec<3> v3 = p4 - p1;

    vol = fabs ((Cross (v1, v2) * v3)) / 6;
    l4 = Dist (p2, p3);
    l5 = Dist (p2, p4);
    l6 = Dist (p3, p4);

    l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

    if (vol <= 1e-8 * l * l * l) return 1e-10;

    return vol/(l*l*l) * 1832.82;    // 6^4 * sqrt(2)
  }




  // static double teterrpow = 2;

  double CalcTetBadness (const Point<3> & p1, const Point<3> & p2,
                         const Point<3> & p3, const Point<3> & p4, double h,
                         const MeshingParameters & mp)
  {
    double vol, l, ll, lll, ll1, ll2, ll3, ll4, ll5, ll6;
    double err;

    Vec<3> v1 (p1, p2);
    Vec<3> v2 (p1, p3);
    Vec<3> v3 (p1, p4);

    vol = Determinant (v1, v2, v3)  * (-0.166666666666666);

    ll1 = v1.Length2();
    ll2 = v2.Length2();
    ll3 = v3.Length2();
    ll4 = Dist2 (p2, p3);
    ll5 = Dist2 (p2, p4);
    ll6 = Dist2 (p3, p4);

    ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
    l = sqrt (ll);
    lll = l * ll;

    if (vol <= 1e-24 * lll)
      return 1e24;

    err = 0.0080187537 * lll / vol;    // sqrt(216) / (6^4 * sqrt(2))

    if (h > 0)
      err += ll / (h * h) + 
        h * h * ( 1 / ll1 + 1 / ll2 + 1 / ll3 + 
                  1 / ll4 + 1 / ll5 + 1 / ll6 ) - 12;
    
    double teterrpow = mp.opterrpow;
    if(teterrpow < 1) teterrpow = 1;
    
    if (teterrpow == 1) return err;
    if (teterrpow == 2) return err*err;
    return pow (err, teterrpow);
  }


  double CalcTetBadnessGrad (const Point<3> & p1, const Point<3> & p2,
                             const Point<3> & p3, const Point<3> & p4, double h,
                             int pi, Vec<3> & grad,
                             const MeshingParameters & mp)
  {
    double vol, l, ll, lll;
    double err;

    const Point<3> *pp1, *pp2, *pp3, *pp4;

    pp1 = &p1;
    pp2 = &p2;
    pp3 = &p3;
    pp4 = &p4;
  
    switch (pi)
      {
      case 2:
        {
          swap (pp1, pp2);
          swap (pp3, pp4);
          break;
        }
      case 3:
        {
          swap (pp1, pp3);
          swap (pp2, pp4);
          break;
        }
      case 4:
        {
          swap (pp1, pp4);
          swap (pp3, pp2);
          break;
        }
      }
  

    Vec<3> v1 (*pp1, *pp2);
    Vec<3> v2 (*pp1, *pp3);
    Vec<3> v3 (*pp1, *pp4);

    Vec<3> v4 (*pp2, *pp3);
    Vec<3> v5 (*pp2, *pp4);
    Vec<3> v6 (*pp3, *pp4);

    vol = Determinant (v1, v2, v3) * (-0.166666666666666);

    Vec<3> gradvol;
    Cross (v5, v4, gradvol);
    gradvol *= (-1.0/6.0);


    double ll1 = v1.Length2();
    double ll2 = v2.Length2();
    double ll3 = v3.Length2();
    double ll4 = v4.Length2();
    double ll5 = v5.Length2();
    double ll6 = v6.Length2();

    ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
    l = sqrt (ll);
    lll = l * ll;

    if (vol <= 1e-24 * lll)
      { 
        grad = Vec<3> (0, 0, 0);
        return 1e24;
      }



    Vec<3> gradll1 (*pp2, *pp1);
    Vec<3> gradll2 (*pp3, *pp1);
    Vec<3> gradll3 (*pp4, *pp1);
    gradll1 *= 2;
    gradll2 *= 2;
    gradll3 *= 2;

    Vec<3> gradll (gradll1);
    gradll += gradll2;
    gradll += gradll3;

    /*
    Vec<3> gradll;
    gradll = v1+v2+v3;
    gradll *= -2;
    */

    err = 0.0080187537 * lll / vol; 


    gradll *= (0.0080187537 * 1.5 * l / vol);
    Vec<3> graderr(gradll);
    gradvol *= ( -0.0080187537 * lll / (vol * vol) );
    graderr += gradvol;
  
    if (h > 0)
      {
        /*
        Vec<3> gradll1 (*pp2, *pp1);
        Vec<3> gradll2 (*pp3, *pp1);
        Vec<3> gradll3 (*pp4, *pp1);
        gradll1 *= 2;
        gradll2 *= 2;
        gradll3 *= 2;
        */
        err += ll / (h*h) + 
          h*h * ( 1 / ll1 + 1 / ll2 + 1 / ll3 + 
                  1 / ll4 + 1 / ll5 + 1 / ll6 ) - 12;

        graderr += (1/(h*h) - h*h/(ll1*ll1)) * gradll1;
        graderr += (1/(h*h) - h*h/(ll2*ll2)) * gradll2;
        graderr += (1/(h*h) - h*h/(ll3*ll3)) * gradll3;
      }

    double errpow;

    double teterrpow = mp.opterrpow;
    if(teterrpow < 1) teterrpow = 1;

    if (teterrpow == 1)
      {
        errpow = err;
        grad = graderr;
      }
    else if (teterrpow == 2)
      {
        errpow = err*err;   
        grad = (2 * err) * graderr;
      }
    else 
      {
        errpow = pow (err, teterrpow);
        grad = (teterrpow * errpow / err) * graderr;
      }
    return errpow;
  }
  




  /*

  double CalcTetBadness (const Point<3> & p1, const Point<3> & p2,
  const Point<3> & p3, const Point<3> & p4, double h)
  {
  double vol, l;
  double err;


  Vec<3> v1 (p1, p2);
  Vec<3> v2 (p1, p3);
  Vec<3> v3 (p1, p4);

  vol = -Determinant (v1, v2, v3) / 6;

  double l1 = v1.Length();
  double l2 = v2.Length();
  double l3 = v3.Length();
  double l4 = Dist (p2, p3);
  double l5 = Dist (p2, p4);
  double l6 = Dist (p3, p4);

  l = l1 + l2 + l3 + l4 + l5 + l6;

  // just for timing
  // l += 1e-40 * CalcTetBadnessNew (p1, p2, p3, p4, h);

  if (vol <= 1e-24 * l * l * l)
  { 
  return 1e24;
  }

  err = (l*l*l) / (1832.82 * vol);    // 6^4 * sqrt(2)
  
  if (h > 0)
  err += l / h + 
  h * (1 / l1 + 1/l2 + 1/l3 + 1/l4 + 1/l5 + 1/l6) - 12;

  return pow (err, teterrpow);
  }


  
  double CalcTetBadnessGrad (const Point<3> & p1, const Point<3> & p2,
  const Point<3> & p3, const Point<3> & p4, double h,
  int pi, Vec<3> & grad)
  {
  double vol, l;
  double err;

  const Point<3> *pp1, *pp2, *pp3, *pp4;

  pp1 = &p1;
  pp2 = &p2;
  pp3 = &p3;
  pp4 = &p4;
  
  switch (pi)
  {
  case 2:
  {
  swap (pp1, pp2);
  swap (pp3, pp4);
  break;
  }
  case 3:
  {
  swap (pp1, pp3);
  swap (pp2, pp4);
  break;
  }
  case 4:
  {
  swap (pp1, pp4);
  swap (pp3, pp2);
  break;
  }
  }
  

  Vec<3> v1 (*pp1, *pp2);
  Vec<3> v2 (*pp1, *pp3);
  Vec<3> v3 (*pp1, *pp4);

  Vec<3> v4 (*pp2, *pp3);
  Vec<3> v5 (*pp2, *pp4);
  Vec<3> v6 (*pp3, *pp4);


  //   Vec<3> n;
  //   Cross (v1, v2, n);
  //   vol = - (n * v3) / 6;


  vol = -Determinant (v1, v2, v3) / 6;  

  Vec<3> gradvol;
  Cross (v5, v4, gradvol);
  gradvol *= (-1.0/6.0);


  double l1 = v1.Length();
  double l2 = v2.Length();
  double l3 = v3.Length();
  double l4 = v4.Length();
  double l5 = v5.Length();
  double l6 = v6.Length();

  l = l1 + l2 + l3 +l4 + l5 + l6;

  Vec<3> gradl1 (*pp2, *pp1);
  Vec<3> gradl2 (*pp3, *pp1);
  Vec<3> gradl3 (*pp4, *pp1);
  gradl1 /= l1;
  gradl2 /= l2;
  gradl3 /= l3;

  Vec<3> gradl (gradl1);
  gradl += gradl2;
  gradl += gradl3;


  if (vol <= 1e-24 * l * l * l)
  { 
  grad = Vec<3> (0, 0, 0);
  return 1e24;
  }


  double c1 = 1.0 / 1832.82;      // 6^4 * sqrt(2)
  err = c1 * (l*l*l) / vol; 


  gradl *= (c1 * 3 * l * l / vol);
  Vec<3> graderr(gradl);
  gradvol *= ( -c1 * l * l * l / (vol * vol) );
  graderr+= gradvol;
  
  if (h > 0)
  {
  err += l / h + 
  h * ( 1 / l1 + 1 / l2 + 1 / l3 + 
  1 / l4 + 1 / l5 + 1 / l6 ) - 12;

  graderr += (1/h - h/(l1*l1)) * gradl1;
  graderr += (1/h - h/(l2*l2)) * gradl2;
  graderr += (1/h - h/(l3*l3)) * gradl3;
  cout << "?";
  }

  double errpow = pow (err, teterrpow);
  grad = (teterrpow * errpow / err) * graderr;
  
  return errpow;
  }
  
  */




  
  /*
    double CalcVolume (const Array<Point<3>> & points,
    const ElementRef & el)
    {
    Vec<3> v1 = points.Get(el[1]) - 
    points.Get(el[0]);
    Vec<3> v2 = points.Get(el[2]) - 
    points.Get(el[0]);
    Vec<3> v3 = points.Get(el[3]) - 
    points.Get(el[0]); 
         
    return -(Cross (v1, v2) * v3) / 6;   
    }  
  */

  double CalcVolume (FlatArray<Point<3>, PointIndex> points, 
                     const Array<Element> & elements)
  {
    double vol;
    Vec<3> v1, v2, v3;
  
    vol = 0;
    for (int i = 0; i < elements.Size(); i++)
      {
        v1 = points[elements[i][1]] - points[elements[i][0]];
        v2 = points[elements[i][2]] - points[elements[i][0]];
        v3 = points[elements[i][3]] - points[elements[i][0]];
        vol -= (Cross (v1, v2) * v3) / 6;        
      }
    return vol;
  }

  
  

  void MeshQuality3d (const Mesh & mesh, Array<int> * inclass)
  { 
    int ncl = 20;
    Array<int> incl(ncl);
    double sum = 0;
    int nontet  = 0;

    for (int i = 0; i < incl.Size(); i++)
      incl[i] = 0;

    for (ElementIndex ei : mesh.VolumeElements().Range())
      {
        if (mesh[ei].GetType() != TET)
          {
            nontet++;
            continue;
          }

        double qual = TetElementQuality (mesh.Point(mesh[ei][0]),
                                         mesh.Point(mesh[ei][1]),
                                         mesh.Point(mesh[ei][2]),
                                         mesh.Point(mesh[ei][3]));

        if (qual > 1) qual = 1;
        signed int cl = int (ncl * qual ) + 1;
     
        if (cl < 1) cl = 1; 
        if (cl > ncl) cl = ncl;

        incl[cl-1]++;
        if (inclass) (*inclass)[ei.Nr0()] = cl;
        sum += 1/qual;
      }

    (*testout) << endl << endl;
    (*testout) << "Points:           " << mesh.GetNP() << endl;
    (*testout) << "Volume Elements:  " << mesh.GetNE() << endl;
    if (nontet)
      (*testout) << nontet << " non tetrahedral elements" << endl;
    (*testout) << endl;

    (*testout) << "Volume elements in qualityclasses:" << endl;
    (*testout) << setprecision(2);
    for (int i = 1; i <= ncl; i++)
      {
        (*testout) << setw(4) << double (i-1)/ncl << " - "
                   << setw(4) << double (i) / ncl << ": "
                   << incl[i-1] << endl;
      }
    (*testout) << "total error: " << sum << endl;
  }


  void SaveEdges (const Mesh & mesh, const char * geomfile, double h, char * filename)
  {
    ofstream of (filename);
  
    of << "edges" << endl;
    of << geomfile << endl;
    of << h << endl;

    of << mesh.GetNP() << endl;
    for (PointIndex pi : mesh.Points().Range())
      of << mesh[pi](0) << " "
         << mesh[pi](1) << " "
         << mesh[pi](2) << "\n";
    
    of << 2 * mesh.GetNSeg() << endl;
    for (auto & seg2 : mesh.LineSegments())
      {
        const Segment * seg = &seg2;

        int seg_face = mesh.HasEdgeDescriptor(*seg) ? mesh.GetEdgeDescriptor(*seg).GetIndex().Nr1() : -1;
        of << (*seg)[1] << " " << (*seg)[0] << " " << seg_face << "\n";
      }
   
  }


  void SaveSurfaceMesh (const Mesh & mesh,
                        double h,
                        char * filename)

  {
    ofstream outfile(filename);

    outfile << "surfacemesh" << endl;
    outfile << h << endl;

    outfile << mesh.GetNP() << endl;
    for (PointIndex pi : mesh.Points().Range())
      outfile << mesh[pi](0) << " "
              << mesh[pi](1) << " "
              << mesh[pi](2) << endl;

  

    outfile << mesh.GetNSE() << endl;
    for (auto & el : mesh.SurfaceElements())
      {

        if (mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0)
          outfile << el[0] << " "
                  << el[1] << " "
                  << el[2] << endl;
        if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0)
          outfile << el[0] << " "
                  << el[2] << " "
                  << el[1] << endl;
      }
  }


#ifdef OLD
  void Save2DMesh (
                   const Mesh & mesh2d,
                   const Array<SplineSegment *> * splines,
                   ostream & outfile)

  {
    outfile.precision (6);
  
    outfile << "areamesh2" << endl;


    outfile << endl;
    outfile << mesh2d.GetNSeg() << endl;
    for (int i = 1; i <= mesh2d.GetNSeg(); i++)
      outfile << mesh2d.LineSegment(i).GetIndex() << "        "
              << mesh2d.LineSegment(i)[0] << " "
              << mesh2d.LineSegment(i)[1] << "  " << endl;
  

    outfile << mesh2d.GetNSE() << endl;
    for (int i = 1; i <= mesh2d.GetNSE(); i++)
      {
        outfile << mesh2d.SurfaceElement(i).GetIndex() << "         ";
        outfile << mesh2d.SurfaceElement(i).GetNP() << " ";
        for (int j = 0; j < mesh2d.SurfaceElement(i).GetNP(); j++)
          outfile << mesh2d.SurfaceElement(i)[j] << " ";
        outfile << endl;
      }

    outfile << mesh2d.GetNP() << endl;
    for (int i = 1; i <= mesh2d.GetNP(); i++)
      outfile << mesh2d.Point(i).X() << " "
              << mesh2d.Point(i).Y() << endl;

    if (splines)
      {
        outfile << splines->Size() << endl;
        for (int i = 1; i <= splines->Size(); i++)
          splines->Get(i) -> PrintCoeff (outfile);
      }
    else
      outfile << "0" << endl;
  }
#endif








  void SaveVolumeMesh (const Mesh & mesh, 
                       const NetgenGeometry & geometry,
                       char * filename)
  {
    ofstream outfile(filename);
    outfile << "volumemesh" << endl;

    outfile << mesh.GetNSE() << endl;
    for (auto & sel : mesh.SurfaceElements())
      {
        if (sel.GetIndex().IsValid())
          outfile << mesh.GetFaceDescriptor(sel.GetIndex ()).SurfNr()
                  << "\t";
        else
          outfile << "0" << "\t";
        outfile << sel[0] << " "
                << sel[1] << " "
                << sel[2] << endl;
      }
    outfile << mesh.GetNE() << endl;
    for (ElementIndex ei : mesh.VolumeElements().Range())
      outfile << mesh[ei].GetIndex() << "\t"
              << mesh[ei][0] << " " << mesh[ei][1] << " "
              << mesh[ei][2] << " " << mesh[ei][3] << endl;

    outfile << mesh.GetNP() << endl;
    for (PointIndex pi : mesh.Points().Range())
      outfile << mesh[pi](0) << " "
              << mesh[pi](1) << " "
              << mesh[pi](2) << endl;

#ifdef SOLIDGEOM
    outfile << geometry.GetNSurf() << endl;
    for (int i = 1; i <= geometry.GetNSurf(); i++)
      geometry.GetSurface(i) -> Print (outfile);
#endif
  }




  int CheckCode ()
  {
    return 1;

    /*
      char st[100];
      ifstream ist("pw");

      if (!ist.good()) return 0;
      ist >> st;
      if (strcmp (st, "JKULinz") == 0) return 1;
      return 0;
    */
  }



  /* ******************** CheckMesh ******************************* */

  /// Checks, whether mesh contains a valid 3d mesh
  int CheckMesh3D (const Mesh & mesh)
  {
    ClosedHashTable<SortedPointIndices<3>, int> faceused(mesh.GetNE()/3);
    int ok = 1;

    for (auto & el : mesh.SurfaceElements())
      {
      
        if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0 ||
            mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0)
          {
            faceused.Set ( { el[0], el[1], el[2] }, 1);
          }
      }
  
    for (auto el : mesh.VolumeElements())
      {

        for (int j = 1; j <= 4; j++)
          {
            PointIndex fp[3];
            int l = 0;
            for (int k = 1; k <= 4; k++)
              if (j != k)
                fp[l++] = el.PNum(k);

            SortedPointIndices<3> i3(fp[0], fp[1], fp[2]);
            if (faceused.Used(i3))
              faceused.Set(i3, faceused.Get(i3)+1);
            else
              faceused.Set (i3, 1);
          }
      }


    for (SurfaceElementIndex i : mesh.SurfaceElements().Range())
      {
        const Element2d & el = mesh[i];

        SortedPointIndices<3> i3(el[0], el[1], el[2]);
        int nel = faceused.Used(i3) ? faceused.Get(i3) : 0;
        if (nel != 2)
          {
            ok = 0;
            (*testout) << "face " << i.Nr1() << " with points " 
                       << i3[0] << "-" << i3[1] << "-" << i3[2] 
                       << " has " << nel << " elements" << endl;
          }
      }
  
    for (ElementIndex ei : mesh.VolumeElements().Range())
      {
        auto el = mesh[ei];

        for (int j = 1; j <= 4; j++)
          {
            PointIndex fp[3];
            int l = 0;
            for (int k = 1; k <= 4; k++)
              if (j != k)
                fp[l++] = el.PNum(k);

            SortedPointIndices<3> i3(fp[0], fp[1], fp[2]);
            int nel = faceused.Used(i3) ? faceused.Get(i3) : 0;
            if (nel != 2)
              {
                ok = 0;
                (*testout) << "element " << ei << " with face " 
                           << i3[0] << "-" << i3[1] << "-"
                           << i3[2] 
                           << " has " << nel << " elements" << endl;
              }
          }
      }





    /*
      for (i = 1; i <= faceused.GetNBags(); i++)
      for (j = 1; j <= faceused.GetBagSize(i); j++)
      {
      faceused.GetData(i, j, i3, k);
      if (k != 2)
      {
      (*testout) << "Face: " << i3.I1() << "-" 
      << i3.I2() << "-" << i3.I3() << " has " 
      << k << " Faces " << endl;
      cerr << "Face Error" << endl;
      ok = 0;
      }
      }
    */


    if (!ok)
      {
        (*testout) << "surfelements: " << endl;
        for (SurfaceElementIndex i : mesh.SurfaceElements().Range())
          {
            const Element2d & el = mesh[i];
            (*testout) << setw(5) << i.Nr1() << ":" 
                       << setw(6) << el.GetIndex() 
                       << setw(6) << el[0] 
                       << setw(4) << el[1] 
                       << setw(4) << el[2]  << endl;
          }
        (*testout) << "volelements: " << endl;
        for (ElementIndex ei : mesh.VolumeElements().Range())
          {
            auto el = mesh[ei];
            (*testout) << setw(5) << ei << ":" 
                       << setw(6) << el.GetIndex() 
                       << setw(6) << el[0] << setw(4) << el[1]
                       << setw(4) << el[2] << setw(4) << el[3] << endl;
          }
      }


    return ok;
  }



  void RemoveProblem (Mesh & mesh, int domainnr)
  {
    mesh.FindOpenElements(domainnr);
    int np = mesh.GetNP();

    Array<bool, PointIndex> ppoints(np);
  
    // int ndom = mesh.GetNDomains();

    PrintMessage (3, "Elements before Remove: ", mesh.GetNE());
    // for (k = 1; k <= ndom; k++)
    int k = domainnr;
      {
        ppoints = false;
      
        for (int i = 1; i <= mesh.GetNOpenElements(); i++)
          {
            const Element2d & sel = mesh.OpenElement(i);
            if (sel.GetIndex().Nr1() == k)
              {
                for (int j = 0; j < sel.GetNP(); j++)
                  ppoints[sel[j]] = true;
              }
          }

        for (auto el : mesh.VolumeElements())
          {
            if (el.GetIndex() == k)
              {
                int todel = 0;
                for (int j = 0; j < el.GetNP(); j++)
                  if (ppoints[el[j]])
                    todel = 1;
              
                if (el.GetNP() != 4)
                  todel = 0;
              
                if (todel)
                  {
                    el.Delete();
                    // ei--;
                  }
              }
          }
      }
  
    mesh.Compress();
    PrintMessage (3, "Elements after Remove: ", mesh.GetNE());
  }

}

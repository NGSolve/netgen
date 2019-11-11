#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>


namespace netgen
{



  int CSGeometry :: changeval = 0;



  TopLevelObject ::  
  TopLevelObject (Solid * asolid,
		  Surface * asurface)
  {
    solid = asolid;
    surface = asurface;

    SetRGB (0, 0, 1);
    SetTransparent (0);
    SetVisible (1); 
    SetLayer (1);

    if (!surface)
      maxh = solid->GetMaxH();
    else
      maxh = surface->GetMaxH();

    SetBCProp (-1);

    bcname = "default";
  }

  void TopLevelObject :: GetData (ostream & ost)
  {
    ost << red << " " << green << " " << blue << " " 
	<< transp << " " << visible << " ";
  }

  void TopLevelObject :: SetData (istream & ist)
  {
    ist >> red >> green >> blue >> transp >> visible;
  }


 
  Box<3> CSGeometry::default_boundingbox (Point<3> (-1000, -1000, -1000),
					  Point<3> ( 1000,  1000,  1000));


  CSGeometry :: CSGeometry ()
    : boundingbox (default_boundingbox),
      identicsurfaces (100), ideps(1e-9), filename("")
  {
    ;
  }

  CSGeometry :: CSGeometry (const string & afilename)
    : boundingbox (default_boundingbox),
      identicsurfaces (100), ideps(1e-9), filename(afilename)
  {
    changeval++;
  }

  CSGeometry :: ~CSGeometry ()
  {
    Clean();
  }

  PointGeomInfo CSGeometry :: ProjectPoint(int surfind, Point<3> & p) const
  {
    Point<3> hp = p;
    GetSurface(surfind)->Project (hp);
    p = hp;
    return PointGeomInfo();
  }

  bool CSGeometry :: ProjectPointGI(int surfind, Point<3> & p, PointGeomInfo & gi) const
  {
    GetSurface(surfind)->Project (p);
    return true;
  }

  void CSGeometry :: ProjectPointEdge(int surfind, INDEX surfind2,
                                      Point<3> & p, EdgePointGeomInfo* /*unused*/) const
  {
    Point<3> hp = p;
    ProjectToEdge (GetSurface(surfind),
                   GetSurface(surfind2), hp);
    p = hp;
  }


  Vec<3> CSGeometry :: GetNormal(int surfind, const Point<3> & p,
                                 const PointGeomInfo* /*unused*/) const
  {
    Vec<3> hn;
    GetSurface(surfind)->CalcGradient(p, hn);
    hn.Normalize();
    return hn;
  }

  void CSGeometry ::
  PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint,
               int surfi,
               const PointGeomInfo & gi1,
               const PointGeomInfo & gi2,
               Point<3> & newp, PointGeomInfo & newgi) const
  {
    Point<3> hnewp;
    hnewp = p1+secpoint*(p2-p1);
    if (surfi != -1)
      {
        GetSurface (surfi) -> Project (hnewp);
        newgi.trignum = 1;
      }

    newp = hnewp;
  }

  void CSGeometry :: PointBetweenEdge(const Point<3> & p1, const Point<3> & p2, double secpoint,
               int surfi1, int surfi2,
               const EdgePointGeomInfo & ap1,
               const EdgePointGeomInfo & ap2,
               Point<3> & newp, EdgePointGeomInfo & newgi) const
  {
    Point<3> hnewp = p1+secpoint*(p2-p1);
    //(*testout) << "hnewp " << hnewp << " s1 " << surfi1 << " s2 " << surfi2 << endl;
    if (surfi1 != -1 && surfi2 != -1 && surfi1 != surfi2)
      {
        netgen::ProjectToEdge (GetSurface(surfi1),
                               GetSurface(surfi2),
                               hnewp);
        // (*testout) << "Pointbetween, newp = " << hnewp << endl
        // << ", err = " << sqrt (sqr (hnewp(0))+ sqr(hnewp(1)) + sqr (hnewp(2))) - 1 << endl;
        newgi.edgenr = 1;
        //(*testout) << "hnewp (a1) " << hnewp << endl;
      }
    else if (surfi1 != -1)
      {
        GetSurface (surfi1) -> Project (hnewp);
        //(*testout) << "hnewp (a2) " << hnewp << endl;
      }

    newp = hnewp;
  };

  Vec<3> CSGeometry :: GetTangent(const Point<3> & p, int surfi1, int surfi2,
                                  const EdgePointGeomInfo & ap1) const
  {
    Vec<3> n1 = GetSurface (surfi1)->GetNormalVector (p);
    Vec<3> n2 = GetSurface (surfi2)->GetNormalVector (p);
    Vec<3> tau = Cross (n1, n2).Normalize();
    return tau;
  }

  void CSGeometry :: Clean ()
  {
    NgArray< Solid* > to_delete;
    
    for (int i = 0; i < solids.Size(); i++)
      if(!to_delete.Contains(solids[i]->S1()))
	to_delete.Append(solids[i]->S1());
    for (int i = 0; i < solids.Size(); i++)
      if(!to_delete.Contains(solids[i]))
	to_delete.Append(solids[i]);

    for(int i = 0; i < to_delete.Size(); i++)
      delete to_delete[i];    
    
    solids.DeleteAll ();

    for (int i = 0; i < splinecurves2d.Size(); i++)
      delete splinecurves2d[i];
    splinecurves2d.DeleteAll();
    
    /*
    for (int i = 0; i < surfaces.Size(); i++)
      delete surfaces[i];
    surfaces.DeleteAll ();
    */
    for(int i = 0; i<delete_them.Size(); i++)
      delete delete_them[i];
    delete_them.DeleteAll();
    surfaces.DeleteAll();

    for (int i = 0; i < toplevelobjects.Size(); i++)
      delete toplevelobjects[i];
    toplevelobjects.DeleteAll ();
    for (int i = 0; i < triapprox.Size(); i++)
      delete triapprox[i];
    triapprox.DeleteAll();

    for(int i = 0; i < identifications.Size(); i++)
      delete identifications[i];
    identifications.DeleteAll();

    for (int i = 0; i < singfaces.Size(); i++)
      delete singfaces[i];
    singfaces.DeleteAll();
    for (int i = 0; i < singedges.Size(); i++)
      delete singedges[i];
    singedges.DeleteAll();
    for (int i = 0; i < singpoints.Size(); i++)
      delete singpoints[i];
    singpoints.DeleteAll();

    changeval++;
  }


  extern int CSGGenerateMesh (CSGeometry & geom, 
			      shared_ptr<Mesh> & mesh, MeshingParameters & mparam);


  int CSGeometry :: GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam)
  {
    return CSGGenerateMesh (*this, mesh, mparam);
  }
  
  class WritePrimitivesIt : public SolidIterator
  {
    ostream & ost;
  public:
    WritePrimitivesIt (ostream & aost) : ost(aost) { ; }
    virtual ~WritePrimitivesIt () { ; }

    virtual void Do (Solid * sol);
  };

  void WritePrimitivesIt :: Do (Solid * sol) 
  {
    Primitive * prim = sol->GetPrimitive();
    if (prim)
      {
	const char * classname;
	NgArray<double> coeffs;

	prim -> GetPrimitiveData (classname, coeffs);

	if (sol->Name())
	  ost << "primitive " 
	      << sol->Name() << " "
	      << classname << "  " << coeffs.Size();
	for (int i = 0; i < coeffs.Size(); i++)
	  ost << " " << coeffs[i];
	ost << endl;
      }
  }


  void CSGeometry :: Save (string filename) const
  {
    ofstream ost (filename.c_str());
    Save (ost);
  }
  
  void CSGeometry :: Save (ostream & ost) const
  {
    ost << "boundingbox "
	<< boundingbox.PMin()(0) << " "
	<< boundingbox.PMin()(1) << " "
	<< boundingbox.PMin()(2) << " "
	<< boundingbox.PMax()(0) << " "
	<< boundingbox.PMax()(1) << " "
	<< boundingbox.PMax()(2) << endl;


    WritePrimitivesIt wpi(ost);
    IterateAllSolids (wpi, 1);

    for (int i = 0; i < solids.Size(); i++)
      {
	if (!solids[i]->GetPrimitive())
	  {
	    ost << "solid " << solids.GetName(i) << " ";
	    solids[i] -> GetSolidData (ost);
	    ost << endl;
	  }
      }

    for (int i = 0; i < GetNTopLevelObjects(); i++)
      {
	TopLevelObject * tlo = GetTopLevelObject (i);
	ost << "toplevel ";
	if (tlo -> GetSurface())
	  ost << "surface " << tlo->GetSolid()->Name() << " "
	      << tlo->GetSurface()->Name() << " ";
	else
	  ost << "solid " << tlo->GetSolid()->Name() << " ";
	tlo->GetData(ost);
	ost << endl;
      }

    for (int i = 0; i < identifications.Size(); i++)
      {
	ost << "identify ";
	identifications[i] -> GetData (ost);
	ost << endl;
      }

    ost << "end" << endl;
  }

 
  void CSGeometry :: Load (istream & ist)
  {
    //  CSGeometry * geo = new CSGeometry;
  
    char key[100], name[100], classname[100], sname[100];
    int ncoeff, i, j;
    NgArray<double> coeff;

    while (ist.good())
      {
	ist >> key;
	if (strcmp (key, "boundingbox") == 0)
	  {
	    Point<3> pmin, pmax;
	    ist >> pmin(0) >> pmin(1) >> pmin(2);
	    ist >> pmax(0) >> pmax(1) >> pmax(2);
	    SetBoundingBox (Box<3> (pmin, pmax));
	  }
	if (strcmp (key, "primitive") == 0)
	  {
	    ist >> name >> classname >> ncoeff;
	    coeff.SetSize (ncoeff);
	    for (i = 0; i < ncoeff; i++)
	      ist >> coeff[i];

	    Primitive * nprim = Primitive::CreatePrimitive (classname);
	    nprim -> SetPrimitiveData (coeff);
	    Solid * nsol = new Solid (nprim);

	    for (j = 0; j < nprim->GetNSurfaces(); j++)
	      {
		sprintf (sname, "%s,%d", name, j);
		AddSurface (sname, &nprim->GetSurface(j));
		nprim -> SetSurfaceId (j, GetNSurf());
	      }
	    SetSolid (name, nsol);
	  }
	else if (strcmp (key, "solid") == 0)
	  {
	    ist >> name;
	    Solid * nsol = Solid::CreateSolid (ist, solids);

	    cout << " I have found solid " << name << " = ";
	    nsol -> GetSolidData (cout);
	    cout << endl;

	    SetSolid (name, nsol);
	  }
	else if (strcmp (key, "toplevel") == 0)
	  {
	    char type[20], solname[50], surfname[50];
	    const Solid * sol = NULL;
	    const Surface * surf = NULL;
	    int nr;

	    ist >> type;
	    if (strcmp (type, "solid") == 0)
	      {
		ist >> solname;
		sol = GetSolid (solname);
	      }
	    if (strcmp (type, "surface") == 0)
	      {
		ist >> solname >> surfname;
		sol = GetSolid (solname);
		surf = GetSurface (surfname);
	      }
	    nr = SetTopLevelObject ((Solid*)sol, (Surface*)surf);
	    GetTopLevelObject (nr) -> SetData (ist);
	  }
	else if (strcmp (key, "identify") == 0)
	  {
	    char type[10], surfname1[50], surfname2[50];
	    const Surface * surf1;
	    const Surface * surf2;


	    ist >> type >> surfname1 >> surfname2;
	    surf1 = GetSurface(surfname1);
	    surf2 = GetSurface(surfname2);
	  
	    AddIdentification (new PeriodicIdentification 
			       (GetNIdentifications(),
				*this, surf1, surf2));
	  }
	else if (strcmp (key, "end") == 0)
	  break;
      }

    changeval++;
  }


  void CSGeometry :: DoArchive(Archive& archive)
  {
    archive & surfaces & solids & toplevelobjects & userpoints & userpoints_ref_factor
      & identpoints & boundingbox & isidenticto & ideps
      & filename & spline_surfaces & splinecurves2d & splinecurves3d & surf2prim;
    if(archive.Input())
      FindIdenticSurfaces(1e-6);
  }

  void CSGeometry :: SaveSurfaces (ostream & out) const
  {
    if(singfaces.Size() > 0 || singedges.Size() > 0 || singpoints.Size() > 0)
      {
	PrintMessage(3,"Singular faces/edges/points => no csg-information in .vol file");
	return;
      }


    
    NgArray<double> coeffs;
    const char * classname;

    out << "csgsurfaces " << GetNSurf() << "\n";
    for(int i=0; i<GetNSurf(); i++)
      {
	const OneSurfacePrimitive * sp = dynamic_cast< const OneSurfacePrimitive * > (GetSurface(i));
	const ExtrusionFace * ef = dynamic_cast< const ExtrusionFace * > (GetSurface(i));
	const RevolutionFace * rf = dynamic_cast< const RevolutionFace * > (GetSurface(i));
	const DummySurface * dummyf = dynamic_cast< const DummySurface * > (GetSurface(i));
        const SplineSurface * splines = dynamic_cast<const SplineSurface *> (GetSurface(i));

        if (splines)
          {
            splines->GetBase()->GetPrimitiveData(classname,coeffs);
            out << classname << " " << coeffs.Size() << "\n";
            for (int j=0; j<coeffs.Size(); j++)
              out << coeffs[j] << " ";
            out << "\n";
            for (auto cut : *(splines->GetCuts()))
              {
                cut->GetPrimitiveData(classname,coeffs);
                out << classname << " " << coeffs.Size() << "\n";
                for (int j=0; j<coeffs.Size(); j++)
                  out << coeffs[j] << " ";
                out << "\n";
              }
            return;
          }

	if(sp)
	  {
	    sp->GetPrimitiveData(classname,coeffs);
	
	    out << classname << " ";
	  }
	else if(ef)
	  {
	    out << "extrusionface ";
	    ef->GetRawData(coeffs);
	  }
	else if(rf)
	  {
	    out << "revolutionface ";
	    rf->GetRawData(coeffs);
	  }
	else if(dummyf)
	  {
	    out << "dummy ";
            coeffs.SetSize(0);
	  }
	else
	  throw NgException ("Cannot write csg surface. Please, contact developers!");
      
	
	out << coeffs.Size() << "\n";
	for(int j=0; j<coeffs.Size(); j++)
	  out << coeffs[j] << " ";
	    
	out << "\n";
      }
  }

  void CSGeometry :: LoadSurfaces (istream & in)
  {
    NgArray<double> coeffs;
    string classname;
    int nsurfaces,size;

    in >> classname;
    
    if(classname == "csgsurfaces")
      in >> nsurfaces;
    else
      nsurfaces = atoi(classname.c_str());
    
    Point<3> dummypoint(0,0,0);
    Vec<3> dummyvec(0,0,0);
    double dummydouble(0.1);

    for(int i=0; i<nsurfaces; i++)
      {
	in >> classname;
	in >> size;

	coeffs.SetSize(size);

	for(int j=0; j<size; j++)
	  in >> coeffs[j];

	if(classname == "plane")
	  {
	    Plane * plane = new Plane(dummypoint,dummyvec);
	    plane->SetPrimitiveData(coeffs);

	    AddSurface(plane);
	    delete_them.Append(plane);
	  }

	else if(classname == "sphere")
	  {
	    Sphere * sphere = new Sphere(dummypoint,dummydouble);
	    sphere->SetPrimitiveData(coeffs);

	    AddSurface(sphere);
	    delete_them.Append(sphere);
	  }

	else if(classname == "cylinder")
	  {
	    Cylinder * cylinder = new Cylinder(coeffs);

	    AddSurface(cylinder);
	    delete_them.Append(cylinder);
	  }

	else if(classname == "ellipticcylinder")
	  {
	    EllipticCylinder * cylinder = new EllipticCylinder(coeffs);
	    AddSurface(cylinder);
	    delete_them.Append(cylinder);
	  }


	else if(classname == "torus")
	  {
	    Torus * torus = new Torus(dummypoint,dummyvec,dummydouble, dummydouble);
	    torus->SetPrimitiveData(coeffs);
	    AddSurface(torus);
	    delete_them.Append(torus);
	  }


	else if(classname == "cone")
	  {
	    Cone * cone = new Cone(dummypoint,dummypoint,dummydouble,dummydouble);
	    cone->SetPrimitiveData(coeffs);

	    AddSurface(cone);
	    delete_them.Append(cone);
	  }

       else if(classname == "ellipticcone")
         {
	    EllipticCone * ellipticcone = new EllipticCone(dummypoint,dummyvec,dummyvec,dummydouble,dummydouble);

	    ellipticcone->SetPrimitiveData(coeffs);

	    AddSurface(ellipticcone);
	    delete_them.Append(ellipticcone);
	  }
	else if(classname == "extrusionface")
	  {
	    ExtrusionFace * ef =
	      new ExtrusionFace(coeffs);

	    AddSurface(ef);
	    delete_them.Append(ef);
	  }

	else if(classname == "revolutionface")
	  {
	    RevolutionFace * rf =
	      new RevolutionFace(coeffs);

	    AddSurface(rf);
	    delete_them.Append(rf);
	  }

	else if(classname == "dummy")
	  {
	    Surface * surf = new DummySurface();
            
	    AddSurface(surf);
	    delete_them.Append(surf);
	  }

      }    
  }
    

  void CSGeometry :: SaveToMeshFile (ostream & ost) const
  {
    SaveSurfaces (ost);
  }




  void CSGeometry :: AddSurface (Surface * surf)
  {
    static int cntsurfs = 0;
    cntsurfs++;
    char name[15];
    sprintf (name, "nnsurf%d", cntsurfs);
    AddSurface (name, surf);
  }
 
  void CSGeometry :: AddSurface (char * name, Surface * surf)
  { 
    (*testout) << "Adding surface " << name << ": " << (*surf) << endl;
    surfaces.Set (name, surf); 
    surf->SetName (name);
    changeval++; 
  }

  void CSGeometry :: AddSurfaces (Primitive * prim)
  {
    for (int i = 0; i < prim->GetNSurfaces(); i++)
      {
	AddSurface (&prim->GetSurface(i));
	prim->SetSurfaceId (i, GetNSurf()-1);
	surf2prim.Append (prim);
      }
  }

  const Surface * CSGeometry :: GetSurface (const char * name) const
  {
    if (surfaces.Used(name))
      return surfaces[name];
    else
      return NULL;
  }

  /*
  const Surface * CSGeometry :: GetSurface (int i) const
  {
    if (i >= 0 && i < surfaces.Size()) 
      return surfaces[i];
    else
      throw NgException ("CSGeometry::GetSurface out of range");
  }
  */




  void CSGeometry :: SetSolid (const char * name, Solid * sol)
  {
    Solid * oldsol = NULL;

    if (solids.Used (name))
      oldsol = solids[name];

    solids.Set (name, sol);
    sol->SetName (name);

    if (oldsol)
      {
	if (oldsol->op != Solid::ROOT ||
	    sol->op != Solid::ROOT)
	  {
	    cerr << "Setsolid: old or new no root" << endl;
	  }
	oldsol -> s1 = sol -> s1;
      }
    changeval++;
  }

  const Solid * CSGeometry :: GetSolid (const char * name) const
  {
    if (solids.Used(name))
      return solids[name];
    else
      return NULL;
  }

  



  const Solid * CSGeometry :: GetSolid (const string & name) const
  {
    if (solids.Used(name))
      return solids[name];
    else
      return NULL;
  }




  void CSGeometry :: SetSplineCurve (const char * name, SplineGeometry<2> * spl)
  {
    splinecurves2d.Set(name,spl);
  }
  void CSGeometry :: SetSplineCurve (const char * name, SplineGeometry<3> * spl)
  {
    splinecurves3d.Set(name,spl);
  }


  const SplineGeometry<2> * CSGeometry :: GetSplineCurve2d (const string & name) const
  {
    if (splinecurves2d.Used(name))
      return splinecurves2d[name];
    else
      return NULL;
  }
  const SplineGeometry<3> * CSGeometry :: GetSplineCurve3d (const string & name) const
  {
    if (splinecurves3d.Used(name))
      return splinecurves3d[name];
    else
      return NULL;
  }



    /*
  class RemoveDummyIterator : public SolidIterator
  {
  public:
  
    RemoveDummyIterator() { ; }
    virtual ~RemoveDummyIterator() { ; }
    virtual void Do(Solid * sol);
  };

  void RemoveDummyIterator :: Do(Solid * sol)
  {
    cerr << "remove dummy iterator is obsolete" << endl;

    if ( (sol->op == Solid::SUB || sol->op == Solid::SECTION || 
	  sol->op == Solid::UNION)
	 && sol->s1->op == Solid::DUMMY)
      sol->s1 = sol->s1->s1;
    if ( (sol->op == Solid::SECTION || sol->op == Solid::UNION)
	 && sol->s2->op == Solid::DUMMY)
      sol->s2 = sol->s2->s1;
  }
    */






  int CSGeometry :: SetTopLevelObject (Solid * sol, Surface * surf)
  {
    toplevelobjects.Append (new TopLevelObject (sol, surf));
    return toplevelobjects.Size()-1;
  }

  TopLevelObject * CSGeometry :: 
  GetTopLevelObject (const Solid * sol, const Surface * surf)
  {
    for (int i = 0; i < toplevelobjects.Size(); i++)
      {
	if (toplevelobjects[i]->GetSolid() == sol &&
	    toplevelobjects[i]->GetSurface() == surf)
	  return (toplevelobjects[i]);
      }
    return NULL;
  }

  void CSGeometry :: RemoveTopLevelObject (Solid * sol, Surface * surf)
  {
    for (int i = 0; i < toplevelobjects.Size(); i++)
      {
	if (toplevelobjects[i]->GetSolid() == sol &&
	    toplevelobjects[i]->GetSurface() == surf)
	  {
	    delete toplevelobjects[i];
	    toplevelobjects.DeleteElement (i+1);
	    changeval++;
	    break;
	  }
      }
  }

  void CSGeometry :: AddIdentification (Identification * ident)
  {
    identifications.Append (ident);
  }

  void CSGeometry :: SetFlags (const char * solidname, const Flags & flags)
  {
    Solid * solid = solids[solidname];
    NgArray<int> surfind;

    int i;
    double maxh = flags.GetNumFlag ("maxh", -1);
    if (maxh > 0 && solid)
      {
	solid->GetSurfaceIndices (surfind);

	for (i = 0; i < surfind.Size(); i++)
	  {
	    if (surfaces[surfind[i]]->GetMaxH() > maxh)
	      surfaces[surfind[i]] -> SetMaxH (maxh);
	  }

	solid->SetMaxH (maxh);
      }

    if ( flags.StringFlagDefined ("bcname") )
      {
	solid->GetSurfaceIndices (surfind);
	string bcn = flags.GetStringFlag("bcname", "default");
	for (i = 0; i < surfind.Size(); i++)
	  {
	    if(surfaces[surfind[i]]->GetBCName() == "default")
	      surfaces[surfind[i]]->SetBCName(bcn);
	  }
      }

    if (flags.StringListFlagDefined ("bcname"))
      {
	auto& bcname = flags.GetStringListFlag("bcname");

	Polyhedra * polyh;
	if(solid->S1())
	  polyh = dynamic_cast<Polyhedra *>(solid->S1()->GetPrimitive());
	else
	  polyh = dynamic_cast<Polyhedra *>(solid->GetPrimitive());

	if(polyh)
	  {
	    NgArray < NgArray<int> * > polysurfs;
	    polyh->GetPolySurfs(polysurfs);
	    if(bcname.Size() != polysurfs.Size())
	      cerr << "WARNING: solid \"" << solidname << "\" has " << polysurfs.Size()
		   << " surfaces and should get " << bcname.Size() << " bc-names!" << endl;
	    
	    for ( i = 0; i < min2(polysurfs.Size(),bcname.Size()); i++)
	      {
		for (int j = 0; j < polysurfs[i]->Size(); j++)
		  {
		    if(surfaces[(*polysurfs[i])[j]]->GetBCName() == "default")
		      surfaces[(*polysurfs[i])[j]]->SetBCName(bcname[i]);
		  }
		delete polysurfs[i];
	      }
	  }
	else
	  {
	    solid->GetSurfaceIndices (surfind);
	    if(bcname.Size() != surfind.Size())
	      cerr << "WARNING: solid \"" << solidname << "\" has " << surfind.Size()
		   << " surfaces and should get " << bcname.Size() << " bc-names!" << endl;
	    
	    for (i = 0; i < min2(surfind.Size(),bcname.Size()); i++)
	      {
		if(surfaces[surfind[i]]->GetBCName() == "default")
		  surfaces[surfind[i]]->SetBCName(bcname[i]);
	      }
	  }
      }

    if (flags.NumFlagDefined ("bc"))
      {
	solid->GetSurfaceIndices (surfind);
	int bc = int (flags.GetNumFlag("bc", -1));
	for (i = 0; i < surfind.Size(); i++)
	  {
	    if (surfaces[surfind[i]]->GetBCProperty() == -1)
	      surfaces[surfind[i]]->SetBCProperty(bc);
	  }
      }
   
    if (flags.NumListFlagDefined ("bc"))
      {
	const auto& bcnum = flags.GetNumListFlag("bc");

	Polyhedra * polyh;
	if(solid->S1())
	  polyh = dynamic_cast<Polyhedra *>(solid->S1()->GetPrimitive());
	else
	  polyh = dynamic_cast<Polyhedra *>(solid->GetPrimitive());

	if(polyh)
	  {
	    NgArray < NgArray<int> * > polysurfs;
	    polyh->GetPolySurfs(polysurfs);
	    if(bcnum.Size() != polysurfs.Size())
	      cerr << "WARNING: solid \"" << solidname << "\" has " << polysurfs.Size()
		   << " surfaces and should get " << bcnum.Size() << " bc-numbers!" << endl;
	    
	    for ( i = 0; i < min2(polysurfs.Size(),bcnum.Size()); i++)
	      {
		for (int j = 0; j < polysurfs[i]->Size(); j++)
		  {
		    if ( surfaces[(*polysurfs[i])[j]]->GetBCProperty() == -1 )
		      surfaces[(*polysurfs[i])[j]]->SetBCProperty(int(bcnum[i]));
		  }
		delete polysurfs[i];
	      }
	  }
	else
	  {
	    solid->GetSurfaceIndices (surfind);
	    if(bcnum.Size() != surfind.Size())
	      cerr << "WARNING: solid \"" << solidname << "\" has " << surfind.Size()
		   << " surfaces and should get " << bcnum.Size() << " bc-numbers!" << endl;
	    
	    for (i = 0; i < min2(surfind.Size(),bcnum.Size()); i++)
	      {
		if (surfaces[surfind[i]]->GetBCProperty() == -1)
		  surfaces[surfind[i]]->SetBCProperty(int(bcnum[i]));
	      }
	  }
      }

  }

  void CSGeometry :: FindIdenticSurfaces (double eps)
  {
    int inv;
    int nsurf = GetNSurf();


    isidenticto.SetSize(nsurf);
    for (int i = 0; i < nsurf; i++)
      isidenticto[i] = i;
  
    for (int i = 0; i < nsurf; i++)
      for (int j = i+1; j < nsurf; j++)
	{
	  if (GetSurface(j) -> IsIdentic (*GetSurface(i), inv, eps))
	    {
	      INDEX_2 i2(i, j);
	      identicsurfaces.Set (i2, inv);
	      isidenticto[j] = isidenticto[i];
	    }
	}

    (*testout) << "identicmap:" << endl;
    for (int i = 0; i < isidenticto.Size(); i++)
      (*testout) << i << " -> " << isidenticto[i] << endl;
  }
  

  
  void CSGeometry ::
  GetSurfaceIndices (const Solid * sol, 
		     const BoxSphere<3> & box, 
		     NgArray<int> & locsurf) const
  {
    ReducePrimitiveIterator rpi(box);
    UnReducePrimitiveIterator urpi;

    const_cast<Solid*> (sol) -> IterateSolid (rpi);
    sol -> GetSurfaceIndices (locsurf);
    const_cast<Solid*> (sol) -> IterateSolid (urpi);

    for (int i = locsurf.Size()-1; i >= 0; i--)
      {
	bool indep = 1;
	for (int j = 0; j < i; j++)
	  if (locsurf[i] == locsurf[j])
	    {
	      indep = 0;
	      break;
	    }

	if (!indep) locsurf.Delete(i);
      }
  }



  
  void CSGeometry ::
  GetIndependentSurfaceIndices (const Solid * sol, 
				const BoxSphere<3> & box, 
				NgArray<int> & locsurf) const
  {
    ReducePrimitiveIterator rpi(box);
    UnReducePrimitiveIterator urpi;

    ((Solid*)sol) -> IterateSolid (rpi);
    sol -> GetSurfaceIndices (locsurf);
    ((Solid*)sol) -> IterateSolid (urpi);

    for (int i = 0; i < locsurf.Size(); i++)
      locsurf[i] = isidenticto[locsurf[i]];

    for (int i = locsurf.Size()-1; i >= 0; i--)
      {
	bool indep = 1;
	for (int j = 0; j < i; j++)
	  if (locsurf[i] == locsurf[j])
	    {
	      indep = 0;
	      break;
	    }

	if (!indep) locsurf.Delete(i);
      }


    /*
    // delete identified
    for (int i = locsurf.Size()-1; i >= 0; i--)
      {
	bool indep = 1;
	for (int j = 0; j < i; j++)
	  {
	    if (identicsurfaces.Used (INDEX_2::Sort (locsurf[i], locsurf[j])) !=
		(isidenticto[locsurf[i]] == isidenticto[locsurf[j]]))
	      {
		cerr << "different result" << endl;
		exit(1);
	      }

	    if (isidenticto[locsurf[i]] == isidenticto[locsurf[j]])
	      {
		indep = 0;
		break;
	      }
	  }
	if (!indep)
	  locsurf.Delete(i);
      }

    for (int i = 0; i < locsurf.Size(); i++)
      locsurf[i] = isidenticto[locsurf[i]];
    */
  }

  /*
  void CSGeometry ::
  GetIndependentSurfaceIndices (const Solid * sol, 
				const Point<3> & p, Vec<3> & v,
				NgArray<int> & locsurf) const
  {
    cout << "very dangerous" << endl;
    Point<3> p2 = p + 1e-2 * v;
    BoxSphere<3> box (p2, p2);
    box.Increase (1e-3);
    box.CalcDiamCenter();
    GetIndependentSurfaceIndices (sol, box, locsurf);
  }
  */

  void CSGeometry ::
  GetIndependentSurfaceIndices (NgArray<int> & locsurf) const
  {
    for (int i = 0; i < locsurf.Size(); i++)
      locsurf[i] = isidenticto[locsurf[i]];

    for (int i = locsurf.Size()-1; i >= 0; i--)
      {
	bool indep = 1;
	for (int j = 0; j < i; j++)
	  if (locsurf[i] == locsurf[j])
	    {
	      indep = 0;
	      break;
	    }

	if (!indep) locsurf.Delete(i);
      }
  }









  void CSGeometry :: 
  CalcTriangleApproximation(double detail, double facets)
  {
    PrintMessage (1, "Calc Triangle Approximation");

    try
      {
    //    FindIdenticSurfaces (1e-6);
  
    int ntlo = GetNTopLevelObjects();

    for (int i = 0; i < triapprox.Size(); i++)
      delete triapprox[i];
    triapprox.SetSize (ntlo);

    NgArray<int> surfind;
    IndexSet iset(GetNSurf());

    for (int i = 0; i < ntlo; i++)
      {
	Solid * sol;
	Surface * surf;
	GetTopLevelObject (i, sol, surf);

	sol -> CalcSurfaceInverse ();

	TriangleApproximation * tams = new TriangleApproximation();
	triapprox[i] = tams;

	// sol -> GetSurfaceIndices (surfind);
	for (int j = 0; j < GetNSurf(); j++)
	  // for (int jj = 0; jj < surfind.Size(); jj++)
	  {
	    // int j = surfind[jj];

	    PrintMessageCR (3, "Surface ", j, "/", GetNSurf());
	    // PrintMessageCR (3, "Surface ", j, "/", surfind.Size());

	    if (surf && surf != GetSurface(j))
	      continue;

	    TriangleApproximation tas;
	    GetSurface (j) -> GetTriangleApproximation (tas, boundingbox, facets);

	    int oldnp = tams -> GetNP();

	    if (!tas.GetNP())
	      continue;

	    for (int k = 0; k < tas.GetNP(); k++)
	      {
		tams -> AddPoint (tas.GetPoint(k));
                Vec<3> n = GetSurface(j) -> GetNormalVector (tas.GetPoint(k)); 
		n.Normalize();
		if (GetSurface(j)->Inverse()) n *= -1;
		tams -> AddNormal (n);
	      }
	  
	    BoxSphere<3> surfbox;

	    if (tas.GetNP())
	      surfbox.Set (tas.GetPoint(0));
	    for (int k = 1; k < tas.GetNP(); k++)
	      surfbox.Add (tas.GetPoint(k));
	    surfbox.Increase (1e-6);
	    surfbox.CalcDiamCenter();

	    Solid * surflocsol = sol -> GetReducedSolid (surfbox);
	    if (!surflocsol)
	      continue;

	    for (int k = 0; k < tas.GetNT(); k++)
	      {
		const TATriangle & tri = tas.GetTriangle (k);

		// check triangle
		BoxSphere<3> box;
		box.Set (tas.GetPoint (tri[0]));
		box.Add (tas.GetPoint (tri[1]));
		box.Add (tas.GetPoint (tri[2]));
		box.Increase (1e-6);
		box.CalcDiamCenter();


		Solid * locsol = surflocsol -> GetReducedSolid (box);
		
		if (locsol)
		  {
		    TATriangle tria(j, 
				    tri[0] + oldnp,
				    tri[1] + oldnp,
				    tri[2] + oldnp);
                    
                    // tams -> AddTriangle (tria);

		    RefineTriangleApprox (locsol, j, box, detail, 
					  tria, *tams, iset, 1);

		    delete locsol;
		  }
	      }
	  }

	tams->RemoveUnusedPoints ();
	PrintMessage (2, "Object ", i, " has ", tams->GetNT(), " triangles");
      }
      }
    catch (exception)
      {
	cerr << "*************************************************************" << endl
	     << "****   out of memory problem in CSG visualization        ****" << endl
	     << "****   Restart netgen, and disable                       ****" << endl
	     << "****   'Draw Geometry' in Geometry -> CSG Options        ****" << endl
	     << "****   before loading the geometry                       ****" << endl
	     << "****   meshing will still work !                         ****" << endl
	     << "*************************************************************" << endl;
	exit(1);
      }
    Change();
  }



  void CSGeometry ::
  RefineTriangleApprox (Solid * locsol, 
			int surfind,
			const BoxSphere<3> & box, 
			double detail,
			const TATriangle & tria, 
			TriangleApproximation & tams,
			IndexSet & iset,
                        int level)
  {
    // if (level > 10) return;

    //tams.AddTriangle (tria);
    //(*testout) << "tria " << tams.GetPoint(tria[0]) << " - " << tams.GetPoint(tria[1]) << " - " << tams.GetPoint(tria[2]) 
    //       << " ( " << tria[0] << " " << tria[1] << " " << tria[2] << ")" <<endl;
    //return;

    int pinds[6];
    NgArrayMem<int,500> surfused(GetNSurf());
  
    ReducePrimitiveIterator rpi(box);
    UnReducePrimitiveIterator urpi;

    locsol -> IterateSolid (rpi);
    //  locsol -> GetSurfaceIndices (lsurfi);


    //    IndexSet iset(GetNSurf());
    locsol -> GetSurfaceIndices (iset);
    const NgArray<int> & lsurfi = iset.GetArray();

    locsol -> IterateSolid (urpi);

    int surfii = -1;
    for (int i = 0; i < lsurfi.Size(); i++)
      if (lsurfi[i] == surfind)
	{
	  surfii = i;
	  break;
	}

    if (surfii == -1)
      return;

    int cntindep = 0;

    for (int i = 0; i < lsurfi.Size(); i++)
      {
	int linkto = isidenticto[lsurfi[i]];
	surfused[linkto] = 0;
      }

    for (int i = 0; i < lsurfi.Size(); i++)
      {
	int linkto = isidenticto[lsurfi[i]];
	if (!surfused[linkto])
	  {
	    surfused[linkto] = 1;
	    cntindep++;
	  }
      }

    int inverse = surfaces[surfind]->Inverse();

    if (cntindep == 1)
      {
	tams.AddTriangle (tria);
	//(*testout) << "pos1 " << tams.GetPoint(tria[0]) << " - " << tams.GetPoint(tria[1]) << " - " << tams.GetPoint(tria[2]) << endl;
	return;
      }

    if (cntindep == 2)
      {
	// just 2 surfaces:
	// if smooth, project inner points to edge and finish

	int otherind = -1;

	for (int i = 0; i < lsurfi.Size(); i++)
	  {
	    INDEX_2 i2 (lsurfi[i], surfind);
	    i2.Sort();
	  
	    if (i != surfii && !identicsurfaces.Used(i2))
	      otherind = lsurfi[i];
	  }

	double kappa = GetSurface(otherind)-> MaxCurvature ();

	if (kappa * box.Diam() < 0.1)
	  {
	    int pnums[6];
	    static int between[3][3] =
	      { { 1, 2, 3 },
		{ 0, 2, 4 },
		{ 0, 1, 5 } };
	    int onsurface[3];

	    for (int j = 0; j < 3; j++)
	      {
		int pi = tria[j];
		pnums[j] = pi;


		onsurface[j] =  
		  !locsol->IsStrictIn (tams.GetPoint (pi), 1e-6) &&
		  locsol->IsIn (tams.GetPoint (pi), 1e-6);
		
		//
		/*
		static int nos=0;
		if(!onsurface[j])
		  {
		    nos++;
		    cout << "NOT ON SURFACE!! "<< nos << endl;
		  }
		*/
	      }
	  
	    for (int j = 0; j < 3; j++)
	      {
		int lpi1 = between[j][0];
		int lpi2 = between[j][1];
		int lpin = between[j][2];
		if (onsurface[lpi1] == onsurface[lpi2])
		  pnums[lpin] = -1;
		else
		  {
		    const Point<3> & p1 = tams.GetPoint (pnums[lpi1]);
		    const Point<3> & p2 = tams.GetPoint (pnums[lpi2]);
		    double f1 = GetSurface(otherind)->CalcFunctionValue (p1);
		    double f2 = GetSurface(otherind)->CalcFunctionValue (p2);

		    Point<3> pn;

		    double l2(100),l1(100);
		    if ( fabs (f1-f2) > 1e-20 )
		      {
			l2 = -f1/(f2-f1);
			l1 = f2/(f2-f1);
			pn = Point<3>(l1 * p1(0) + l2 * p2(0),
				      l1 * p1(1) + l2 * p2(1),
				      l1 * p1(2) + l2 * p2(2));
		      }
		    else
		      pn = p1;

// 		    if(fabs(pn(0)) > 4 || fabs(pn(1)) > 4 || fabs(pn(2)) > 4)
// 		      {
// 			cout << "p1 " << p1 << " p2 " << p2 
// 			     << " f1 " << f1 << " f2 " << f2
// 			     << " l1 " << l1 << " l2 " << l2 
// 			     << " pn " << pn << endl;

// 		      }

		    
		    //GetSurface (surfind)->Project (pn);
		    
		    pnums[lpin] = tams.AddPoint (pn);

		    GetSurface (surfind)->Project (pn);
		    
		    Vec<3> n;
		    n = GetSurface (surfind)->GetNormalVector (pn);
		    if (inverse) n *= -1;
		    tams.AddNormal(n);
		  }
	      }
	  
	    int vcase = 0;
	    if (onsurface[0]) vcase++;
	    if (onsurface[1]) vcase+=2;
	    if (onsurface[2]) vcase+=4;
	  
	    static int trias[8][6] =
	      { { 0, 0, 0,   0, 0, 0 },
		{ 1, 6, 5,   0, 0, 0 },
		{ 2, 4, 6,   0, 0, 0 },
		{ 1, 2, 4,   1, 4, 5 },
		{ 3, 5, 4,   0, 0, 0 },
		{ 1, 6, 4,   1, 4, 3 },
		{ 2, 3, 6,   3, 5, 6 },
		{ 1, 2, 3,   0, 0, 0 } };
	    static int ntrias[4] =
	      { 0, 1, 2, 1 };

	    int nvis = 0;
	    for (int j = 0; j < 3; j++)
	      if (onsurface[j])
		nvis++;

	    for (int j = 0; j < ntrias[nvis]; j++)
	      {
		TATriangle ntria(tria.SurfaceIndex(),
				 pnums[trias[vcase][3*j]-1],
				 pnums[trias[vcase][3*j+1]-1],
				 pnums[trias[vcase][3*j+2]-1]);
		//(*testout) << "pos2 " << tams.GetPoint(ntria[0]) << " - " << tams.GetPoint(ntria[1]) << " - " << tams.GetPoint(ntria[2]) << endl
		//	   << "( " << ntria[0] << " - " << ntria[1] << " - " << ntria[2] << ")" << endl;
		tams.AddTriangle (ntria);
	      }

	    /* saturn changes:

	    int pvis[3];
	    for (j = 0; j < 3; j++)
	    pvis[j] = !locsol->IsStrictIn (tams.GetPoint (j+1), 1e-6) &&
	    locsol->IsIn (tams.GetPoint (j+1), 1e-6);
	  
	    int newpi[3];
	    for (j = 0; j < 3; j++)
	    {
	    int pi1 = j;
	    int pi2 = (j+1) % 3;
	    int pic = j;

	    if (pvis[pi1] != pvis[pi2])
	    {
	    Point<3> hp = Center (tams.GetPoint (tria.PNum (pi1+1)),
	    tams.GetPoint (tria.PNum (pi2+1)));

	    newpi[j] = tams.AddPoint (hp);
	    Vec<3> n = tams.GetNormal (pi1);
	    tams.AddNormal (n);
	    }
	    else
	    newpi[j] = 0;
	    }

	    int nvis = 0;
	    for (j = 0; j <= nvis; j++)
	    if (pvis[j]) nvis++;

	    int si = tria.SurfaceIndex();
	    switch (nvis)
	    {
	    case 0:
	    break;
	    case 1:
	    {
	    int visj;
	    for (j = 0; j < 3; j++)
	    if (pvis[j]) visj = j;
	    int pivis = tria.PNum (visj+1);
	    int pic1 = newpi[(visj+1)%3];
	    int pic2 = newpi[(visj+2)%3];
		
	    cout << pivis << "," << pic1 << "," << pic2 << endl;
		
	    tams.AddTriangle (TATriangle (si, pivis, pic1,pic2));
	    break;
	    }
	    case 2:
	    {
	    int nvisj;
	    for (j = 0; j < 3; j++)
	    if (!pvis[j]) nvisj = j;

	    int pivis1 = tria.PNum ((nvisj+1)%3+1);
	    int pivis2 = tria.PNum ((nvisj+2)%3+1);
	    int pic1 = newpi[nvisj];
	    int pic2 = newpi[(nvisj+2)%3];

	    tams.AddTriangle (TATriangle (si, pivis1, pic1,pic2));
	    tams.AddTriangle (TATriangle (si, pivis1, pic1,pivis2));
	    break;
	    }
	    case 3:
	    {
	    tams.AddTriangle (tria);
	    break;
	    }
	    }

	    */
	    return;
	  }
      }

    // bisection
    if (box.Diam() < detail)
      {
	//cout << "returning" << endl;
	return;
      }

    for (int i = 0; i < 3; i++)
      pinds[i] = tria[i];
  
    static int between[3][3] =
      { { 0, 1, 5 },
	{ 0, 2, 4 },
	{ 1, 2, 3 } };
  
    for (int i = 0; i < 3; i++)
      {
	// int pi1 = tria[between[i][0]];

	Point<3> newp = Center (tams.GetPoint (tria[between[i][0]]),
				tams.GetPoint (tria[between[i][1]]));
	Vec<3> n;
	
	GetSurface(surfind)->Project (newp);

	n = GetSurface(surfind)->GetNormalVector (newp);
      
	pinds[between[i][2]] = tams.AddPoint (newp);
	if (inverse) n *= -1;
	tams.AddNormal (n);
      }
  
    static int trias[4][4] =
      { { 0, 5, 4 },
	{ 5, 1, 3 },
	{ 4, 3, 2 },
	{ 3, 4, 5 } };
 
    for (int i = 0; i < 4; i++)
      {
	TATriangle ntri(surfind,
			pinds[trias[i][0]],
			pinds[trias[i][1]],
			pinds[trias[i][2]]);

	// check triangle
	BoxSphere<3> nbox;
	nbox.Set (tams.GetPoint (ntri[0]));
	nbox.Add (tams.GetPoint (ntri[1]));
	nbox.Add (tams.GetPoint (ntri[2]));
	nbox.Increase (1e-8);
	nbox.CalcDiamCenter();

	Solid * nsol = locsol -> GetReducedSolid (nbox);

	if (nsol)
	  {
	    RefineTriangleApprox (nsol, surfind, nbox, 
				  detail, ntri, tams, iset, level+1);
	  
	    delete nsol;
	  }
      }
  }




  class ClearVisitedIt : public SolidIterator
  {
  public:
    ClearVisitedIt () { ; }
    virtual ~ClearVisitedIt () { ; }

    virtual void Do (Solid * sol)
    { 
      sol -> visited = 0;
    }
  };


  void CSGeometry :: 
  IterateAllSolids (SolidIterator & it, bool only_once) const
  {
    if (only_once)
      {
	ClearVisitedIt clit;
	for (int i = 0; i < solids.Size(); i++)
	  solids[i] -> IterateSolid (clit, 0);
      }

    for (int i = 0; i < solids.Size(); i++)
      solids[i] -> IterateSolid (it, only_once);
  }


  double CSGeometry ::  MaxSize () const
  {
    double maxs, mins;
    maxs = max3 (boundingbox.PMax()(0), 
		 boundingbox.PMax()(1), 
		 boundingbox.PMax()(2));
    mins = min3 (boundingbox.PMin()(0), 
		 boundingbox.PMin()(1), 
		 boundingbox.PMin()(2));
    return max2 (maxs, -mins) * 1.1;
  }



  class CSGeometryRegister : public GeometryRegister
  {
  public:
    virtual NetgenGeometry * Load (string filename) const;
    virtual NetgenGeometry * LoadFromMeshFile (istream & ist) const;
    // virtual VisualScene * GetVisualScene (const NetgenGeometry * geom) const;
  };

  extern CSGeometry * ParseCSG (istream & istr, CSGeometry *instance=nullptr);

  NetgenGeometry *  CSGeometryRegister :: Load (string filename) const
  {
    const char * cfilename = filename.c_str();
    if (strcmp (&cfilename[strlen(cfilename)-3], "geo") == 0)
      {
	PrintMessage (1, "Load CSG geometry file ", cfilename);

	ifstream infile(cfilename);

	CSGeometry * hgeom = ParseCSG (infile);
	if (!hgeom)
	  throw NgException ("geo-file should start with 'algebraic3d'");

	hgeom -> FindIdenticSurfaces(1e-8 * hgeom->MaxSize()); 
	return hgeom;
      }

    if (strcmp (&cfilename[strlen(cfilename)-3], "ngg") == 0)
      {
	PrintMessage (1, "Load new CSG geometry file ", cfilename);

	ifstream infile(cfilename);
	CSGeometry * hgeom = new CSGeometry("");
	hgeom -> Load (infile);

	return hgeom;
      }


    
    return NULL;
  }

  NetgenGeometry * CSGeometryRegister :: LoadFromMeshFile (istream & ist) const 
  {
    string auxstring;
    if (ist.good())
      {
	ist >> auxstring;
	if (auxstring == "csgsurfaces")
	  {
	    CSGeometry * geometry = new CSGeometry ("");
	    geometry -> LoadSurfaces(ist);
	    return geometry;
	  }
	// else
	// ist.putback (auxstring);
      }
    return NULL;
  }


  class CSGInit
  {
  public:
    CSGInit()
    {
      geometryregister.Append (new CSGeometryRegister);
    }
  };

  CSGInit csginit;
  static RegisterClassForArchive<CSGeometry, NetgenGeometry> regcsg;
}

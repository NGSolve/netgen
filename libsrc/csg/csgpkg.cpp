#include <mystdlib.h>
#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>


#include <inctcl.hpp>
#include <visual.hpp>


#include "vscsg.hpp"


extern "C" int Ng_CSG_Init (Tcl_Interp * interp);



namespace netgen
{
  extern DLL_HEADER shared_ptr<NetgenGeometry> ng_geometry;
  extern DLL_HEADER shared_ptr<Mesh> mesh;

  static VisualSceneGeometry vsgeom;
 
  char * err_needscsgeometry = (char*) "This operation needs an CSG geometry";
  char * err_needsmesh = (char*) "This operation needs a mesh";
  char * err_jobrunning = (char*) "Meshing Job already running";

 


  int Ng_ParseGeometry (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    CSGeometry * csgeom = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (csgeom)
      {
	double detail = atof (Tcl_GetVar (interp, "::geooptions.detail", 0));
	double facets = atof (Tcl_GetVar (interp, "::geooptions.facets", 0));
      
	if (atoi (Tcl_GetVar (interp, "::geooptions.drawcsg", 0)))
	  csgeom->CalcTriangleApproximation(detail, facets);
      }
    return TCL_OK;
  }




  int Ng_GeometryOptions (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());


    const char * command = argv[1];

    if (strcmp (command, "get") == 0)
      {
	if (geometry)
	  {
	    char buf[20];
	    Point3d pmin = geometry->BoundingBox ().PMin();
	    Point3d pmax = geometry->BoundingBox ().PMax();
	    
	    snprintf (buf, size(buf), "%5.1lf", pmin.X());
	    Tcl_SetVar (interp, "::geooptions.minx", buf, 0);
	    snprintf (buf, size(buf), "%5.1lf", pmin.Y());
	    Tcl_SetVar (interp, "::geooptions.miny", buf, 0);
	    snprintf (buf, size(buf), "%5.1lf", pmin.Z());
	    Tcl_SetVar (interp, "::geooptions.minz", buf, 0);
	    
	    snprintf (buf, size(buf), "%5.1lf", pmax.X());
	    Tcl_SetVar (interp, "::geooptions.maxx", buf, 0);
	    snprintf (buf, size(buf), "%5.1lf", pmax.Y());
	    Tcl_SetVar (interp, "::geooptions.maxy", buf, 0);
	    snprintf (buf, size(buf), "%5.1lf", pmax.Z());
	    Tcl_SetVar (interp, "::geooptions.maxz", buf, 0);
	  }
      }
    else if (strcmp (command, "set") == 0)
      {
        Point<3> pmin (atof (Tcl_GetVar (interp, "::geooptions.minx", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.miny", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.minz", 0)));
        Point<3> pmax (atof (Tcl_GetVar (interp, "::geooptions.maxx", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.maxy", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.maxz", 0)));
	Box<3> box (pmin, pmax);
	if (geometry)
	  geometry -> SetBoundingBox (box);
	CSGeometry::SetDefaultBoundingBox (box);
      }

    return TCL_OK;
  }





  // attempt of a simple modeller

  int Ng_CreatePrimitive (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {/*
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }


    tcl_const char * classname = argv[1];
    tcl_const char * name = argv[2];

    cout << "Create primitive, class = " << classname
	 << ", name = " << name << endl;

    Primitive * nprim = Primitive::CreatePrimitive (classname);
    Solid * nsol = new Solid (nprim);

    char sname[100];
    for (int j = 1; j <= nprim->GetNSurfaces(); j++)
      {
	sprintf (sname, "%s,%d", name, j);
	geometry -> AddSurface (sname, &nprim->GetSurface(j));
	nprim -> SetSurfaceId (j, geometry->GetNSurf());
      }

    geometry->SetSolid (name, nsol);
	*/
    return TCL_OK;
  }


  int Ng_SetPrimitiveData (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }


    tcl_const char * name = argv[1];
    tcl_const char * value = argv[2];

    NgArray<double> coeffs;


    cout << "Set primitive data, name = " << name
	 << ", value = " << value  << endl;


    istringstream vst (value);
    double val;
    while (!vst.eof())
      {
	vst >> val;
	coeffs.Append (val);
      }

    ((Primitive*)
     geometry->GetSolid (name)->GetPrimitive())->SetPrimitiveData (coeffs);

    return TCL_OK;
  }



  int Ng_SetSolidData (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {/*
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }


    tcl_const char * name = argv[1];
    tcl_const char * val = argv[2];

    cout << "Set Solid Data, name = " << name
	 << ", value = " << val << endl;

    istringstream vst (val);

    Solid * nsol = Solid::CreateSolid (vst, geometry->GetSolids());
    geometry->SetSolid (name, nsol);
	*/
    return TCL_OK;
  }


  int Ng_GetPrimitiveData (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }


    tcl_const char * name = argv[1];
    tcl_const char * classnamevar = argv[2];
    tcl_const char * valuevar = argv[3];

    const char * classname;

    NgArray<double> coeffs;

    geometry->GetSolid (name)->GetPrimitive()->GetPrimitiveData (classname, coeffs);

    ostringstream vst;
    for (int i = 1; i <= coeffs.Size(); i++)
      vst << coeffs.Get(i) << " ";

    cout << "GetPrimitiveData, name = " << name
	 << ", classnamevar = " << classnamevar
	 << ", classname = " << classname << endl
	 << " valuevar = " << valuevar
	 << ", values = " << vst.str() << endl;

    Tcl_SetVar  (interp, classnamevar, (char*)classname, 0);
    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }

  int Ng_GetSolidData (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {/*
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }

    tcl_const char * name = argv[1];
    tcl_const char * valuevar = argv[2];

    ostringstream vst;

    const Solid * sol = geometry->GetSolid (name);
    sol->GetSolidData (vst);

    cout << "GetSolidData, name = " << name << ", data = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);
	*/
    return TCL_OK;
  }


  int Ng_GetPrimitiveList (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }


    tcl_const char * valuevar = argv[1];
    int i;

    stringstream vst;

    for (i = 1; i <= geometry->GetNSolids(); i++)
      {
	const Solid * sol = geometry->GetSolid(i);
	if (sol->GetPrimitive())
	  vst << sol->Name() << " ";
      }

    cout << "primnames = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }



  int Ng_GetSurfaceList (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }


    tcl_const char * valuevar = argv[1];
    int i;

    stringstream vst;

    for (i = 1; i <= geometry->GetNSurf(); i++)
      {
	const Surface * surf = geometry->GetSurface(i);
	vst << surf->Name() << " ";
      }

    cout << "surfnames = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }


  int Ng_GetSolidList (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }

    tcl_const char * valuevar = argv[1];
    int i;

    stringstream vst;

    for (i = 1; i <= geometry->GetNSolids(); i++)
      {
	const Solid * sol = geometry->GetSolid(i);
	if (!sol->GetPrimitive())
	  vst << sol->Name() << " ";
      }

    cout << "solnames = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }


  int Ng_TopLevel (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }


    int i;
    /*
      for (i = 0; i < argc; i++)
      cout << argv[i] << ", ";
      cout << endl;
    */

    if (strcmp (argv[1], "getlist") == 0)
      {
	stringstream vst;

	for (i = 0; i < geometry->GetNTopLevelObjects(); i++)
	  {
	    const Solid * sol;
	    const Surface * surf;
	    geometry->GetTopLevelObject (i, sol, surf);

	    if (!surf)
	      vst << "{ " << sol->Name() << " } ";
	    else
	      vst << "{ " << sol->Name() << " " << surf->Name() << " } ";
	  }

	tcl_const char * valuevar = argv[2];
	Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);
      }

    if (strcmp (argv[1], "set") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	geometry->SetTopLevelObject (sol, surf);
      }

    if (strcmp (argv[1], "remove") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	geometry->RemoveTopLevelObject (sol, surf);
      }

    if (strcmp (argv[1], "setprop") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	tcl_const char * propvar = argv[4];
	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	TopLevelObject * tlo = geometry->GetTopLevelObject (sol, surf);

	if (!tlo) return TCL_OK;

	char varname[50];
	snprintf (varname, size(varname), "%s(red)", propvar);
	double red = atof (Tcl_GetVar (interp, varname, 0));
	snprintf (varname, size(varname), "%s(blue)", propvar);
	double blue = atof (Tcl_GetVar (interp, varname, 0));
	snprintf (varname, size(varname), "%s(green)", propvar);
	double green = atof (Tcl_GetVar (interp, varname, 0));
	tlo -> SetRGB (red, green, blue);

	snprintf (varname, size(varname), "%s(visible)", propvar);
	tlo -> SetVisible (bool(atoi (Tcl_GetVar (interp, varname, 0))));
	snprintf (varname, size(varname), "%s(transp)", propvar);
	tlo -> SetTransparent (bool(atoi (Tcl_GetVar (interp, varname, 0))));
      }

    if (strcmp (argv[1], "getprop") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	tcl_const char * propvar = argv[4];

	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	TopLevelObject * tlo = geometry->GetTopLevelObject (sol, surf);

	if (!tlo) return TCL_OK;

	char varname[50], varval[10];

	snprintf (varname, size(varname), "%s(red)", propvar);
	snprintf (varval, size(varval), "%lf", tlo->GetRed());
	Tcl_SetVar (interp, varname, varval, 0);

	snprintf (varname, size(varname), "%s(green)", propvar);
	snprintf (varval, size(varval), "%lf", tlo->GetGreen());
	Tcl_SetVar (interp, varname, varval, 0);

	snprintf (varname, size(varname), "%s(blue)", propvar);
	snprintf (varval, size(varval), "%lf", tlo->GetBlue());
	Tcl_SetVar (interp, varname, varval, 0);

	snprintf (varname, size(varname), "%s(visible)", propvar);
	snprintf (varval, size(varval), "%d", tlo->GetVisible());
	Tcl_SetVar (interp, varname, varval, 0);

	snprintf (varname, size(varname), "%s(transp)", propvar);
	snprintf (varval, size(varval), "%d", tlo->GetTransparent());
	Tcl_SetVar (interp, varname, varval, 0);
      }


    return TCL_OK;
  }




  int Ng_SingularEdgeMS (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }

    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    // double globh = mparam.maxh;
    for (int i = 1; i <= geometry->singedges.Size(); i++)
      geometry->singedges.Get(i)->SetMeshSize (*mesh, 1e99 /* globh*/);
    return TCL_OK;
  }


  int Ng_SingularPointMS (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    if (!geometry)
      {
	Tcl_SetResult (interp, err_needscsgeometry, TCL_STATIC);
	return TCL_ERROR;
      }

    // double globh = mparam.maxh;
    for (int i = 1; i <= geometry->singpoints.Size(); i++)
      geometry->singpoints.Get(i)->SetMeshSize (*mesh, 1e99 /* globh */ );
    return TCL_OK;
  }



  int Ng_SelectSurface (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    int surfnr = atoi (argv[1]);
    vsgeom.SelectSurface (surfnr);
    return TCL_OK;
  }


  class CSGeometryVisRegister : public GeometryRegister
  {
  public:
    virtual NetgenGeometry * Load (const filesystem::path & filename) const { return NULL; }
    virtual VisualScene * GetVisualScene (const NetgenGeometry * geom) const;
  };

  VisualScene * CSGeometryVisRegister :: GetVisualScene (const NetgenGeometry * geom) const
  {
    const CSGeometry * geometry = dynamic_cast<const CSGeometry*> (geom);
    if (geometry)
      {
	vsgeom.SetGeometry (const_cast<CSGeometry*>(geometry));
	return &vsgeom;
      }
    return NULL;
  }
}


using namespace netgen;

int Ng_CSG_Init (Tcl_Interp * interp)
{
  geometryregister.Append (new CSGeometryVisRegister);
  if (interp == NULL) return TCL_OK;
  

  Tcl_CreateCommand (interp, "Ng_ParseGeometry", Ng_ParseGeometry,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  // geometry
  Tcl_CreateCommand (interp, "Ng_CreatePrimitive", Ng_CreatePrimitive,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_SetPrimitiveData", Ng_SetPrimitiveData,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_GetPrimitiveData", Ng_GetPrimitiveData,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_GetPrimitiveList", Ng_GetPrimitiveList,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);


  Tcl_CreateCommand (interp, "Ng_GetSurfaceList", Ng_GetSurfaceList,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);



  Tcl_CreateCommand (interp, "Ng_SetSolidData", Ng_SetSolidData,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_GetSolidData", Ng_GetSolidData,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_GetSolidList", Ng_GetSolidList,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);


  Tcl_CreateCommand (interp, "Ng_TopLevel", Ng_TopLevel,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_GeometryOptions", Ng_GeometryOptions,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_SingularEdgeMS", Ng_SingularEdgeMS,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);
  
  Tcl_CreateCommand (interp, "Ng_SingularPointMS", Ng_SingularPointMS,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SelectSurface", Ng_SelectSurface,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


  return TCL_OK;
}




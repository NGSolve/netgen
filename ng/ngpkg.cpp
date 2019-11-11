/*

The interface between the GUI and the netgen library

*/

#include <mystdlib.h>
#include <myadt.hpp>
#include <linalg.hpp>

#include <meshing.hpp>


#include <inctcl.hpp>
#include <visual.hpp>


#include <csg.hpp>

#ifdef SOCKETS
#include "../libsrc/sockets/sockets.hpp"
#include "../libsrc/sockets/socketmanager.hpp"
#endif


// to be sure to include the 'right' togl-version
#include "Togl2.1/togl.h"
// EXTERN int Togl_PixelScale (const Togl * togl);

#include "fonts.hpp"

extern bool nodisplay;

#include <nginterface.h>


#include "../libsrc/interface/writeuser.hpp"

namespace netgen
{
  DLL_HEADER extern MeshingParameters mparam;
  DLL_HEADER extern void ImportSolution2(const char * filename);
#include "demoview.hpp"
}


#ifdef ACIS
#include "ng_acis.hpp"
#endif


#ifdef JPEGLIB
#include <jpeglib.h>
#endif

#ifdef FFMPEG
#include "encoding.hpp"
#endif

#ifdef NGSOLVE
extern "C" void NGSolve_Exit();
#endif

// extern void * ngsolve_handle;


namespace netgen
{
  extern Flags parameters;

  /*
  NetgenOutStream operator<< ( ostream & ost, Imp  imp )
  {
    return ( NetgenOutStream ( &ost, imp ) );
  }

  NetgenOutStream operator<< ( ostream & ost, Proc proc )
  {
    return ( NetgenOutStream ( &ost, proc ) );
  }


  NetgenOutStream operator<< ( ostream & ost, Procs & procs )
  {
    return ( NetgenOutStream ( &ost, procs ) );
  }
  */

  DLL_HEADER extern std::shared_ptr<NetgenGeometry> ng_geometry;
  DLL_HEADER extern std::shared_ptr<Mesh> mesh;
  Tcl_Interp * tcl_interp;


#ifdef SOCKETS
  AutoPtr<ClientSocket> clientsocket;
  ServerSocketManager serversocketmanager;
  //NgArray< AutoPtr < ServerInfo > > servers;
  NgArray< ServerInfo* > servers;
  AutoPtr<ServerSocketUserNetgen> serversocketusernetgen;
#endif



  // visualization scenes, pointer vs selects which one is drawn:

  static VisualScene vscross;
  DLL_HEADER extern VisualSceneSurfaceMeshing vssurfacemeshing;
  DLL_HEADER extern VisualSceneMeshDoctor vsmeshdoc;

  static VisualSceneSpecPoints vsspecpoints;



  VisualScene *vs = &vscross;



  extern char * err_needsmesh;// = (char*) "This operation needs a mesh";
  extern char * err_jobrunning;// = (char*) "Meshing Job already running";





#ifndef SMALLLIB
//  // Destination for messages, errors, ...
#ifndef WIN32
  DLL_HEADER void Ng_PrintDest(const char * s)
  {
    /*
#ifdef PARALLEL
    int id, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
#else
    int id = 0; int ntasks = 1;
#endif
    */

    if (id == 0)
      (*mycout) << s << flush;

    /*
    if ( ntasks == 1 )
      (*mycout) << s << flush;
    else
      (*mycout) << "p" << id << ": " << s << flush ;
    */
  }
#endif
  void MyError2(const char * ch)
  {
    cout << ch;
    (*testout) << "Error !!! " << ch << endl << flush;
  }
#endif

  static clock_t starttimea;
  void ResetTime2 ()
  {
    starttimea = clock();
  }

#ifndef SMALLLIB
  double GetTime2 ()
  {
    return double(clock() - starttimea) / CLOCKS_PER_SEC;
  }
#endif





  // file handling ..
  int Ng_New (ClientData clientData,
	      Tcl_Interp * interp,
	      int argc, tcl_const char *argv[])
  {
    if (strcmp (argv[1], "mesh") == 0)
      mesh.reset();

    if (strcmp (argv[1], "geom") == 0)
      {
        /*
	delete ng_geometry;
	ng_geometry = new NetgenGeometry;
        */
        ng_geometry = make_shared<NetgenGeometry>();
      }

    return TCL_OK;
  }




  int Ng_ImportMesh (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[]);

  int Ng_LoadMesh (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    string filename (argv[1]);

    if (filename.find(".vol") == string::npos) 
      {
	return Ng_ImportMesh(clientData,interp,argc,argv);
      }

    PrintMessage (1, "load mesh from file ", filename);

    mesh = make_shared<Mesh>();
    try
      {
        istream * infile;
        // if (filename.substr (filename.length()-3, 3) == ".gz")
        if (filename.find(".vol.gz") != string::npos)
          infile = new igzstream (filename.c_str());
        else
          infile = new ifstream (filename.c_str());

	// ifstream infile(filename.c_str());
	mesh -> Load(*infile);
        // vsmesh.SetMesh (mesh);
        SetGlobalMesh (mesh);

#ifdef PARALLEL
	MyMPI_SendCmd ("mesh");
	mesh -> Distribute();
#endif
	for (int i = 0; i < geometryregister.Size(); i++)
	  {
	    NetgenGeometry * hgeom = geometryregister[i]->LoadFromMeshFile (*infile);
	    if (hgeom)
	      {
                ng_geometry = shared_ptr<NetgenGeometry>(hgeom);
		break;
	      }
	  }
        delete infile;

	/*
	string auxstring;
	if(infile.good())
	  {
	    infile >> auxstring;
	    if(auxstring == "csgsurfaces")
	      {
		CSGeometry * geometry = new CSGeometry ("");
		geometry -> LoadSurfaces(infile);

		delete ng_geometry;
		ng_geometry = geometry;
	      }
	  }
	*/
      }
    catch (NgException e)
      {
	PrintMessage (3, e.What());
	return TCL_ERROR;
      }

    PrintMessage (2,  mesh->GetNP(), " Points, ",
		  mesh->GetNE(), " Elements.");

    return TCL_OK;
  }





  int Ng_SaveMesh (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    string filename (argv[1]);
    PrintMessage (1, "Save mesh to file ", filename, ".... Please Wait!");
    
    ostream * outfile;
    if (filename.substr (filename.length()-3, 3) == ".gz")
      outfile = new ogzstream (filename.c_str());
    else
      outfile = new ofstream (filename.c_str());

    mesh -> Save (*outfile);
    // *outfile << endl << endl << "endmesh" << endl << endl;

    if (ng_geometry && !mesh->GetGeometry())
      ng_geometry -> SaveToMeshFile (*outfile);

    delete outfile;
    PrintMessage (1, "Save mesh to file .... DONE!");
    return TCL_OK;
  }





  int Ng_MergeMesh (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    string filename (argv[1]);

    PrintMessage (1, "merge with mesh from file ", filename);

    try
      {
	CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry.get());
    
	//mesh -> Merge (filename);
	ifstream infile(filename.c_str());
	const int offset = (geometry) ? geometry->GetNSurf() : 0;
	mesh -> Merge(infile,offset);

	string auxstring;
	if(infile.good())
	  {
	    infile >> auxstring;
	    if(geometry && auxstring == "csgsurfaces")
	      geometry -> LoadSurfaces(infile);
	  }
      }
    catch (NgException e)
      {
	PrintMessage (3, e.What());
	return TCL_ERROR;
      }

    PrintMessage (2,  mesh->GetNP(), " Points, ",
		  mesh->GetNSE(), " Surface Elements.");

    return TCL_OK;
  }


  int Ng_GetExportFormats (ClientData clientData,
                           Tcl_Interp * interp,
                           int argc, tcl_const char *argv[])
  {
    NgArray<const char*> userformats;
    NgArray<const char*> extensions;
    RegisterUserFormats (userformats, extensions);
    
    ostringstream fstr;
    for (int i = 1; i <= userformats.Size(); i++)
      fstr << "{ {" << userformats.Get(i) << "} {" << extensions.Get(i) << "} }\n";
    
    Tcl_SetResult (interp, const_cast<char*>(fstr.str().c_str()), TCL_VOLATILE);
    return TCL_OK;
  }


  int Ng_ExportMesh (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    string filename (argv[1]);
    string filetype (argv[2]);
    PrintMessage (1, "Export mesh to file ", filename, ".... Please Wait!");

    // CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry);
    if (WriteUserFormat (filetype, *mesh, /* *ng_geometry, */ filename))
      {
	ostringstream ost;
	ost << "Sorry, nothing known about file format " << filetype << endl;
	Tcl_SetResult (interp, (char*)ost.str().c_str(), TCL_VOLATILE);
	return TCL_ERROR;
      }

    PrintMessage (1, "Export mesh to file .... DONE!");
    return TCL_OK;
  }



  int Ng_ImportMesh (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    const string filename (argv[1]);
    PrintMessage (1, "import mesh from ", filename);

    mesh = make_shared<Mesh>();

    ReadFile (*mesh, filename);
    PrintMessage (2, mesh->GetNP(), " Points, ",
		  mesh->GetNE(), " Elements.");

    SetGlobalMesh (mesh);
    mesh->SetGlobalH (mparam.maxh);
    mesh->CalcLocalH(mparam.grading);

    return TCL_OK;
  }



  int Ng_ImportSolution (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    const char * filename = argv[1];
    PrintMessage (1, "Import solution from file ", filename);

    ImportSolution2 (filename);
    return TCL_OK;
  }



  static DemoView * demoview = 0;
  int Ng_ShowDemo (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    const char * filename = argv[1];
    PrintMessage (1, "Show demo ", filename);
    demoview = new DemoView (filename);

    return TCL_OK;
  }



  int Ng_DemoSetTime (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
    cout << "demosettime, time = " << argv[1] << endl;
    int result = -1;

    static char strminusone[] = "-1";
    static char str0[] = "0";

    if (demoview)
      result = demoview->SetTime (atof (argv[1]));

    if (result == -1)
      Tcl_SetResult (interp, strminusone, TCL_STATIC);
    else
      Tcl_SetResult (interp, str0, TCL_STATIC);

    return TCL_OK;
  }







  int Ng_SaveSolution (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    const char * filename = argv[1];
    PrintMessage (1, "Save solution to file ", filename);

    netgen::GetVSSolution().SaveSolutionData (filename);
    return TCL_OK;
  }




  int Ng_SetNextTimeStamp  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argqc, tcl_const char *argv[])
  {
    if (mesh)
      mesh -> SetNextTimeStamp();
    return TCL_OK;
  }




  int Ng_LoadGeometry (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    tcl_const char * lgfilename = argv[1];

#ifdef LOG_STREAM
    (*logout) << "Load geometry file: " << lgfilename << endl;
#endif

#ifdef STAT_STREAM
    (*statout) << lgfilename << " & " << endl;
#endif


    try
      {
	for (int i = 0; i < geometryregister.Size(); i++)
	  {
	    NetgenGeometry * hgeom = geometryregister[i]->Load (lgfilename);
	    if (hgeom)
	      {
                // delete ng_geometry;
		// ng_geometry = hgeom;
                ng_geometry = shared_ptr<NetgenGeometry> (hgeom);
                geometryregister[i]->SetParameters(interp);
		
		mesh.reset();
		return TCL_OK;
	      }
	  }


	ifstream infile(lgfilename);

	if (strlen(lgfilename) < 4)
	  {
	    cout << "ERROR: cannot recognise file format!" << endl;
	  }
	else
	  {
	    if ((strcmp (&lgfilename[strlen(lgfilename)-4], "iges") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "igs") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "IGS") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "IGES") == 0))
	      {
		Tcl_SetResult (interp, (char*)"IGES import requires the OpenCascade geometry kernel. "
			       "Please install OpenCascade as described in the Netgen-website",
			       TCL_STATIC);
		return TCL_ERROR;
	      }

	    else if (strcmp (&lgfilename[strlen(lgfilename)-3], "sat") == 0)
	      {
#ifdef ACIS
		PrintMessage (1, "Load ACIS geometry file ", lgfilename);
		acisgeometry = netgen::LoadACIS_SAT (lgfilename);
#endif
	      }
	    else if ((strcmp (&lgfilename[strlen(lgfilename)-4], "step") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "stp") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "STP") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "STEP") == 0))
	      {
#ifdef ACISxxx
		PrintMessage (1, "Load STEP geometry file ", lgfilename);
		acisgeometry = netgen::LoadACIS_STEP (lgfilename);
#else
		Tcl_SetResult (interp, (char*)"IGES import requires the OpenCascade geometry kernel. "
			       "Please install OpenCascade as described in the Netgen-website",
			       TCL_STATIC);
		return TCL_ERROR;
#endif
	      }
	    else if ((strcmp (&lgfilename[strlen(lgfilename)-4], "brep") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "Brep") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "BREP") == 0))
	      {
		Tcl_SetResult (interp, (char*)"BREP import requires the OpenCascade geometry kernel. "
			       "Please install OpenCascade as described in the Netgen-website",
			       TCL_STATIC);
		return TCL_ERROR;
	      }
	  }
      }

    catch (NgException e)
      {
	Tcl_SetResult (interp, const_cast<char*> (e.What().c_str()), TCL_VOLATILE);
	return TCL_ERROR;
      }

    mesh.reset();
    return TCL_OK;
  }










  int Ng_SaveGeometry (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (argc == 2)
      {
	const char * cfilename = argv[1];

	try
	  {
	    ng_geometry -> Save (string (cfilename));
	  }
	catch (NgException e)
	  {
	    Tcl_SetResult (interp, const_cast<char*> (e.What().c_str()), TCL_VOLATILE);
	    return TCL_ERROR;
	  }

	PrintMessage (1, "Save geometry to file ", cfilename);

	if (strlen(cfilename) < 4) {cout << "ERROR: can not recognise file format!!!" << endl;}
	else
	  {
#ifdef ACIS
	    if (acisgeometry)
	      {
		char * filename = const_cast<char*> (argv[1]);
		if (strcmp (&filename[strlen(filename)-3], "sat") == 0)
		  {
		    acisgeometry -> SaveSATFile (filename);
		  }
	      }
#endif
	    /*
	    if (strcmp (&cfilename[strlen(cfilename)-3], "ngg") == 0)
	      {
		CSGeometry * geometry = dynamic_cast<CSGeometry*> (ng_geometry);
		if (geometry)
		  {
		    ofstream of(cfilename);
		    geometry->Save (of);
		  }
	      }
	    */
	  }
      }

    return TCL_OK;
  }









  int Ng_ReadStatus (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    char buf[20], lstring[200];
    static int prev_np = -1;
    static int prev_ne = -1;
    static int prev_nse = -1;
    
    if (mesh)
      {
        if (prev_np != mesh->GetNP())
          {
            sprintf (buf, "%u", unsigned(mesh->GetNP()));
            Tcl_SetVar  (interp, "::status_np", buf, 0);
            prev_np = mesh->GetNP();
          }

        if (prev_ne != mesh->GetNE())
          {
            sprintf (buf, "%u", unsigned(mesh->GetNE()));
            Tcl_SetVar  (interp, "::status_ne", buf, 0);
            prev_ne = mesh->GetNE();
          }

        if (prev_nse != mesh->GetNSE())
          {
            sprintf (buf, "%u", unsigned(mesh->GetNSE()));
            Tcl_SetVar  (interp, "::status_nse", buf, 0);
            prev_nse = mesh->GetNSE();
          }
        
        auto tets_in_qualclass = mesh->GetQualityHistogram();
        lstring[0] = 0;
        for (int i = 0; i < tets_in_qualclass.Size(); i++)
          {
            sprintf (buf, " %d", tets_in_qualclass[i]);
            strcat (lstring, buf);
          }
        for (int i = tets_in_qualclass.Size(); i < 20; i++)
            strcat (lstring, " 0");
        Tcl_SetVar  (interp, "::status_tetqualclasses", lstring, 0);
      }
    else
      {
        if (prev_np != 0)
          {
            Tcl_SetVar  (interp, "::status_np", "0", 0);
            prev_np = 0;
          }

        if (prev_ne != 0)
          {
            Tcl_SetVar  (interp, "::status_ne", "0", 0);
            prev_ne = 0;
          }

        if (prev_nse != 0)
          {
            Tcl_SetVar  (interp, "::status_nse", "0", 0);
            prev_nse = 0;
          }
        Tcl_SetVar  (interp, "::status_tetqualclasses", "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0", 0);
      }

    static string prev_working;
    string working = multithread.running ? "working" : "       ";
    if (working != prev_working)
      {
        Tcl_SetVar (interp, "::status_working", working.c_str(), 0);
        prev_working = working;
      }

    /*
    if (multithread.running)
      Tcl_SetVar (interp, "::status_working", "working", 0);
    else
      Tcl_SetVar (interp, "::status_working", "       ", 0);
    */
    
    static string prev_task;
    if (prev_task != string(multithread.task))
      {
        prev_task = multithread.task;
        Tcl_SetVar (interp, "::status_task", prev_task.c_str(), 0);
      }

    static double prev_percent = -1;
    if (prev_percent != multithread.percent)
      {
        prev_percent = multithread.percent;
        sprintf (buf, "%lf", prev_percent);
        Tcl_SetVar  (interp, "::status_percent", buf, 0);
      }
        

    {
      lock_guard<mutex> guard(tcl_todo_mutex);
      if (multithread.tcl_todo->length())
        {
          Tcl_Eval (interp, multithread.tcl_todo->c_str());
          *multithread.tcl_todo = "";
        }
    }

    return TCL_OK;
  }


  int Ng_MemInfo (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {/*
    if (argc < 2) return TCL_ERROR;

    if (strcmp (argv[1], "usedmb") == 0)
      { // returns string of 512 '0' or '1'

	static char usedmb[513];
	for (int i = 0; i < 512; i++)
	  usedmb[i] = (i % 7 == 0) ? '1' : '0';

	usedmb[512] = 0;
	BaseDynamicMem::GetUsed (512, usedmb);
	Tcl_SetResult (interp, usedmb, TCL_STATIC);
	return TCL_OK;
      }
	*/
    return TCL_ERROR;
  }



  int Ng_BCProp (ClientData clientData,
		 Tcl_Interp * interp,
		 int argc, tcl_const char *argv[])
  {
    static char buf[100];

    if (argc < 2)
      {
	Tcl_SetResult (interp, (char*)"Ng_BCProp needs arguments", TCL_STATIC);
	return TCL_ERROR;
      }

    if (strcmp (argv[1], "setbc") == 0)
      {
	int facenr = atoi (argv[2]);
	int bcnr = atoi (argv[3]);
	if (mesh && facenr >= 1 && facenr <= mesh->GetNFD())
	  mesh->GetFaceDescriptor (facenr).SetBCProperty (bcnr);
      }

    if (strcmp (argv[1], "setall") == 0)
      {
	int bcnr = atoi (argv[2]);
	if (mesh)
	  {
	    int nfd = mesh->GetNFD();
	    for (int i = 1; i <= nfd; i++)
	      mesh->GetFaceDescriptor (i).SetBCProperty (bcnr);
	  }
      }

    if (strcmp (argv[1], "getbc") == 0)
      {
	int facenr = atoi (argv[2]);
	if (mesh && facenr >= 1 && facenr <= mesh->GetNFD())
	  {
	    sprintf (buf, "%d", mesh->GetFaceDescriptor(facenr).BCProperty());
	  }
	else
	  {
	    strcpy (buf, "0");
	  }
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }

    if (strcmp (argv[1], "getbcname") == 0)
      {
	int facenr = atoi (argv[2]);
	if (mesh && facenr >= 1 && facenr <= mesh->GetNFD())
	  {
	    sprintf (buf, "%s", mesh->GetFaceDescriptor(facenr).GetBCName().c_str());
	  }
	else
	  {
	    strcpy (buf, "-");
	  }
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }


    if (strcmp (argv[1], "getactive") == 0)
      {
	sprintf (buf, "%d", vsmesh.SelectedFace());
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }

    if (strcmp (argv[1], "setactive") == 0)
      {
	int facenr = atoi (argv[2]);
	if (mesh && facenr >= 1 && facenr <= mesh->GetNFD())
	  {
	    vsmesh.SetSelectedFace (facenr);
	  }
      }

    if (strcmp (argv[1], "getnfd") == 0)
      {
	if (mesh)
	  sprintf (buf, "%d", mesh->GetNFD());
	else
	  sprintf (buf, "0");
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }

    return TCL_OK;
  }





  int Ng_Refine  (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {
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

#ifdef ACIS
    if (acisgeometry)
      {
	ACISRefinementSurfaces ref (*acisgeometry);
	ACISMeshOptimize2dSurfaces opt(*acisgeometry);
	ref.Set2dOptimizer(&opt);
	ref.Refine (*mesh);
      }
    else
#endif
      {
	// ng_geometry -> GetRefinement().Refine(*mesh);
        mesh->GetGeometry()->GetRefinement().Refine(*mesh);
      }

//redo second order refinement if desired
    if (mparam.secondorder)
      const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement()).MakeSecondOrder(*mesh);

    return TCL_OK;
  }

  int Ng_SecondOrder  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
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
    
    const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement()).MakeSecondOrder (*mesh);

    return TCL_OK;
  }


  void * HighOrderDummy (void *)
  {
    //  mparam.elementorder = atoi (Tcl_GetVar (interp, "options.elementorder", 0));
    const char * savetask = multithread.task;

    Refinement & ref = const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement());
    mesh -> GetCurvedElements().BuildCurvedElements (&ref, mparam.elementorder);

    multithread.task = savetask;
    multithread.running = 0;
    multithread.terminate = 1;

    mesh -> SetNextMajorTimeStamp();
    return 0;
  }

  int Ng_HighOrder  (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
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

    multithread.running = 1;
    multithread.terminate = 0;

    mparam.elementorder = atoi(argv[1]);


    HighOrderDummy(NULL);

    return TCL_OK;
  }



  void * ValidateDummy (void *)
  {
    Refinement & ref = const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement());
    ref.ValidateSecondOrder (*mesh);

    multithread.running = 0;
    return NULL;
  }



  int Ng_ValidateSecondOrder  (ClientData clientData,
			       Tcl_Interp * interp,
			       int argc, tcl_const char *argv[])
  {
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

    multithread.running = 1;
    RunParallel (ValidateDummy, NULL);

    return TCL_OK;
  }


  int Ng_ZRefinement  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
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

    ZRefinementOptions opt;
    opt.minref = 5;

    if (argc >= 2) opt.minref = atoi (argv[1]);

    ZRefinement (*mesh, ng_geometry.get(), opt);

    return TCL_OK;
  }

  int Ng_HPRefinement  (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
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

    int levels = atoi (argv[1]);


    Refinement & ref = const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement());
    HPRefinement (*mesh, &ref, levels);
    return TCL_OK;
  }


  int Ng_LoadMeshSize  (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
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

    mesh->LoadLocalMeshSize(argv[1]);
    return TCL_OK;
  }


  int Ng_MeshSizeFromSurfaceMesh  (ClientData clientData,
				   Tcl_Interp * interp,
				   int argc, tcl_const char *argv[])
  {
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

    mesh->SetGlobalH (mparam.maxh);
    mesh->CalcLocalH(mparam.grading);

    return TCL_OK;
  }





  // Philippose Rajan - 13 June 2009
  // Added a new TCL function call for the generation 
  // of prismatic boundary layers
  int Ng_GenerateBoundaryLayer (ClientData clientData,
           Tcl_Interp * interp,
           int argc, tcl_const char *argv[])
  {
     if (!mesh)
     {
        Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
        return TCL_ERROR;
     }

     if(multithread.running)
     {
        Tcl_SetResult(interp, err_jobrunning, TCL_STATIC);
        return TCL_ERROR;
     }




     
     cout << "Generate Prismatic Boundary Layers (Experimental)...." << endl;
     
     // Use an array to support creation of boundary 
     // layers for multiple surfaces in the future...
     NgArray<int> surfid;
     int surfinp = 0;
     int prismlayers = 1;
     double hfirst = 0.01;
     double growthfactor = 1.0;

     
     while(surfinp >= 0)
       {
         cout << "Enter Surface ID (-1 to end list): ";
         cin >> surfinp;
         if(surfinp >= 0) surfid.Append(surfinp);
      }

     cout << "Number of surfaces entered = " << surfid.Size() << endl; 
     cout << "Selected surfaces are:" << endl;
     
     for(int i = 1; i <= surfid.Size(); i++)
       cout << "Surface " << i << ": " << surfid.Elem(i) << endl;
     
     cout << endl << "Enter number of prism layers: ";
     cin >> prismlayers;
     if(prismlayers < 1) prismlayers = 1;
     
     cout << "Enter height of first layer: ";
     cin >> hfirst;
     if(hfirst <= 0.0) hfirst = 0.01;
     
     cout << "Enter layer growth / shrink factor: ";
     cin >> growthfactor;
     if(growthfactor <= 0.0) growthfactor = 0.5;
     
     BoundaryLayerParameters blp;
     blp.surfid = surfid;
     blp.prismlayers = prismlayers;
     blp.hfirst = blp.hfirst;
     blp.growthfactor = growthfactor;
     GenerateBoundaryLayer (*mesh, blp);
     return TCL_OK;
  }


  int Ng_InsertVirtualBL (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
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

    InsertVirtualBoundaryLayer (*mesh);
    return TCL_OK;
  }

  int Ng_CutOffAndCombine (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    Mesh othermesh;
    othermesh.Load (argv[1]);
    othermesh.SetGlobalH (mparam.maxh);
    othermesh.CalcLocalH(mparam.grading);

    CutOffAndCombine (*mesh, othermesh);
    return TCL_OK;
  }


  int Ng_HelmholtzMesh (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    HelmholtzMesh (*mesh);
    return TCL_OK;
  }




  int Ng_SetMeshingParameters  (ClientData clientData,
				Tcl_Interp * interp,
				int argc, tcl_const char *argv[])
  {
    mparam.maxh = atof (Tcl_GetVar (interp, "::options.meshsize", 0));
    mparam.minh = atof (Tcl_GetVar (interp, "::options.minmeshsize", 0));

    mparam.meshsizefilename = Tcl_GetVar (interp, "::options.meshsizefilename", 0);
    // if (!strlen (mparam.meshsizefilename)) mparam.meshsizefilename = NULL;

    mparam.curvaturesafety = atof (Tcl_GetVar (interp, "::options.curvaturesafety", 0));
    mparam.segmentsperedge = atof (Tcl_GetVar (interp, "::options.segmentsperedge", 0));
    mparam.badellimit = atof (Tcl_GetVar (interp, "::options.badellimit", 0));
    mparam.secondorder = atoi (Tcl_GetVar (interp, "::options.secondorder", 0));
    mparam.elementorder = atoi (Tcl_GetVar (interp, "::options.elementorder", 0));
    mparam.quad = atoi (Tcl_GetVar (interp, "::options.quad", 0));
    mparam.try_hexes = atoi (Tcl_GetVar (interp, "::options.try_hexes", 0));

    mparam.inverttets = atoi (Tcl_GetVar (interp, "::options.inverttets", 0));
    mparam.inverttrigs = atoi (Tcl_GetVar (interp, "::options.inverttrigs", 0));
    mparam.uselocalh = atoi (Tcl_GetVar (interp, "::options.localh", 0));
    mparam.grading = atof (Tcl_GetVar (interp, "::options.grading", 0));
    mparam.delaunay = atoi (Tcl_GetVar (interp, "::options.delaunay", 0));
    mparam.checkoverlap = atoi (Tcl_GetVar (interp, "::options.checkoverlap", 0));
    mparam.checkoverlappingboundary = atoi (Tcl_GetVar (interp, "::options.checkoverlappingboundary", 0));
    mparam.checkchartboundary = atoi (Tcl_GetVar (interp, "::options.checkchartboundary", 0));
    mparam.optsteps3d = atoi (Tcl_GetVar (interp, "::options.optsteps3d", 0));
    mparam.optsteps2d = atoi (Tcl_GetVar (interp, "::options.optsteps2d", 0));
    mparam.opterrpow = atof (Tcl_GetVar (interp, "::options.opterrpow", 0));

    mparam.parthread = atoi (Tcl_GetVar (interp, "::options.parthread", 0));
    mparam.elsizeweight = atof (Tcl_GetVar (interp, "::options.elsizeweight", 0));

    mparam.autozrefine = atoi (Tcl_GetVar (interp, "::options.autozrefine", 0));

    extern int printmessage_importance;
    extern int printdots;
    printmessage_importance = atoi (Tcl_GetVar (interp, "::options.printmsg", 0));
    printdots = (printmessage_importance >= 4);

    mparam.parallel_meshing = atoi (Tcl_GetVar (interp, "::options.parallel_meshing", 0));
    mparam.nthreads = atoi (Tcl_GetVar (interp, "::options.nthreads", 0));
    if(atoi(Tcl_GetVar (interp, "::stloptions.resthcloseedgeenable", 0)))
      mparam.closeedgefac = atof(Tcl_GetVar (interp, "::stloptions.resthcloseedgefac", 0));
    else
      mparam.closeedgefac = {};

    //BaseMoveableMem::totalsize = 0;
    // 1048576 * atoi (Tcl_GetVar (interp, "::options.memory", 0));
    if (mesh)
      {
	mesh->SetGlobalH (mparam.maxh);
	mesh->SetMinimalH (mparam.minh);
      }

#ifdef PARALLELGL
    MyMPI_SendCmd ("bcastparthread");
    MyMPI_Bcast (mparam.parthread, MPI_COMM_WORLD);
#endif

    return TCL_OK;
  }



  int Ng_SetDebugParameters  (ClientData clientData,
			      Tcl_Interp * interp,
			      int argc, tcl_const char *argv[])
  {
    debugparam.slowchecks = atoi (Tcl_GetVar (interp, "::debug.slowchecks", 0));
    debugparam.debugoutput = atoi (Tcl_GetVar (interp, "::debug.debugoutput", 0));
    debugparam.haltexistingline = atoi (Tcl_GetVar (interp, "::debug.haltexistingline", 0));
    debugparam.haltoverlap = atoi (Tcl_GetVar (interp, "::debug.haltoverlap", 0));
    debugparam.haltsuccess = atoi (Tcl_GetVar (interp, "::debug.haltsuccess", 0));
    debugparam.haltnosuccess = atoi (Tcl_GetVar (interp, "::debug.haltnosuccess", 0));
    debugparam.haltlargequalclass = atoi (Tcl_GetVar (interp, "::debug.haltlargequalclass", 0));
    debugparam.haltsegment = atoi (Tcl_GetVar (interp, "::debug.haltsegment", 0));
    debugparam.haltnode = atoi (Tcl_GetVar (interp, "::debug.haltnode", 0));
    debugparam.haltface = atoi (Tcl_GetVar (interp, "::debug.haltface", 0));
    debugparam.haltsegmentp1 = atoi (Tcl_GetVar (interp, "::debug.haltsegmentp1", 0));
    debugparam.haltsegmentp2 = atoi (Tcl_GetVar (interp, "::debug.haltsegmentp2", 0));
    debugparam.haltfacenr = atoi (Tcl_GetVar (interp, "::debug.haltfacenr", 0));
    return TCL_OK;
  }



  int Ng_SetCommandLineParameter  (ClientData clientData,
				   Tcl_Interp * interp,
				   int argc, tcl_const char *argv[])
  {
    if (argc != 2)
      {
	Tcl_SetResult (interp, (char*)"Ng_SetCommandLineParameter needs 1 parameter",
                       TCL_STATIC);
	return TCL_ERROR;
      }

    if (argv[1][0] == '-')
      parameters.SetCommandLineFlag (argv[1]);
    else
      {
        if (strstr(argv[1], ".py"))
          parameters.SetFlag ("py", argv[1]);
        else
          parameters.SetFlag ("geofile", argv[1]);
      }
    return TCL_OK;
  }


  int Ng_GetCommandLineParameter  (ClientData clientData,
				   Tcl_Interp * interp,
				   int argc, tcl_const char *argv[])
  {
    if (argc != 2)
      {
	Tcl_SetResult (interp, (char*)"Ng_GetCommandLineParameter needs 1 parameter",
                       TCL_STATIC);
	return TCL_ERROR;
      }

    static char buf[10];

    if (parameters.StringFlagDefined (argv[1]))
        Tcl_SetResult (interp,
                       const_cast<char*>(parameters.GetStringFlag (argv[1], NULL).c_str()), TCL_VOLATILE);
    else if (parameters.NumFlagDefined (argv[1]))
      {
	sprintf (buf, "%lf", parameters.GetNumFlag (argv[1], 0));
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }
    else if (parameters.GetDefineFlag (argv[1]))
      Tcl_SetResult (interp, (char*)"defined", TCL_STATIC);
    else
      Tcl_SetResult (interp, (char*)"undefined", TCL_STATIC);

    return TCL_OK;
  }


  static int perfstepsstart;
  static int perfstepsend;

  static char* optstring = NULL;
  static char* optstringcsg = NULL;

  void * MeshingDummy (void *)
  {

    const char * savetask = multithread.task;
    multithread.task = "Generate Mesh";

    ResetTime();


    try
      {

#ifdef LOG_STREAM
	(*logout) << "Start meshing" << endl;
	(*logout) << "Meshing parameters:" << endl;
	mparam.Print (*logout);
#endif

#ifdef ACIS
	if (acisgeometry)
	  {
	    ACISGenerateMesh(*acisgeometry, mesh.Ptr(), perfstepsstart, perfstepsend, optstring);
	  }
	else
#endif
          if (ng_geometry)
	    {
              if (perfstepsstart == 1)
                {
                  mesh = make_shared<Mesh> ();
                  // vsmesh.SetMesh (mesh);
                  SetGlobalMesh (mesh);
                  mesh -> SetGeometry(ng_geometry);
                }
              if(!mesh)
                throw Exception("Need existing global mesh");
              mparam.perfstepsstart = perfstepsstart;
	      mparam.perfstepsend = perfstepsend;
              if(optstring)
                mparam.optimize3d = *optstring;
              int res = ng_geometry -> GenerateMesh (mesh, mparam);

	      if (res != MESHING3_OK) 
		{
		  multithread.task = savetask;
		  multithread.running = 0;
		  return 0;
		}
	    }
          else if (mesh)
            {
              if(perfstepsstart > 1 && perfstepsstart < 5)
                throw Exception("Need geometry for surface mesh operations!");
              MeshVolume(mparam, *mesh);
              OptimizeVolume(mparam, *mesh);
              return 0;
            }
          else // no ng_geometry
            {
              multithread.task = savetask;
              multithread.running = 0;
              return 0;
            }



	if (mparam.autozrefine)
	  {
	    ZRefinementOptions opt;
	    opt.minref = 5;
	    ZRefinement (*mesh, ng_geometry.get(), opt);
	    mesh -> SetNextMajorTimeStamp();
	  }
	
	if (mparam.secondorder)
	  {
	    const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement()).MakeSecondOrder (*mesh);
	    mesh -> SetNextMajorTimeStamp();
	  }

	if (mparam.elementorder > 1)
	  {
	    mesh -> GetCurvedElements().BuildCurvedElements (&const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement()),
							     mparam.elementorder);

	    mesh -> SetNextMajorTimeStamp();
	  }


	PrintMessage (1, "Meshing done, time = ", GetTime(), " sec");
      }

    catch (NgException e)
      {
	cout << e.What() << endl;
      }

    multithread.task = savetask;
    multithread.running = 0;

#ifdef OCCGEOMETRYorig
    // currently not active
    OCCGeometry * occgeometry = dynamic_cast<OCCGeometry*> (ng_geometry);
    if (occgeometry && occgeometry->ErrorInSurfaceMeshing())
      {
	char script[] = "rebuildoccdialog";
	Tcl_GlobalEval (tcl_interp, script);
      }
#endif
    return NULL;
  }

  
  int MeshingVal(tcl_const char* str)
  {
    if (strcmp(str, "ag") == 0) {return MESHCONST_ANALYSE;}
    if (strcmp(str, "me") == 0) {return MESHCONST_MESHEDGES;}
    if (strcmp(str, "ms") == 0) {return MESHCONST_MESHSURFACE;}
    if (strcmp(str, "os") == 0) {return MESHCONST_OPTSURFACE;}
    if (strcmp(str, "mv") == 0) {return MESHCONST_MESHVOLUME;}
    if (strcmp(str, "ov") == 0) {return MESHCONST_OPTVOLUME;}

    cout << "TCL TK ERROR, wrong meshing value, return='" << str << "'" << endl;
    return 0;
  }



  int Ng_GenerateMesh  (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    multithread.running = 1;
    multithread.terminate = 0;
    
    extern void Render(bool blocking);
    mparam.render_function = &Render;

    for (int i = 0; i < geometryregister.Size(); i++)
      geometryregister[i] -> SetParameters (interp);


    Ng_SetMeshingParameters (clientData, interp, 0, argv);

    perfstepsstart = 1;
    perfstepsend = 6;

    if (optstringcsg) delete optstringcsg;
    optstringcsg = NULL;
    if (optstring) delete optstring;
    optstring = NULL;

    if (argc == 2)
      {
	perfstepsstart = 1;
	perfstepsend = MeshingVal(argv[1]);
      }
    else if (argc == 3)
      {
	perfstepsstart = MeshingVal(argv[1]);
	perfstepsend = MeshingVal(argv[2]);
      }
    else if (argc == 4)
      {
	perfstepsstart = MeshingVal(argv[1]);
	perfstepsend = MeshingVal(argv[2]);
	optstring = new char[strlen(argv[3])+1];
	strcpy(optstring, argv[3]);
	optstringcsg = new char[strlen(argv[3])+1];
	strcpy(optstringcsg, argv[3]);
      }

    RunParallel (MeshingDummy, NULL);

    return TCL_OK;
  }


  int Ng_StopMeshing  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    multithread.terminate = 1;
    return TCL_OK;
  }

  int Ng_MeshInfo  (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    ostringstream str;

    if (argc >= 2 && strcmp (argv[1], "dim") == 0)
      str << mesh->GetDimension();
    else if (argc >= 2 && strcmp (argv[1], "np") == 0)
      str << mesh->GetNP();
    else if (argc >= 2 && strcmp (argv[1], "ne") == 0)
      str << mesh->GetNE();
    else if (argc >= 2 && strcmp (argv[1], "nse") == 0)
      str << mesh->GetNSE();
    else if (argc >= 2 && strcmp (argv[1], "nseg") == 0)
      str << mesh->GetNSeg();
    else if (argc >= 2 && strcmp (argv[1], "bbox") == 0)
      {
	Point3d pmin, pmax;
	mesh->GetBox (pmin, pmax);
	str << pmin.X() << " " << pmax.X() << " "
	    << pmin.Y() << " " << pmax.Y() << " "
	    << pmin.Z() << " " << pmax.Z() << endl;
      }
    else
      {
	cout << "argv[1] = " << argv[1] << endl;
	Tcl_SetResult (interp, (char*)"Ng_MeshInfo requires an argument out of \n dim np ne", TCL_STATIC);
	return TCL_ERROR;
      }

    Tcl_SetResult  (interp, (char*)str.str().c_str(), TCL_VOLATILE);
    return TCL_OK;
  }

  int Ng_MeshQuality  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
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

    double angles[4];
    char buf[10];

    if (mesh)
      mesh->CalcMinMaxAngle(mparam.badellimit, angles);
    else
      {
	angles[0] = angles[1] = angles[2] = angles[3] = 0;
      }
    sprintf (buf, "%5.1lf", angles[0]);
    Tcl_SetVar (interp, argv[1], buf, 0);
    sprintf (buf, "%5.1lf", angles[1]);
    Tcl_SetVar (interp, argv[2], buf, 0);
    sprintf (buf, "%5.1lf", angles[2]);
    Tcl_SetVar (interp, argv[3], buf, 0);
    sprintf (buf, "%5.1lf", angles[3]);
    Tcl_SetVar (interp, argv[4], buf, 0);

    return TCL_OK;
  }

  int Ng_CheckSurfaceMesh  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
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

    mesh->FindOpenElements();
    if (mesh->CheckConsistentBoundary())
      {
	PrintMessage (1, "surface mesh not consistent, trying orientation");
	mesh->SurfaceMeshOrientation();
      }
    else
      {
	PrintMessage (1, "surface mesh consistent");
      }

    mesh->CheckOverlappingBoundary();
    return TCL_OK;
  }

  int Ng_CheckVolumeMesh  (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
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

    mesh->CheckVolumeMesh();
    return TCL_OK;
  }


  int Ng_DeleteVolMesh  (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    if (mesh)
      mesh->ClearVolumeElements();

    return TCL_OK;
  }


  int Ng_SplitSeparatedFaces (ClientData clientData,
			      Tcl_Interp * interp,
			      int argc, tcl_const char *argv[])
  {
    if (mesh)
      mesh->SplitSeparatedFaces ();
    return TCL_OK;
  }



  int Ng_RestrictH  (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
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

    if (argc != 3)
      return TCL_OK;
    if (!mesh)
      return TCL_OK;


    double loch = atof (argv[2]);
    if (strcmp (argv[1], "face") == 0)
      {
	cout << "Restrict h at face to " << loch << endl;
	mesh -> RestrictLocalH  (RESTRICTH_FACE, vsmesh.SelectedFace(), loch);
      }
    if (strcmp (argv[1], "edge") == 0)
      {
	cout << "Restrict h at edge to " << loch << endl;
	mesh -> RestrictLocalH  (RESTRICTH_EDGE, vsmesh.SelectedEdge(), loch);
      }
    if (strcmp (argv[1], "point") == 0)
      {
	cout << "Restrict h at point to " << loch << endl;
	mesh -> RestrictLocalH  (RESTRICTH_POINT, vsmesh.SelectedPoint(), loch);
      }

    return TCL_OK;
  }





  int Ng_Anisotropy  (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
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

    if (argc != 2)
      return TCL_OK;
    if (!mesh)
      return TCL_OK;

    if (strcmp (argv[1], "edge") == 0)
      {
	int edgenr = vsmesh.SelectedEdge();
	for (int i = 1; i <= mesh->GetNSeg(); i++)
	  {
	    Segment & seg = mesh->LineSegment(i);
	    if (seg.edgenr == edgenr)
	      {
		seg.singedge_left = 1 - seg.singedge_left;
		seg.singedge_right = 1 - seg.singedge_right;
	      }
	  }
      }

    return TCL_OK;
  }




  BisectionOptions biopt;

  void * BisectDummy (void *)
  {
    const Refinement & ref = mesh->GetGeometry()->GetRefinement();
    MeshOptimize2d * opt = NULL;

    /*
#ifdef ACIS
    if (acisgeometry)
      {
	// ref = new ACISRefinementSurfaces(*acisgeometry);
	opt = new ACISMeshOptimize2dSurfaces(*acisgeometry);
	ref->Set2dOptimizer(opt);
      }
#endif
    else
      {
	ref = new RefinementSurfaces(*geometry);
	opt = new MeshOptimize2dSurfaces(*geometry);
	ref->Set2dOptimizer(opt);
      }
    */

    if(!mesh->LocalHFunctionGenerated())
      mesh->CalcLocalH(mparam.grading);

    mesh->LocalHFunction().SetGrading (mparam.grading);
    ref . Bisect (*mesh, biopt);
    mesh -> UpdateTopology();
    mesh -> GetCurvedElements().BuildCurvedElements (&ref, mparam.elementorder);

    multithread.running = 0;

    delete opt;
    return NULL;
  }


  int Ng_Bisect  (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {
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


    if (multithread.running)
      {
	cout << "Thread alrad running" << endl;
	return TCL_OK;
      }
    multithread.running = 1;

    biopt.outfilename = NULL; // "ngfepp.vol";
    biopt.femcode = "fepp";
    biopt.refinementfilename = NULL;

    if (argc >= 2)
      biopt.refinementfilename = argv[1];


    BisectDummy (0);

    /*
      extern void BisectTets (Mesh &, const CSGeometry *);
      BisectTets (*mesh, geometry);
    */
    return TCL_OK;
  }



  //   int Ng_BisectCopyMesh  (ClientData clientData,
  // 			  Tcl_Interp * interp,
  // 			  int argc, tcl_const char *argv[])
  //   {
  //     if (!mesh)
  //       {
  // 	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
  // 	return TCL_ERROR;
  //       }
  //     if (multithread.running)
  //       {
  // 	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
  // 	return TCL_ERROR;
  //       }

  //     BisectTetsCopyMesh (*mesh, geometry.Ptr(), biopt);
  //     return TCL_OK;
  //   }






  int Ng_Split2Tets  (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
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

    mesh->Split2Tets ();

    return TCL_OK;
  }








  extern int Ng_MeshDoctor (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[]);







  SymbolTable<VisualScene*> & GetVisualizationScenes ()
  {
    static SymbolTable<VisualScene*> vss;
    return vss;
  }

  void AddVisualizationScene (const string & name,
			      VisualScene * avs)
  {
    GetVisualizationScenes().Set (name.c_str(), avs);
  }


  void SetVisualScene (Tcl_Interp * interp)
  {
    const char * vismode = vispar.selectvisual;
    // Tcl_GetVar (interp, "selectvisual", 0);
    vs = &vscross;
    if (GetVisualizationScenes().Used(vismode))
      {
	vs = GetVisualizationScenes()[vismode];
      }
    else if (vismode)
      {
	if (strcmp (vismode, "geometry") == 0)
	  {
	    for (int i = 0; i < geometryregister.Size(); i++)
	      {
		VisualScene * hvs = geometryregister[i]->GetVisualScene (ng_geometry.get());
		if (hvs)
		  {
		    vs = hvs;
		    return;
		  }
	      }

#ifdef ACIS
	    else if (acisgeometry)
	      vs = &vsacisgeom;
#endif // ACIS
	  }
	
	if (strcmp (vismode, "mesh") == 0)
	  {
	    if (!meshdoctor.active)
	      vs = &vsmesh;
	    else
	      vs = &vsmeshdoc;
	  }

	if (strcmp (vismode, "surfmeshing") == 0) vs = &vssurfacemeshing;
	if (strcmp (vismode, "specpoints") == 0) vs = &vsspecpoints;
        if (strcmp (vismode, "solution") == 0) vs = &netgen::GetVSSolution();
      }
  }



  Font * font = nullptr;
  Togl * togl = NULL;

  void MyOpenGLText_GUI (const char * text)
  {
    glListBase (font->getDisplayListsBase());
    glCallLists (GLsizei(strlen(text)), GL_UNSIGNED_BYTE, text);
  }

  static int
  Ng_ToglVersion(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    Tcl_SetResult (interp,  (char*)"2", TCL_STATIC);
    return TCL_OK;
  }
  
  static int
  init(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    // cout << "call init" << endl;


    if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) 
      return TCL_ERROR;

    // possible values: 12,14,16,18,20,22,24,28,32
    font = selectFont(18);

    LoadOpenGLFunctionPointers();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);

    SetVisualScene (Togl_Interp(togl));
    vs->DrawScene();
    Set_OpenGLText_Callback (&MyOpenGLText_GUI);
    return TCL_OK;
  }

  static int
  zap(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    return TCL_OK;
  }


  static int
  draw(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    SetVisualScene (interp);

    glPushMatrix();

    glLoadIdentity();
    vs->DrawScene();
    Togl_SwapBuffers(togl);

    glPopMatrix();

    return TCL_OK;
  }

  static int
  reshape(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    int w = Togl_Width (togl);
    int h = Togl_Height (togl);

    // glViewport(0, 0, w, h);
    int res[4];
    glGetIntegerv(GL_VIEWPORT, res);
    // cout << "w = " << w << " h = " << h << endl;
    w = res[2];
    h = res[3];
    /*
    cout << "viewport: "
         << res[0] << " "
         << res[1] << " "
         << res[2] << " "
         << res[3] << endl;
    */
    // change font size according to window width
    font = selectFont(w/80);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // OpenGL near and far clipping planes
    double pnear = 0.1;
    double pfar = 10;

    gluPerspective(20.0f, double(w) / h, pnear, pfar);
    glMatrixMode(GL_MODELVIEW);

    return TCL_OK;
  }

  static int Ng_SnapShot(ClientData clientData, Tcl_Interp *interp, int argc, Tcl_Obj *const *argv)
  {
    struct Togl *togl;
    if (Togl_GetToglFromObj(interp, argv[1], &togl) != TCL_OK) 
      return TCL_ERROR;
    const char * filename = Tcl_GetString(argv[2]);

    int len = strlen(filename);
    int w = Togl_PixelScale(togl)*Togl_Width (togl);
    int h = Togl_PixelScale(togl)*Togl_Height (togl);

    NgArray<unsigned char> buffer(w*h*3);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadPixels (0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &buffer[0]);

#ifdef JPEGLIB
    if (strcmp ("jpg", filename+len-3) == 0)
      {
        cout << "Snapshot to file '" << filename << "'" << endl;

        struct jpeg_compress_struct cinfo;
        struct jpeg_error_mgr jerr;
        FILE *outfile = fopen(filename,"wb");
        JSAMPROW row_pointer[1];
        int row_stride, quality = 100; // 1...100

        cinfo.err = jpeg_std_error( &jerr );
        jpeg_create_compress( &cinfo );
        jpeg_stdio_dest( &cinfo, outfile );


        cinfo.image_width = w;
        cinfo.image_height = h;
        cinfo.input_components = 3;
        cinfo.in_color_space = JCS_RGB;

        jpeg_set_defaults( &cinfo );
        jpeg_set_quality( &cinfo, quality, FALSE );   // TRUE
        jpeg_start_compress( &cinfo, TRUE );

        row_stride = 3*w;
        while( cinfo.next_scanline < cinfo.image_height ) {
          row_pointer[0] = &buffer[ (h-1-cinfo.next_scanline) * row_stride ];
          (void)jpeg_write_scanlines( &cinfo, row_pointer, 1 );
        }

        jpeg_finish_compress( &cinfo );
        fclose( outfile );

        jpeg_destroy_compress( &cinfo );
        fprintf( stdout, "done [ok]\n" );
        fflush( stdout );

        return TCL_OK;
      }
#endif // JPEGLIB
    {
        string command;
        string filename2;

        filename2 = filename;

        if(filename2.substr(len-3) != ".ppm")
            filename2 += ".ppm";

        cout << "Snapshot to file '" << filename << endl;

        // int w = Togl_Width (togl);
        // int h = Togl_Height (togl);

        ofstream outfile(filename2);
        outfile << "P6" << endl
          << "# CREATOR: Netgen" << endl
          << w << " " << h << endl
          << "255" << endl;
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j++)
                for (int k = 0; k < 3; k++)
                    outfile.put (buffer[k+3*j+3*w*(h-i-1)]);
        outfile << flush;

        if (filename2 == string(filename))
            return TCL_OK;
        else
        {
          // convert image file (Unix/Linux only):
          command = string("convert -quality 100 ") + filename2 + " " + filename;

          int err = system(command.c_str());
          if (err != 0)
          {
              Tcl_SetResult (Togl_Interp(togl), (char*)"Cannot convert image file, stored as .ppm", TCL_VOLATILE);
              return TCL_ERROR;
          }

          command  = string("rm ") + filename2;
          err = system(command.c_str());

          if (err != 0)
          {
              Tcl_SetResult (Togl_Interp(togl), (char*)"Cannot delete temporary file", TCL_VOLATILE);
              return TCL_ERROR;
          }
          return TCL_OK;
        }
    }
  }



#ifdef FFMPEG
  static int Ng_VideoClip(ClientData clientData, Tcl_Interp *interp, int argc, Tcl_Obj *const *argv)
    {
      static Mpeg mpeg;
      struct Togl *togl;
    if (Togl_GetToglFromObj(interp, argv[1], &togl) != TCL_OK) 
      return TCL_ERROR;

    if (strcmp (Tcl_GetString(argv[2]), "init") == 0)
      {
        // Can't initialize when running:
        //-------------------------------
        if( mpeg.IsStarted() )  {
          cout << "cannot initialize: already running" << endl;
          return TCL_ERROR;
        }

        const char * filename = Tcl_GetString(argv[3]);
        mpeg.Start(filename);

        return TCL_OK;
      }



    else if (strcmp (Tcl_GetString(argv[2]), "addframe") == 0)
      {
        if(mpeg.AddFrame())
          return TCL_ERROR;
      }



    else if (strcmp (Tcl_GetString(argv[2]), "finalize") == 0)
      {
        mpeg.Stop();
      }
    return TCL_OK;
  }

#else // FFMPEG
  static int Ng_VideoClip(ClientData clientData, Tcl_Interp *interp, int argc, Tcl_Obj *const *argv)
  {
    Tcl_SetResult (Togl_Interp(togl), (char*)"Video not available, Netgen was not compiled with FFMPEG library", TCL_STATIC);
    return TCL_ERROR;
  }
#endif // FFMPEG







  int Ng_MouseMove (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    int oldx, oldy;
    int newx, newy;

    oldx = atoi (argv[1]);
    oldy = atoi (argv[2]);
    newx = atoi (argv[3]);
    newy = atoi (argv[4]);

    SetVisualScene(interp);
    vs->MouseMove (oldx, oldy, newx, newy, argv[5][0]);

    return TCL_OK;
  }


  int Ng_MouseDblClick (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    int px = Togl_PixelScale(togl)*atoi (argv[1]);
    int py = Togl_PixelScale(togl)*atoi (argv[2]);

    SetVisualScene(interp);
    vs->MouseDblClick (px, py);

    return TCL_OK;
  }


  int Ng_ZoomAll (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    vs->BuildScene (1);

    return TCL_OK;
  }


  int Ng_Center (ClientData clientData,
		 Tcl_Interp * interp,
		 int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    vs->BuildScene (2);

    return TCL_OK;
  }


  int Ng_StandardRotation (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    vs->StandardRotation (argv[1]);

    return TCL_OK;
  }

  int Ng_ArbitraryRotation (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    NgArray<double> alpha;
    NgArray<Vec3d> vec;

    for(int i=1; i<argc; i+=4)
      {
	alpha.Append(atof(argv[i]));
	vec.Append(Vec3d(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3])));
      }

    vs->ArbitraryRotation (alpha,vec);

    return TCL_OK;
  }



  int Ng_Metis (ClientData clientData,
		Tcl_Interp * interp,
		int argc, tcl_const char *argv[])
  {
#ifdef PARALLEL
    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    int nparts = atoi (argv[1]);
    ntasks = nparts+1;
    cout << "calling metis ... " << flush;
    mesh->ParallelMetis();
    cout << "done" << endl;
    ntasks = 1;

    // for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
    // (*mesh)[ei].SetIndex ( (*mesh)[ei].GetPartition() );

    return TCL_OK;

#else
    Tcl_SetResult (interp, (char*)"metis not available", TCL_STATIC);
    return TCL_ERROR;
    
#endif



#ifdef OLDOLD
    // METIS Partitioning

    if (mesh->GetDimension() == 3)
      {
	using namespace metis;

	int ne = mesh->GetNE();
	if (ne < 3)
	  {
	    Tcl_SetResult (interp, "This operation needs a volume mesh", TCL_STATIC);
	    return TCL_ERROR;
	  }

	int nn = mesh->GetNP();

	ELEMENT_TYPE elementtype = mesh->VolumeElement(1).GetType();
	int npe = mesh->VolumeElement(1).GetNP();

	for (int i = 2; i<=ne; i++)
	  if (mesh->VolumeElement(i).GetType() != elementtype)
	    {
	      Tcl_SetResult (interp, "Works in 3D only uniformal tet or hex meshes", TCL_STATIC);
	      return TCL_ERROR;
	    }

	idxtype *elmnts;
	elmnts = new idxtype[ne*npe];

	int etype;
	if (elementtype == TET)
	  etype = 2;
	else if (elementtype == HEX)
	  etype = 3;
	else
	  {
	    Tcl_SetResult (interp, "Works in 3D only uniformal tet or hex meshes", TCL_STATIC);
	    return TCL_ERROR;
	  }

	for (int i=1; i<=ne; i++)
	  for (int j=1; j<=npe; j++)
	    elmnts[(i-1)*npe+(j-1)] = mesh->VolumeElement(i).PNum(j)-1;


	int numflag = 0;
	int nparts = atoi (argv[1]);
	int edgecut;
	idxtype *epart, *npart;
	epart = new idxtype[ne];
	npart = new idxtype[nn];

	cout << "Starting Metis (" << ne << " Elements, " << nn << " Nodes, " << nparts << " Partitions) ... " << flush;

	METIS_PartMeshNodal (&ne, &nn, elmnts, &etype, &numflag, &nparts,
			     &edgecut, epart, npart);

	cout << "done" << endl;

	cout << "edge-cut: " << edgecut << ", balance: " << ComputeElementBalance(ne, nparts, epart) << endl;

	for (int i=1; i<=ne; i++)
	  mesh->VolumeElement(i).SetPartition(epart[i-1]);

	mesh->SetNextTimeStamp();
      }


#endif
    return TCL_OK;
  }





  void SelectFaceInOCCDialogTree (int facenr)
  {
    char script[50];
    sprintf (script, "selectentity {Face %i}", facenr);
    Tcl_GlobalEval (tcl_interp, script);
  }





#ifndef ACIS
  int Ng_ACISCommand (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
    if (argc >= 2)
      {
	if (strcmp (argv[1], "isACISavailable") == 0)
	  {
            Tcl_SetResult (interp, (char*)"no", TCL_STATIC);
	    return TCL_OK;
	  }
      }
    Tcl_SetResult (interp, (char*)"undefined ACiS command", TCL_STATIC);
    return TCL_ERROR;
  }
#endif


  // from ng_interface
  void Ng_SetVisualizationParameter (const char * name, const char * value)
  {
    // #ifdef OPENGL
    // #ifndef NOTCL
    char buf[100];
    sprintf (buf, "visoptions.%s", name);
    if (printmessage_importance>0)
      {
	cout << "name = " << name << ", value = " << value << endl;
	cout << "set tcl-variable " << buf << " to " << value << endl;
      }
    Tcl_SetVar (tcl_interp, buf, const_cast<char*> (value), 0);
    Tcl_Eval (tcl_interp, "Ng_Vis_Set parameters;");
    // #endif
    // #endif
  }
}


using namespace netgen;

void Ng_SetMouseEventHandler (netgen::MouseEventHandler * handler)
{
  vsmesh.SetMouseEventHandler (handler);
}

void Ng_SetUserVisualizationObject (netgen::UserVisualizationObject * vis)
{
  netgen::GetVSSolution().AddUserVisualizationObject (vis);
}




namespace netgen
{


int firsttime = 1;
int animcnt = 0;
void PlayAnimFile(const char* name, int speed, int maxcnt)
{
  //extern Mesh * mesh;

  /*
    if (mesh) mesh->DeleteMesh();
    if (!mesh) mesh = new Mesh();
  */
  mesh = make_shared<Mesh>();

  int ne, np, i;

  char str[80];
  char str2[80];

  //int tend = 5000;
  //  for (ti = 1; ti <= tend; ti++)
  //{
  int rti = (animcnt%(maxcnt-1)) + 1;
  animcnt+=speed;
  
  sprintf(str2,"%05i.sol",rti);
  strcpy(str,"mbssol/");
  strcat(str,name);
  strcat(str,str2);

  if (printmessage_importance>0)
    cout << "read file '" << str << "'" << endl;
  
  ifstream infile(str);
  infile >> ne;
  for (i = 1; i <= ne; i++)
    {
      int j;
      Element2d tri(TRIG);
      tri.SetIndex(1); //faceind
      
      for (j = 1; j <= 3; j++)
	infile >> tri.PNum(j);

      infile >> np;
      for (i = 1; i <= np; i++)
	{
	  Point3d p;
	  infile >> p.X() >> p.Y() >> p.Z();
	  if (firsttime)
	    mesh->AddPoint (p);
	  else
	    mesh->Point(i) = Point<3> (p);
	}

      //firsttime = 0;
      Ng_Redraw();
    }
}




  int Ng_SetVisParameters  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    if (!Tcl_GetVar (interp, "::viewoptions.light.amb", TCL_GLOBAL_ONLY))
      return TCL_ERROR;

    vispar.lightamb = atof (Tcl_GetVar (interp, "::viewoptions.light.amb", TCL_GLOBAL_ONLY));
    vispar.lightdiff = atof (Tcl_GetVar (interp, "::viewoptions.light.diff", TCL_GLOBAL_ONLY));
    vispar.lightspec = atof (Tcl_GetVar (interp, "::viewoptions.light.spec", TCL_GLOBAL_ONLY));
    vispar.shininess = atof (Tcl_GetVar (interp, "::viewoptions.mat.shininess", TCL_GLOBAL_ONLY));
    vispar.locviewer = atoi (Tcl_GetVar (interp, "::viewoptions.light.locviewer", TCL_GLOBAL_ONLY));
    vispar.transp = atof (Tcl_GetVar (interp, "::viewoptions.mat.transp", TCL_GLOBAL_ONLY));
    
    VisualizationParameters::Clipping hclip;
    hclip.normal.X() = atof (Tcl_GetVar (interp, "::viewoptions.clipping.nx", TCL_GLOBAL_ONLY));
    hclip.normal.Y() = atof (Tcl_GetVar (interp, "::viewoptions.clipping.ny", TCL_GLOBAL_ONLY));
    hclip.normal.Z() = atof (Tcl_GetVar (interp, "::viewoptions.clipping.nz", TCL_GLOBAL_ONLY));
    hclip.dist = atof (Tcl_GetVar (interp, "::viewoptions.clipping.dist", TCL_GLOBAL_ONLY));
    hclip.dist2 = atof (Tcl_GetVar (interp, "::viewoptions.clipping.dist2", TCL_GLOBAL_ONLY));
    hclip.enable = atoi (Tcl_GetVar (interp, "::viewoptions.clipping.enable", TCL_GLOBAL_ONLY));
    vispar.clipdomain =
      atoi (Tcl_GetVar (interp, "::viewoptions.clipping.onlydomain", TCL_GLOBAL_ONLY));
    vispar.donotclipdomain =
      atoi (Tcl_GetVar (interp, "::viewoptions.clipping.notdomain", TCL_GLOBAL_ONLY));

    if ( ! (hclip == vispar.clipping) )
      {
	vispar.clipping = hclip;
	vispar.clipping.timestamp = NextTimeStamp();
      }

    vispar.whitebackground = atoi (Tcl_GetVar (interp, "::viewoptions.whitebackground", TCL_GLOBAL_ONLY));
    vispar.drawcoordinatecross = atoi (Tcl_GetVar (interp, "::viewoptions.drawcoordinatecross", TCL_GLOBAL_ONLY));
    vispar.drawcolorbar = atoi (Tcl_GetVar (interp, "::viewoptions.drawcolorbar", TCL_GLOBAL_ONLY));
    vispar.drawnetgenlogo = atoi (Tcl_GetVar (interp, "::viewoptions.drawnetgenlogo", TCL_GLOBAL_ONLY));
    vispar.stereo = atoi (Tcl_GetVar (interp, "::viewoptions.stereo", TCL_GLOBAL_ONLY));
    vispar.colormeshsize = atoi (Tcl_GetVar (interp, "::viewoptions.colormeshsize", TCL_GLOBAL_ONLY));
    VisualScene :: SetBackGroundColor (vispar.whitebackground ? 1 : 0);

    strcpy (vispar.selectvisual, Tcl_GetVar (interp, "::selectvisual", TCL_GLOBAL_ONLY));

    //  vispar.showstltrias = atoi (Tcl_GetVar (interp, "::viewoptions.stl.showtrias", TCL_GLOBAL_ONLY));

    vispar.stlshowtrias =
      atoi (Tcl_GetVar (interp, "::stloptions.showtrias", TCL_GLOBAL_ONLY));
    vispar.stlshowfilledtrias =
      atoi (Tcl_GetVar (interp, "::stloptions.showfilledtrias", TCL_GLOBAL_ONLY));
    vispar.stlshowedges =
      atoi (Tcl_GetVar (interp, "::stloptions.showedges", TCL_GLOBAL_ONLY));
    vispar.stlshowmarktrias =
      atoi (Tcl_GetVar (interp, "::stloptions.showmarktrias", TCL_GLOBAL_ONLY));
    vispar.stlshowactivechart =
      atoi (Tcl_GetVar (interp, "::stloptions.showactivechart", TCL_GLOBAL_ONLY));
    vispar.stlchartnumber =
      atoi (Tcl_GetVar (interp, "::stloptions.chartnumber", TCL_GLOBAL_ONLY));
    vispar.stlchartnumberoffset =
      atoi (Tcl_GetVar (interp, "::stloptions.chartnumberoffset", TCL_GLOBAL_ONLY));

    vispar.occshowsurfaces =
      atoi (Tcl_GetVar (interp, "::occoptions.showsurfaces", TCL_GLOBAL_ONLY));
    vispar.occshowedges =
      atoi (Tcl_GetVar (interp, "::occoptions.showedges", TCL_GLOBAL_ONLY));

    vispar.drawoutline =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawoutline", TCL_GLOBAL_ONLY));
    vispar.drawfilledtrigs =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawfilledtrigs", TCL_GLOBAL_ONLY));
    vispar.subdivisions =
      atoi (Tcl_GetVar (interp, "::visoptions.subdivisions", TCL_GLOBAL_ONLY));
    vispar.drawbadels =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawbadels", TCL_GLOBAL_ONLY));
    vispar.drawedges =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedges", TCL_GLOBAL_ONLY));

    vispar.drawtetsdomain =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawtetsdomain", TCL_GLOBAL_ONLY));
    vispar.drawtets =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawtets", TCL_GLOBAL_ONLY));
    vispar.drawprisms =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawprisms", TCL_GLOBAL_ONLY));
    vispar.drawpyramids =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawpyramids", TCL_GLOBAL_ONLY));
    vispar.drawhexes =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawhexes", TCL_GLOBAL_ONLY));
    /*
    vispar.shrink =
      atof (Tcl_GetVar (interp, "::viewoptions.shrink", TCL_GLOBAL_ONLY));
    */
    double hshrink = atof (Tcl_GetVar (interp, "::viewoptions.shrink", TCL_GLOBAL_ONLY));
    if (hshrink != vispar.shrink)
      { vispar.shrink = hshrink; vispar.clipping.timestamp = NextTimeStamp();}
    vispar.drawidentified =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawidentified", TCL_GLOBAL_ONLY));
    vispar.drawpointnumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawpointnumbers", TCL_GLOBAL_ONLY));
    vispar.drawedgenumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedgenumbers", TCL_GLOBAL_ONLY));
    vispar.drawfacenumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawfacenumbers", TCL_GLOBAL_ONLY));
    vispar.drawelementnumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawelementnumbers", TCL_GLOBAL_ONLY));
    vispar.drawdomainsurf =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawdomainsurf", TCL_GLOBAL_ONLY));

    vispar.drawededges =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawededges", TCL_GLOBAL_ONLY));
    vispar.drawedpoints =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedpoints", TCL_GLOBAL_ONLY));
    vispar.drawedpointnrs =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedpointnrs", TCL_GLOBAL_ONLY));
    vispar.drawedtangents =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedtangents", TCL_GLOBAL_ONLY));
    vispar.drawededgenrs =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawededgenrs", TCL_GLOBAL_ONLY));

    vispar.drawcurveproj =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawcurveproj", TCL_GLOBAL_ONLY));
    vispar.drawcurveprojedge =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawcurveprojedge", TCL_GLOBAL_ONLY));

    vispar.centerpoint =
      atoi (Tcl_GetVar (interp, "::viewoptions.centerpoint", TCL_GLOBAL_ONLY));
    vispar.use_center_coords =
      atoi (Tcl_GetVar (interp, "::viewoptions.usecentercoords", TCL_GLOBAL_ONLY)) > 0;
    vispar.centerx =
      atof (Tcl_GetVar (interp, "::viewoptions.centerx", TCL_GLOBAL_ONLY));
    vispar.centery =
      atof (Tcl_GetVar (interp, "::viewoptions.centery", TCL_GLOBAL_ONLY));
    vispar.centerz =
      atof (Tcl_GetVar (interp, "::viewoptions.centerz", TCL_GLOBAL_ONLY));
    vispar.drawelement =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawelement", TCL_GLOBAL_ONLY));
    vispar.drawmetispartition =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawmetispartition", TCL_GLOBAL_ONLY));

    vispar.drawspecpoint =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawspecpoint", TCL_GLOBAL_ONLY));
    vispar.specpointx =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointx", TCL_GLOBAL_ONLY));
    vispar.specpointy =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointy", TCL_GLOBAL_ONLY));
    vispar.specpointz =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointz", TCL_GLOBAL_ONLY));


    vsspecpoints.len =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointvlen", TCL_GLOBAL_ONLY));

    vispar.occdeflection = pow(10.0,-1-atof (Tcl_GetVar (interp, "::occoptions.deflection", TCL_GLOBAL_ONLY)));



#ifdef PARALLELGL
    vsmesh.Broadcast ();
#endif

    return TCL_OK;
  }



  int Ng_BuildFieldLines (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    netgen::GetVSSolution().BuildFieldLinesPlot();
    return TCL_OK;
  }



  int Ng_Exit (ClientData clientData,
	       Tcl_Interp * interp,
	       int argc, tcl_const char *argv[])
  {
    /*
#ifdef PARALLEL
    int id, rc, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if ( id != 0 )
      return TCL_OK;
#endif
    */

    /*
    if (ngsolve_handle)
      {
        void (*ngs_exit)();
        ngs_exit = ( void (*)() ) dlsym (ngsolve_handle, "NGSolve_Exit");
        if (ngs_exit) (*ngs_exit)();
      }
    */

#ifdef NGSOLVE
    NGSolve_Exit ();
#endif



#ifdef ACIS
    outcome res;
    res = api_terminate_faceter();
    if(!res.ok())
      cerr << "problem with terminating acis faceter" << endl;
    res = api_terminate_constructors();
    if(!res.ok())
      cerr << "problem with terminating acis constructors" << endl;
    res = api_terminate_kernel();
    if(!res.ok())
      cerr << "problem with terminating acis kernel" << endl;
    res = api_stop_modeller();
    if(!res.ok())
      cerr << "problem with terminating acis modeller" << endl;

    //cout << "stopped acis, outcome = " << res.ok() << endl;
#endif

#ifdef PARALLELGL
    if (id == 0) MyMPI_SendCmd ("end");
    MPI_Finalize();
#endif

    mesh.reset();
    ng_geometry.reset();
    
    if (testout != &cout)
      delete testout;

    return TCL_OK;
  }


#ifdef SOCKETS
  void * ServerSocketManagerRunDummy ( void * nix )
  {
    serversocketmanager.Run();
    return NULL;
  }

  extern "C" int Ng_ServerSocketManagerRun( void );
  int Ng_ServerSocketManagerRun( void )
  {
    if(mparam.parthread)
      RunParallel(ServerSocketManagerRunDummy,NULL);
    else
      serversocketmanager.Run();

    return TCL_OK;
  }

  extern "C" int Ng_ServerSocketManagerInit(int port);
  int Ng_ServerSocketManagerInit(int port)
  {
    serversocketmanager.Init(port);
    return TCL_OK;
  }
#endif //SOCKETS


  extern "C" int Ng_Init (Tcl_Interp * interp);
  extern "C" int Ng_CSG_Init (Tcl_Interp * interp);

  extern "C" int Ng_stl_Init (Tcl_Interp * interp);
  extern "C" int Ng_geom2d_Init (Tcl_Interp * interp);
  #ifdef OCCGEOMETRY
    extern "C" int Ng_occ_Init (Tcl_Interp * interp);
#endif


  // extern "C" int Ng_Geom2d_Init (Tcl_Interp * interp); 

  //   int main_Eero (ClientData clientData,
  // 	       Tcl_Interp * interp,
  // 		 int argc, tcl_const char *argv[]);


  int Ng_Init (Tcl_Interp * interp)
  {
#ifdef SOCKETS
    if(serversocketmanager.Good())
      serversocketusernetgen.Reset(new ServerSocketUserNetgen (serversocketmanager, mesh, geometry));
#endif

    Ng_CSG_Init(interp);
    Ng_stl_Init(interp);
    Ng_geom2d_Init (interp);
#ifdef OCCGEOMETRY
    Ng_occ_Init (interp);
#endif


    // Ng_Geom2d_Init(interp);

    tcl_interp = interp;

    //     Tcl_CreateCommand (interp, "Ng_Eero", main_Eero,
    // 		       (ClientData)NULL,
    // 		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_New", Ng_New,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    //     Tcl_CreateCommand (interp, "Ng_Lock", Ng_Lock,
    // 		       (ClientData)NULL,
    // 		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_LoadGeometry", Ng_LoadGeometry,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SaveGeometry", Ng_SaveGeometry,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_LoadMesh", Ng_LoadMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SaveMesh", Ng_SaveMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MergeMesh", Ng_MergeMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetExportFormats", Ng_GetExportFormats,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ExportMesh", Ng_ExportMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ImportMesh", Ng_ImportMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ImportSolution", Ng_ImportSolution,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ShowDemo", Ng_ShowDemo,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_DemoSetTime", Ng_DemoSetTime,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SaveSolution", Ng_SaveSolution,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);




    

    // meshing
    Tcl_CreateCommand (interp, "Ng_GenerateMesh", Ng_GenerateMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_StopMeshing", Ng_StopMeshing,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MeshInfo", Ng_MeshInfo,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MeshQuality", Ng_MeshQuality,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_CheckSurfaceMesh", Ng_CheckSurfaceMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_CheckVolumeMesh", Ng_CheckVolumeMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_DeleteVolMesh", Ng_DeleteVolMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SplitSeparatedFaces", Ng_SplitSeparatedFaces,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetNextTimeStamp", Ng_SetNextTimeStamp,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Refine", Ng_Refine,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SecondOrder", Ng_SecondOrder,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_HighOrder", Ng_HighOrder,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ValidateSecondOrder", Ng_ValidateSecondOrder,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_RestrictH", Ng_RestrictH,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Anisotropy", Ng_Anisotropy,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Bisect", Ng_Bisect,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    //     Tcl_CreateCommand (interp, "Ng_BisectCopyMesh", Ng_BisectCopyMesh,
    // 		       (ClientData)NULL,
    // 		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Split2Tets", Ng_Split2Tets,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ZRefinement", Ng_ZRefinement,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_HPRefinement", Ng_HPRefinement,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_LoadMeshSize", Ng_LoadMeshSize,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MeshSizeFromSurfaceMesh", Ng_MeshSizeFromSurfaceMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_GenerateBoundaryLayer", Ng_GenerateBoundaryLayer,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_InsertVirtualBL", Ng_InsertVirtualBL,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_CutOffAndCombine", Ng_CutOffAndCombine,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_HelmholtzMesh", Ng_HelmholtzMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ReadStatus", Ng_ReadStatus,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MemInfo", Ng_MemInfo,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_MeshDoctor", Ng_MeshDoctor,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_BCProp", Ng_BCProp,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);




    Tcl_CreateCommand (interp, "Ng_ACISCommand",
		       Ng_ACISCommand,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_MouseMove", Ng_MouseMove,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MouseDblClick", Ng_MouseDblClick,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_ZoomAll", Ng_ZoomAll,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Center", Ng_Center,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_StandardRotation", Ng_StandardRotation,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ArbitraryRotation", Ng_ArbitraryRotation,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_SetVisParameters", Ng_SetVisParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetMeshingParameters", Ng_SetMeshingParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetDebugParameters", Ng_SetDebugParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetCommandLineParameter",
		       Ng_SetCommandLineParameter,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetCommandLineParameter",
		       Ng_GetCommandLineParameter,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Exit",
		       Ng_Exit,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Metis",
		       Ng_Metis,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_BuildFieldLines",
		       Ng_BuildFieldLines,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    /*
     * Specify the C callback functions for widget creation, display,
     * and reshape.
     */
    Tcl_CreateObjCommand(interp, "Ng_GetToglVersion", Ng_ToglVersion, NULL, NULL);

    if (!nodisplay)
      {
	if (Togl_Init(interp) == TCL_ERROR) 
	  return TCL_ERROR;
	
	
	Tcl_CreateObjCommand(interp, "init", init, NULL, NULL);
	Tcl_CreateObjCommand(interp, "zap", zap, NULL, NULL);
	Tcl_CreateObjCommand(interp, "draw", draw, NULL, NULL);
	Tcl_CreateObjCommand(interp, "reshape", reshape, NULL, NULL);
	
	//   Togl_TimerFunc(  idle );
	Tcl_CreateObjCommand(interp, "Ng_SnapShot", Ng_SnapShot, NULL, NULL);
        Tcl_CreateObjCommand(interp, "Ng_VideoClip", Ng_VideoClip, NULL, NULL);
      }


    multithread.pause = 0;
    multithread.testmode = 0;
    multithread.redraw = 0;
    multithread.drawing = 1;
    multithread.terminate = 0;
    multithread.running = 0;

    multithread.task = "";
    multithread.percent = 20;


    Tcl_LinkVar (interp, "multithread_pause",
		 (char*)&multithread.pause, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_testmode",
		 (char*)&multithread.testmode, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_redraw",
		 (char*)&multithread.redraw, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_drawing",
		 (char*)&multithread.drawing, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_terminate",
		 (char*)&multithread.terminate, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_running",
		 (char*)&multithread.running, TCL_LINK_INT);


    //testout->setstate(ios_base::badbit);

    myerr = &cerr;
    extern ostream * mycout;
    mycout = &cout;

    testmode = 0;



#ifdef ACIS
    outcome res;
    res = api_start_modeller (0);
    if(!res.ok())
      cerr << "problem with starting acis modeller" << endl;

#ifdef ACIS_R17
    unlock_spatial_products_661();
#endif
    res = api_initialize_kernel();
    if(!res.ok())
      cerr << "problem with starting acis kernel" << endl;
    res = api_initialize_constructors();
    if(!res.ok())
      cerr << "problem with starting acis constructors" << endl;
    res = api_initialize_faceter();
    if(!res.ok())
      cerr << "problem with starting acis faceter" << endl;
#endif


    return TCL_OK;
  }


}



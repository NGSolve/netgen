#include <mystdlib.h> 
#include <inctcl.hpp>
#include <meshing.hpp>

#include "../libsrc/interface/writeuser.hpp"
#include <pybind11/pybind11.h>

namespace netgen
{
  DLL_HEADER extern Flags parameters;
  DLL_HEADER extern bool netgen_executable_started;
  DLL_HEADER extern int printmessage_importance;
}

using netgen::parameters;
using netgen::ngdir;
using netgen::verbose;
using netgen::Array;
using netgen::RegisterUserFormats;

extern "C" int Ng_ServerSocketManagerInit (int port);
extern "C" int Ng_ServerSocketManagerRun (void);
 
int Tcl_AppInit(Tcl_Interp * interp);

extern bool nodisplay;
extern bool shellmode;

void main_gui()
{
  if(netgen::netgen_executable_started)
      return;

  if ( netgen::id == 0 )
    {
      cout << "NETGEN-" << PACKAGE_VERSION << endl;
      
      cout << "Developed by Joachim Schoeberl at" << endl
	   << "2010-xxxx Vienna University of Technology" << endl
	   << "2006-2010 RWTH Aachen University" << endl
           << "1996-2006 Johannes Kepler University Linz" << endl;
      
#ifdef OCCGEOMETRY
      cout << "Including OpenCascade geometry kernel" << endl;
#endif
      
#ifdef ACIS
      cout << "Including ACIS geometry kernel" << endl;
#endif

#ifdef DEBUG
      cout << "You are running the debug version !" << endl;
#endif


    }

//   netgen::h_argc = argc;
//   netgen::h_argv = argv;

  if (getenv ("NETGENDIR") && strlen (getenv ("NETGENDIR")))
    ngdir = getenv ("NETGENDIR");
  else
    ngdir = ".";
  
  verbose = parameters.GetDefineFlag ("V");

  if (verbose)
    cout << "NETGENDIR = " << ngdir << endl;
  

  if ( netgen::id == 0 )
    {
      if (parameters.StringFlagDefined ("testout"))      
        netgen::testout = new ofstream (parameters.GetStringFlag ("testout", "test.out"));


      if(parameters.GetDefineFlag("batchmode"))
        nodisplay = true;

    
      if(parameters.GetDefineFlag("shellmode"))
        {
          nodisplay = true;
          shellmode = true;
        }

      Tcl_FindExecutable(NULL);

      // initialize application
      Tcl_Interp * myinterp = Tcl_CreateInterp ();
      if (Tcl_AppInit (myinterp) == TCL_ERROR)
        {
          cerr << "Exit Netgen due to initialization problem" << endl;
          exit (1);
        }



      // parse tcl-script
      int errcode;

      bool internaltcl = false;
      if (shellmode)
        internaltcl = false;
  
      if (verbose)
        {
          cout << "Tcl header version = " << TCL_PATCH_LEVEL << endl;
          Tcl_Eval (myinterp, "puts \"Tcl runtime version = [info patchlevel] \";");
        }

      if (parameters.GetDefineFlag ("internaltcl"))
        internaltcl=true;
      if (parameters.GetDefineFlag ("externaltcl"))
        internaltcl=false;

 

      if (internaltcl)
        {
          if (verbose)
            cout << "using internal Tcl-script" << endl;
      
          // connect to one string 
          extern const char * ngscript[];
          const char ** hcp = ngscript;
          int len = 0;
          while (*hcp)
            len += strlen (*hcp++); 

          char * tr1 = new char[len+1];
          *tr1 = 0;
          hcp = ngscript;
      
          char * tt1 = tr1;
          while (*hcp)
            {
              strcat (tt1, *hcp); 
              tt1 += strlen (*hcp++);
            }
      
          errcode = Tcl_Eval (myinterp, tr1);
          delete [] tr1;
        }

      else

        {
          string startfile = ngdir + "/ng.tcl";
      
          if (verbose)
            cout << "Load Tcl-script from " << startfile << endl;
      
          errcode = Tcl_EvalFile (myinterp, (char*)startfile.c_str());
        }

      if (errcode)
        {
          cout << "Error in Tcl-Script:" << endl;
          // cout << "result = " << myinterp->result << endl;
          cout << "result = " << Tcl_GetStringResult (myinterp) << endl;
          // cout << "in line " << myinterp->errorLine << endl;

          cout << "\nMake sure to set environment variable NETGENDIR to directory containing ng.tcl" << endl;
          exit (1);
        }


      // lookup user file formats and insert into format list:
      Array<const char*> userformats;
      Array<const char*> extensions;
      RegisterUserFormats (userformats, extensions);

      ostringstream fstr;
      
      tcl_const char * exportft = Tcl_GetVar (myinterp, "exportfiletype", 0);
      for (int i = 1; i <= userformats.Size(); i++)
	{
	  fstr << ".ngmenu.file.filetype add radio -label \"" 
	       << userformats.Get(i) << "\" -variable exportfiletype\n";
	  fstr << "lappend meshexportformats { {" << userformats.Get(i) << "} {" << extensions.Get(i) << "} }\n";
	}

      Tcl_Eval (myinterp, (char*)fstr.str().c_str());
      Tcl_SetVar (myinterp, "exportfiletype", exportft, 0);


      // start event-loop
      Tk_MainLoop();
      Tcl_DeleteInterp (myinterp); 
//       Tcl_Exit(0);
    }
}



extern "C" int Ng_Init (Tcl_Interp * interp);
extern "C" int Ng_Vis_Init (Tcl_Interp * interp);

extern "C" void Ng_TclCmd(string);

struct GuiThread {
    std::thread t;
    ~GuiThread() {
        Ng_TclCmd(".ngmenu.file invoke \"Quit\";");
        if(netgen::printmessage_importance>2)
          cout << "waiting for GUI to finish..." << endl;
        t.join();
    }
};

static GuiThread gui_thread;

PYBIND11_PLUGIN(gui) {
  pybind11::module m("gui", "pybind gui");
    gui_thread.t = std::thread([]() 
                {
                  main_gui();
                });
  return m.ptr();
}

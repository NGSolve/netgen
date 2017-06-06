#include <mystdlib.h> 
#include <inctcl.hpp>
#include <meshing.hpp>

#ifdef WIN32
  #define DLL_HEADER_IMPORT __declspec(dllimport)
  #define DLL_HEADER_EXPORT __declspec(dllexport)
#else
  #define DLL_HEADER_IMPORT
  #define DLL_HEADER_EXPORT
#endif


namespace netgen
{
  DLL_HEADER_EXPORT Flags parameters;
}

DLL_HEADER_EXPORT bool nodisplay = false;

extern "C" int Ng_Init (Tcl_Interp * interp);
extern "C" int Ng_Vis_Init (Tcl_Interp * interp);
extern "C" void Ng_TclCmd(string);

// tcl package dynamic load
extern "C" int DLL_HEADER_EXPORT Gui_Init (Tcl_Interp * interp)
{
  if (Ng_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Ng_Init: " << endl;
    cout << "result = " << Tcl_GetStringResult (interp) << endl;
    return TCL_ERROR;
  }
 
  if (!nodisplay && Ng_Vis_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Ng_Vis_Init: " << endl;
    cout << "result = " << Tcl_GetStringResult (interp) << endl;
    return TCL_ERROR;
  }

  return TCL_OK;
}

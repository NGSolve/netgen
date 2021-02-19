#include <mystdlib.h> 
#include <inctcl.hpp>
#include <meshing.hpp>
#include <core/ngcore_api.hpp>

namespace netgen
{
  NGCORE_API_EXPORT Flags parameters;
}

NGCORE_API_EXPORT bool nodisplay = false;

extern "C" int Ng_Init (Tcl_Interp * interp);
extern "C" int Ng_Vis_Init (Tcl_Interp * interp);
extern "C" void Ng_TclCmd(string);

// tcl package dynamic load
extern "C" int NGCORE_API_EXPORT Gui_Init (Tcl_Interp * interp)
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

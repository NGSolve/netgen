#include <meshing.hpp>
#include <inctcl.hpp>
#include "../libsrc/meshing/visual_interface.hpp"


static void Impl_Ng_Tcl_SetResult(Tcl_Interp *interp, char *result, Tcl_FreeProc *freeProc)
{
    Tcl_SetResult(interp, result, freeProc);
}

static void Impl_Ng_Tcl_CreateCommand(Tcl_Interp *interp, const char *cmdName, Tcl_CmdProc *proc)
{
    Tcl_CreateCommand(interp, cmdName, proc, nullptr, nullptr);
}

static bool dummy_init_pointers = [](){
    netgen::Ptr_Ng_Tcl_SetResult = Impl_Ng_Tcl_SetResult;
    netgen::Ptr_Ng_Tcl_CreateCommand = Impl_Ng_Tcl_CreateCommand;
    return true;
}();

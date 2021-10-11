#include "ngtcl.hpp"

#include <inctcl.hpp>

namespace netgen
{
    void Ng_Tcl_SetResult(Tcl_Interp *interp, char *result, const int freeProc)
    {
      Tcl_SetResult(interp, result, (Tcl_FreeProc*)freeProc);
    }

    void Ng_Tcl_CreateCommand(Tcl_Interp *interp, const char *cmdName, Tcl_CmdProc *proc)
    {
      Tcl_CreateCommand(interp, cmdName, proc, nullptr, nullptr);
    }
}

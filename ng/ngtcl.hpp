#ifndef FILE_NG_TCL_HPP
#define FILE_NG_TCL_HPP

#include <myadt.hpp>

class Tcl_Interp;
class Tcl_cmdProc;

namespace netgen
{
    typedef int (Tcl_CmdProc) (void * clientData, Tcl_Interp *interp,
            int argc, const char *argv[]);

    inline constexpr int NG_TCL_VOLATILE = 1;
    inline constexpr int NG_TCL_STATIC   = 0;
    inline constexpr int NG_TCL_DYNAMIC  = 3;

    inline constexpr int NG_TCL_OK       = 0;
    inline constexpr int NG_TCL_ERROR    = 1;
    inline constexpr int NG_TCL_RETURN   = 2;
    inline constexpr int NG_TCL_BREAK    = 3;
    inline constexpr int NG_TCL_CONTINUE = 4;

    DLL_HEADER void Ng_Tcl_SetResult(Tcl_Interp *interp, char *result, const int freeProc);
    DLL_HEADER void Ng_Tcl_CreateCommand(Tcl_Interp *interp,
				const char *cmdName, Tcl_CmdProc *proc);

}

#endif // FILE_NG_TCL_HPP

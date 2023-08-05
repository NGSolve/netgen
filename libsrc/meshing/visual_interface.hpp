#ifndef VISUAL_INTERFACE_HPP_INCLUDED
#define VISUAL_INTERFACE_HPP_INCLUDED

#include <mystdlib.h>
#include <meshing.hpp>
#include <myadt.hpp>

struct Ng_SolutionData;

// Function pointers for visualization purposed, all set to nullptr by default and initialized correctly when the GUI library is loaded

DLL_HEADER extern void (*Ptr_Ng_ClearSolutionData) ();
DLL_HEADER extern void (*Ptr_Ng_InitSolutionData) (Ng_SolutionData * soldata);
DLL_HEADER extern void (*Ptr_Ng_SetSolutionData) (Ng_SolutionData * soldata);
DLL_HEADER extern void (*Ptr_Ng_Redraw) (bool blocking);

// Tcl wrapper functions
struct Tcl_Interp;
typedef int (Tcl_CmdProc) (void * clientData, Tcl_Interp *interp,
        int argc, const char *argv[]);
typedef void (Tcl_FreeProc) (char *blockPtr);

namespace netgen {
  /*
  inline constexpr int NG_TCL_VOLATILE = 1;
  inline constexpr int NG_TCL_STATIC   = 0;
  inline constexpr int NG_TCL_DYNAMIC  = 3;
  */

#define NG_TCL_VOLATILE		((Tcl_FreeProc *) 1)
#define NG_TCL_STATIC		((Tcl_FreeProc *) 0)
#define NG_TCL_DYNAMIC		((Tcl_FreeProc *) 3)

    inline constexpr int NG_TCL_OK       = 0;
    inline constexpr int NG_TCL_ERROR    = 1;
    inline constexpr int NG_TCL_RETURN   = 2;
    inline constexpr int NG_TCL_BREAK    = 3;
    inline constexpr int NG_TCL_CONTINUE = 4;
    DLL_HEADER extern void (*Ptr_Ng_Tcl_SetResult)(Tcl_Interp *interp, char *result, Tcl_FreeProc *freeProc);
    DLL_HEADER extern void (*Ptr_Ng_Tcl_CreateCommand)(Tcl_Interp *interp,
                                    const char *cmdName, Tcl_CmdProc *proc);

    DLL_HEADER extern void (*Ptr_Render)(bool);
    DLL_HEADER extern void (*Ptr_UpdateVisSurfaceMeshData)(int,
            shared_ptr<NgArray<Point<3>>>,
            shared_ptr<NgArray<INDEX_2>>,
            shared_ptr<NgArray<Point<2>>>
            );

    inline void Render(bool blocking = false) { if(Ptr_Render) Ptr_Render(blocking); }
    inline void UpdateVisSurfaceMeshData(int oldnl,
            shared_ptr<NgArray<Point<3>>> locpointsptr = nullptr,
            shared_ptr<NgArray<INDEX_2>> loclinesptr = nullptr,
            shared_ptr<NgArray<Point<2>>> plainpointsptr = nullptr
            ) {
        if(Ptr_UpdateVisSurfaceMeshData) Ptr_UpdateVisSurfaceMeshData(oldnl, locpointsptr, loclinesptr, plainpointsptr);
    }

    inline void Ng_Tcl_SetResult(Tcl_Interp *interp, char *result, Tcl_FreeProc *freeProc)
    {
        if(Ptr_Ng_Tcl_SetResult)
            Ptr_Ng_Tcl_SetResult(interp, result, freeProc);
    }
    inline void Ng_Tcl_CreateCommand(Tcl_Interp *interp, const char *cmdName, Tcl_CmdProc *proc)
    {
        if(Ptr_Ng_Tcl_CreateCommand)
            Ptr_Ng_Tcl_CreateCommand(interp, cmdName, proc);
    }
}

#endif // VISUAL_INTERFACE_HPP_INCLUDED

#include "visual_interface.hpp"
#include "../include/nginterface.h"

void (*Ptr_Ng_ClearSolutionData) () = nullptr;
void (*Ptr_Ng_InitSolutionData) (Ng_SolutionData*)  = nullptr;
void (*Ptr_Ng_SetSolutionData) (Ng_SolutionData*) = nullptr;
void (*Ptr_Ng_Redraw) (bool blocking) = nullptr;

void Ng_ClearSolutionData () { if(Ptr_Ng_ClearSolutionData) Ptr_Ng_ClearSolutionData(); }
void Ng_InitSolutionData (Ng_SolutionData * soldata) { if(Ptr_Ng_InitSolutionData) Ptr_Ng_InitSolutionData(soldata); }
void Ng_SetSolutionData (Ng_SolutionData * soldata) { if(Ptr_Ng_SetSolutionData) Ptr_Ng_SetSolutionData(soldata); }
void Ng_Redraw (bool blocking) { if(Ptr_Ng_Redraw) Ptr_Ng_Redraw(blocking); }

namespace netgen
{
    void (*Ptr_Ng_Tcl_SetResult)(Tcl_Interp *interp, char *result, Tcl_FreeProc *freeProc) = nullptr;
    void (*Ptr_Ng_Tcl_CreateCommand)(Tcl_Interp *interp,
                                    const char *cmdName, Tcl_CmdProc *proc) = nullptr;
    void (*Ptr_Render)(bool) = nullptr;
    void (*Ptr_UpdateVisSurfaceMeshData)(int,
            shared_ptr<NgArray<Point<3>>>,
            shared_ptr<NgArray<INDEX_2>>,
            shared_ptr<NgArray<Point<2>>>
            ) = nullptr;
} // namespace netgen


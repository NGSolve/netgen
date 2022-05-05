#ifndef VISUAL_INTERFACE_HPP_INCLUDED
#define VISUAL_INTERFACE_HPP_INCLUDED

#include <mystdlib.h>
#include <meshing.hpp>

class Ng_SolutionData;

// Function pointers for visualization purposed, all set to nullptr by default and initialized correctly when the GUI library is loaded

DLL_HEADER extern void (*Ptr_Ng_ClearSolutionData) ();
DLL_HEADER extern void (*Ptr_Ng_InitSolutionData) (Ng_SolutionData * soldata);
DLL_HEADER extern void (*Ptr_Ng_SetSolutionData) (Ng_SolutionData * soldata);
DLL_HEADER extern void (*Ptr_Ng_Redraw) (bool blocking);

namespace netgen {
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
}

#endif // VISUAL_INTERFACE_HPP_INCLUDED

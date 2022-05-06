#include "../libsrc/meshing/visual_interface.hpp"

#include <inctcl.hpp>

static bool dummy_init_pointers = [](){
    Ptr_Ng_Tcl_SetResult = Tcl_SetResult;
    Ptr_Ng_Tcl_CreateCommand = Tcl_CreateCommand;
    return true;
}();

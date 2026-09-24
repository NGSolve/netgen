#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <tcl.h>
#include <tk.h>

#if TCL_MAJOR_VERSION >= 9
// Tcl 9 implements Tcl_SetResult as a macro, which lets temporaries like
// str.str().c_str() die before the string is copied
#undef Tcl_SetResult
inline void Tcl_SetResult(Tcl_Interp *interp, const char *result, Tcl_FreeProc *freeProc)
{
  Tcl_SetObjResult(interp, Tcl_NewStringObj(result, -1));
  if (result != NULL && freeProc != NULL && freeProc != TCL_VOLATILE)
    {
      if (freeProc == TCL_DYNAMIC)
        Tcl_Free((void *)result);
      else
        (*freeProc)((void *)result);
    }
}
#endif

#if TK_MAJOR_VERSION>8 || (TK_MAJOR_VERSION==8 && TK_MINOR_VERSION>=4)
#define tcl_const const
#else
#define tcl_const
#endif

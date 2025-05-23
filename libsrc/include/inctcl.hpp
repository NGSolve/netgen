#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <tcl.h>
#include <tk.h>

#if TK_MAJOR_VERSION==8 && TK_MINOR_VERSION>=4
#define tcl_const const
#else
#define tcl_const
#endif

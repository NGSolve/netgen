#ifndef FILE_VISUAL
#define FILE_VISUAL

/* *************************************************************************/
/* File:   visual.hpp                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   02. Dec. 01                                                     */
/* *************************************************************************/

/* 

Visualization

*/

// #ifdef PARALLEL
// #define PARALLELGL
// #endif

/*** Windows headers ***/
#ifdef _MSC_VER
# define WIN32_LEAN_AND_MEAN
# ifndef NO_PARALLEL_THREADS
#  ifdef MSVC_EXPRESS
#  else
#   include <afxwin.h>
#   include <afxmt.h>
#  endif // MSVC_EXPRESS
# endif
# include <windows.h>
# undef WIN32_LEAN_AND_MEAN
# include <winnt.h>

#else // Not using MC VC++


#endif



#include "visual_api.hpp"
#include "../include/incopengl.hpp"

#include "../meshing/visual_interface.hpp"
#include "../meshing/soldata.hpp"
#include "vispar.hpp"
#include "mvdraw.hpp"

#include <complex>
#include "vssolution.hpp"
#include "meshdoc.hpp"

#endif

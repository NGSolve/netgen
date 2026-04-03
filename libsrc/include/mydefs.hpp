#ifndef FILE_MYDEFS
#define FILE_MYDEFS

/**************************************************************************/
/* File:   mydefs.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Mar. 98                                                    */
/**************************************************************************/

/*
  defines for graphics, testmodes, ...
*/

#include <core/ngcore.hpp>
#define PACKAGE_VERSION "6.2-dev"

// #define DEBUG

#if defined(nglib_EXPORTS)
   #define DLL_HEADER   NGCORE_API_EXPORT
#else
   #define DLL_HEADER   NGCORE_API_IMPORT
#endif




#ifndef __assume
#ifdef __GNUC__
#define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
#else
#define __assume(cond)
#endif
#endif


#ifndef NG_INLINE
#ifdef __INTEL_COMPILER
#ifdef WIN32
#define NG_INLINE __forceinline inline
#else
#define NG_INLINE __forceinline inline
#endif
#else
#ifdef __GNUC__
#define NG_INLINE __attribute__ ((__always_inline__)) inline
#define VLA
#else
#define NG_INLINE inline
#endif
#endif
#endif


// #define BASE0
// #define DEBUG


#define noDEMOVERSION
#define noDEVELOP
#define noSTEP
#define noSOLIDGEOM

#define noDEMOAPP
#define noMODELLER

#define noSTAT_STREAM
#define noLOG_STREAM

#endif

#ifndef NETGEN_CORE_NGCORE_API_HPP
#define NETGEN_CORE_NGCORE_API_HPP

#ifdef WIN32
        #define NGCORE_API_EXPORT __declspec(dllexport)
        #define NGCORE_API_IMPORT __declspec(dllimport)
#else
        #define NGCORE_API_EXPORT __attribute__((visibility("default")))
        #define NGCORE_API_IMPORT __attribute__((visibility("default")))
#endif

#ifdef NGCORE_EXPORTS
        #define NGCORE_API NGCORE_API_EXPORT
#else
        #define NGCORE_API NGCORE_API_IMPORT
#endif

#ifdef __INTEL_COMPILER
  #ifdef WIN32
    #define NETGEN_INLINE __forceinline inline
    #define NETGEN_LAMBDA_INLINE
  #else
    #define NETGEN_INLINE __forceinline inline
    #define NETGEN_LAMBDA_INLINE __attribute__ ((__always_inline__))
  #endif
#else
  #ifdef __GNUC__
    #define NETGEN_INLINE __attribute__ ((__always_inline__)) inline
    #define NETGEN_LAMBDA_INLINE __attribute__ ((__always_inline__))
    #define NETGEN_VLA
  #else
    #define NETGEN_INLINE inline
    #define NETGEN_LAMBDA_INLINE
  #endif
#endif

#endif // NETGEN_CORE_NGCORE_API_HPP

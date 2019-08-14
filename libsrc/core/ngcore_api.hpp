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

#ifdef __MAC_OS_X_VERSION_MIN_REQUIRED
#if __MAC_OS_X_VERSION_MIN_REQUIRED < 101400
// The c++ standard library on MacOS 10.13 and earlier has no aligned new operator,
// thus implement it here globally
#include <mm_malloc.h>
#ifdef __clang__
#pragma clang diagnostic ignored "-Winline-new-delete"
#endif
inline void * operator new (size_t s, std::align_val_t al)
{
  if (int(al) > __STDCPP_DEFAULT_NEW_ALIGNMENT__)
    return _mm_malloc(s, int(al));
  else
    return new char[s];
}

inline void * operator new[] (size_t s, std::align_val_t al)
{
  if (int(al) > __STDCPP_DEFAULT_NEW_ALIGNMENT__)
    return _mm_malloc(s, int(al));
  else
    return new char[s];
}

inline void operator delete  ( void* ptr, std::align_val_t al ) noexcept
{
  if (int(al) > __STDCPP_DEFAULT_NEW_ALIGNMENT__)
     _mm_free(ptr);
  else
    delete (char*)ptr;
}

inline void operator delete[]( void* ptr, std::align_val_t al ) noexcept
{
  if (int(al) > __STDCPP_DEFAULT_NEW_ALIGNMENT__)
     _mm_free(ptr);
  else
    delete[] (char*)ptr;
}
#endif // __MAC_OS_X_VERSION_MIN_REQUIRED
#endif // __MAC_OS_X_VERSION_MIN_REQUIRED < 101300

#endif // NETGEN_CORE_NGCORE_API_HPP

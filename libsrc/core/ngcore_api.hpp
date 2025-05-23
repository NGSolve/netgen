#ifndef NETGEN_CORE_NGCORE_API_HPP
#define NETGEN_CORE_NGCORE_API_HPP

#include "netgen_config.hpp"

#ifdef WIN32

// This function or variable may be unsafe. Consider using _ftime64_s instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. See online help for details.
#pragma warning(disable:4244)
#pragma warning(disable:4996)

// multiple inheritance via dominance
#pragma warning(disable:4250)

// needs to have dll-interface to be used by clients of class
#pragma warning(disable:4251)

// size_t to int conversion:
#pragma warning(disable:4267)

// non dll-interface class 'std::exception' used as base for dll-interface class
#pragma warning(disable:4275)

// C++ exception specification ignored except to indicate a function is not __declspec(nothrow)
#pragma warning(disable:4290)

// no suitable definition provided for explicit template instantiation request
#pragma warning(disable:4661)

// bool-int conversion
#pragma warning(disable:4800)

// '__declspec(dllexport)' and 'extern' are incompatible on an explicit instantiation
#pragma warning(disable:4910)

#endif // WIN32


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

// Set __host__ __device__ for all inline functions
#ifdef __CUDACC__
  #define NETGEN_HD __host__ __device__
#else // __CUDACC__
  #define NETGEN_HD
#endif // __CUDACC__

#ifdef __INTEL_COMPILER
  #define NETGEN_ALWAYS_INLINE __forceinline
  #define NETGEN_INLINE __forceinline inline
  #ifdef WIN32
    #define NETGEN_LAMBDA_INLINE
  #else
    #define NETGEN_LAMBDA_INLINE __attribute__ ((__always_inline__))
  #endif
#else
  #ifdef __GNUC__
    #define NETGEN_ALWAYS_INLINE __attribute__ ((__always_inline__))
    #define NETGEN_INLINE __attribute__ ((__always_inline__)) inline NETGEN_HD
    #define NETGEN_LAMBDA_INLINE __attribute__ ((__always_inline__)) NETGEN_HD
    #define NETGEN_VLA
  #else
    #define NETGEN_ALWAYS_INLINE
    #define NETGEN_INLINE inline
    #define NETGEN_LAMBDA_INLINE
  #endif
#endif

#if defined(__amd64__) || defined(_M_AMD64)
#define NETGEN_ARCH_AMD64
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
#define NETGEN_ARCH_ARM64
#endif

#if defined(__arm__) || defined(_M_ARM)
#define NETGEN_ARCH_ARM
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

inline void operator delete  ( void* ptr, std::size_t sz, std::align_val_t al ) noexcept
{
  if (int(al) > __STDCPP_DEFAULT_NEW_ALIGNMENT__)
     _mm_free(ptr);
  else
    delete (char*)ptr;
}

inline void operator delete[]( void* ptr, std::size_t sz, std::align_val_t al ) noexcept
{
  if (int(al) > __STDCPP_DEFAULT_NEW_ALIGNMENT__)
     _mm_free(ptr);
  else
    delete[] (char*)ptr;
}

#endif // __MAC_OS_X_VERSION_MIN_REQUIRED
#endif // __MAC_OS_X_VERSION_MIN_REQUIRED < 101300

#endif // NETGEN_CORE_NGCORE_API_HPP

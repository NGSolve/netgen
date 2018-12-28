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

namespace ngcore
{
#if defined(__GNUC__)
  inline bool likely (bool x) { return bool(__builtin_expect(long(x), 1L)); }
  inline bool unlikely (bool x) { return bool(__builtin_expect(long(x), 0L)); }
#else
  inline bool likely (bool x) { return x; }
  inline bool unlikely (bool x) { return x; }
#endif
} // namespace ngcore

#endif // NETGEN_CORE_NGCORE_API_HPP

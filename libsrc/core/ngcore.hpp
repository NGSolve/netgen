#ifndef NG_CORE_HPP
#define NG_CORE_HPP

// std includes
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <type_traits>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstring>
#include <complex>

#ifdef WIN32
        #define NGCORE_API_EXPORT __declspec(dllexport)
        #define NGCORE_API_IMPORT __declspec(dllimport)
#else
        #define NGCORE_API_EXPORT
        #define NGCORE_API_IMPORT
#endif

#ifdef NGCORE_EXPORTS
        #define NGCORE_API NGCORE_API_EXPORT
#else
        #define NGCORE_API NGCORE_API_IMPORT
#endif

namespace ngcore
{
#if defined(__GNUC__)
  inline bool likely (bool x) { return __builtin_expect((x), true); }
  inline bool unlikely (bool x) { return __builtin_expect((x), false); }
#else
  inline bool likely (bool x) { return x; }
  inline bool unlikely (bool x) { return x; }
#endif
}

// own includes
#include "type_traits.hpp"
#include "basearchive.hpp"
#include "version.hpp"
#include "archive.hpp"

#endif // NG_CORE_HPP

#ifndef NETGEN_CORE_UTILS_HPP
#define NETGEN_CORE_UTILS_HPP

#include <string>
#include <sstream>

#include "ngcore_api.hpp"       // for NGCORE_API

namespace ngcore
{
  NGCORE_API std::string Demangle(const char* typeinfo);

#if defined(__GNUC__)
  inline bool likely (bool x) { return bool(__builtin_expect(long(x), 1L)); }
  inline bool unlikely (bool x) { return bool(__builtin_expect(long(x), 0L)); }
#else
  inline bool likely (bool x) { return x; }
  inline bool unlikely (bool x) { return x; }
#endif

} // namespace ngcore

#endif // NETGEN_CORE_UTILS_HPP

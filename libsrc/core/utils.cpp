#include "utils.hpp"

#ifndef WIN32
#include <cxxabi.h>
#endif

namespace ngcore
{
#ifdef WIN32
  // windows does demangling in typeid(T).name()
  NGCORE_API std::string Demangle(const char* typeinfo) { return typeinfo; }
#else
  NGCORE_API std::string Demangle(const char* typeinfo) { int status; return abi::__cxa_demangle(typeinfo,
                                                                                      nullptr,
                                                                                      nullptr,
                                                                                      &status); }
} // namespace ngcore

#endif


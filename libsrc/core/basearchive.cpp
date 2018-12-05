
#include "ngcore.hpp"

#ifndef WIN
#include <cxxabi.h>
#endif

namespace ngcore
{
  std::map<std::string, VersionInfo>& GetLibraryVersions()
  {
    static std::map<std::string, VersionInfo> library_versions;
    return library_versions;
  }
#ifdef WIN
  // windows does demangling in typeid(T).name()
  std::string demangle(const char* typeinfo) { return typeinfo; }
#else
  std::string demangle(const char* typeinfo) { int status; return abi::__cxa_demangle(typeinfo, 0, 0, &status); }
#endif

  std::map<std::string, ClassArchiveInfo>& GetArchiveRegister()
  {
    static std::map<std::string, ClassArchiveInfo> type_register;
    return type_register;
  }
}

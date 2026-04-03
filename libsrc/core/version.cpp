#include <map>

#include <netgen_version.hpp>
#include "exception.hpp"
#include "version.hpp"

namespace ngcore
{
  // clang-tidy should ignore this static object
  static std::map<std::string, VersionInfo> library_versions;  // NOLINT

  const VersionInfo& GetLibraryVersion(const std::string& library)
  { return library_versions[library]; }

  const std::map<std::string, VersionInfo>& GetLibraryVersions()
  { return library_versions; }

  void SetLibraryVersion(const std::string& library, const VersionInfo& version)
  {
    if(library_versions.count(library) && (library_versions[library] != version))
      throw Exception("Failed to set library version for " + library + " to " + version.to_string() + ": version already set to " + library_versions[library].to_string());
    library_versions[library] = version;
  }

  static bool dummy = [](){
    SetLibraryVersion("netgen", NETGEN_VERSION);
    return true;
  }();
} // namespace ngcore

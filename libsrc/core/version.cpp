
#include "archive.hpp"
#include "version.hpp"

namespace ngcore
{
  void VersionInfo :: DoArchive(Archive& ar)
  {
    ar & mayor_ & minor_ & release & patch & git_hash;
  }
} // namespace ngcore


#include "archive.hpp"
#include "register_archive.hpp"
#include "version.hpp"

#ifndef WIN32
#include <cxxabi.h>
#endif

namespace ngcore
{
  std::map<std::string, detail::ClassArchiveInfo> & GetTypeRegister()
  {
    static std::map<std::string, detail::ClassArchiveInfo> type_register;
    return type_register;
  }
  
  const detail::ClassArchiveInfo& Archive :: GetArchiveRegister(const std::string& classname)
  {
    // if(type_register == nullptr) type_register =
    // std::make_unique<std::map<std::string, detail::ClassArchiveInfo>>();
    return GetTypeRegister()[classname];
  }
  void Archive :: SetArchiveRegister(const std::string& classname, const detail::ClassArchiveInfo& info)
  {
    // if(type_register == nullptr) type_register =
    // std::make_unique<std::map<std::string, detail::ClassArchiveInfo>>();
    GetTypeRegister()[classname] = info;
  }
  bool Archive :: IsRegistered(const std::string& classname)
  {
    // if(type_register == nullptr) type_register =
    // std::make_unique<std::map<std::string, detail::ClassArchiveInfo>>();
    return GetTypeRegister().count(classname) != 0;
  }

  size_t Archive :: BytesLeft ()
  {
    if (auto bin = dynamic_cast<BinaryInArchive*> (this))
      return bin->StreamBytesLeft();
    return size_t(-1);
  }

  void Archive :: CheckSize (size_t n, size_t elsize)
  {
    if (elsize == 0) elsize = 1;
    if (is_output || n <= 65536/elsize) return;   // small ones are caught at the end of the archive
    size_t left = BytesLeft();
    if (left == size_t(-1) || n <= left/elsize) return;
    throw Exception("Archive: block of " + std::to_string(n) + " entries of " + std::to_string(elsize) +
                    " bytes, but only " + std::to_string(left) + " bytes are left. The archive was written by netgen " +
                    GetVersion("netgen").to_string() + ", is its format supported?");
  }

  void Archive :: ThrowEnd ()
  {
    throw Exception("Archive: read beyond the end, or invalid size. The archive was written by netgen " +
                    GetVersion("netgen").to_string() + ", is its format supported?");
  }

  Archive& Archive::Shallow(std::any& val)
  {
    if (shallow_to_python)
      {
        if (is_output)
          ShallowOutAny(val);
        else
          ShallowInAny(val);
      }
    return *this;
  }

} // namespace ngcore

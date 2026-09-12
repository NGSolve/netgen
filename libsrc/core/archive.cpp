
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

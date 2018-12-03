
#include "ngcore.hpp"

namespace ngcore
{
  std::map<std::string, ClassArchiveInfo>& GetArchiveRegister()
  {
    static std::map<std::string, ClassArchiveInfo> type_register;
    return type_register;
  }
}


#include "ngcore.hpp"

namespace ngcore
{
  std::map<std::string, std::function<void*()>>& GetArchiveRegister()
  {
    static std::map<std::string, std::function<void*()>> type_register = {};
    return type_register;
  }
}

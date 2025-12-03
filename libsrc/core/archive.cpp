
#include "archive.hpp"
#include "register_archive.hpp"
#include "version.hpp"

#ifndef WIN32
#include <cxxabi.h>
#endif

namespace ngcore
{
  // clang-tidy should ignore this static object
  // static std::map<std::string, detail::ClassArchiveInfo> type_register;  // NOLINT

  auto& GetTypeRegister()
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

#ifdef NETGEN_PYTHON
  pybind11::object CastAnyToPy(const std::any& a)
  {
    auto info = Archive::GetArchiveRegister(Demangle(a.type().name()));
    return info.anyToPyCaster(a);
  }

  std::any CastPyToAny(pybind11::object& obj)
  {
    auto name = Demangle(pybind11::detail::get_type_info((PyTypeObject*) pybind11::type::of(obj).ptr())->cpptype->name());
    auto info = Archive::GetArchiveRegister(name);
    if(!info.pyToAnyCaster)
      throw Exception("Need to register class " + name + " for Archive using std::any");
    return info.pyToAnyCaster(obj);
  }
#endif // NETGEN_PYTHON

} // namespace ngcore

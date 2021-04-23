#include "utils.hpp"
#include "logging.hpp"

#ifndef WIN32
#include <cxxabi.h>
#endif
#include <iostream>

namespace ngcore
{
  // parallel netgen
  int id = 0, ntasks = 1;

#ifdef WIN32
  // windows does demangling in typeid(T).name()
  NGCORE_API std::string Demangle(const char* typeinfo) {
      std::string name = typeinfo;
      // remove "class " and "struct " at beginning of type names to be consistent with demangled names of gcc/clang
      if(name.find("class ") == 0)
          name.erase(0,6);
      if(name.find("struct ") == 0)
          name.erase(0,7);
      return name;
  }
#else
  NGCORE_API std::string Demangle(const char* typeinfo)
  {
    int status=0;
    try
      {
        char *s = abi::__cxa_demangle(typeinfo, nullptr, nullptr, &status);
        std::string result{s};
        free(s);
        return result;
      }
    catch( const std::exception & e )
      {
        GetLogger("utils")->warn("{}:{} cannot demangle {}, status: {}, error:{}", __FILE__, __LINE__, typeinfo, status, e.what());
      }
    return typeinfo;
  }
#endif

  double seconds_per_tick = [] () noexcept
  {
      auto tick_start = GetTimeCounter();
      double tstart = WallTime();
      double tend = WallTime()+0.001;

      // wait for 1ms and compare wall time with time counter
      while(WallTime()<tend);

      auto tick_end = GetTimeCounter();
      tend = WallTime();

      return (tend-tstart)/static_cast<double>(tick_end-tick_start);
  }();

  const std::chrono::time_point<TClock> wall_time_start = TClock::now();

} // namespace ngcore


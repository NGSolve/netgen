#include "utils.hpp"
#include "logging.hpp"

#ifndef WIN32
#include <cxxabi.h>
#endif
#include <iostream>

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

  double ticks_per_second = [] () noexcept
  {
      auto tick_start = GetTimeCounter();
      double tstart = WallTime();
      double tend = WallTime()+0.001;

      // wait for 1ms and compare wall time with time counter
      while(WallTime()<tend);

      auto tick_end = GetTimeCounter();
      tend = WallTime();

      return (tick_end-tick_start)/(tend-tstart);
  }();

  const std::chrono::time_point<TClock> wall_time_start = TClock::now();

} // namespace ngcore

#endif


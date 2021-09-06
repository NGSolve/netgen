#include "utils.hpp"
#include "logging.hpp"

#ifndef WIN32
#include <cxxabi.h>
#endif
#include <array>
#include <iostream>
#include <regex>

#include "ngstream.hpp"

namespace ngcore
{
    namespace detail
    {
        // see https://github.com/RobotLocomotion/drake/blob/master/common/nice_type_name.cc
        static const auto demangle_regexes =
            std::array<std::pair<std::regex, std::string>, 8>{
                // Remove unwanted keywords and following space. (\b is word boundary.)
                std::make_pair(std::regex("\\b(class|struct|enum|union) "), ""),
                // Tidy up anonymous namespace.
                {std::regex("[`(]anonymous namespace[')]"), "(anonymous)"},
                // Replace Microsoft __int64 with long long.
                {std::regex("\\b__int64\\b"), "long long"},
                // Temporarily replace spaces we want to keep with "!". (\w is
                // alphanumeric or underscore.)
                {std::regex("(\\w) (\\w)"), "$1!$2"},
                {std::regex(" "), ""},  // Delete unwanted spaces.
                // Some compilers throw in extra namespaces like "__1" or "__cxx11".
                // Delete them.
                {std::regex("\\b__[[:alnum:]_]+::"), ""},
                {std::regex("!"), " "},  // Restore wanted spaces.

                // Recognize std::string's full name and abbreviate.
                {std::regex("\\bstd::basic_string<char,std::char_traits<char>,"
                        "std::allocator<char>>"), "std::string"}
            };
        std::string CleanupDemangledName( std::string s )
        {
            for(const auto & [r, sub] : demangle_regexes)
                s = std::regex_replace (s,r,sub);

            return s;
        }
    } // namespace detail

  // parallel netgen
  int id = 0, ntasks = 1;

#ifdef WIN32
  // windows does demangling in typeid(T).name()
  NGCORE_API std::string Demangle(const char* typeinfo) {
      std::string name = typeinfo;
      return detail::CleanupDemangledName(name);
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
        result = detail::CleanupDemangledName(result);
        return result;
      }
    catch( const std::exception & e )
      {
        GetLogger("utils")->warn("{}:{} cannot demangle {}, status: {}, error:{}", __FILE__, __LINE__, typeinfo, status, e.what());
      }
    std::string name = typeinfo;
    return detail::CleanupDemangledName(name);
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

  int printmessage_importance = 0;
  bool NGSOStream :: glob_active = true;

} // namespace ngcore


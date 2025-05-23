#ifndef NETGEN_CORE_LOGGING_HPP
#define NETGEN_CORE_LOGGING_HPP

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "exception.hpp"
#include "ngcore_api.hpp"
#include "utils.hpp"

#ifndef NETGEN_DEBUG_LOG
#define NETGEN_DEBUG_LOG(logger, ...)
#endif  // NETGEN_DEBUG_LOG

namespace spdlog
{
  class logger;
} // namespace spdlog

namespace ngcore
{
  NGCORE_API extern std::ostream* testout; // NOLINT
  
  namespace level
  {
    enum level_enum
    {
      trace = 0,
      debug = 1,
      info = 2,
      warn = 3,
      err = 4,
      critical = 5,
      off = 6
    };
  } // namespace level

  class Logger
  {
    static NGCORE_API level::level_enum global_level;

  public:
    static auto SetGlobalLoggingLevel( level::level_enum level )
    {
      auto oldval = global_level;
      global_level = level;
      return oldval;
    }

    std::shared_ptr<spdlog::logger> logger;

    Logger(std::shared_ptr<spdlog::logger> l) : logger(std::move(l)) {}

    void NGCORE_API log( level::level_enum level, std::string && s);

    template<typename T>
    std::string replace(std::string s, const T & t)
    {
      auto p0 = s.find_first_of('{');
      auto p1 = s.find_first_of('}', p0);
      if(p0==std::string::npos || p1==std::string::npos)
        throw Exception("invalid format string");
      s.replace(p0, p1-p0+1, ToString(t));
      return s;
    }

    std::string log_helper(std::string s)
    {
      return s;
    }

    template<typename T>
    std::string log_helper(std::string s,  const T &t)
    {
      return replace(s,t);
    }

    template<typename T, typename ... Args>
    std::string log_helper( std::string s, const T &t, Args ... args)
    {
      return log_helper(replace(s,t), args...);
    }

    template<typename ... Args>
    void log( level::level_enum level, const char* str, Args ... args)
    {
      log(level, log_helper(std::string(str), args...));
    }

    template<typename ... Args>
    void trace( const char* str, Args ... args) { log(level::level_enum::trace, str, args...); }
    template<typename ... Args>
    void debug( const char* str, Args ... args) { log(level::level_enum::debug, str, args...); }
    template<typename ... Args>
    void info( const char* str, Args ... args) { log(level::level_enum::info, str, args...); }
    template<typename ... Args>
    void warn( const char* str, Args ... args) { log(level::level_enum::warn, str, args...); }
    template<typename ... Args>
    void error( const char* str, Args ... args) { log(level::level_enum::err, str, args...); }
    template<typename ... Args>
    void critical( const char* str, Args ... args) { log(level::level_enum::critical, str, args...); }
  };




  NGCORE_API std::shared_ptr<Logger> GetLogger(const std::string& name);
  NGCORE_API void SetLoggingLevel(level::level_enum level, const std::string& name);
  NGCORE_API void AddFileSink(const std::string& filename, level::level_enum level, const std::string& logger);
  NGCORE_API void AddConsoleSink(level::level_enum level, const std::string& logger);
  NGCORE_API void ClearLoggingSinks(const std::string& logger);
  NGCORE_API void FlushOnLoggingLevel(level::level_enum level, const std::string& logger);
} // namespace ngcore

#endif // NETGEN_CORE_LOGGING_HPP

#ifndef NETGEN_CORE_LOGGING_HPP
#define NETGEN_CORE_LOGGING_HPP

#include <memory>
#include <string>
#include <vector>

#include "ngcore_api.hpp"

#ifdef NETGEN_USE_SPDLOG

#ifdef NETGEN_LOG_DEBUG
#define SPDLOG_DEBUG_ON
#endif // NETGEN_LOG_DEBUG

#include <spdlog/logger.h>
#include <spdlog/sinks/sink.h>
#include <spdlog/spdlog.h>

#define NETGEN_DEBUG_LOG(logger, ...) SPDLOG_DEBUG(logger, __VA_ARGS__)

#else // NETGEN_USE_SPDLOG

#include <iostream>

namespace spdlog
{
  // Dummys if Netgen is compiled with USE_SPDLOG=OFF.
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

  class logger
  {
  public:
    void log_helper() {}

    template<typename T>
    void log_helper( T t) { std::clog << t; }

    template<typename T, typename ... Args>
    void log_helper( T t, Args ... args)
    {
        std::clog << t;
        log_helper(args...);
        std::clog << ", ";
    }

    template<typename ... Args>
    void log( level::level_enum level, const char* fmt, Args ... args)
    {
        std::clog << level << ": " << fmt << "\t Arguments: ";
        log_helper(args...);
        std::clog << "\n";
    }

    template<typename ... Args>
    void trace( const char* fmt, Args ... args) { log(level::level_enum::trace, fmt, args...); }
    template<typename ... Args>
    void debug( const char* fmt, Args ... args) { log(level::level_enum::debug, fmt, args...); }
    template<typename ... Args>
    void info( const char* fmt, Args ... args) { log(level::level_enum::info, fmt, args...); }
    template<typename ... Args>
    void warn( const char* fmt, Args ... args) { log(level::level_enum::warn, fmt, args...); }
    template<typename ... Args>
    void error( const char* fmt, Args ... args) { log(level::level_enum::err, fmt, args...); }
    template<typename ... Args>
    void critical( const char* fmt, Args ... args) { log(level::level_enum::critical, fmt, args...); }
  };

  namespace sinks
  {
    class sink {};
  } // namespace sinks

} //namespace spdlog

#define NETGEN_DEBUG_LOG(logger, ...)

#endif // NETGEN_USE_SPDLOG

namespace ngcore
{
  namespace detail
  {
    std::vector<std::shared_ptr<spdlog::sinks::sink>>& GetDefaultSinks();
    inline std::shared_ptr<spdlog::logger> CreateDefaultLogger(const std::string& name);
  } //namespace detail

  NGCORE_API std::shared_ptr<spdlog::logger> GetLogger(const std::string& name);
  NGCORE_API void SetLoggingLevel(spdlog::level::level_enum level, const std::string& name);
  NGCORE_API void AddFileSink(const std::string& filename, spdlog::level::level_enum level, const std::string& logger);
  NGCORE_API void AddConsoleSink(spdlog::level::level_enum level, const std::string& logger);
  NGCORE_API void ClearLoggingSinks(const std::string& logger);
  NGCORE_API void FlushOnLoggingLevel(spdlog::level::level_enum level, const std::string& logger);
} // namespace ngcore

#endif // NETGEN_CORE_LOGGING_HPP

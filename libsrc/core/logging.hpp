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
  class logger
  {
  public:
    template<typename T>
    void trace(const T& /*unused*/) {}
    template<typename T>
    void debug(const T& /*unused*/) {}
    template<typename T>
    void info(const T& text) { std::cout << text << std::endl; }
    template<typename T>
    void warn(const T& text) { std::cout << text << std::endl; }
    template<typename T>
    void error(const T& text) { std::cout << text << std::endl; }
    template<typename T>
    void critical(const T& text) { std::cout << text << std::endl; }
  };

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

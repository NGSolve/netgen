#include "logging.hpp"

#include <iostream>

namespace ngcore
{
  std::ostream* testout = new std::ostream(nullptr); // NOLINT

  level::level_enum Logger::global_level = level::warn;

  void Logger::log(level::level_enum level, std::string && s)
  {
    if(level>=global_level)
      std::clog << s << '\n';
  }

} //namespace ngcore

namespace spdlog
{
  class logger
  {
    public:
    logger() = default;
  };
} // namespace spdlog

namespace ngcore
{

  // Dummy functions if no spdlog is available

  std::shared_ptr<Logger> GetLogger(const std::string& /*unused*/)
  {
    return std::make_shared<Logger>(std::make_shared<spdlog::logger>());
  }

  void SetLoggingLevel(level::level_enum level, const std::string& /*unused*/)
  {
    Logger::SetGlobalLoggingLevel(level);
  }

  void AddFileSink(const std::string& /*unused*/, level::level_enum /*unused*/,
      const std::string& /*unused*/)
  {}
  void AddConsoleSink(level::level_enum /*unused*/, const std::string& /*unused*/) {}
  void ClearLoggingSinks(const std::string& /*unused*/) {}
  void FlushOnLoggingLevel(level::level_enum /*unused*/, const std::string& /*unused*/) {}
} //namespace ngcore

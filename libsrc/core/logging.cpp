#include "logging.hpp"

#ifdef NETGEN_USE_SPDLOG

#include <spdlog/spdlog.h>
#include <spdlog/sinks/ansicolor_sink.h>
#include <spdlog/sinks/basic_file_sink.h>

#else // NETGEN_USE_SPDLOG
#include <iostream>
#endif // NETGEN_USE_SPDLOG


namespace ngcore
{
  std::ostream* testout = new std::ostream(nullptr); // NOLINT

  level::level_enum Logger::global_level = level::warn;

  void Logger::log(level::level_enum level, std::string && s)
  {
#ifdef NETGEN_USE_SPDLOG
    logger->log(spdlog::level::level_enum(level), s);
#else // NETGEN_USE_SPDLOG
    if(level>=global_level)
      std::clog << s << '\n';
#endif // NETGEN_USE_SPDLOG
  }

#ifdef NETGEN_USE_SPDLOG
  namespace detail
  {
    std::vector<std::shared_ptr<spdlog::sinks::sink>>& GetDefaultSinks()
    {
      static std::vector<std::shared_ptr<spdlog::sinks::sink>> sinks =
        { std::make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>() };
      return sinks;
    }
    std::shared_ptr<spdlog::logger> CreateDefaultLogger(const std::string& name)
    {
      auto& default_sinks = GetDefaultSinks();
      auto logger = std::make_shared<spdlog::logger>(name, default_sinks.begin(), default_sinks.end());
      spdlog::details::registry::instance().register_and_init(logger);
      return logger;
    }
  } // namespace detail

  std::shared_ptr<Logger> GetLogger(const std::string& name)
  {
    auto logger = spdlog::get(name);
    if(!logger)
      logger = detail::CreateDefaultLogger(name);
    return std::make_shared<Logger>(logger);
  }

  void SetLoggingLevel(spdlog::level::level_enum level, const std::string& name)
  {
    if(!name.empty())
      spdlog::get(name)->set_level(level);
    else
      spdlog::set_level(level);
  }

  void AddFileSink(const std::string& filename, spdlog::level::level_enum level, const std::string& logger)
  {
    auto sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename);
    sink->set_level(level);
    if(!logger.empty())
      GetLogger(logger)->logger->sinks().push_back(sink);
    else
      {
        detail::GetDefaultSinks().push_back(sink);
        spdlog::details::registry::instance().apply_all([sink](auto logger) { logger->sinks().push_back(sink); });
      }
  }

  void AddConsoleSink(spdlog::level::level_enum level, const std::string& logger)
  {
    auto sink = std::make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>();
    sink->set_level(level);
    if(!logger.empty())
      GetLogger(logger)->logger->sinks().push_back(sink);
    else
      {
        detail::GetDefaultSinks().push_back(sink);
        spdlog::details::registry::instance().apply_all([sink](auto logger) { logger->sinks().push_back(sink); });
      }
  }

  void ClearLoggingSinks(const std::string& logger)
  {
    if(!logger.empty())
      GetLogger(logger)->logger->sinks().clear();
    else
      {
        detail::GetDefaultSinks().clear();
        spdlog::details::registry::instance().apply_all([](auto logger) { logger->sinks().clear(); });
      }
  }

  void FlushOnLoggingLevel(spdlog::level::level_enum level, const std::string& logger)
  {
    if(!logger.empty())
      GetLogger(logger)->logger->flush_on(level);
    else
      spdlog::flush_on(level);
  }

#else // NETGEN_USE_SPDLOG
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

#endif // NETGEN_USE_SPDLOG

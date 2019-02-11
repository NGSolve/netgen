
#include <pybind11/pybind11.h>

#include "logging.hpp"

namespace py = pybind11;
using namespace ngcore;

PYBIND11_MODULE(pyngcore, m) // NOLINT
{
  py::enum_<level::level_enum>(m, "LOG_LEVEL", "Logging level")
    .value("Trace", level::trace)
    .value("Debug", level::debug)
    .value("Info", level::info)
    .value("Warn", level::warn)
    .value("Error", level::err)
    .value("Critical", level::critical)
    .value("Off", level::off);

  m.def("SetLoggingLevel", &SetLoggingLevel, py::arg("level"), py::arg("logger")="",
        "Set logging level, if name is given only to the specific logger, else set the global logging level");
  m.def("AddFileSink", &AddFileSink, py::arg("filename"), py::arg("level"), py::arg("logger")="",
        "Add File sink, either only to logger specified or globally to all loggers");
  m.def("AddConsoleSink", &AddConsoleSink, py::arg("level"), py::arg("logger")="",
        "Add console output for specific logger or all if none given");
  m.def("ClearLoggingSinks", &ClearLoggingSinks, py::arg("logger")="",
        "Clear sinks of specific logger, or all if none given");
  m.def("FlushOnLoggingLevel", &FlushOnLoggingLevel, py::arg("level"), py::arg("logger")="",
        "Flush every message with level at least `level` for specific logger or all loggers if none given.");
}

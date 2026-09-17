#ifndef NETGEN_CORE_MEMTRACER_HPP
#define NETGEN_CORE_MEMTRACER_HPP

#include <iostream>
#include <string>

#include "ngcore_api.hpp"

namespace ngcore
{
  NGCORE_API void PrintMemoryUsage(const char * file = nullptr, int line = 0, std::string msg = "", std::ostream & out = std::cout);
  NGCORE_API size_t GetRSSMemory();
  NGCORE_API size_t GetPageSize();
} // namespace ngcore

#endif // NETGEN_CORE_MEMTRACER_HPP

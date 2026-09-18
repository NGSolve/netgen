#ifndef NETGEN_CORE_MEMTRACE_HPP
#define NETGEN_CORE_MEMTRACE_HPP

#include <cstddef>

#include "ngcore_api.hpp"

namespace ngcore
{
  extern NGCORE_API class PajeTrace * memtrace;

  NGCORE_API void MemTraceRecord (const void * p, size_t bytes, unsigned char kind);

  NETGEN_INLINE void MemTraceAlloc (const void * p, size_t bytes)
  { if(memtrace) MemTraceRecord(p, bytes, 0); }
  NETGEN_INLINE void MemTraceFree (const void * p, size_t bytes)
  { if(memtrace) MemTraceRecord(p, bytes, 1); }
  NETGEN_INLINE void MemTraceDeviceAlloc (const void * p, size_t bytes)
  { if(memtrace) MemTraceRecord(p, bytes, 2); }
  NETGEN_INLINE void MemTraceDeviceFree (const void * p, size_t bytes)
  { if(memtrace) MemTraceRecord(p, bytes, 3); }
} // namespace ngcore

#endif // NETGEN_CORE_MEMTRACE_HPP

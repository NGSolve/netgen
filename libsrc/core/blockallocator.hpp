#ifndef NETGEN_CORE_BLOCKALLOCATOR_HPP
#define NETGEN_CORE_BLOCKALLOCATOR_HPP

#include <mutex>
#include "array.hpp"
#include "ngcore_api.hpp"

namespace ngcore
{
  /// allocates many fixed-size objects at once, keeps a free list of released ones (thread-safe)
  class BlockAllocator
  {
    size_t size;               // bytes per element
    size_t blocks;             // elements per block
    void * freelist = nullptr;
    Array<char*> bablocks;
    size_t nels = 0;
    std::mutex mut;
  public:
    NGCORE_API BlockAllocator (size_t asize, size_t ablocks = 100);
    NGCORE_API ~BlockAllocator ();
    NGCORE_API void * Alloc ();
    NGCORE_API void Free (void * p);
    size_t NumElements () const { return nels; }
  };
}

NETGEN_INLINE void * operator new (size_t /* size */, ngcore::BlockAllocator & ball)
{
  return ball.Alloc();
}

NETGEN_INLINE void operator delete (void * p, ngcore::BlockAllocator & ball)
{
  ball.Free (p);
}

#endif // NETGEN_CORE_BLOCKALLOCATOR_HPP

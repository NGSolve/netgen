#include "blockallocator.hpp"
#include "memtrace.hpp"

namespace ngcore
{
  BlockAllocator :: BlockAllocator (size_t asize, size_t ablocks)
    : size(std::max(asize, sizeof(void*))), blocks(ablocks)
  { }

  BlockAllocator :: ~BlockAllocator ()
  {
    std::lock_guard<std::mutex> guard(mut);
    for (char * block : bablocks)
      {
        MemTraceFree(block, size * blocks);
        delete [] block;
      }
    bablocks.SetSize(0);
  }

  void * BlockAllocator :: Alloc ()
  {
    std::lock_guard<std::mutex> guard(mut);
    if (!freelist)
      {
        char * hcp = new char [size * blocks];
        MemTraceAlloc(hcp, size * blocks);
        bablocks.Append (hcp);
        for (size_t i = 0; i < blocks-1; i++)
          *(void**)&(hcp[i * size]) = &(hcp[(i+1) * size]);
        *(void**)&(hcp[(blocks-1)*size]) = nullptr;
        freelist = hcp;
      }
    void * p = freelist;
    freelist = *(void**)freelist;
    nels++;
    return p;
  }

  void BlockAllocator :: Free (void * p)
  {
    std::lock_guard<std::mutex> guard(mut);
    if (bablocks.Size())
      {
        *(void**)p = freelist;
        freelist = p;
        nels--;
      }
  }
}

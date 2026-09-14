#ifndef FILE_OPTMEM
#define FILE_OPTMEM

/**************************************************************************/
/* File:   optmem.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

#include <mydefs.hpp>
#include <core/array.hpp>


namespace netgen
{
  using namespace ngcore;

/** 
    Optimized Memory allocation classes
*/

class BlockAllocator
{
private:
  ///
  unsigned size, blocks;
  ///
  void * freelist;
  ///
  Array<char*> bablocks;
  mutex block_allocator_mutex;
public:
  ///
  DLL_HEADER BlockAllocator (unsigned asize, unsigned ablocks = 100);
  ///
  DLL_HEADER ~BlockAllocator ();
  ///
  DLL_HEADER void * Alloc ();
  ///
  DLL_HEADER void Free (void * p);
};

}

#endif

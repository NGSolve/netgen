/**************************************************************************/
/* File:   localheap.cpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2002                                                  */
/**************************************************************************/

#include <exception>
#include <string>

#include "localheap.hpp"
#include "taskmanager.hpp"

namespace ngcore
{

  LocalHeap :: LocalHeap (size_t asize, const char * aname, bool mult_by_threads)
  {
    if (mult_by_threads)
      asize *= TaskManager::GetMaxThreads();
    totsize = asize;
    try
      {
        data = new char[asize];
      }
    catch (std::exception & e)
      {
        throw Exception (ToString ("Could not allocate localheap, heapsize = ") + ToString(asize));
      }

    next = data + totsize;
    p = data;
    owner = true;
    name = aname;
    CleanUp();   // align pointer
  }

  LocalHeap LocalHeap :: Split() const
  {
    int pieces = TaskManager::GetNumThreads();
    int i = TaskManager::GetThreadId();
    size_t freemem = totsize - (p - data);
    size_t size_of_piece = freemem / pieces;
    return LocalHeap (p + i * size_of_piece, size_of_piece, name);
  }

  void LocalHeap :: ThrowException() // throw (LocalHeapOverflow)
  {
    /*
    cout << "allocated: " << (p-data) << endl;
    cout << "throw LocalHeapOverflow, totsize = "<< totsize << endl;
    cout << "heap name = " << name << endl;
    */
    throw LocalHeapOverflow(totsize, name);
  }

  size_t tl_heap_size = 50*1000*1000;
  thread_local LocalHeap tl_heap(tl_heap_size, "tlheap");

  LocalHeap& TLHeap()
  {
    return tl_heap;
  }

  
  void SetTLHeapSize(size_t s)
  {
    tl_heap_size = s;
    // if (tl_heap.Available() == tl_heap.Size()) // check if not in used
    tl_heap = LocalHeap(tl_heap_size, "tlheap");
  }

  

  LocalHeapOverflow :: LocalHeapOverflow (size_t size, const char *name)
    : Exception("Local Heap overflow\n")
  {
    std::stringstream str;
    if(name)
        str << "\tName: " << name << '\n';
    str << "\tSize: " << size << '\n';
    Append (str.str());
    // Append ("please use 'define constant heapsize = xxx' with larger value\n");
  }
  
  LocalHeapOverflow :: ~LocalHeapOverflow ()
  {
    ; 
  }

}


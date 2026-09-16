#ifndef FILE_PARTHREADS
#define FILE_PARTHREADS

/**************************************************************************/
/* File:   parthreads.hh                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   22. Nov. 2000                                                  */
/**************************************************************************/

/*
  Parallel thread
*/
#include <functional>

namespace netgen
{



  
  typedef void (*NgTaskManager)(std::function<void(int,int)>);
  typedef void (*NgTracer)(std::string, bool);  // false .. start, true .. stop

  inline void DummyTaskManager (std::function<void(int,int)> func)
  {
    func(0,2);
    func(1,2);
  }

  inline void DummyTracer (std::string, bool) { ; }
  
  template <typename FUNC>
  inline void ParallelFor (NgTaskManager tm, size_t n, FUNC func)
  {
    (*tm) ([n,func] (size_t nr, size_t nums)
           {
             size_t begin = nr*n / nums;
             size_t end = (nr+1)*n / nums;

             for (size_t i = begin; i < end; i++)
               func(i);
           });
  }
  
  template <typename FUNC>
  inline void ParallelForRange (NgTaskManager tm, size_t n, FUNC func)
  {
    (*tm) ([n,func] (size_t nr, size_t nums)
           {
             size_t begin = nr*n / nums;
             size_t end = (nr+1)*n / nums;
             func(begin, end);
           });
  }
                    

  
}

#endif

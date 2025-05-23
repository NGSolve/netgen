#ifndef FILE_PARTHREADS
#define FILE_PARTHREADS

/**************************************************************************/
/* File:   parthreads.hh                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   22. Nov. 2000                                                  */
/**************************************************************************/

/*
  Parallel thread, Mutex,
*/
#include <functional>

namespace netgen
{

#ifdef NO_PARALLEL_THREADS

class NgMutex { };

class NgLock
{
public:
  NgLock (NgMutex & mut, bool lock = 0) { ; }
  void Lock () { ; }
  void UnLock () { ; }
};

#else

typedef std::mutex NgMutex;

class NgLock
{
  NgMutex & mut;
  bool locked;
public:
  NgLock (NgMutex & ngmut, bool lock = false)
    : mut (ngmut)
  {
    if (lock)
      mut.lock();

    locked = lock;
  };

  ~NgLock()
  {
    if (locked)
      mut.unlock();
  }

  void Lock ()
  {
    mut.lock();
    locked = true;
  }
  void UnLock ()
  {
    mut.unlock();
    locked = false;
  }
  /*
  int TryLock ()
  {
    return mut.try_lock();
  }
  */
};


#endif


// Simple ParallelFor function to replace OpenMP
template<typename TFunc>
void ParallelFor( int first, int next, const TFunc & f )
{
  int nthreads = std::thread::hardware_concurrency();
  std::thread * threads = new std::thread[nthreads];
  for (int i=0; i<nthreads; i++)
    {
      int myfirst = first + (next-first)*i/nthreads;
      int mynext = first + (next-first)*(i+1)/nthreads;
      threads[i] = std::thread( [myfirst,mynext,&f] ()
        {
          f(myfirst, mynext);
        });
    }

  for (int i=0; i<nthreads; i++)
    threads[i].join();
  delete [] threads;
}


  
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

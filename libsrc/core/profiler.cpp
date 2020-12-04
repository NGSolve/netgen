#include <mutex>

#include "profiler.hpp"

namespace ngcore
{
  std::vector<NgProfiler::TimerVal> NgProfiler::timers(NgProfiler::SIZE); // NOLINT

  std::string NgProfiler::filename;

  std::array<size_t,NgProfiler::SIZE> NgProfiler::dummy_thread_times;
  size_t * NgProfiler::thread_times = NgProfiler::dummy_thread_times.data(); // NOLINT
  std::array<size_t,NgProfiler::SIZE> NgProfiler::dummy_thread_flops;
  size_t * NgProfiler::thread_flops = NgProfiler::dummy_thread_flops.data(); // NOLINT

  std::shared_ptr<Logger> NgProfiler::logger = GetLogger("Profiler"); // NOLINT

  NgProfiler :: NgProfiler()
  {
    for (auto & t : timers)
    {
        t.tottime = 0.0;
        t.usedcounter = 0;
        t.flops = 0.0;
    }
  }

  NgProfiler :: ~NgProfiler()
  {
    if (filename.length())
      {
        logger->debug( "write profile to file {}", filename );
        FILE *prof = fopen(filename.c_str(),"w"); // NOLINT
        Print (prof);
        fclose(prof); // NOLINT
      }

    if (getenv ("NGPROFILE"))
      {
       std::string filename = "netgen.prof";
#ifdef PARALLEL
       filename += "."+ToString(id);
#endif
       if (id == 0) logger->info( "write profile to file {}", filename );
       FILE *prof = fopen(filename.c_str(),"w"); // NOLINT
       Print (prof);
       fclose(prof); // NOLINT
      }
  }

  void NgProfiler :: Print (FILE * prof)
  {
    int i = 0;
    for (auto & t : timers)
    {
      if (t.count != 0 || t.usedcounter != 0)
        {
          fprintf(prof,"job %3i calls %8li, time %6.4f sec",i,t.count,t.tottime); // NOLINT
          if(t.flops != 0.0)
            fprintf(prof,", MFlops = %6.2f",t.flops / (t.tottime) * 1e-6); // NOLINT
          if(t.loads != 0.0)
            fprintf(prof,", MLoads = %6.2f",t.loads / (t.tottime) * 1e-6); // NOLINT
          if(t.stores != 0.0)
            fprintf(prof,", MStores = %6.2f",t.stores / (t.tottime) * 1e-6); // NOLINT
          if(t.usedcounter)
            fprintf(prof," %s",t.name.c_str()); // NOLINT
          fprintf(prof,"\n"); // NOLINT
        }
      i++;
    }
  }


  int NgProfiler :: CreateTimer (const std::string & name)
  {
    static std::mutex createtimer_mutex;
    int nr = -1;
    {
      std::lock_guard<std::mutex> guard(createtimer_mutex);
      for (int i = SIZE-1; i > 0; i--)
      {
        auto & t = timers[i];
        if (!t.usedcounter)
          {
            t.usedcounter = 1;
            t.name = name;
            nr = i;
            break;
          }
      }
    }
    if (nr > -1) return nr;
    static bool first_overflow = true;
    if (first_overflow)
      {
        first_overflow = false;
        NgProfiler::logger->warn("no more timer available, reusing last one");
      }
    return 0;
  }

  void NgProfiler :: Reset ()
  {
      for(auto & t : timers)
      {
          t.tottime = 0.0;
          t.count = 0;
          t.flops = 0.0;
          t.loads = 0;
          t.stores = 0;
      }
  }

  NgProfiler prof; // NOLINT

#ifdef NETGEN_TRACE_MEMORY
  std::vector<std::string> MemoryTracer::names{"all"};
  std::vector<int> MemoryTracer::parents{-1};
#endif // NETGEN_TRACE_MEMORY

} // namespace ngcore

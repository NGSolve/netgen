#ifndef NETGEN_CORE_PROFILER_HPP
#define NETGEN_CORE_PROFILER_HPP

#include <array>
#include <chrono>
#include <string>

#include "logging.hpp"
#include "paje_trace.hpp"
#include "utils.hpp"

namespace ngcore
{
  class NgProfiler
  {
  public:
    /// maximal number of timers
    enum { SIZE = 8*1024 };

    struct TimerVal
    {
        TimerVal() = default;

        double tottime = 0.0;
        TTimePoint starttime=0;
        double flops = 0.0;
        double loads = 0.0;
        double stores = 0.0;
        long count = 0;
        std::string name = "";
        int usedcounter = 0;
    };

    NGCORE_API static std::vector<TimerVal> timers;

    NGCORE_API static TTimePoint * thread_times;
    NGCORE_API static TTimePoint * thread_flops;
    NGCORE_API static std::shared_ptr<Logger> logger;
    NGCORE_API static std::array<size_t, NgProfiler::SIZE> dummy_thread_times;
    NGCORE_API static std::array<size_t, NgProfiler::SIZE> dummy_thread_flops;
  private:

    NGCORE_API static std::string filename;
  public:
    NgProfiler();
    ~NgProfiler();

    NgProfiler(const NgProfiler &) = delete;
    NgProfiler(NgProfiler &&) = delete;
    void operator=(const NgProfiler &) = delete;
    void operator=(NgProfiler &&) = delete;

    static void SetFileName (const std::string & afilename) { filename = afilename; }

    /// create new timer, use integer index
    NGCORE_API static int CreateTimer (const std::string & name);

    NGCORE_API static void Reset ();


    /// start timer of index nr
    static void StartTimer (int nr)
    {
      timers[nr].starttime = GetTimeCounter(); timers[nr].count++;
    }

    /// stop timer of index nr
    static void StopTimer (int nr)
    {
      timers[nr].tottime += (GetTimeCounter()-timers[nr].starttime)*seconds_per_tick;
    }

    static void StartThreadTimer (size_t nr, size_t tid)
    {
      thread_times[tid*SIZE+nr] -= GetTimeCounter(); // NOLINT
    }

    static void StopThreadTimer (size_t nr, size_t tid)
    {
      thread_times[tid*SIZE+nr] += GetTimeCounter(); // NOLINT
    }

    static void AddThreadFlops (size_t nr, size_t tid, size_t flops)
    {
      thread_flops[tid*SIZE+nr] += flops; // NOLINT
    }

    /// if you know number of flops, provide them to obtain the MFlop - rate
    static void AddFlops (int nr, double aflops) { timers[nr].flops += aflops; }
    static void AddLoads (int nr, double aloads) { timers[nr].loads += aloads; }
    static void AddStores (int nr, double astores) { timers[nr].stores += astores; }

    static int GetNr (const std::string & name)
    {
      for (int i = SIZE-1; i >= 0; i--)
        if (timers[i].name == name)
          return i;
      return -1;
    }

    static double GetTime (int nr)
    {
      return timers[nr].tottime;
    }

    static double GetTime (const std::string & name)
    {
      for (int i = SIZE-1; i >= 0; i--)
        if (timers[i].name == name)
          return GetTime (i);
      return 0;
    }

    static long int GetCounts (int nr)
    {
      return timers[nr].count;
    }

    static double GetFlops (int nr)
    {
      return timers[nr].flops;
    }

    /// change name
    static void SetName (int nr, const std::string & name) { timers[nr].name = name; }
    static std::string GetName (int nr) { return timers[nr].name; }
    /// print profile
    NGCORE_API static void Print (FILE * prof);

    class RegionTimer
    {
      int nr;
    public:
      /// start timer
      RegionTimer (int anr) : nr(anr) { NgProfiler::StartTimer(nr); }
      /// stop timer
      ~RegionTimer () { NgProfiler::StopTimer(nr); }

      RegionTimer() = delete;
      RegionTimer(const RegionTimer &) = delete;
      RegionTimer(RegionTimer &&) = delete;
      void operator=(const RegionTimer &) = delete;
      void operator=(RegionTimer &&) = delete;
    };
  };



  class NGCORE_API Timer
  {
    int timernr;
    int priority;
  public:
    Timer (const std::string & name, int apriority = 1)
      : priority(apriority)
    {
      timernr = NgProfiler::CreateTimer (name);
    }
    void SetName (const std::string & name)
    {
      NgProfiler::SetName (timernr, name);
    }
    void Start ()
    {
      if (priority <= 2)
	NgProfiler::StartTimer (timernr);
      if (priority <= 1)
        if(trace) trace->StartTimer(timernr);
    }
    void Stop ()
    {
      if (priority <= 2)
	NgProfiler::StopTimer (timernr);
      if (priority <= 1)
        if(trace) trace->StopTimer(timernr);
    }
    void AddFlops (double aflops)
    {
      if (priority <= 2)
	NgProfiler::AddFlops (timernr, aflops);
    }

    double GetTime () { return NgProfiler::GetTime(timernr); }
    long int GetCounts () { return NgProfiler::GetCounts(timernr); }
    double GetMFlops ()
    { return NgProfiler::GetFlops(timernr)
        / NgProfiler::GetTime(timernr) * 1e-6; }
    operator int () { return timernr; }
  };


  /**
     Timer object.
       Start / stop timer at constructor / destructor.
  */
  class RegionTimer
  {
    Timer & timer;
  public:
    /// start timer
    RegionTimer (Timer & atimer) : timer(atimer) { timer.Start(); }
    /// stop timer
    ~RegionTimer () { timer.Stop(); }

    RegionTimer() = delete;
    RegionTimer(const RegionTimer &) = delete;
    RegionTimer(RegionTimer &&) = delete;
    void operator=(const RegionTimer &) = delete;
    void operator=(RegionTimer &&) = delete;
  };

  class ThreadRegionTimer
  {
    size_t nr;
    size_t tid;
  public:
    /// start timer
    ThreadRegionTimer (size_t _nr, size_t _tid) : nr(_nr), tid(_tid)
    { NgProfiler::StartThreadTimer(nr, tid); }
    /// stop timer
    ~ThreadRegionTimer ()
    { NgProfiler::StopThreadTimer(nr, tid); }

    ThreadRegionTimer() = delete;
    ThreadRegionTimer(ThreadRegionTimer &&) = delete;
    ThreadRegionTimer(const ThreadRegionTimer &) = delete;
    void operator=(const ThreadRegionTimer &) = delete;
    void operator=(ThreadRegionTimer &&) = delete;
  };

  class RegionTracer
    {
      int nr;
      int thread_id;
    public:
      static constexpr int ID_JOB = PajeTrace::Task::ID_JOB;
      static constexpr int ID_NONE = PajeTrace::Task::ID_NONE;
      static constexpr int ID_TIMER = PajeTrace::Task::ID_TIMER;

      RegionTracer() = delete;
      RegionTracer(RegionTracer &&) = delete;
      RegionTracer(const RegionTracer &) = delete;
      void operator=(const RegionTracer &) = delete;
      void operator=(RegionTracer &&) = delete;

      /// start trace
      RegionTracer (int athread_id, int region_id, int id_type = ID_NONE, int additional_value = -1 )
        : thread_id(athread_id)
        {
	  if (trace)
          nr = trace->StartTask (athread_id, region_id, id_type, additional_value);
        }
      /// start trace with timer
      RegionTracer (int athread_id, Timer & timer, int additional_value = -1 )
        : thread_id(athread_id)
        {
	  if (trace)
          nr = trace->StartTask (athread_id, static_cast<int>(timer), ID_TIMER, additional_value);
        }

      /// set user defined value
      void SetValue( int additional_value )
      {
	  if (trace)
        trace->SetTask( thread_id, nr, additional_value );
      }

      /// stop trace
      ~RegionTracer ()
        {
	  if (trace)
          trace->StopTask (thread_id, nr);
        }
    };


  // Helper function for timings
  // Run f() at least min_iterations times until max_time seconds elapsed
  // returns minimum runtime for a call of f()
  template<typename TFunc>
  double RunTiming( TFunc f, double max_time = 0.5, int min_iterations = 10 )
  {
      // Make sure the whole test run does not exceed maxtime
      double tend = WallTime()+max_time;

      // warmup
      f();

      double tres = std::numeric_limits<double>::max();
      int iteration = 0;
      while(WallTime()<tend || iteration++ < min_iterations)
      {
          double t = -WallTime();
          f();
          t += WallTime();
          tres = std::min(tres, t);
      }

      return tres;
  }

} // namespace ngcore


#endif // NETGEN_CORE_PROFILER_HPP

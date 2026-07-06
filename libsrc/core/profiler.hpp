#ifndef NETGEN_CORE_PROFILER_HPP
#define NETGEN_CORE_PROFILER_HPP

#include <array>
#include <chrono>
#include <functional>
#include <string>

#include "array.hpp"
#include "logging.hpp"
#include "paje_trace.hpp"
#include "taskmanager.hpp"
#include "utils.hpp"

namespace ngcore
{
  class NgProfiler
  {
  public:
    static constexpr int SIZE = WorkerData::MAX_TIMERS;
    /// maximal number of timers

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

    NGCORE_API static std::shared_ptr<Logger> logger;
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
      StartTimer(nr, TaskManager::GetTimerThreadId());
    }

    /// stop timer of index nr
    static void StopTimer (int nr)
    {
      StopTimer(nr, TaskManager::GetTimerThreadId());
    }

    static void StartTimer (size_t nr, int tid)
    {
      if(tid == -1)
          return;
      else if(tid==-2)
      {
          timers[nr].starttime = GetTimeCounter();
          timers[nr].count++;
      }
      else
          TaskManager::GetWorkerData()->times[nr] -= GetTimeCounter();
    }

    static void StopTimer (size_t nr, int tid)
    {
      if(tid == -1)
          return;
      else if(tid==-2)
          timers[nr].tottime += (GetTimeCounter()-timers[nr].starttime)*seconds_per_tick;
      else
          TaskManager::GetWorkerData()->times[nr] += GetTimeCounter();
    }

    static void AddFlops (size_t nr, size_t flops, int tid)
    {
      if(tid == -1)
          return;
      else if(tid==-2)
          timers[nr].flops += flops;
      else
          TaskManager::GetWorkerData()->flops[nr] += flops;
    }

    /// if you know number of flops, provide them to obtain the MFlop - rate
    static void AddFlops (int nr, double aflops) { AddFlops(nr, aflops, TaskManager::GetTimerThreadId()); }

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
  };

  
  struct TNoTracing{ static constexpr bool do_tracing=false; };
  struct TTracing{ static constexpr bool do_tracing=true; };

  struct TNoTiming{ static constexpr bool do_timing=false; };
  struct TTiming{ static constexpr bool do_timing=true; };

  namespace detail {

      template<typename T>
      constexpr bool is_tracing_type_v = std::is_same_v<T, TNoTracing> || std::is_same_v<T, TTracing>;

      template<typename T>
      constexpr bool is_timing_type_v = std::is_same_v<T, TNoTiming> || std::is_same_v<T, TTiming>;
  }

  [[maybe_unused]] static TNoTracing NoTracing;
  [[maybe_unused]] static TNoTiming NoTiming;

  template<typename TTracing=TTracing, typename TTiming=TTiming>
  class Timer
  {
    int timernr;
    int Init( const std::string & name )
    {
      return NgProfiler::CreateTimer (name);
    }
  public:
    static constexpr bool do_tracing = TTracing::do_tracing;
    static constexpr bool do_timing = TTiming::do_timing;

    Timer (const std::string & name) : timernr(Init(name)) { }

    template<std::enable_if_t< detail::is_tracing_type_v<TTracing>, bool> = false>
    Timer( const std::string & name, TTracing ) : timernr(Init(name)) { }

    template<std::enable_if_t< detail::is_timing_type_v<TTiming>, bool> = false>
    Timer( const std::string & name, TTiming ) : timernr(Init(name)) { }

    Timer( const std::string & name, TTracing, TTiming ) : timernr(Init(name)) { }

    void SetName (const std::string & name)
    {
      NgProfiler::SetName (timernr, name);
    }
    void Start () const
    {
      Start(TaskManager::GetTimerThreadId());
    }
    void Stop () const
    {
      Stop(TaskManager::GetTimerThreadId());
    }
    void Start (int tid, int trace_value = -1) const
    {
      if constexpr(do_timing)
          NgProfiler::StartTimer(timernr, tid);
      if constexpr(do_tracing)
          if(trace) trace->StartTask (tid, timernr, PajeTrace::Task::ID_TIMER, trace_value);
    }
    void Stop (int tid) const
    {
      if constexpr(do_timing)
        NgProfiler::StopTimer(timernr, tid);
      if constexpr(do_tracing)
        if(trace) trace->StopTask (tid, timernr, PajeTrace::Task::ID_TIMER);
    }
    
    void AddFlops (double aflops) const
    {
      if constexpr(do_timing)
	NgProfiler::AddFlops (timernr, aflops);
    }
    
    void AddFlops (double aflops, int tid) const
    {
      if constexpr(do_timing)
	NgProfiler::AddFlops (timernr, aflops, tid);
    }

    double GetTime () { return NgProfiler::GetTime(timernr); }
    long int GetCounts () { return NgProfiler::GetCounts(timernr); }
    double GetMFlops ()
    { return NgProfiler::GetFlops(timernr)
        / NgProfiler::GetTime(timernr) * 1e-6; }
    operator int () const { return timernr; }
  };


  /**
     Timer object.
       Start / stop timer at constructor / destructor.
  */
  template<typename TTimer>
  struct RegionTimer
  {
    static_assert(!std::is_same_v<TTimer, int>, "RegionTimer should be used with Timer objects, not with timer indices");
    const TTimer & timer;
    const int tid;
    /// start timer
    RegionTimer (const TTimer & atimer, int trace_value = -1) 
        : timer(atimer), tid(atimer.do_timing ? TaskManager::GetTimerThreadId() : -1)
    {
      timer.Start(tid, trace_value);
    }

    /// stop timer
    ~RegionTimer () { timer.Stop(tid); }

    RegionTimer() = delete;
    RegionTimer(const RegionTimer &) = delete;
    RegionTimer(RegionTimer &&) = delete;
    void operator=(const RegionTimer &) = delete;
    void operator=(RegionTimer &&) = delete;
    void AddFlops (double aflops) { timer.AddFlops(aflops, tid); }
  };

  class RegionTracer
    {
      int nr;
      int thread_id;
      int type;
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
          trace->StartTask (athread_id, region_id, id_type, additional_value);
          type = id_type;
          nr = region_id;
        }
      /// start trace with timer
      template<typename TTimer>
      RegionTracer (int athread_id, TTimer & timer, int additional_value = -1 )
        : thread_id(athread_id)
        {
          nr = timer;
          type = ID_TIMER;
	  if (trace)
            trace->StartTask (athread_id, nr, type, additional_value);
        }

      /// stop trace
      ~RegionTracer ()
        {
	  if (trace)
            trace->StopTask (thread_id, nr, type);
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

// Helper macro to easily add multiple timers in a function for profiling
// Usage: NETGEN_TIMER_FROM_HERE("my_timer_name")
// Effect: define static Timer and RegionTimer with given name and line number
#define NETGEN_TOKEN_CONCAT(x, y) x ## y
#define NETGEN_TOKEN_CONCAT2(x, y) NETGEN_TOKEN_CONCAT(x, y)
#define NETGEN_TIMER_FROM_HERE(name) \
  static Timer NETGEN_TOKEN_CONCAT2(timer_, __LINE__)( string(name)+"_"+ToString(__LINE__)); \
  RegionTimer NETGEN_TOKEN_CONCAT2(rt_,__LINE__)(NETGEN_TOKEN_CONCAT2(timer_,__LINE__));


#endif // NETGEN_CORE_PROFILER_HPP

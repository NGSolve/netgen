#ifndef NETGEN_CORE_PAJE_TRACE_HPP
#define NETGEN_CORE_PAJE_TRACE_HPP

#include <limits>
#include <vector>

#include "logging.hpp"       // for logger
#include "ngcore_api.hpp"    // for NGCORE_API
#include "utils.hpp"

namespace ngcore
{

  extern NGCORE_API class PajeTrace *trace;
  class PajeTrace
    {
    public:
      using TClock = std::chrono::system_clock;

    protected:
      std::shared_ptr<Logger> logger = GetLogger("PajeTrace");
    private:
      NGCORE_API static size_t max_tracefile_size;
      NGCORE_API static bool trace_thread_counter;
      NGCORE_API static bool trace_threads;

      bool tracing_enabled;
      TTimePoint start_time;
      int nthreads;

    public:
      void WriteSunburstHTML();

      // Approximate number of events to trace. Tracing will
      // be stopped if any thread reaches this number of events
      unsigned int max_num_events_per_thread;

      static void SetTraceThreads( bool atrace_threads )
        {
          trace_threads = atrace_threads;
        }

      static void SetTraceThreadCounter( bool trace_threads )
        {
          trace_thread_counter = trace_threads;
        }

      static void SetMaxTracefileSize( size_t max_size )
        {
          max_tracefile_size = max_size;
        }

      std::string tracefile_name;

      struct Job
        {
          int job_id;
          const std::type_info *type;
          TTimePoint start_time;
          TTimePoint stop_time;
        };

      struct Task
        {
          int thread_id;

          int id;
          int id_type;

          int additional_value;

          TTimePoint start_time;
          TTimePoint stop_time;

          static constexpr int ID_NONE = -1;
          static constexpr int ID_JOB = 1;
          static constexpr int ID_TIMER = 2;
        };

      struct TimerEvent
        {
          int timer_id;
          TTimePoint time;
          bool is_start;
          int thread_id;

          bool operator < (const TimerEvent & other) const { return time < other.time; }
        };

      struct ThreadLink
        {
          int thread_id;
          int key;
          TTimePoint time;
          bool is_start;
          bool operator < (const ThreadLink & other) const { return time < other.time; }
        };

      std::vector<std::vector<Task> > tasks;
      std::vector<Job> jobs;
      std::vector<TimerEvent> timer_events;
      std::vector<std::vector<ThreadLink> > links;

    public:
      NGCORE_API void StopTracing();

      PajeTrace() = delete;
      PajeTrace(const PajeTrace &) = delete;
      PajeTrace(PajeTrace &&) = delete;
      NGCORE_API PajeTrace(int anthreads, std::string aname = "");
      NGCORE_API ~PajeTrace();

      void operator=(const PajeTrace &) = delete;
      void operator=(PajeTrace &&) = delete;

      void StartTimer(int timer_id)
        {
          if(!tracing_enabled) return;
          if(unlikely(timer_events.size() == max_num_events_per_thread))
            StopTracing();
          timer_events.push_back(TimerEvent{timer_id, GetTimeCounter(), true});
        }

      void StopTimer(int timer_id)
        {
          if(!tracing_enabled) return;
          if(unlikely(timer_events.size() == max_num_events_per_thread))
            StopTracing();
          timer_events.push_back(TimerEvent{timer_id, GetTimeCounter(), false});
        }

      NETGEN_INLINE int StartTask(int thread_id, int id, int id_type = Task::ID_NONE, int additional_value = -1)
        {
          if(!tracing_enabled) return -1;
          if(!trace_threads && !trace_thread_counter) return -1;
	  if(unlikely(tasks[thread_id].size() == max_num_events_per_thread))
            StopTracing();
          int task_num = tasks[thread_id].size();
          tasks[thread_id].push_back( Task{thread_id, id, id_type, additional_value, GetTimeCounter()} );
          return task_num;
        }

      void StopTask(int thread_id, int task_num)
        {
          if(!trace_threads && !trace_thread_counter) return;
          if(task_num>=0)
            tasks[thread_id][task_num].stop_time = GetTimeCounter();
        }

      void SetTask(int thread_id, int task_num, int additional_value) {
          if(!trace_threads && !trace_thread_counter) return;
          if(task_num>=0)
            tasks[thread_id][task_num].additional_value = additional_value;
      }

      void StartJob(int job_id, const std::type_info & type)
        {
          if(!tracing_enabled) return;
          if(jobs.size() == max_num_events_per_thread)
            StopTracing();
          jobs.push_back( Job{job_id, &type, GetTimeCounter()} );
        }

      void StopJob()
        {
          if(tracing_enabled)
            jobs.back().stop_time = GetTimeCounter();
        }

      void StartLink(int thread_id, int key)
        {
          if(!tracing_enabled) return;
          if(links[thread_id].size() == max_num_events_per_thread)
            StopTracing();
          links[thread_id].push_back( ThreadLink{thread_id, key, GetTimeCounter(), true} );
        }

      void StopLink(int thread_id, int key)
        {
          if(!tracing_enabled) return;
          if(links[thread_id].size() == max_num_events_per_thread)
            StopTracing();
          links[thread_id].push_back( ThreadLink{thread_id, key, GetTimeCounter(), false} );
        }

      void Write( const std::string & filename );

    };
} // namespace ngcore

#endif // NETGEN_CORE_PAJE_TRACE_HPP

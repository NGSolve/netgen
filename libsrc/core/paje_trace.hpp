#ifndef NETGEN_CORE_PAJE_TRACE_HPP
#define NETGEN_CORE_PAJE_TRACE_HPP

#include <algorithm>
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
      NGCORE_API static bool mem_tracing_enabled;

      bool tracing_enabled;
      TTimePoint start_time;
      int nthreads;
      size_t n_memory_events_at_start;

    public:
      NGCORE_API void WriteTimingChart();
#ifdef NETGEN_TRACE_MEMORY
      NGCORE_API void WriteMemoryChart( std::string fname );
#endif // NETGEN_TRACE_MEMORY

      // Approximate number of events to trace. Tracing will
      // be stopped if any thread reaches this number of events
      unsigned int max_num_events_per_thread;

      static void SetTraceMemory( bool trace_memory )
        {
          mem_tracing_enabled = trace_memory;
        }

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

          TTimePoint time;
          bool is_start;

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

      struct UserEvent
        {
          TTimePoint t_start = 0, t_end = 0;
          std::string data = "";
          int container = 0;
          int id = 0;

          bool operator < (const UserEvent & other) const { return t_start < other.t_start; }
        };

      struct ThreadLink
        {
          int thread_id;
          int key;
          TTimePoint time;
          bool is_start;
          bool operator < (const ThreadLink & other) const { return time < other.time; }
        };

      struct MemoryEvent
        {
          TTimePoint time;
          size_t size;
          int id;
          bool is_alloc;

          bool operator < (const MemoryEvent & other) const { return time < other.time; }
        };

      std::vector<std::vector<Task> > tasks;
      std::vector<Job> jobs;
      std::vector<TimerEvent> timer_events;
      std::vector<UserEvent> user_events;
      std::vector<std::tuple<std::string, int>> user_containers;
      std::vector<TimerEvent> gpu_events;
      std::vector<std::vector<ThreadLink> > links;
      NGCORE_API static std::vector<MemoryEvent> memory_events;

    public:
      NGCORE_API void StopTracing();

      PajeTrace() = delete;
      PajeTrace(const PajeTrace &) = delete;
      PajeTrace(PajeTrace &&) = delete;
      NGCORE_API PajeTrace(int anthreads, std::string aname = "");
      NGCORE_API ~PajeTrace();

      void operator=(const PajeTrace &) = delete;
      void operator=(PajeTrace &&) = delete;

      int AddUserContainer(std::string name, int parent=-1)
      {
        if(auto pos = std::find(user_containers.begin(), user_containers.end(), std::tuple{name,parent}); pos != user_containers.end())
          return pos - user_containers.begin();
        int id = user_containers.size();
        user_containers.push_back({name, parent});
        return id;
      }

      void AddUserEvent(UserEvent ue)
      {
          if(!tracing_enabled) return;
          user_events.push_back(ue);
      }
      void StartGPU(int timer_id = 0)
        {
          if(!tracing_enabled) return;
          if(unlikely(gpu_events.size() == max_num_events_per_thread))
            StopTracing();
          gpu_events.push_back(TimerEvent{timer_id, GetTimeCounter(), true});
        }

      void StopGPU(int timer_id)
        {
          if(!tracing_enabled) return;
          if(unlikely(gpu_events.size() == max_num_events_per_thread))
            StopTracing();
          gpu_events.push_back(TimerEvent{timer_id, GetTimeCounter(), false});
        }

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

      void AllocMemory(int id, size_t size)
        {
          if(!mem_tracing_enabled) return;
          memory_events.push_back(MemoryEvent{GetTimeCounter(), size, id, true});
        }

      void FreeMemory(int id, size_t size)
        {
          if(!mem_tracing_enabled) return;
          memory_events.push_back(MemoryEvent{GetTimeCounter(), size, id, false});
        }

      void ChangeMemory(int id, long long size)
        {
          if(size>0)
            AllocMemory(id, size);
          if(size<0)
            FreeMemory(id, -size);
        }


      int StartTask(int thread_id, int id, int id_type = Task::ID_NONE, int additional_value = -1)
        {
          if(!tracing_enabled) return -1;
          if(!trace_threads && !trace_thread_counter) return -1;
	  if(unlikely(tasks[thread_id].size() == max_num_events_per_thread))
            StopTracing();
          int task_num = tasks[thread_id].size();
          tasks[thread_id].push_back( Task{thread_id, id, id_type, additional_value, GetTimeCounter(), true} );
          return task_num;
        }

      void StopTask(int thread_id, int id, int id_type = Task::ID_NONE)
        {
          if(!trace_threads && !trace_thread_counter) return;
          tasks[thread_id].push_back( Task{thread_id, id, id_type, 0, GetTimeCounter(), false} );
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

      void SendData(); // MPI parallel data reduction

    };
} // namespace ngcore

#endif // NETGEN_CORE_PAJE_TRACE_HPP

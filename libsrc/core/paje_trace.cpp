#include <algorithm>
#include <atomic>
#include <iostream>
#include <map>
#include <set>
#include <thread>

#include "archive.hpp"           // for Demangle
#include "paje_trace.hpp"
#include "profiler.hpp"
#include "mpi_wrapper.hpp"

extern const char *header;

constexpr int MPI_PAJE_WRITER = 1;

namespace ngcore
{
  static std::string GetTimerName( int id )
  {
#ifndef PARALLEL
    return NgProfiler::GetName(id);
#else // PARALLEL
    if(id<NgProfiler::SIZE)
      return NgProfiler::GetName(id);

    NgMPI_Comm comm(MPI_COMM_WORLD);
    return NgProfiler::GetName(id-NgProfiler::SIZE*comm.Rank());
#endif // PARALLEL
  }

  std::vector<PajeTrace::MemoryEvent> PajeTrace::memory_events;

  // Produce no traces by default
  size_t PajeTrace::max_tracefile_size = 0;

  // If true, produce variable counting active threads
  // increases trace by a factor of two
  bool PajeTrace::trace_thread_counter = false;
  bool PajeTrace::trace_threads = true;
  bool PajeTrace::mem_tracing_enabled = true;

  PajeTrace :: PajeTrace(int anthreads, std::string aname)
  {

    nthreads = anthreads;
    tracefile_name = std::move(aname);

    int bytes_per_event=33;
    max_num_events_per_thread = std::min( static_cast<size_t>(std::numeric_limits<int>::max()), max_tracefile_size/bytes_per_event/(nthreads+1+trace_thread_counter*nthreads)*10/7);
    if(max_num_events_per_thread>0)
    {
      logger->info( "Tracefile size = {}MB", max_tracefile_size/1024/1024);
      logger->info( "Tracing {} events per thread", max_num_events_per_thread);
    }

    tasks.resize(nthreads);
    int reserve_size = std::min(1000000U, max_num_events_per_thread);
    for(auto & t : tasks)
        t.reserve(reserve_size);

    links.resize(nthreads);
    for(auto & l : links)
      l.reserve(reserve_size);

    jobs.reserve(reserve_size);
    timer_events.reserve(reserve_size);
    memory_events.reserve(1024*1024);

    // sync start time when running in parallel
#ifdef PARALLEL
    NgMPI_Comm comm(MPI_COMM_WORLD);
    for(auto i : Range(5))
        comm.Barrier();
#endif // PARALLEL

    start_time = GetTimeCounter();
    tracing_enabled = true;
    mem_tracing_enabled = true;
    n_memory_events_at_start = memory_events.size();
  }

  PajeTrace :: ~PajeTrace()
  {
    for(auto & ltask : tasks)
        for(auto & task : ltask)
          {
            task.start_time -= start_time;
            task.stop_time -= start_time;
          }
    for(auto & job : jobs)
      {
        job.start_time -= start_time;
        job.stop_time -= start_time;
      }
    for(auto & event : timer_events)
        event.time -= start_time;

    for(auto & llink : links)
        for(auto & link : llink)
            link.time -= start_time;

    for(auto i : IntRange(n_memory_events_at_start, memory_events.size()))
      memory_events[i].time -= start_time;

    NgMPI_Comm comm(MPI_COMM_WORLD);

    if(comm.Size()==1)
    {
      Write(tracefile_name);
    }
    else
    {
      // make sure the timer id is unique across all ranks
      for(auto & event : timer_events)
        event.timer_id += NgProfiler::SIZE*comm.Rank();

      if(comm.Rank() == MPI_PAJE_WRITER)
        Write(tracefile_name);
      else
        SendData();
    }
  }


  void PajeTrace::StopTracing()
    {
      if(tracing_enabled && max_num_events_per_thread>0)
        {
          logger->warn("Maximum number of traces reached, tracing is stopped now.");
        }
      tracing_enabled = false;
    }

  class PajeFile
    {
    public:
      static void Hue2RGB ( double x, double &r, double &g, double &b )
        {
          double d = 1.0/6.0;
          if(x<d)
            r=1, g=6*x,b=0;
          else if (x<2*d)
            r=1.0-6*(x-d),g=1,b=0;
          else if (x<3*d)
            r=0, g=1,b=6*(x-2*d);
          else if (x<4*d)
            r=0, g=1-6*(x-3*d),b=1;
          else if (x<5*d)
            r=6*(x-4*d), g=0,b=1;
          else
            r=1, g=0,b=1-5*(x-d);
        };

      int alias_counter;

      FILE * ctrace_stream;
      std::shared_ptr<Logger> logger = GetLogger("PajeTrace");


      double ConvertTime(TTimePoint t) {
          // return time in milliseconds as double
        // return std::chrono::duration<double>(t-start_time).count()*1000.0;
        // return std::chrono::duration<double>(t-start_time).count() / 2.7e3;
        return 1000.0*static_cast<double>(t) * seconds_per_tick;
      }

      enum PType
        {
          SET_VARIABLE=1,
          ADD_VARIABLE,
          SUB_VARIABLE,
          PUSH_STATE,
          POP_STATE,
          START_LINK,
          STOP_LINK
        };

      struct PajeEvent
        {
          PajeEvent( int aevent_type, double atime, int atype, int acontainer, double avar_value )
            : time(atime), var_value(avar_value), event_type(aevent_type), type(atype), container(acontainer)
            { }

          PajeEvent( int aevent_type, double atime, int atype, int acontainer, int avalue = 0, int aid = 0, bool avalue_is_alias = true )
            : time(atime), event_type(aevent_type), type(atype), container(acontainer), value(avalue), id(aid), value_is_alias(avalue_is_alias)
            { }

          PajeEvent( int aevent_type, double atime, int atype, int acontainer, int avalue, int astart_container, int akey )
            : time(atime), event_type(aevent_type), type(atype), container(acontainer), value(avalue), start_container(astart_container), id(akey)
            { }

          double time;
          double var_value = 0.0;
          int event_type;
          int type;
          int container;
          int value = 0;
          int start_container = 0;
          int id = 0;
          bool value_is_alias = true;

          bool operator < (const PajeEvent & other) const {
              // Same start and stop times can occur for very small tasks -> take "starting" events first (eg. PajePushState before PajePopState)
              if(time == other.time)
                return event_type < other.event_type;
              return (time < other.time);
          }

          int write(FILE *stream)
            {
              const int &key = id;
              const int &end_container = start_container;
              switch(event_type)
                {
                case PajeSetVariable:
                  return fprintf( stream, "%d\t%.15g\ta%d\ta%d\t%.15g\n", PajeSetVariable, time, type, container, var_value ); // NOLINT
                case PajeAddVariable:
                  return fprintf( stream, "%d\t%.15g\ta%d\ta%d\t%.15g\n", PajeAddVariable, time, type, container, var_value ); // NOLINT
                case PajeSubVariable:
                  return fprintf( stream, "%d\t%.15g\ta%d\ta%d\t%.15g\n", PajeSubVariable, time, type, container, var_value ); // NOLINT
                case PajePushState:
                  if(value_is_alias)
                    return fprintf( stream, "%d\t%.15g\ta%d\ta%d\ta%d\t%d\n", PajePushState, time, type, container, value, id); // NOLINT
                  else
                    return fprintf( stream, "%d\t%.15g\ta%d\ta%d\t%d\t%d\n", PajePushState, time, type, container, value, id); // NOLINT
                case PajePopState:
                  return fprintf( stream, "%d\t%.15g\ta%d\ta%d\n", PajePopState, time, type, container ); // NOLINT
                case PajeStartLink:
                  return fprintf( stream, "%d\t%.15g\ta%d\ta%d\t%d\ta%d\t%d\n", PajeStartLink, time, type, container, value, start_container, key ); // NOLINT
                case PajeEndLink:
                  return fprintf( stream, "%d\t%.15g\ta%d\ta%d\t%d\ta%d\t%d\n", PajeEndLink, time, type, container, value, end_container, key ); // NOLINT
                }
              return 0;
            }
        };

      std::vector<PajeEvent> events;

    public:
      PajeFile() = delete;
      PajeFile(const PajeFile &) = delete;
      PajeFile(PajeFile &&) = delete;
      void operator=(const PajeFile &) = delete;
      void operator=(PajeFile &&) = delete;

      PajeFile( const std::string & filename)
        {
          std::string fname = filename + ".trace";
          ctrace_stream = fopen (fname.c_str(),"w"); // NOLINT
          fprintf(ctrace_stream, "%s", header ); // NOLINT
          alias_counter = 0;
        }

      ~PajeFile()
        {
          fclose (ctrace_stream); // NOLINT
        }

      int DefineContainerType ( int parent_type, const std::string & name )
        {
          int alias = ++alias_counter;
          if(parent_type!=0)
            fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\n", PajeDefineContainerType, alias, parent_type, name.c_str() ); // NOLINT
          else
            fprintf( ctrace_stream, "%d\ta%d\t%d\t\"%s\"\n", PajeDefineContainerType, alias, parent_type, name.c_str() ); // NOLINT
          return alias;
        }

      int DefineVariableType ( int container_type, const std::string & name )
        {
          int alias = ++alias_counter;
          fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\t\"1.0 1.0 1.0\"\n", PajeDefineVariableType, alias, container_type, name.c_str() ); // NOLINT
          return alias;
        }

      int DefineStateType ( int type, const std::string & name )
        {
          int alias = ++alias_counter;
          fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\n", PajeDefineStateType, alias, type, name.c_str() ); // NOLINT
          return alias;
        }

      //       int DefineEventType ()
      //         {
      //           Write("event not implemented");
      //         }

      int DefineLinkType (int parent_container_type, int start_container_type, int stop_container_type, const std::string & name)
        {
          int alias = ++alias_counter;
          fprintf( ctrace_stream, "%d\ta%d\ta%d\ta%d\ta%d\t\"%s\"\n", PajeDefineLinkType, alias, parent_container_type, start_container_type, stop_container_type, name.c_str() ); // NOLINT
          return alias;
        }

      int DefineEntityValue (int type, const std::string & name, double hue = -1)
        {
          if(hue==-1)
            {
              std::hash<std::string> shash;
              size_t h = shash(name);
              h ^= h>>32U;
              h = static_cast<uint32_t>(h);
              hue = h*1.0/std::numeric_limits<uint32_t>::max();
            }

          int alias = ++alias_counter;
          double r;
          double g;
          double b;
          Hue2RGB( hue, r, g, b );
          fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\t\"%.15g %.15g %.15g\"\n", PajeDefineEntityValue, alias, type, name.c_str(), r,g,b ); // NOLINT
          return alias;
        }

      int CreateContainer ( int type, int parent, const std::string & name )
        {
          int alias = ++alias_counter;
          if(parent!=0)
            fprintf( ctrace_stream, "%d\t0\ta%d\ta%d\ta%d\t\"%s\"\n", PajeCreateContainer, alias, type, parent, name.c_str() ); // NOLINT
          else
            fprintf( ctrace_stream, "%d\t0\ta%d\ta%d\t%d\t\"%s\"\n", PajeCreateContainer, alias, type, parent, name.c_str() ); // NOLINT
          return alias;
        }
      void DestroyContainer ()
        {}

      void SetVariable (TTimePoint time, int type, int container, double value )
        {
          events.emplace_back( PajeEvent( PajeSetVariable, ConvertTime(time), type, container, value ) );
        }

      void AddVariable (TTimePoint time, int type, int container, double value )
        {
          events.emplace_back( PajeEvent( PajeAddVariable, ConvertTime(time), type, container, value ) );
        }

      void SubVariable (TTimePoint time, int type, int container, double value )
        {
          events.emplace_back( PajeEvent( PajeSubVariable, ConvertTime(time), type, container, value ) );
        }

      void SetState ()
        {}

      void PushState ( TTimePoint time, int type, int container, int value, int id = 0, bool value_is_alias = true )
        {
          events.emplace_back( PajeEvent( PajePushState, ConvertTime(time), type, container, value, id, value_is_alias) );
        }

      void PopState ( TTimePoint time, int type, int container )
        {
          events.emplace_back( PajeEvent( PajePopState, ConvertTime(time), type, container ) );
        }

      void ResetState ()
        {}

      void StartLink ( TTimePoint time, int type, int container, int value, int start_container, int key )
        {
          events.emplace_back( PajeEvent( PajeStartLink, ConvertTime(time), type, container, value, start_container, key ) );
        }

      void EndLink ( TTimePoint time, int type, int container, int value, int end_container, int key )
        {
          events.emplace_back( PajeEvent(  PajeEndLink, ConvertTime(time), type, container, value, end_container, key ) );
        }

      void NewEvent ()
        {}

      void WriteEvents()
        {
          logger->info("Sorting traces...");
          std::sort (events.begin(), events.end());

          logger->info("Writing traces... ");
          for (auto & event : events)
          {
              event.write( ctrace_stream );
//               fprintf( ctrace_stream, "%s", buf ); // NOLINT
          }
          logger->info("Done");
        }

    private:
      enum
        {
          PajeDefineContainerType = 0,
          PajeDefineVariableType = 1,
          PajeDefineStateType = 2,
          PajeDefineEventType = 3,
          PajeDefineLinkType = 4,
          PajeDefineEntityValue = 5,
          PajeCreateContainer = 6,
          PajeDestroyContainer = 7,
          PajeSetVariable = 8,
          PajeAddVariable = 9,
          PajeSubVariable = 10,
          PajeSetState = 11,
          PajePushState = 12,
          PajePopState = 13,
          PajeResetState = 14,
          PajeStartLink = 15,
          PajeEndLink = 16,
          PajeNewEvent = 17
        };

    };

  NGCORE_API PajeTrace *trace;

  void PajeTrace::Write( const std::string & filename )
    {
      auto n_events = jobs.size() + timer_events.size();
      for(auto & vtasks : tasks)
        n_events += vtasks.size();

      logger->info("{} events traced",  n_events);

      if(n_events==0)
        {
          logger->info("No data traced, skip writing trace file");
          return;
        }

      if(!tracing_enabled)
        {
            logger->warn("Tracing stopped during computation due to tracefile size limit of {} megabytes.", max_tracefile_size/1024/1024);
        }

      PajeFile paje(filename);

      const int container_type_task_manager = paje.DefineContainerType( 0, "Task Manager" );
      const int container_type_node = paje.DefineContainerType( container_type_task_manager, "Node");
      const int container_type_thread = paje.DefineContainerType( container_type_task_manager, "Thread");
      const int container_type_timer = container_type_thread; //paje.DefineContainerType( container_type_task_manager, "Timers");
      const int container_type_jobs = paje.DefineContainerType( container_type_task_manager, "Jobs");
      const int container_type_memory = paje.DefineContainerType( container_type_task_manager, "Memory usage");

      const int state_type_job = paje.DefineStateType( container_type_jobs, "Job" );
      const int state_type_task = paje.DefineStateType( container_type_thread, "Task" );
      const int state_type_timer = paje.DefineStateType( container_type_timer, "Timer state" );

      int variable_type_active_threads = 0;
      if(trace_thread_counter)
          variable_type_active_threads = paje.DefineVariableType( container_type_jobs, "Active threads" );

      const int container_task_manager = paje.CreateContainer( container_type_task_manager, 0, "The task manager" );
      const int container_jobs = paje.CreateContainer( container_type_jobs, container_task_manager, "Jobs" );

      int variable_type_memory = 0;
      const int container_memory = paje.CreateContainer( container_type_memory, container_task_manager, "Memory" );
      if(mem_tracing_enabled)
      {
        variable_type_memory = paje.DefineVariableType( container_type_task_manager, "Memory [MB]" );
      }


      int num_nodes = 1; //task_manager ? task_manager->GetNumNodes() : 1;
      std::vector <int> thread_aliases;
      std::vector<int> container_nodes;

#ifdef PARALLEL
      // Hostnames
      NgMPI_Comm comm(MPI_COMM_WORLD);
      auto rank = comm.Rank();
      auto nranks = comm.Size();
      if(nranks>1)
      {
        nthreads = nranks;
        thread_aliases.reserve(nthreads);

        std::array<char, MPI_MAX_PROCESSOR_NAME+1> ahostname;
        int len;
        MPI_Get_processor_name(ahostname.data(), &len);
        std::string hostname = ahostname.data();

        std::map<std::string, int> host_map;

        std::string name;
        for(auto i : IntRange(0, nranks))
        {
          if(i!=MPI_PAJE_WRITER)
            comm.Recv(name, i, 0);
          else
            name = hostname;
          if(host_map.count(name)==0)
          {
            host_map[name] = container_nodes.size();
            container_nodes.emplace_back( paje.CreateContainer( container_type_node, container_task_manager, name) );
          }
          thread_aliases.emplace_back( paje.CreateContainer( container_type_thread, container_nodes[host_map[name]], "Rank " + ToString(i) ) );
        }
      }
      else
#endif // PARALLEL
      {
        container_nodes.reserve(num_nodes);
        for(int i=0; i<num_nodes; i++)
          container_nodes.emplace_back( paje.CreateContainer( container_type_node, container_task_manager, "Node " + ToString(i)) );

        thread_aliases.reserve(nthreads);
        if(trace_threads)
          for (int i=0; i<nthreads; i++)
          {
            auto name = "Thread " + ToString(i);
            thread_aliases.emplace_back( paje.CreateContainer( container_type_thread, container_nodes[i*num_nodes/nthreads], name ) );
          }
      }

      std::map<const std::type_info *, int> job_map;
      std::map<const std::type_info *, int> job_task_map;

      for(Job & j : jobs)
        if(job_map.find(j.type) == job_map.end())
          {
            std::string name = Demangle(j.type->name());
            job_map[j.type] = paje.DefineEntityValue( state_type_job, name, -1 );
            job_task_map[j.type] = paje.DefineEntityValue( state_type_task, name, -1 );
          }

      for(Job & j : jobs)
        {
          paje.PushState( j.start_time, state_type_job, container_jobs, job_map[j.type] );
          paje.PopState( j.stop_time, state_type_job, container_jobs );
        }

      size_t memory_at_start = 0;

      for(const auto & i : IntRange(0, n_memory_events_at_start))
      {
        if(memory_events[i].is_alloc)
            memory_at_start += memory_events[i].size;
        else
            memory_at_start -= memory_events[i].size;
      }

      paje.SetVariable( 0, variable_type_memory, container_memory, 1.0*memory_at_start/(1024*1024));

      for(const auto & i : IntRange(n_memory_events_at_start, memory_events.size()))
      {
        auto & m = memory_events[i];
        if(m.size==0)
            continue;
        double size = 1.0*m.size/(1024*1024);
        if(m.is_alloc)
          paje.AddVariable( m.time, variable_type_memory, container_memory, size);
        else
          paje.SubVariable( m.time, variable_type_memory, container_memory, size);
      }

      std::set<int> timer_ids;
      std::map<int,int> timer_aliases;
      std::map<int,std::string> timer_names;

      for(auto & event : timer_events)
          timer_ids.insert(event.timer_id);


      // Timer names
      for(auto & vtasks : tasks)
          for (Task & t : vtasks)
              if(t.id_type==Task::ID_TIMER)
                  timer_ids.insert(t.id);

      for(auto id : timer_ids)
          timer_names[id] = GetTimerName(id);

#ifdef PARALLEL
      if(nranks>1)
      {
        for(auto src : IntRange(0, nranks))
        {
          if(src==MPI_PAJE_WRITER)
            continue;

          size_t n_timers;
          comm.Recv (n_timers, src, 0);

          int id;
          std::string name;
          for(auto i : IntRange(n_timers))
          {
            comm.Recv (id, src, 0);
            comm.Recv (name, src, 0);
            timer_ids.insert(id);
            timer_names[id] = name;
          }
        }
      }
#endif // PARALLEL

      for(auto id : timer_ids)
          timer_aliases[id] = paje.DefineEntityValue( state_type_timer, timer_names[id], -1 );

      int timerdepth = 0;
      int maxdepth = 0;
      for(auto & event : timer_events)
        {
          if(event.is_start)
            {
              timerdepth++;
              maxdepth = timerdepth>maxdepth ? timerdepth : maxdepth;
            }
          else
            timerdepth--;
        }

      std::vector<int> timer_container_aliases;
      timer_container_aliases.resize(maxdepth);
      for(int i=0; i<maxdepth; i++)
        {
          auto name = "Timer level " + ToString(i);
          timer_container_aliases[i] =  paje.CreateContainer( container_type_timer, container_task_manager, name );
        }

      timerdepth = 0;
      for(auto & event : timer_events)
        {
          if(event.is_start)
            paje.PushState( event.time, state_type_timer, timer_container_aliases[timerdepth++], timer_aliases[event.timer_id] );
          else
            paje.PopState( event.time, state_type_timer, timer_container_aliases[--timerdepth] );
        }

      for(auto & vtasks : tasks)
        {
          for (Task & t : vtasks) {
              int value_id = t.id;

              switch(t.id_type)
                {
                case Task::ID_JOB:
                  value_id = job_task_map[jobs[t.id-1].type];
                  if(trace_thread_counter)
                    {
                      paje.AddVariable( t.start_time, variable_type_active_threads, container_jobs, 1.0 );
                      paje.SubVariable( t.stop_time, variable_type_active_threads, container_jobs, 1.0 );
                    }
                  if(trace_threads)
                    {
                      paje.PushState( t.start_time, state_type_task, thread_aliases[t.thread_id], value_id, t.additional_value, true );
                      paje.PopState( t.stop_time, state_type_task, thread_aliases[t.thread_id] );
                    }
                  break;
                case Task::ID_TIMER:
                  value_id = timer_aliases[t.id];
                  paje.PushState( t.start_time, state_type_timer, thread_aliases[t.thread_id], value_id, t.additional_value, true );
                  paje.PopState( t.stop_time, state_type_timer, thread_aliases[t.thread_id] );
                  break;
                default:
                  paje.PushState( t.start_time, state_type_task, thread_aliases[t.thread_id], value_id, t.additional_value, false );
                  paje.PopState( t.stop_time, state_type_task, thread_aliases[t.thread_id] );
                  break;
                }
          }
        }

#ifdef PARALLEL
      if(nranks>1)
      {
        for(auto & event : timer_events)
        {
          if(event.is_start)
            paje.PushState( event.time, state_type_timer, thread_aliases[MPI_PAJE_WRITER], timer_aliases[event.timer_id] );
          else
            paje.PopState( event.time, state_type_timer, thread_aliases[MPI_PAJE_WRITER] );
        }

        // Timer events
        Array<int> timer_id;
        Array<TTimePoint> time;
        Array<bool> is_start;
        Array<int> thread_id;

        for(auto src : IntRange(0, nranks))
        {
          if(src==MPI_PAJE_WRITER)
            continue;

          comm.Recv (timer_id, src, 0);
          comm.Recv (time, src, 0);
          comm.Recv (is_start, src, 0);
          comm.Recv (thread_id, src, 0);

          for(auto i : Range(timer_id.Size()))
          {
            TimerEvent event;
            event.timer_id = timer_id[i];
            event.time = time[i];
            event.is_start = is_start[i];
            event.thread_id = thread_id[i];

            if(event.is_start)
              paje.PushState( event.time, state_type_timer, thread_aliases[src], timer_aliases[event.timer_id] );
            else
              paje.PopState( event.time, state_type_timer, thread_aliases[src] );
          }
        }
      }
#endif // PARALLEL

      // Merge link event
      int nlinks = 0;
      for( auto & l : links)
        nlinks += l.size();

      std::vector<ThreadLink> links_merged;
      links_merged.reserve(nlinks);
      std::vector<unsigned int> pos(nthreads);

      int nlinks_merged = 0;
      while(nlinks_merged < nlinks)
        {
          int minpos = -1;
          TTimePoint mintime = -1;
          for (int t = 0; t<nthreads; t++)
            {
              if(pos[t] < links[t].size() && (minpos==-1 || links[t][pos[t]].time < mintime))
                {
                  minpos = t;
                  mintime = links[t][pos[t]].time;
                }
            }
          links_merged.push_back( links[minpos][pos[minpos]] );
          pos[minpos]++;
          nlinks_merged++;
        }

      std::vector<ThreadLink> started_links;

      int link_type = paje.DefineLinkType(container_type_node, container_type_thread, container_type_thread, "links");

      // match links
      for ( auto & l : links_merged )
        {
          if(l.is_start)
            {
              started_links.push_back(l);
            }
          else
            {
              unsigned int i = 0;
              while(i<started_links.size())
                {
                  while(i<started_links.size() && started_links[i].key == l.key)
                    {
                      ThreadLink & sl = started_links[i];
                      // Avoid links on same thread
                      if(sl.thread_id != l.thread_id)
                        {
                          paje.StartLink( sl.time, link_type, container_nodes[sl.thread_id*num_nodes/nthreads], l.key, thread_aliases[sl.thread_id], l.key);
                          paje.EndLink(    l.time, link_type, container_nodes[l.thread_id*num_nodes/nthreads], l.key, thread_aliases[l.thread_id], l.key);
                        }
                      started_links.erase(started_links.begin()+i);
                    }
                  i++;
                }
            }
        }
      WriteTimingChart();
#ifdef NETGEN_TRACE_MEMORY
      WriteMemoryChart("");
#endif // NETGEN_TRACE_MEMORY
      paje.WriteEvents();
    }

  void PajeTrace::SendData( )
    {
#ifdef PARALLEL
      // Hostname
      NgMPI_Comm comm(MPI_COMM_WORLD);
      auto rank = comm.Rank();
      auto nranks = comm.Size();

      std::string hostname;
        {
          std::array<char, MPI_MAX_PROCESSOR_NAME+1> ahostname;
          int len;
          MPI_Get_processor_name(ahostname.data(), &len);
          hostname = ahostname.data();
        }

      comm.Send(hostname, MPI_PAJE_WRITER, 0);

      // Timer names
      std::set<int> timer_ids;
      std::map<int,std::string> timer_names;

      for(auto & event : timer_events)
          timer_ids.insert(event.timer_id);

      for(auto id : timer_ids)
          timer_names[id] = GetTimerName(id);
      size_t size = timer_ids.size();
      comm.Send(size, MPI_PAJE_WRITER, 0);
      for(auto id : timer_ids)
        {
          comm.Send(id, MPI_PAJE_WRITER, 0);
          comm.Send(timer_names[id], MPI_PAJE_WRITER, 0);
        }


      // Timer events
      Array<int> timer_id;
      Array<TTimePoint> time;
      Array<bool> is_start;
      Array<int> thread_id;

      for(auto & event : timer_events)
        {
          timer_id.Append(event.timer_id);
          time.Append(event.time);
          is_start.Append(event.is_start);
          thread_id.Append(event.thread_id);
        }

      comm.Send (timer_id, MPI_PAJE_WRITER, 0);
      comm.Send (time, MPI_PAJE_WRITER, 0);
      comm.Send (is_start, MPI_PAJE_WRITER, 0);
      comm.Send (thread_id, MPI_PAJE_WRITER, 0);
#endif // PARALLEL
    }

  ///////////////////////////////////////////////////////////////////
  // Write HTML file drawing a sunburst chart with cumulated timings
  struct TreeNode
  {
      int id = 0;
      std::map<int, TreeNode> children;
      double chart_size = 0.0; // time without children (the chart lib accumulates children sizes again)
      double size = 0.0;
      double min_size = 1e99;
      double max_size = 0.0;
      std::string name;

      size_t calls = 0;
      TTimePoint start_time = 0;
  };

  void PrintNode (const TreeNode &n, std::ofstream & f)
  {
      f << "{ name: \"" + n.name + "\"";
      f << ", calls: " << n.calls;
      f << ", size: " << n.chart_size;
      f << ", value: " << n.size;
      f << ", min: " << n.min_size;
      f << ", max: " << n.max_size;
      if(n.calls)
        f << ", avg: " << n.size/n.calls;
      int size = n.children.size();
      if(size>0)
      {
          int i = 0;
          f << ", children: [";
          for(auto & c : n.children)
          {
              PrintNode(c.second, f);
              if(++i<size)
                  f << " , ";
          }
          f << ']';
      }
      f << '}';
  }

  void WriteSunburstHTML( TreeNode & root, std::string filename, bool time_or_memory )
  {
    std::ofstream f(filename+".html");
    f.precision(4);
    f << R"CODE_(
<head>
  <script src="https://d3js.org/d3.v5.min.js"></script>
  <script src="https://unpkg.com/sunburst-chart"></script>

  <style>body { margin: 0 }</style>
)CODE_";
    if(!time_or_memory)
      f << "<title>Maximum Memory Consumption</title>\n";
    f << R"CODE_(
</head>
<body>
  <div id="chart"></div>

  <script>
    const data = 
)CODE_";
      PrintNode(root, f);
      f << ";\n\n";
      if(time_or_memory)
        f << "const chart_type = 'time';\n";
      else
        f << "const chart_type = 'memory';\n";
      f << R"CODE_(
    const color = d3.scaleOrdinal(d3.schemePaired);

    let getTime = (t) =>
    {
       if(t>=1000)  return (t/1000).toPrecision(4) + '  s';
       if(t>=0.1)   return t.toPrecision(4) + ' ms';
       if(t>=1e-4)  return (t*1e3).toPrecision(4) + ' us';

       return (t/1e6).toPrecision(4) + ' ns';
    };

    const KB_ = 1024;
    const MB_ = KB_*1024;
    const GB_ = MB_*1024;
    let getMemory = (m) =>
    {
       if(m>=GB_)  return (m/GB_).toPrecision(4) + ' GB';
       if(m>=MB_)  return (m/MB_).toPrecision(4) + ' MB';
       if(m>=KB_)  return (m/KB_).toPrecision(4) + ' KB';
       return m.toPrecision(4) + ' B';
    };

    Sunburst()
      .data(data)
      .size('size')
      .color(d => color(d.name))
      .tooltipTitle((d, node) => { return node.parent ? node.parent.data.name + " &rarr; " + d.name : d.name; })
      .tooltipContent((d, node) => {
        if(chart_type=="memory")
        {
          return `Total Memory: <i>${getMemory(d.value)}</i> <br>`
               + `Memory: <i>${getMemory(d.size)}</i>`
        }
        else
        {
          return `Time: <i>${getTime(d.value)}</i> <br>`
               + `calls: <i>${d.calls}</i> <br>`
               + `min: <i>${getTime(d.min)}</i> <br>`
               + `max: <i>${getTime(d.max)}</i> <br>`
               + `avg: <i>${getTime(d.avg)}</i>`
        }
      })
      (document.getElementById('chart'));

      // Line breaks in tooltip
      var all = document.getElementsByClassName('sunbirst-tooltip');
      for (var i = 0; i < all.length; i++) {
          all[i].white_space = "";
      }
  </script>
</body>
)CODE_" << std::endl;


  }

#ifdef NETGEN_TRACE_MEMORY
  void PajeTrace::WriteMemoryChart( std::string fname )
  {
    if(fname=="")
      fname = tracefile_name + "_memory";
    size_t mem_allocated = 0;
    size_t max_mem_allocated = 0;
    size_t imax_mem_allocated = 0;

    const auto & names = MemoryTracer::GetNames();
    const auto & parents = MemoryTracer::GetParents();
    size_t N = names.size();

    Array<size_t> mem_allocated_id;
    mem_allocated_id.SetSize(N);
    mem_allocated_id = 0;

    // Find point with maximum memory allocation, check for missing allocs/frees
    for(auto i : IntRange(memory_events.size()))
    {
      const auto & ev = memory_events[i];

      if(ev.is_alloc)
      {
        mem_allocated += ev.size;
        mem_allocated_id[ev.id] += ev.size;
        if(mem_allocated > max_mem_allocated && i>=n_memory_events_at_start)
        {
          imax_mem_allocated = i;
          max_mem_allocated = mem_allocated;
        }
      }
      else
      {
        if(ev.size > mem_allocated)
          {
            std::cerr << "Error in memory tracer: have total allocated memory < 0" << std::endl;
            mem_allocated = 0;
          }
        else
          mem_allocated -= ev.size;
        if(ev.size > mem_allocated_id[ev.id])
          {
            std::cerr << "Error in memory tracer: have allocated memory < 0 in tracer " << names[ev.id] << std::endl;
            mem_allocated_id[ev.id] = 0;
          }
        else
          mem_allocated_id[ev.id] -= ev.size;
      }
    }

    // reconstruct again the memory consumption after event imax_mem_allocated
    mem_allocated_id = 0;
    for(auto i : IntRange(imax_mem_allocated+1))
    {
      const auto & ev = memory_events[i];

      if(ev.is_alloc)
        mem_allocated_id[ev.id] += ev.size;
      else
        {
          if(ev.size > mem_allocated_id[ev.id])
            mem_allocated_id[ev.id] = 0;
          else
            mem_allocated_id[ev.id] -= ev.size;
        }
    }

    TreeNode root;
    root.name="all";

    Array<TreeNode*> nodes;
    nodes.SetSize(N);
    nodes = nullptr;
    nodes[0] = &root;
    Array<Array<int>> children(N);

    Array<size_t> sorting; // topological sorting (parents before children)
    sorting.SetAllocSize(N);

    for(auto i : IntRange(1, N))
        children[parents[i]].Append(i);

    ArrayMem<size_t, 100> stack;
    sorting.Append(0);
    stack.Append(0);

    while(stack.Size())
    {
      auto current = stack.Last();
      stack.DeleteLast();

      for(const auto child : children[current])
      {
        sorting.Append(child);
        if(children[child].Size())
          stack.Append(child);
      }
    }

    for(auto i : sorting)
    {
      if(i==0)
          continue;

      TreeNode * parent = nodes[parents[i]];

      auto & node = parent->children[i];
      nodes[i] = &node;
      node.id = i;
      node.chart_size = mem_allocated_id[i];
      node.size = mem_allocated_id[i];
      node.name = names[i];
    }

    for(auto i_ : Range(sorting))
    {
      // reverse topological order to accumulate total memory usage of all children
      auto i = sorting[sorting.Size()-1-i_];
      if(i==0)
          continue;
      nodes[parents[i]]->size += nodes[i]->size;
    }

    WriteSunburstHTML( root, fname, false );

  }
#endif // NETGEN_TRACE_MEMORY

  void PajeTrace::WriteTimingChart( )
  {
      std::vector<TimerEvent> events;

      TreeNode root;
      root.name="all";
      TreeNode *current = &root;

      std::vector<TreeNode*> node_stack;

      node_stack.push_back(&root);

      TTimePoint stop_time = 0;

      for(auto & event : timer_events)
      {
          events.push_back(event);
          stop_time = std::max(event.time, stop_time);
      }

      std::map<std::string, int> jobs_map;
      std::vector<std::string> job_names;
      for(auto & job : jobs)
      {
          auto name = Demangle(job.type->name());
          int id = job_names.size();
          if(jobs_map.count(name)==0)
          {
              jobs_map[name] = id;
              job_names.push_back(name);
          }
          else
              id = jobs_map[name];

          events.push_back(TimerEvent{-1, job.start_time, true, id});
          events.push_back(TimerEvent{-1, job.stop_time, false, id});
          stop_time = std::max(job.stop_time, stop_time);
      }

      std::sort (events.begin(), events.end());

      root.size = 1000.0*static_cast<double>(stop_time) * seconds_per_tick;
      root.calls = 1;
      root.min_size = root.size;
      root.max_size = root.size;

      for(auto & event : events)
      {
          bool is_timer_event = event.timer_id != -1;
          int id = is_timer_event ? event.timer_id : event.thread_id;

          if(event.is_start)
          {
              bool need_init = !current->children.count(id);

              node_stack.push_back(current);
              current = &current->children[id];

              if(need_init)
              {
                  current->name = is_timer_event ? GetTimerName(id) : job_names[id];
                  current->size = 0.0;
                  current->id = id;
              }

              current->start_time = event.time;
          }
          else
          {
              if(node_stack.size()==0) {
                std::cout << "node stack empty!" << std::endl;
                break;
              }
              double size = 1000.0*static_cast<double>(event.time-current->start_time) * seconds_per_tick;
              current->size += size;
              current->chart_size += size;
              current->min_size = std::min(current->min_size, size);
              current->max_size = std::max(current->max_size, size);
              current->calls++;

              current = node_stack.back();
              current->chart_size -= size;
              node_stack.pop_back();
          }
      }

      root.chart_size = 0.0;

      ngcore::WriteSunburstHTML( root, tracefile_name, true );
  }

} // namespace ngcore

const char *header =
        "%EventDef PajeDefineContainerType 0 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineVariableType 1 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%       Color color \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineStateType 2 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineEventType 3 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%       Color color \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineLinkType 4 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       StartContainerType string \n"
        "%       EndContainerType string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineEntityValue 5 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%       Color color \n"
        "%EndEventDef \n"
        "%EventDef PajeCreateContainer 6 \n"
        "%       Time date \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDestroyContainer 7 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeSetVariable 8 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value double \n"
        "%EndEventDef\n"
        "%EventDef PajeAddVariable 9 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value double \n"
        "%EndEventDef\n"
        "%EventDef PajeSubVariable 10 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value double \n"
        "%EndEventDef\n"
        "%EventDef PajeSetState 11 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%EndEventDef\n"
        "%EventDef PajePushState 12 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%       Id string \n"
        "%EndEventDef\n"
        "%EventDef PajePopState 13 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%EndEventDef\n"
        "%EventDef PajeResetState 14 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%EndEventDef\n"
        "%EventDef PajeStartLink 15 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%       StartContainer string \n"
        "%       Key string \n"
        "%EndEventDef\n"
        "%EventDef PajeEndLink 16 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%       EndContainer string \n"
        "%       Key string \n"
        "%EndEventDef\n"
        "%EventDef PajeNewEvent 17 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%EndEventDef\n";

#include <algorithm>
#include <atomic>
#include <iostream>
#include <map>
#include <set>
#include <thread>

#include "archive.hpp"           // for Demangle
#include "paje_trace.hpp"
#include "profiler.hpp"

extern const char *header;

namespace ngcore
{
  // Produce no traces by default
  size_t PajeTrace::max_tracefile_size = 0;

  // If true, produce variable counting active threads
  // increases trace by a factor of two
  bool PajeTrace::trace_thread_counter = true;
  bool PajeTrace::trace_threads = true;

  PajeTrace :: PajeTrace(int anthreads, std::string aname)
  {
    start_time = GetTimeCounter();

    nthreads = anthreads;
    tracefile_name = std::move(aname);

    int bytes_per_event=33;
    max_num_events_per_thread = std::min( static_cast<size_t>(std::numeric_limits<int>::max()), max_tracefile_size/bytes_per_event/(2*nthreads+1)*10/7);
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

    tracing_enabled = true;
  }

  PajeTrace :: ~PajeTrace()
  {
    if(!tracefile_name.empty())
      Write(tracefile_name);
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
      TTimePoint start_time;
      std::shared_ptr<Logger> logger = GetLogger("PajeTrace");


      double ConvertTime(TTimePoint t) {
          // return time in milliseconds as double
        // return std::chrono::duration<double>(t-start_time).count()*1000.0;
        // return std::chrono::duration<double>(t-start_time).count() / 2.7e3;
        return 1000.0*static_cast<double>(t-start_time) / ticks_per_second;
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

      PajeFile( const std::string & filename, TTimePoint astart_time )
        {
          start_time = astart_time;
          ctrace_stream = fopen (filename.c_str(),"w"); // NOLINT
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

      PajeFile paje(filename, start_time);

      const int container_type_task_manager = paje.DefineContainerType( 0, "Task Manager" );
      const int container_type_node = paje.DefineContainerType( container_type_task_manager, "Node");
      const int container_type_thread = paje.DefineContainerType( container_type_task_manager, "Thread");
      const int container_type_timer = container_type_thread; //paje.DefineContainerType( container_type_task_manager, "Timers");
      const int container_type_jobs = paje.DefineContainerType( container_type_task_manager, "Jobs");

      const int state_type_job = paje.DefineStateType( container_type_jobs, "Job" );
      const int state_type_task = paje.DefineStateType( container_type_thread, "Task" );
      const int state_type_timer = paje.DefineStateType( container_type_timer, "Timer state" );

      const int variable_type_active_threads = paje.DefineVariableType( container_type_jobs, "Active threads" );

      const int container_task_manager = paje.CreateContainer( container_type_task_manager, 0, "The task manager" );
      const int container_jobs = paje.CreateContainer( container_type_jobs, container_task_manager, "Jobs" );
      paje.SetVariable( start_time, variable_type_active_threads, container_jobs, 0.0 );

      const int num_nodes = 1; //task_manager ? task_manager->GetNumNodes() : 1;

      std::vector<int> container_nodes;
      container_nodes.reserve(num_nodes);
      for(int i=0; i<num_nodes; i++)
          container_nodes.emplace_back( paje.CreateContainer( container_type_node, container_task_manager, "Node " + ToString(i)) );

      std::vector <int> thread_aliases;
      thread_aliases.reserve(nthreads);
      if(trace_threads)
        for (int i=0; i<nthreads; i++)
          {
            auto name = "Thread " + ToString(i);
            thread_aliases.emplace_back( paje.CreateContainer( container_type_thread, container_nodes[i*num_nodes/nthreads], name ) );
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

      std::set<int> timer_ids;
      std::map<int,int> timer_aliases;

      for(auto & event : timer_events)
        timer_ids.insert(event.timer_id);


      for(auto & vtasks : tasks)
        for (Task & t : vtasks)
          if(t.id_type==Task::ID_TIMER)
            timer_ids.insert(t.id);

      for(auto id : timer_ids)
        timer_aliases[id] = paje.DefineEntityValue( state_type_timer, NgProfiler::GetName(id), -1 );

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
      WriteSunburstHTML();
      paje.WriteEvents();
    }

  ///////////////////////////////////////////////////////////////////
  // Write HTML file drawing a sunburst chart with cumulated timings
  struct TreeNode
  {
      int id = 0;
      std::map<int, TreeNode> children;
      double time = 0.0;
      std::string name;
      TTimePoint start_time = 0;
  };

  void PrintNode (const TreeNode &n, int &level, std::ofstream & f);
  void PrintNode (const TreeNode &n, int &level, std::ofstream & f)
  {
      f << "{ name: \"" + n.name + "\", size: " + ToString(n.time);
      int size = n.children.size();
      if(size>0)
      {
          int i = 0;
          f << ", children: [";
          for(auto & c : n.children)
          {
              PrintNode(c.second, level, f);
              if(++i<size)
                  f << " , ";
          }
          f << ']';
      }
      f << '}';
  }

  void PajeTrace::WriteSunburstHTML( )
  {
      std::vector<TimerEvent> events;

      TreeNode root;
      root.time=0;
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

      root.time = 1000.0*static_cast<double>(stop_time-start_time)/ticks_per_second;

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
                  current->name = is_timer_event ? NgProfiler::GetName(id) : job_names[id];
                  current->time = 0.0;
                  current->id = id;
              }

              current->start_time = event.time;
          }
          else
          {
              double time = 1000.0*static_cast<double>(event.time-current->start_time)/ticks_per_second;
              current->time += time;
              current = node_stack.back();
              current->time -= time;
              node_stack.pop_back();
          }
      }

      int level = 0;
      std::ofstream f(tracefile_name+".html");
      f.precision(4);
      f << R"CODE_(
<head>
  <script src="https://d3js.org/d3.v5.min.js"></script>
  <script src="https://unpkg.com/sunburst-chart"></script>

  <style>body { margin: 0 }</style>
</head>
<body>
  <div id="chart"></div>

  <script>
    const data = 
)CODE_";
      PrintNode(root, level, f);
      f << R"CODE_( ;

    const color = d3.scaleOrdinal(d3.schemePaired);

    Sunburst()
      .data(data)
      .size('size')
      .color(d => color(d.name))
      .tooltipContent((d, node) => `Time: <i>${node.value.toPrecision(6)}ms</i>`)
      (document.getElementById('chart'));
  </script>
</body>
)CODE_" << std::endl;
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

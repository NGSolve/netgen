/********************************************************************/
/* File:   taskmanager.cpp                                          */
/* Author: M. Hochsterger, J. Schoeberl                             */
/* Date:   10. Mar. 2015                                            */
/********************************************************************/

#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>

#include "concurrentqueue.h"
#include "mpi_wrapper.hpp"
#include "paje_trace.hpp"
#include "profiler.hpp"
#include "taskmanager.hpp"

#ifdef USE_MKL
#include <mkl.h>
#endif



namespace ngcore
{
  using std::mutex;
  using std::lock_guard;
  using std::memory_order_release;
  using std::memory_order_relaxed;
  using std::make_tuple;

  struct TNestedTask
  {
    const function<void(TaskInfo&)> * func;
    int mynr;
    int total;
    int producing_thread;
    atomic<int> * endcnt;

    TNestedTask () { ; }
    TNestedTask (const function<void(TaskInfo&)> & _func,
                 int _mynr, int _total,
                 atomic<int> & _endcnt, int prod_tid)
      : func(&_func), mynr(_mynr), total(_total), producing_thread(prod_tid), endcnt(&_endcnt)
    {
      ;
    }
  };

  typedef moodycamel::ConcurrentQueue<TNestedTask> TQueue;
  typedef moodycamel::ProducerToken TPToken;
  typedef moodycamel::ConsumerToken TCToken;

  thread_local TaskManager * task_manager = nullptr;
  bool TaskManager :: use_paje_trace = false;
  int TaskManager :: max_threads = getenv("NGS_NUM_THREADS") ? atoi(getenv("NGS_NUM_THREADS")) : std::thread::hardware_concurrency();

  static std::thread::id main_thread_id = std::this_thread::get_id();
  
  thread_local int TaskManager :: thread_id = std::this_thread::get_id() == main_thread_id ? -2 : -1;
  thread_local WorkerData* TaskManager :: worker_data = nullptr;
  
  #ifdef WIN32
      TaskManager * GetTaskManager() { return task_manager; }
  #endif

  WorkerData :: ~WorkerData()
  {
      delete ex;
      delete static_cast<TPToken*>(produce_token);
      delete static_cast<TCToken*>(consume_token);
  }

  int EnterTaskManager (int nthreads)
  {
    if(nthreads == 0) return 0;
    if(nthreads == -1) nthreads = TaskManager::GetMaxThreads();
    if (task_manager)
      {
        // no task manager started
        return 0;
      }

    task_manager = new TaskManager(nthreads);

    GetLogger("TaskManager")->info("task-based parallelization (C++11 threads) using {} threads", task_manager->GetNumThreads());

#ifdef USE_NUMA
    numa_run_on_node (0);
#endif

#if !defined(WIN32) && !defined(EMSCRIPTEN)
    // master has maximal priority !
    int policy;
    struct sched_param param;
    pthread_getschedparam(pthread_self(), &policy, &param);
    param.sched_priority = sched_get_priority_max(policy);
    pthread_setschedparam(pthread_self(), policy, &param);
#endif // !defined(WIN32) && !defined(EMSCRIPTEN)

    
    task_manager->StartWorkers();

    ParallelFor (Range(100), [&] (int i) { ; });    // startup
    return task_manager->GetNumThreads();
  }


  void ExitTaskManager (int num_threads)
  {
    if(num_threads > 0)
      {
        task_manager->StopWorkers();
        delete task_manager;
        task_manager = nullptr;
      }
  }

  void RunWithTaskManager (function<void()> alg)
  {
    int num_threads = EnterTaskManager();
    alg();
    ExitTaskManager(num_threads);
  }




  void TaskManager :: SetNumThreads(int amax_threads)
    { 
      if(task_manager && task_manager->active_workers>0)
        {
          std::cerr << "Warning: can't change number of threads while TaskManager active!" << std::endl;
          return;
        }
      max_threads = amax_threads;
    }


  TaskManager :: TaskManager(int anthreads)
    {
      taskqueue_ptr = new TQueue;
      num_threads = anthreads;
      // if (MyMPI_GetNTasks() > 1) num_threads = 1;

#ifdef USE_NUMA
      numa_available();
      num_nodes = numa_max_node() + 1;
      if (num_nodes > num_threads) num_nodes = num_threads;

      for (int j = 0; j < num_nodes; j++)
        {
          void * mem = numa_alloc_onnode (sizeof(NodeData), j);
          nodedata[j] = new (mem) NodeData;
	  complete[j] = -1;
          workers_on_node[j] = 0;          
        }
#else
      num_nodes = 1;
      nodedata[0] = new NodeData;
      complete[0] = -1;
      workers_on_node[0] = 0;
#endif

      jobnr = 0;
      done = 0;
      sleep = false;
      sleep_usecs = 1000;
      active_workers = 0;

      static int cnt = 0;
      if (use_paje_trace)
          trace = new PajeTrace(num_threads, "ng" + ToString(cnt++));
    }


  TaskManager :: ~TaskManager ()
  {
    delete static_cast<TQueue*>(taskqueue_ptr);
    if (use_paje_trace)
      {
        delete trace;
        trace = nullptr;
      }
    num_threads = 1;
#ifdef USE_NUMA
      for (int j = 0; j < num_nodes; j++)
          numa_free (nodedata[j], sizeof(NodeData));
#else
      delete nodedata[0];
#endif
    task_manager = nullptr;
  }

#ifdef WIN32
  int TaskManager :: GetThreadId()
  {
    return thread_id;
  }
  WorkerData* TaskManager :: GetWorkerData()
  {
    return worker_data;
  }
#endif
  
  void TaskManager :: StartWorkers()
  {
    done = false;

    thread_id = 0;
    workers.resize(num_threads);
    auto &queue = *static_cast<TQueue*>(taskqueue_ptr);
    for (int i = 0; i < num_threads; i++)
    {
        workers[i] = WorkerData();
        workers[i].produce_token = new TPToken(queue);
        workers[i].consume_token = new TCToken(queue);
        if(i>0) workers[i].thread = std::thread([this,i]() { this->Loop(i); });
    }

    Loop(0); // sets all thread local variables
    
    while (active_workers < num_threads-1)
      ;
  }

  static size_t calibrate_init_tsc = GetTimeCounter();
  typedef std::chrono::system_clock TClock;
  static TClock::time_point calibrate_init_clock = TClock::now();
  
  void TaskManager :: StopWorkers()
  {
    static std::mutex timers_mutex;
    done = true;
    double delta_tsc = GetTimeCounter()-calibrate_init_tsc;
    double delta_sec = std::chrono::duration<double>(TClock::now()-calibrate_init_clock).count();
    double frequ = (delta_sec != 0) ? delta_tsc/delta_sec : 2.7e9;
    
    for(auto i : IntRange(1, workers.size()))
        workers[i].thread.join();

    std::lock_guard<std::mutex> guard(timers_mutex);
    for (size_t i = 0; i < num_threads; i++)
      for (size_t j = NgProfiler::SIZE; j-- > 0; )
        {
          if (!NgProfiler::timers[j].usedcounter) break;
          NgProfiler::timers[j].tottime += 1.0/frequ * workers[i].times[j];
          NgProfiler::timers[j].flops += workers[i].flops[j];
        }
    workers.clear();
    thread_id = -1;
  }

  /////////////////////// NEW: nested tasks using concurrent queue

  void TaskManager :: AddTask (const function<void(TaskInfo&)> & afunc,
                atomic<int> & endcnt)
                
  {
    auto &taskqueue = *static_cast<TQueue*>(taskqueue_ptr);
    int num = endcnt;
    auto tid = TaskManager::GetThreadId();
    for (int i = 0; i < num; i++)
      taskqueue.enqueue (*static_cast<TPToken*>(worker_data->produce_token), { afunc, i, num, endcnt, tid });
  }

  bool TaskManager :: ProcessTask()
  {
    // static Timer t("process task");
    TNestedTask task;
    auto &taskqueue = *static_cast<TQueue*>(taskqueue_ptr);
    auto tid = TaskManager::GetThreadId();
    if (taskqueue.try_dequeue(*static_cast<TCToken*>(worker_data->consume_token), task))
      {
        TaskInfo ti;
        ti.task_nr = task.mynr;
        ti.ntasks = task.total;
        ti.thread_nr = tid;
        ti.nthreads = TaskManager::GetNumThreads();
        /*
        {
          lock_guard<mutex> guard(m);
          cout << "process nested, nr = " << ti.task_nr << "/" << ti.ntasks << endl;
        }
        */
        // if(trace && task.producing_thread != ti.thread_nr)
        // trace->StartTask (ti.thread_nr, t, PajeTrace::Task::ID_TIMER, task.producing_thread);

        (*task.func)(ti);
        --*task.endcnt;

        // if(trace && task.producing_thread != ti.thread_nr)
        // trace->StopTask (ti.thread_nr, t);
        return true;
      }
    return false;
  }


  void TaskManager :: CreateJob (const function<void(TaskInfo&)> & afunc,
                                 int antasks)
  {
    auto *tm = GetTaskManager();
    if (!tm || tm->num_threads == 1 ) //  || func)
      {
        if (tm && tm->startup_function) (*tm->startup_function)();
        if (antasks == -1) antasks = 1;
        
        TaskInfo ti;
        ti.ntasks = antasks;
        ti.thread_nr = 0; ti.nthreads = 1;
        // ti.node_nr = 0; ti.nnodes = 1;
        for (ti.task_nr = 0; ti.task_nr < antasks; ti.task_nr++)
          afunc(ti);

        if (tm && tm->cleanup_function) (*tm->cleanup_function)();
        return;
      }


    if (tm->func)
      { // we are already parallel, use nested tasks
        // startup for inner function not supported ...
        // if (startup_function) (*startup_function)();
        if (antasks == -1) antasks = tm->GetNumThreads();

        if (antasks == 1)
          {
            TaskInfo ti;
            ti.task_nr = 0;
            ti.ntasks = 1;
            ti.thread_nr = 0; ti.nthreads = 1;
            afunc(ti);
            return;
          }
        
        atomic<int> endcnt(antasks);
        tm->AddTask (afunc, endcnt);
        while (endcnt > 0)
          {
            tm->ProcessTask();
          }
        
        // if (cleanup_function) (*cleanup_function)();
        return;
      }


    class StartStop
    {
    public:
      StartStop(const function<void(TaskInfo&)> & afunc)
      {
        if (trace)
          trace->StartJob(GetTaskManager()->jobnr, afunc.target_type());
      }
      ~StartStop()
      {
        if (trace)
          trace->StopJob();
      }
    };
    
    if (antasks == 1)
      {
        StartStop startstop(afunc);
        //if (trace)
        // trace->StartJob(jobnr, afunc.target_type());
        tm->jobnr++;
        if (tm->startup_function) (*tm->startup_function)();
        TaskInfo ti;
        ti.task_nr = 0;
        ti.ntasks = 1;
        ti.thread_nr = 0; ti.nthreads = 1;
        {
          RegionTracer t(ti.thread_nr, tm->jobnr, RegionTracer::ID_JOB, ti.task_nr);
          afunc(ti);
        }
        if (tm->cleanup_function) (*tm->cleanup_function)();
        // if (trace)
        // trace->StopJob();
        return;
      }

    if (antasks == -1) antasks = tm->GetNumThreads();

    StartStop startstop(afunc);    
    // if (trace)
    // trace->StartJob(jobnr, afunc.target_type());

    tm->func = &afunc;

    tm->ntasks.store (antasks); // , memory_order_relaxed);
    tm->ex = nullptr;


    tm->nodedata[0]->start_cnt.store (0, memory_order_relaxed);

    tm->jobnr++;
    
    for (int j = 0; j < tm->num_nodes; j++)
      tm->nodedata[j]->participate |= 1;

    if (tm->startup_function) (*tm->startup_function)();
    
    int thd = 0;
    int thds = tm->GetNumThreads();
    int mynode = tm->num_nodes * thd/thds;

    IntRange mytasks = Range(int(tm->ntasks)).Split (mynode, tm->num_nodes);
    NodeData & mynode_data = *((tm->nodedata)[mynode]);

    TaskInfo ti;
    ti.nthreads = thds;
    ti.thread_nr = thd;
    // ti.nnodes = num_nodes;
    // ti.node_nr = mynode;

    try
      {
        while (1)
          {
            int mytask = mynode_data.start_cnt++;
            if (mytask >= mytasks.Size()) break;
            
            ti.task_nr = mytasks.First()+mytask;
            ti.ntasks = tm->ntasks;

            {
              RegionTracer t(ti.thread_nr, tm->jobnr, RegionTracer::ID_JOB, ti.task_nr);
              (*tm->func)(ti); 
            }
          }

      }
    catch (Exception & e)
      {
        {
          lock_guard<mutex> guard(tm->copyex_mutex);
          delete tm->ex;
          tm->ex = new Exception (e);
          mynode_data.start_cnt = mytasks.Size();
        }
      }

    if (tm->cleanup_function) (*tm->cleanup_function)();
    
    for (int j = 0; j < tm->num_nodes; j++)
      if (tm->workers_on_node[j])
        {
          while (tm->complete[j] != tm->jobnr)
          {
#ifdef NETGEN_ARCH_AMD64
            _mm_pause();
#endif // NETGEN_ARCH_AMD64
          }
        }

    tm->func = nullptr;
    if (tm->ex)
      throw Exception (*tm->ex);

    // if (trace)
    //    trace->StopJob();
  }
    
  void TaskManager :: Loop(int thd)
  {
    /*
    static Timer tADD("add entry counter");
    static Timer tCASready1("spin-CAS ready tick1");
    static Timer tCASready2("spin-CAS ready tick2");
    static Timer tCASyield("spin-CAS yield");
    static Timer tCAS1("spin-CAS wait");
    static Timer texit("exit zone");
    static Timer tdec("decrement");
    */
    task_manager = this;
    thread_id = thd;

    worker_data = &workers[thd];

    if(thd == 0)
        return;

    int thds = num_threads;
    int mynode = num_nodes * thd/thds;
    NodeData & mynode_data = *(nodedata[mynode]);

    TaskInfo ti;
    ti.nthreads = thds;
    ti.thread_nr = thd;
    // ti.nnodes = num_nodes;
    // ti.node_nr = mynode;

      
#ifdef USE_NUMA
    numa_run_on_node (mynode);
#endif
    active_workers++;
    workers_on_node[mynode]++;
    int jobdone = 0;


#ifdef USE_MKL
    auto mkl_max = mkl_get_max_threads();
    mkl_set_num_threads_local(1);
#endif

    
    size_t no_job_counter = 0;
    while (!done)
      {
        if (complete[mynode] > jobdone)
          jobdone = complete[mynode];

        if (jobnr == jobdone)
          {
            no_job_counter++;
            // RegionTracer t(ti.thread_nr, tCASyield, ti.task_nr);
            while (ProcessTask()) no_job_counter = 0; // do the nested tasks
                   
            if(sleep)
              std::this_thread::sleep_for(std::chrono::microseconds(sleep_usecs));
            else if(no_job_counter > 30000)
              std::this_thread::sleep_for(std::chrono::microseconds(1000));
            else if(no_job_counter > 20000)
              std::this_thread::sleep_for(std::chrono::microseconds(100));
            else if(no_job_counter > 10000)
              std::this_thread::sleep_for(std::chrono::microseconds(10));
            else
              {
#ifdef WIN32
                std::this_thread::yield();
#else  // WIN32
                sched_yield();
#endif // WIN32
              }
            continue;
          }

        {
          // RegionTracer t(ti.thread_nr, tADD, ti.task_nr);

          // non-atomic fast check ...
          if ( (mynode_data.participate & 1) == 0) continue;

          int oldval = mynode_data.participate += 2;
          if ( (oldval & 1) == 0)
            { // job not active, going out again
              mynode_data.participate -= 2;
              continue;
            }
        }

        if (startup_function) (*startup_function)();
        
        IntRange mytasks = Range(int(ntasks)).Split (mynode, num_nodes);
          
        try
          {
            while (1)
              {
                if (mynode_data.start_cnt >= mytasks.Size()) break;
		int mytask = mynode_data.start_cnt.fetch_add(1, memory_order_relaxed);
                if (mytask >= mytasks.Size()) break;
                
                ti.task_nr = mytasks.First()+mytask;
                ti.ntasks = ntasks;
                no_job_counter = 0;
                
                {
                  RegionTracer t(ti.thread_nr, jobnr, RegionTracer::ID_JOB, ti.task_nr);
                  (*func)(ti);
                }
              }

          }
        catch (Exception & e)
          {
            {
              // cout << "got exception in TM" << endl; 
              lock_guard<mutex> guard(copyex_mutex);
              delete ex;
              ex = new Exception (e);
              mynode_data.start_cnt = mytasks.Size();
            }
          }

#ifndef __MIC__
        atomic_thread_fence (memory_order_release);     
#endif // __MIC__

        if (cleanup_function) (*cleanup_function)();

        jobdone = jobnr;

        mynode_data.participate-=2;

	{
	  int oldpart = 1;
	  if (mynode_data.participate.compare_exchange_strong (oldpart, 0))
	    {
              if (jobdone < jobnr.load())
                { // reopen gate
                  mynode_data.participate |= 1;                  
                }
              else
                {
                  if (mynode != 0)
                    mynode_data.start_cnt = 0;
                  complete[mynode] = jobnr.load(); 
                }
	    }	      
	}
      }
    

#ifdef USE_MKL
    mkl_set_num_threads_local(mkl_max);
#endif

    workers_on_node[mynode]--;
    active_workers--;
  }


  std::list<std::tuple<std::string,double>> TaskManager :: Timing ()
  {
    /*
    list<tuple<string,double>>timings;
    double time =
      RunTiming
      ( [&] ()
        {
          ParallelJob ( [] (TaskInfo ti) { ; } ,
                        TasksPerThread(1) );
        });
    timings.push_back (make_tuple("parallel job with 1 task per thread", time*1e9));
    
    time =
      RunTiming
      ( [&] ()
        {
          ParallelJob ( [] (TaskInfo ti) { ; } ,
                        TasksPerThread(10) );
        });
    timings.push_back (make_tuple("parallel job with 10 tasks per thread", time*1e9));

    time =
      RunTiming
      ( [&] ()
        {
          ParallelJob ( [] (TaskInfo ti) { ; } ,
                        TasksPerThread(100) );
        });
    timings.push_back (make_tuple("parallel job with 100 tasks per thread", time*1e9));

    return timings;
    */


    
    // this is the old function moved from the py-interface:
    std::list<std::tuple<std::string,double>>timings;           
    double starttime, time;
    double maxtime = 0.5;
    size_t steps;
    
    starttime = WallTime();
    steps = 0;
    do
      {
        for (size_t i = 0; i < 1000; i++)
          ParallelJob ( [] (TaskInfo ti) { ; },
                        TasksPerThread(1));
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("ParallelJob 1 task/thread", time/steps*1e9));


    starttime = WallTime();
    steps = 0;
    do
      {
        for (size_t i = 0; i < 1000; i++)
          ParallelJob ( [] (TaskInfo ti) { ; },
                        TasksPerThread(100));
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("ParallelJob 100 task/thread", time/steps*1e9));

    
    starttime = WallTime();
    steps = 0;
    do
      {
        for (int k = 0; k < 10000; k++)
          {
            SharedLoop2 sl(1000);
            steps += 1;
          }
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop init", time/steps*1e9));
    
    starttime = WallTime();
    steps = 0;
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            SharedLoop sl(5);
            ParallelJob ( [&sl] (TaskInfo ti)
                          {
                            for (auto i : sl)
                              (void)i;  // silence warning
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("short SharedLoop", time/steps*1e9));
    

    starttime = WallTime();
    steps = 0;
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            SharedLoop sl1(5), sl2(5), sl3(5), sl4(5), sl5(5);
            ParallelJob ( [&sl1, &sl2, &sl3, &sl4, &sl5] (TaskInfo ti)
                          {
                            for (auto i : sl1)
                              (void)i;  // silence warning
                            for (auto i : sl2)
                              (void)i;  // silence warning
                            for (auto i : sl3)
                              (void)i;  // silence warning
                            for (auto i : sl4)
                              (void)i;  // silence warning
                            for (auto i : sl5)
                              (void)i;  // silence warning
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("5 short SharedLoops", time/steps*1e9));
    

    starttime = WallTime();
    steps = 0;
    SharedLoop2 sl2(5);
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            sl2.Reset(5);
            ParallelJob ( [&sl2] (TaskInfo ti)
                          {
                            for (auto i : sl2)
                              (void)i;  // silence warning                              
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("short SharedLoop2", time/steps*1e9));

    {
    starttime = WallTime();
    steps = 0;
    SharedLoop2 sl1(5), sl2(5), sl3(5), sl4(5), sl5(5);
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            sl1.Reset(5);
            sl2.Reset(5);
            sl3.Reset(5);
            sl4.Reset(5);
            sl5.Reset(5);
            ParallelJob ( [&sl1,&sl2,&sl3,&sl4,&sl5] (TaskInfo ti)
                          {
                            for (auto i : sl1)
                              (void)i;  // silence warning                              
                            for (auto i : sl2)
                              (void)i;  // silence warning                              
                            for (auto i : sl3)
                              (void)i;  // silence warning                              
                            for (auto i : sl4)
                              (void)i;  // silence warning                              
                            for (auto i : sl5)
                              (void)i;  // silence warning                              
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("5 short SharedLoop2", time/steps*1e9));
    }

    
    starttime = WallTime();
    steps = 0;
    {
    SharedLoop2 sl(1000);
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            sl.Reset(1000);
            ParallelJob ( [&sl] (TaskInfo ti)
                          {
                            for (auto i : sl)
                              (void)i;  // silence warning                               
                          } );
            steps += 1000;
          }
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop2 1000, time per iteration", time/steps*1e9));
    }

    {
    starttime = WallTime();
    steps = 0;
    SharedLoop2 sl(1000000);
    do
      {
        sl.Reset(1000000);
        ParallelJob ( [&sl] (TaskInfo ti)
                      {
                        for (auto i : sl)
                          (void)i;  // silence warning
                      } );
        steps += 1000000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop2 1000000, time per iteration", time/steps*1e9));
    }
    
    return timings;
  }
  
}

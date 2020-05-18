#include <mystdlib.h>
#include "meshing.hpp"
#include <netgen_version.hpp>


namespace netgen
{

  class TraceGlobal
  {
    string name;
  public:
    TraceGlobal(string _name) : name(_name) { cout << "init global " << name << endl; }
    ~TraceGlobal() { cout << "exit global " << name << endl; }
  };
  
  // stringstream emptystr;
  // ostream * testout = &emptystr;
  // testout -> clear(ios::failbit);

  // ostream * testout = &cout;

  // NetgenOutStream * testout = new NetgenOutStream;

  const string netgen_version = NETGEN_VERSION;

  ostream * mycout = &cout;
  ostream * myerr = &cerr;

  // some functions (visualization) still need a global mesh
  // TraceGlobal glob1("global1");
  DLL_HEADER shared_ptr<Mesh> mesh;
  DLL_HEADER shared_ptr<NetgenGeometry> ng_geometry;
  // TraceGlobal glob2("global2");

  // global communicator for netgen
  // DLL_HEADER NgMPI_Comm ng_comm;
  
  weak_ptr<Mesh> global_mesh;
  void SetGlobalMesh (shared_ptr<Mesh> m)
  {
    PrintMessage(5, "set global mesh");
    global_mesh = m;
  }
  
  // true if netgen was started using the netgen executable
  // false if netgen.gui was imported from python
  DLL_HEADER bool netgen_executable_started = false;
  
  //  Flags parameters;
  int silentflag = 0;
  int testmode = 0;

  volatile multithreadt multithread;

  string ngdir = ".";

  void Ng_PrintDest(const char * s)
  {
    if (id == 0)
      (*mycout) << s << flush;
  }

  DLL_HEADER void MyError(const char * ch)
  {
    cout << ch;
    (*testout) << "Error !!! " << ch << endl << flush;
  }

  static double starttimea;
  void ResetTime ()
  {
    starttimea = WallTime();
  }

  double GetTime ()
  {
    return WallTime() - starttimea;
  }



  mutex tcl_todo_mutex;

  int h_argc = 0;
  char ** h_argv = NULL;

  multithreadt :: multithreadt()
  {
    pause =0;
    testmode = 0;
    redraw = 0;
    drawing = 0;
    terminate = 0;
    running = 0;
    percent = 0;
    task = "";
  }

  DebugParameters debugparam;
  bool verbose = 0;

  size_t timestamp = 0;
  /*
  int GetTimeStamp() 
  { 
    return timestamp; 
  }

  int NextTimeStamp()
  {
    timestamp++;
    return timestamp;
  }
  */
}

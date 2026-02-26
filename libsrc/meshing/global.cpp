#include <mystdlib.h>
#include "global.hpp"
#include <netgen_version.hpp>
#include "msghandler.hpp"
#include "meshtype.hpp"

namespace netgen
{
  class NetgenGeometry;
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
  void(*on_set_global_mesh)(shared_ptr<Mesh>) = nullptr;

  void SetGlobalMesh (shared_ptr<Mesh> m)
  {
    if(GetGlobalMesh() == m)
      return;
    PrintMessage(5, "set global mesh");
    global_mesh = m;
    if (on_set_global_mesh)
      on_set_global_mesh(m);
  }

  shared_ptr<Mesh> GetGlobalMesh ()
  {
    try {
      return global_mesh.lock();
    } catch (const bad_weak_ptr & e) {
      return nullptr;
    }
  }
  
  // true if netgen was started using the netgen executable
  // false if netgen.gui was imported from python
  DLL_HEADER bool netgen_executable_started = false;
  
  //  Flags parameters;
  int silentflag = 0;
  int testmode = 0;


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

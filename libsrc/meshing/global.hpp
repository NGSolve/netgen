#ifndef NETGEN_GLOBAL_HPP
#define NETGEN_GLOBAL_HPP


/**************************************************************************/
/* File:   global.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
  global functions and variables
*/

#include <mydefs.hpp>

namespace netgen
{
  using namespace ngcore;
  ///
  DLL_HEADER extern double GetTime ();
  DLL_HEADER extern void ResetTime ();

  ///
  DLL_HEADER extern int testmode;

  /// calling parameters
  // extern Flags parameters;

  // extern DLL_HEADER MeshingParameters mparam;

  DLL_HEADER extern mutex tcl_todo_mutex;

  class DLL_HEADER multithreadt
  {
  public:
    int pause;
    int testmode;
    int redraw;
    int drawing;
    int terminate;
    int running;
    double percent;
    const char * task;
    bool demorunning;
    string * tcl_todo = new string("");  // tcl commands set from parallel thread
    multithreadt();
  };

  DLL_HEADER extern volatile multithreadt multithread;

  class DebugParameters;
  class Mesh;

  DLL_HEADER extern string ngdir;
  DLL_HEADER extern DebugParameters debugparam;
  DLL_HEADER extern bool verbose;

  DLL_HEADER extern int h_argc;
  DLL_HEADER extern char ** h_argv;


  DLL_HEADER extern weak_ptr<Mesh> global_mesh;
  DLL_HEADER void SetGlobalMesh (shared_ptr<Mesh> m);

  // global communicator for netgen (dummy if no MPI)
  // extern DLL_HEADER NgMPI_Comm ng_comm;
  
} // namespace netgen

#endif // NETGEN_GLOBAL_HPP

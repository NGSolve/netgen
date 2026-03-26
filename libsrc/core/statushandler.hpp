#ifndef NETGEN_CORE_STATUSHANDLER
#define NETGEN_CORE_STATUSHANDLER

#include <string>
#include "utils.hpp"

namespace ngcore
{
  
  class NGCORE_API multithreadt
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
    std::string * tcl_todo = new std::string("");  // tcl commands set from parallel thread
    multithreadt();
  };

  NGCORE_API extern volatile multithreadt multithread;

  
  extern NGCORE_API void SetStatMsg(const std::string& s);

  extern NGCORE_API void PushStatus(const std::string& s);
  extern NGCORE_API void PushStatusF(const std::string& s);
  extern NGCORE_API void PopStatus();
  extern NGCORE_API void SetThreadPercent(double percent);
  extern NGCORE_API void GetStatus(std::string & s, double & percentage);
}
#endif

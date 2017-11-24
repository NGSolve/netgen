#ifndef FILE_NG_PROFILER
#define FILE_NG_PROFILER

/**************************************************************************/
/* File:   profiler.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jan. 2005                                                  */
/**************************************************************************/



#ifdef VTRACE
#include "vt_user.h"
#else
  #define VT_USER_START(n)
  #define VT_USER_END(n)
  #define VT_TRACER(n)
#endif


// #define USE_TSC
#ifdef USE_TSC
#include <x86intrin.h>   // for __rdtsc()  CPU time step counter
#endif

namespace netgen
{

class NgProfiler
{
  enum { SIZE = 1000 };

  static long int tottimes[SIZE];
  static long int starttimes[SIZE];
  static long int counts[SIZE];
  static string names[SIZE];
  static int usedcounter[SIZE];

  int total_timer;
public: 
  NgProfiler();
  ~NgProfiler();
  static int CreateTimer (const string & name);
#ifndef USE_TSC
  static void StartTimer (int nr) 
  { 
    starttimes[nr] = clock(); counts[nr]++; 
    // VT_USER_START (const_cast<char*> (names[nr].c_str())); 
    VT_USER_START ( (char * const) (names[nr].c_str())); 
  }
  static void StopTimer (int nr) 
  { 
    tottimes[nr] += clock()-starttimes[nr]; 
    VT_USER_END (const_cast<char*> (names[nr].c_str())); 
  }
#else
  static void StartTimer (int nr) 
  { 
    starttimes[nr] = __rdtsc(); counts[nr]++; 
    // VT_USER_START (const_cast<char*> (names[nr].c_str())); 
    // VT_USER_START ( (char * const) (names[nr].c_str())); 
  }
  static void StopTimer (int nr) 
  { 
    tottimes[nr] += __rdtsc()-starttimes[nr]; 
    VT_USER_END (const_cast<char*> (names[nr].c_str())); 
  }
#endif
  
  
  //static void Print (ostream & ost);
  static void Print (FILE * prof);

  static void ClearTimers ();

  class RegionTimer
  {
    int nr;
  public:
    RegionTimer (int anr) : nr(anr)
      { StartTimer (nr); }
    ~RegionTimer () { StopTimer (nr); }
  };
};

}

#endif

#ifndef FILE_MYSTDLIB
#define FILE_MYSTDLIB

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cctype>
#include <ctime>
#include <cstring>
#include <climits>
#include <algorithm>
#include <memory>
#include <thread>
#include <mutex>
#include <atomic>
#include <optional>
#include <cassert>

#include <new>
#include <string>
#include <typeinfo>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*** Windows headers ***/
#ifdef _MSC_VER
# define WIN32_LEAN_AND_MEAN
# ifndef NO_PARALLEL_THREADS
#  ifdef MSVC_EXPRESS
#  else
// #   include <afxwin.h>
// #   include <afxmt.h>
#  endif // MSVC_EXPRESS
# endif
// # include <windows.h>
# undef WIN32_LEAN_AND_MEAN
// # include <winnt.h>
#else // Not using MC VC++
#endif


// using namespace std;
namespace netgen
{
  using namespace std;  
}


#endif


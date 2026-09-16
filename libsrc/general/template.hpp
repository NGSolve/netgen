#ifndef FILE_TEMPLATE
#define FILE_TEMPLATE

/**************************************************************************/
/* File:   template.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <core/utils.hpp>

namespace netgen 
{
  using namespace ngcore;
/*
   templates, global types, defines and variables
*/

  DLL_HEADER extern const std::string netgen_version;

///     The following value may be adapted to the hardware !
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000
#endif


// #include <iostream>
/** output stream for testing.
  testout is opened by main */


/**
  BOOL is a typedef for boolean variables
  */
// typedef int BOOL;








/*




///
template <class T>
inline T min2 (T a, T b)
{
  ///
  return (a < b) ? a : b;
}
///
template <class T>
inline T max2 (T a, T b)
{
  ///
  return (a > b) ? a : b;
}
///
template <class T>
inline T min3 (T a, T b, T c)
{
  ///
  return (a < b) ? (a < c) ? a : c
    : (b < c) ? b : c;
}
///
template <class T>
inline T max3 (T a, T b, T c)
{
  ///
  return (a > b) ? ((a > c) ? a : c)
    : ((b > c) ? b : c);
}

///


///
template <class T>
inline int sgn (T a)
{
  return (a > 0) ? 1 : (   ( a < 0) ? -1 : 0 );
}

///
template <class T>
inline T sqr (const T a)
{
  return a * a; 
}

///
template <class T>
inline T pow3 (const T a)
{
  return a * a * a; 
}
*/



/*
template <class T>
void BubbleSort (int size, T * data);

template <class T>
void MergeSort (int size, T * data, T * help);
*/



}

namespace netgen
{

inline void SetInvalid (int & i) { i = -1; }
inline bool IsInvalid (int i) { return i == -1; }
inline size_t HashValue (int i, size_t size) { return (113*size_t(i)) % size; }




}

#endif

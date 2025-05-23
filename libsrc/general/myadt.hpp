#ifndef FILE_MYADT
#define FILE_MYADT

/**************************************************************************/
/* File:   myadt.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   include for all abstract data types
*/


#include "../include/mystdlib.h"
#include "../include/mydefs.hpp"


#include <core/ngcore.hpp>
namespace netgen
{
  using namespace ngcore;
  using NgException = Exception;
}
#include "parthreads.hpp"
// #include "moveablemem.hpp"
// #include "dynamicmem.hpp"

#include "template.hpp"
#include "ngarray.hpp"
#include "table.hpp"
#include "hashtabl.hpp"


#include "ngbitarray.hpp"
#include "spbita2d.hpp"

#include "seti.hpp"
#include "optmem.hpp"
// #include "autoptr.hpp"
#include "sort.hpp"
#include "stack.hpp"
#include "mystring.hpp"

// #include "mpi_interface.hpp"
#include "netgenout.hpp"

#endif

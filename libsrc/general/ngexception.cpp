/**************************************************************************/
/* File:   ngexception.cpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 02                                                    */
/**************************************************************************/

#include <myadt.hpp>

namespace netgen
{
  //using namespace netgen;



  NgException :: NgException (const string & s) 
    : m_what(s)
  {
    ; 
  }


  NgException :: ~NgException () 
  {
    ;
  }

  /// append string to description
  void NgException :: Append (const string & s)
  { 
    m_what += s;
  }

}

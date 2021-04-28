#ifndef NETGEN_CORE_XBOOL_HPP
#define NETGEN_CORE_XBOOL_HPP

/**************************************************************************/
/* File:   xbool.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Nov. 07                                                    */
/**************************************************************************/


namespace ngcore
{
  // an extended bool with values false/maybe/true

  enum TMAYBE { maybe };
  
  class xbool
  {
    uint8_t state;

  public:
    xbool (bool b) : state(b ? 2 : 0) { ; }
    xbool (TMAYBE /* x */) : state(1) { ; }
    xbool () = default;
    xbool (const xbool &) = default;
    
    xbool & operator= (bool b) { state = b ? 2 : 0; return *this; }
    xbool & operator= (TMAYBE /* x */) { state = 1; return *this; }
    
    bool IsTrue () const { return state == 2; }
    bool IsMaybe () const { return state == 1; }
    bool IsFalse () const { return state == 0; }
    bool IsMaybeTrue() const { return state >= 1; } 
    bool IsMaybeFalse() const { return state <= 1; } 
    friend ostream & operator<< (ostream & ost, xbool xb);
  };


  static char output[] = "0?1";
  inline ostream & operator<< (ostream & ost, xbool xb)
  {
    return ost << output[xb.state];
  }
  
} // namespace ngcore

#endif // NETGEN_CORE_XBOOL_HPP

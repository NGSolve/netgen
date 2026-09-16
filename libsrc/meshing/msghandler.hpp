#ifndef FILE_MSGHANDLER
#define FILE_MSGHANDLER

/**************************************************************************/
/* File:   msghandler.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

#include <filesystem>
#include <sstream>
#include <string>

#include <core/ngstream.hpp>

namespace netgen
{
  DLL_HEADER extern int printwarnings;
  DLL_HEADER extern int printerrors;
  DLL_HEADER extern int printfnstart;

  DLL_HEADER void Ng_PrintDest(const char * s);

  extern DLL_HEADER void PrintDot(char ch = '.');

  namespace detail
  {
    template <typename T>
    inline void PrintArg (std::ostream & ost, const T & t) { ost << t; }
    inline void PrintArg (std::ostream & ost, const std::filesystem::path & p) { ost << p.string(); }

    template <typename... Args>
    void PrintDestArgs (const char * prefix, const char * suffix, const Args &... args)
    {
      std::ostringstream ost;
      ost << prefix;
      (PrintArg(ost, args), ...);
      ost << suffix;
      Ng_PrintDest(ost.str().c_str());
    }
  }

  //Message Pipeline:

  //importance: importance of message: 1=very important, 3=middle, 5=low, 7=unimportant
  template <typename... Args>
  void PrintMessage (int importance, const Args &... args)
  {
    if (importance <= printmessage_importance)
      detail::PrintDestArgs(" ", "\n", args...);
  }

  // CR without line-feed
  template <typename... Args>
  void PrintMessageCR (int importance, const Args &... args)
  {
    if (importance <= printmessage_importance)
      detail::PrintDestArgs(" ", "\r", args...);
  }

  template <typename... Args>
  void PrintFnStart (const Args &... args)
  {
    if (printfnstart)
      detail::PrintDestArgs(" Start Function: ", "\n", args...);
  }

  template <typename... Args>
  void PrintWarning (const Args &... args)
  {
    if (printwarnings)
      detail::PrintDestArgs(" WARNING: ", "\n", args...);
  }

  template <typename... Args>
  void PrintError (const Args &... args)
  {
    if (printerrors)
      detail::PrintDestArgs(" ERROR: ", "\n", args...);
  }

  template <typename... Args>
  void PrintFileError (const Args &... args)
  {
    if (printerrors)
      detail::PrintDestArgs(" FILE ERROR: ", "\n", args...);
  }

  template <typename... Args>
  void PrintSysError (const Args &... args)
  {
    if (printerrors)
      detail::PrintDestArgs(" SYSTEM ERROR: ", "\n", args...);
  }

  template <typename... Args>
  void PrintUserError (const Args &... args)
  {
    detail::PrintDestArgs(" USER ERROR: ", "\n", args...);
  }

  template <typename... Args>
  void PrintTime (const Args &... args)
  {
    if (printmessage_importance >= 3)
      detail::PrintDestArgs(" Time = ", "\n", args...);
  }


  inline void PushStatusF(const std::string& s)
  {
    PushStatus (s);
    PrintFnStart(s);
  }

}


#endif


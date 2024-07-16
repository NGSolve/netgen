#ifndef NETGEN_CORE_FLAGS_HPP
#define NETGEN_CORE_FLAGS_HPP


/**************************************************************************/
/* File:   flags.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

#include <iostream>
#include <memory>
#include <string>
#include <any>

#include "array.hpp"
#include "symboltable.hpp"
#include "xbool.hpp"

namespace ngcore
{

  /** 
      A storage for command-line flags.
      The flag structure maintains string flags, numerical flags,
      define flags, string list flags, num list flags.
  */
  class NGCORE_API Flags 
  {
    /// string flags
    SymbolTable<std::string> strflags;
    /// numerical flags
    SymbolTable<double> numflags;
    /// define flags
    SymbolTable<bool> defflags;
    /// string list flags
    SymbolTable<std::shared_ptr<Array<std::string>>> strlistflags;
    /// numerical list flags
    SymbolTable<std::shared_ptr<Array<double>>> numlistflags;
    /// flags list flags
    SymbolTable<Flags> flaglistflags;
    /// any object can be stored as a flag
    SymbolTable<std::any> anyflags;
  public:
    /// no flags
    Flags ();
    /// copy flags 
    Flags (const Flags & flags);
    /// steal flags
    Flags (Flags && flags);
    ///
    Flags (std::initializer_list<std::string> list);
    /// 
    Flags (std::string f1, std::string f2 = "", std::string f3 = "", std::string f4 = "", std::string f5 = "");
    /// delete mem
    ~Flags ();

    Flags & operator= (const Flags & f2) = default;
    Flags & operator= (Flags && f2) = default;

    void DoArchive(class Archive& ar);

    void Update(const Flags& other);
  
    /// Deletes all flags
    void DeleteFlags ();

    /// Sets string flag, overwrite if exists
    Flags & SetFlag (const char * name, const std::string & val);
    /// Sets string flag, overwrite if exists
    Flags & SetFlag (const char * name, const char * str)
    { return SetFlag (name, std::string(str)); }
    /// Sets numerical flag, overwrite if exists
    Flags & SetFlag (const char * name, double val) &;
    /// Sets numerical flag, overwrite if exists    
    Flags & SetFlag (const char * name, int val)
    { return SetFlag (name, double(val)); }
    /// Sets boolean flag
    Flags & SetFlag (const char * name, bool b = true) &;
    /// Sets numerical flag, overwrite if exists
    Flags & SetFlag (const char * name, Flags & val) &;

    /// Sets string flag, overwrite if exists
    Flags & SetFlag (const std::string & name, const std::string & val);
    Flags & SetFlag (const std::string & name, const char * str)
    { return SetFlag (name, std::string(str)); }
    /// Sets numerical flag, overwrite if exists
    Flags &  SetFlag (const std::string & name, double val);
    /// Sets numerical flag, overwrite if exists    
    Flags &  SetFlag (const std::string & name, int val)
    { return SetFlag (name, double(val)); }
    /// Sets boolean flag
    Flags &  SetFlag (const std::string & name, bool b = true);
    /// Sets numerical flag, overwrite if exists
    Flags &  SetFlag (const std::string & name, Flags & val);
    /// Sets string array flag
    Flags &  SetFlag (const std::string & name, const Array<std::string> & val);
    /// Sets double array flag
    Flags &  SetFlag (const std::string & name, const Array<double> & val);
    /// Sets any flag
    Flags &  SetFlag(const std::string& name, const std::any& val);


    Flags SetFlag (const char * name, bool b = true) &&;
    Flags SetFlag (const char * name, double val) &&;



    /// Save flags to file
    void SaveFlags (const char * filename) const;
    void SaveFlags (ostream & str) const;
    /// write flags to stream
    void PrintFlags (ostream & ost) const;
    /// Load flags from file
    void LoadFlags (const char * filename, SymbolTable<Flags> * sf = nullptr);
    void LoadFlags (std::istream & str, SymbolTable<Flags> * sf = nullptr);
    /**
       Set command line flag.
       Flag must be in form: -name=hello -val=0.5 -defflag 
       -names=[Joe,Jim] -values=[1,3,4] -solverflags=*abc
    */
    void SetCommandLineFlag (const char * st, SymbolTable<Flags> * sf = nullptr);
    /// Returns string flag, default value if not exists
    std::string GetStringFlag (const std::string & name, const char * def) const;
    /// Returns std::string flag, default value if not exists
    std::string GetStringFlag (const std::string & name, std::string def = "") const;
    /// Returns numerical flag, default value if not exists
    double GetNumFlag (std::string_view name, double def) const;
    /// Returns address of numerical flag, null if not exists
    const double * GetNumFlagPtr (std::string_view name) const;
    /// Returns address of numerical flag, null if not exists
    double * GetNumFlagPtr (const std::string & name);
    /// Returns boolean flag
    // int GetDefineFlag (const char * name) const;
    bool GetDefineFlag (std::string_view name) const  noexcept;
    xbool GetDefineFlagX (std::string_view name) const  noexcept;
    /// Returns string list flag, empty array if not exist
    const Array<std::string> & GetStringListFlag (const std::string & name) const;
    /// Returns num list flag, empty array if not exist
    const Array<double> & GetNumListFlag (const std::string & name) const;
    /// Returns flag list flag, empty flag if not exist
    const Flags & GetFlagsFlag (const std::string & name) const;
    const std::any& GetAnyFlag (const std::string& name) const;


    /// Test, if string flag is defined
    bool StringFlagDefined (std::string_view name) const noexcept;
    /// Test, if num flag is defined
    bool NumFlagDefined (std::string_view name) const noexcept;
    /// Test, if num flag is defined
    bool FlagsFlagDefined (std::string_view name) const noexcept;
    /// Test, if string list flag is defined
    bool StringListFlagDefined (std::string_view name) const noexcept;
    /// Test, if num list flag is defined
    bool NumListFlagDefined (std::string_view name) const noexcept;
    bool AnyFlagDefined (std::string_view name) const noexcept;

    /// number of string flags
    int GetNStringFlags () const { return strflags.Size(); }
    /// number of num flags
    int GetNNumFlags () const { return numflags.Size(); }
    /// number of num flags
    int GetNFlagsFlags () const { return flaglistflags.Size(); }
    /// number of define flags
    int GetNDefineFlags () const { return defflags.Size(); }
    /// number of string-list flags
    int GetNStringListFlags () const { return strlistflags.Size(); }
    /// number of num-list flags
    int GetNNumListFlags () const { return numlistflags.Size(); }
    int GetNAnyFlags() const { return anyflags.Size(); }

    ///
    const std::string & GetStringFlag (int i, std::string & name) const
    { name = strflags.GetName(i); return strflags[i]; }
    double GetNumFlag (int i, std::string & name) const
    { name = numflags.GetName(i); return numflags[i]; }
    bool GetDefineFlag (int i, std::string & name) const
    { name = defflags.GetName(i); return defflags[i]; }
    const std::shared_ptr<Array<double>> GetNumListFlag (int i, std::string & name) const
    { name = numlistflags.GetName(i).c_str(); return numlistflags[i]; }
    const std::shared_ptr<Array<std::string>> GetStringListFlag (int i, std::string & name) const
    { name = strlistflags.GetName(i); return strlistflags[i]; }
    const Flags & GetFlagsFlag (int i, std::string & name) const
    { name = flaglistflags.GetName(i); return flaglistflags[i]; }
    const std::any& GetAnyFlag(int i, std::string& name) const
    { name = anyflags.GetName(i); return anyflags[i]; }
  };

  /// Print flags
  inline std::ostream & operator<< (std::ostream & s, const Flags & flags)
  {
    flags.PrintFlags (s);
    return s;
  }
} // namespace ngcore

  
#endif // NETGEN_CORE_FLAGS_HPP


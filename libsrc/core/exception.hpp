#ifndef NETGEN_CORE_EXCEPTION_HPP
#define NETGEN_CORE_EXCEPTION_HPP

#include <sstream>         // for stringstream
#include <stdexcept>       // for exception
#include <string>          // for string

#include "ngcore_api.hpp"  // for NGCORE_API


namespace ngcore
{

  NGCORE_API std::string GetBackTrace();

  // Exception for code that shouldn't be executed
  class NGCORE_API UnreachableCodeException : public std::exception
  {
    const char* what() const noexcept override
    {
      return "Shouldn't get here, something went wrong!";
    }
  };

  // Default exception class
  class NGCORE_API Exception : public std::exception
  {
    /// a verbal description of the exception
    std::string m_what;
  public:
    Exception() = default;
    Exception(const Exception&) = default;
    Exception(Exception&&) = default;
    Exception(const std::string& s); //  : m_what(s) {}
    Exception(const char* s); //  : m_what(s) {}
    ~Exception() override = default;

    Exception& operator =(const Exception&) = default;
    Exception& operator =(Exception&&) noexcept = default;

    /// append string to description
    Exception & Append (const std::string & s) { m_what += s; return *this; }
  /// append string to description
    Exception & Append (const char * s) { m_what += s; return *this; }

    /// verbal description of exception
    const std::string & What() const { return m_what; }

    /// implement virtual function of std::exception
    const char* what() const noexcept override { return m_what.c_str(); }
  };
  
  NGCORE_API void ThrowException(const std::string & s);
  NGCORE_API void ThrowException(const char * s);
  
  // Out of Range exception
  class NGCORE_API RangeException : public Exception
  {
  public:
    /// where it occurs, index, minimal and maximal indices
    RangeException (const std::string & where,
                    int ind, int imin, int imax) : Exception("")
      {
        std::stringstream str;
        str << where << ": index " << ind << " out of range [" << imin << "," << imax << ")\n";
        Append (str.str());
        Append (GetBackTrace());
      }

    template<typename T>
    RangeException(const std::string& where, const T& value)
    {
      std::stringstream str;
      str << where << " called with wrong value " << value << "\n";
      Append(str.str());
    }
  };

  // Exception used if no simd implementation is available to fall back to standard evaluation
  class NGCORE_API ExceptionNOSIMD : public Exception
  { public: using Exception::Exception; };
} // namespace ngcore

#define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
#define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)

// Convenience macro to append file name and line of exception origin to the string
#define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))

#ifdef NETGEN_ENABLE_CHECK_RANGE
#define NETGEN_CHECK_RANGE(value, min, max_plus_one) \
  { if ((value)<(min) ||  (value)>=(max_plus_one)) \
      throw ngcore::RangeException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", (value), (min), (max_plus_one)); }
#define NETGEN_CHECK_SHAPE(a,b) \
  { if(a.Shape() != b.Shape()) \
      throw ngcore::Exception(__FILE__": shape don't match"); }
#else // NETGEN_ENABLE_CHECK_RANGE
#define NETGEN_CHECK_RANGE(value, min, max)
#define NETGEN_CHECK_SHAPE(a,b)

#endif // NETGEN_ENABLE_CHECK_RANGE


  
#endif // NETGEN_CORE_EXCEPTION_HPP

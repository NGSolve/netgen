#ifndef NETGEN_CORE_EXCEPTION_HPP
#define NETGEN_CORE_EXCEPTION_HPP

#include <cstddef>
#include <sstream>         // for stringstream
#include <stdexcept>       // for exception
#include <string>          // for string

#include "ngcore_api.hpp"  // for NGCORE_API
#include "utils.hpp"       // for ToString


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
    Exception(std::string_view s1, std::string_view s2);
    Exception(std::string_view s1, std::string_view s2, std::string_view s3);
    ~Exception() override = default;

    [[noreturn]] static void Throw (std::string_view s1);
    [[noreturn]] static void Throw (std::string_view s1, std::string_view s2);
    [[noreturn]] static void Throw (std::string_view s1, std::string_view s2, std::string_view s3);
    
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
  
  [[noreturn]] NGCORE_API void ThrowException(const std::string & s);
  [[noreturn]] NGCORE_API void ThrowException(const char * s);
  
  // Out of Range exception
  class NGCORE_API RangeException : public Exception
  {
  public:
    /// where it occurs, index, minimal and maximal indices
    RangeException (// const std::string & where,
                    const char * where,
                    ptrdiff_t ind, ptrdiff_t imin, ptrdiff_t imax);
    /*
    : Exception("")
      {
        std::stringstream str;
        str << where << ": index " << ind << " out of range [" << imin << "," << imax << ")\n";
        Append (str.str());
        Append (GetBackTrace());
      }
    */
    template<typename T>
    RangeException(const std::string& where, const T& value)
    {
      std::stringstream str;
      str << where << " called with wrong value " << value << "\n";
      Append(str.str());
    }
  };

  [[noreturn]] NGCORE_API void ThrowRangeException(const char * s, ptrdiff_t ind, ptrdiff_t imin, ptrdiff_t imax);
  [[noreturn]] NGCORE_API void ThrowNotTheSameException(const char * s, ptrdiff_t a, ptrdiff_t b);
  
  
  // Exception used if no simd implementation is available to fall back to standard evaluation
  class NGCORE_API ExceptionNOSIMD : public Exception
  { public: using Exception::Exception; };

  template <typename T>
  struct IsSafe {
    constexpr operator bool() const { return false; } };

  namespace detail {
    template <typename T, typename Tmin, typename Tmax>
    inline static constexpr void CheckRange(const char * s, const T& n, Tmin first, Tmax next)
    {
      if constexpr (!IsSafe<decltype(n)>())
        if (n<first || n>=next)
          ThrowRangeException(s, ptrdiff_t(n), ptrdiff_t(first), ptrdiff_t(next));
    }

    template <typename Ta, typename Tb>
    inline static constexpr void CheckSame(const char * s, const Ta& a, const Tb& b)
    {
     if constexpr (!IsSafe<decltype(a)>() || !IsSafe<decltype(b)>())
      if(a != b)
      {
        if constexpr(std::is_integral_v<decltype(a)> && std::is_same_v<decltype(a),decltype(b)>)
          ThrowNotTheSameException(s, long(a), long(b)); \
        else
          throw Exception(std::string(s) + "\t: not the same, a="+ToString(a) + ", b="+ngcore::ToString(b) + GetBackTrace());
      }
    }
  } // namespace detail
} // namespace ngcore

#define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
#define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)

// Convenience macro to append file name and line of exception origin to the string
#define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))

#if defined(NETGEN_ENABLE_CHECK_RANGE) && !defined(__CUDA_ARCH__)
#define NETGEN_CHECK_RANGE(value, min, max_plus_one) ngcore::detail::CheckRange(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", value, min, max_plus_one);
#define NETGEN_CHECK_SAME(a,b) ngcore::detail::CheckSame(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", a, b);

#define NETGEN_NOEXCEPT 
#else // defined(NETGEN_ENABLE_CHECK_RANGE) && !defined(__CUDA_ARCH__)
#define NETGEN_CHECK_RANGE(value, min, max)
#define NETGEN_CHECK_SAME(a,b)
// #define NETGEN_CHECK_SHAPE(a,b)
#define NETGEN_NOEXCEPT noexcept
#endif // defined(NETGEN_ENABLE_CHECK_RANGE) && !defined(__CUDA_ARCH__)


  
#endif // NETGEN_CORE_EXCEPTION_HPP

#ifndef NETGEN_CORE_UTILS_HPP
#define NETGEN_CORE_UTILS_HPP

#include <chrono>
#include <map>
#include <ostream>
#include <sstream>
#include <string>

#ifdef WIN32
#include <intrin.h>   // for __rdtsc()  CPU time step counter
#else
#include <x86intrin.h>   // for __rdtsc()  CPU time step counter
#endif // WIN32

#include "ngcore_api.hpp"       // for NGCORE_API

namespace ngcore
{
  // MPI rank, nranks TODO: Rename
  extern NGCORE_API int id, ntasks;
  
  NGCORE_API std::string Demangle(const char* typeinfo);

#if defined(__GNUC__)
  inline bool likely (bool x) { return bool(__builtin_expect(long(x), 1L)); }
  inline bool unlikely (bool x) { return bool(__builtin_expect(long(x), 0L)); }
#else
  inline bool likely (bool x) { return x; }
  inline bool unlikely (bool x) { return x; }
#endif

  using TClock = std::chrono::system_clock;
  extern NGCORE_API const std::chrono::time_point<TClock> wall_time_start;

  // Time in seconds since program start
  inline double WallTime () noexcept
  {
      std::chrono::time_point<TClock> now = TClock::now();
      std::chrono::duration<double> elapsed_seconds = now-wall_time_start;
      return elapsed_seconds.count();
  }

  // High precision clock counter register
  using TTimePoint = size_t;
  extern NGCORE_API double ticks_per_second;

  inline TTimePoint GetTimeCounter() noexcept
  {
      return TTimePoint(__rdtsc());
  }

  template <class T>
  inline std::string ToString (const T& t)
  {
      std::stringstream ss;
      ss << t;
      return ss.str();
  }

  template<typename T1, typename T2>
  std::ostream& operator << (std::ostream& ost, const std::map<T1,T2>& map)
  {
    for(auto& val : map)
      ost << "\n" << val.first << ": " << val.second;
    return ost;
  }

  template <class T>
  NETGEN_INLINE void Swap (T & a, T & b)
  {
      T temp = std::move(a);
      a = std::move(b);
      b = std::move(temp);
  }

  template <typename T>
  class AlignedAlloc
  {
    protected:
      static void * aligned_malloc(size_t s)
      {
        // Assume 16 byte alignment of standard library
        if(alignof(T)<=16)
            return malloc(s);
        else
            return  _mm_malloc(s, alignof(T));
      }

      static void aligned_free(void *p)
      {
        if(alignof(T)<=16)
            free(p);
        else
            _mm_free(p);
      }

  public:
    void * operator new (size_t s, void *p) { return p; }
    void * operator new (size_t s) { return aligned_malloc(s); }
    void * operator new[] (size_t s) { return aligned_malloc(s); }
    void operator delete (void * p) { aligned_free(p); }
    void operator delete[] (void * p) { aligned_free(p); }
  };

} // namespace ngcore

#endif // NETGEN_CORE_UTILS_HPP

#ifndef NETGEN_CORE_UTILS_HPP
#define NETGEN_CORE_UTILS_HPP

#include <atomic>
#include <chrono>
#include <filesystem>
#include <map>
#include <ostream>
#include <sstream>
#include <string>

#include "ngcore_api.hpp"       // for NGCORE_API and CPU arch macros

#if defined(__APPLE__) && defined(NETGEN_ARCH_ARM64)
#include <mach/mach_time.h>
#endif

#ifdef NETGEN_ARCH_AMD64
#ifdef WIN32
#include <intrin.h>   // for __rdtsc()  CPU time step counter
#else
#include <x86intrin.h>   // for __rdtsc()  CPU time step counter
#endif // WIN32
#endif // NETGEN_ARCH_AMD64

namespace ngcore
{
  // MPI rank, nranks TODO: Rename
  // [[deprecated("don't use global id/ntasks")]]       
  extern NGCORE_API int id;
  // [[deprecated("don't use global id/ntasks")]]         
  extern NGCORE_API int ntasks;
  
  NGCORE_API std::string Demangle(const char* typeinfo);

  template<typename T>
  std::string GetName(const T& obj)
  { return Demangle(typeid(obj).name()); }

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
  extern NGCORE_API double seconds_per_tick;

  inline TTimePoint GetTimeCounter() noexcept
  {
#if defined(__APPLE__) && defined(NETGEN_ARCH_ARM64)
    return mach_absolute_time();
#elif defined(NETGEN_ARCH_AMD64)
    return __rdtsc();
#elif defined(NETGEN_ARCH_ARM64) && defined(__GNUC__)
    // __GNUC__ is also defined by CLANG. Use inline asm to read Generic Timer
    unsigned long long tics;
    __asm __volatile("mrs %0, CNTVCT_EL0" : "=&r" (tics));
    return tics;
#else
#warning "Unsupported CPU architecture"
    return 0;
#endif
  }

  template <class T>
  inline std::string ToString (const T& t)
  {
      std::stringstream ss;
      ss << t;
      return ss.str();
  }

  inline std::string ToLower( const std::string & s )
  {
    std::string res;
    res.reserve(s.size());

    for(auto & c : s)
        res.push_back(tolower(c));

    return res;
  }

  inline std::string ToLower( const std::filesystem::path & p )
  {
    return ToLower(p.string());
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

  // checks if string starts with sequence
  inline bool StartsWith(const std::string& str, const std::string& start)
  {
    if(start.size() > str.size())
      return false;
    return std::equal(start.begin(), start.end(), str.begin());
  }

  // checks if string ends with sequence
  inline bool EndsWith(const std::string& str, const std::string& end)
  {
    if(end.size() > str.size())
      return false;
    return std::equal(end.rbegin(), end.rend(), str.rbegin());
  }

  template<typename T>
  NETGEN_INLINE std::atomic<T> & AsAtomic (T & d)
  {
    return reinterpret_cast<std::atomic<T>&> (d);
  }

  NETGEN_INLINE double AtomicAdd( double & sum, double val )
  {
      std::atomic<double> & asum = AsAtomic(sum);
      double current = asum.load();
      while (!asum.compare_exchange_weak(current, current + val))
          ;
      return current;
  }

  template<typename T>
  NETGEN_INLINE T AtomicMin( T & minval, T val )
  {
      std::atomic<T> & aminval = AsAtomic(minval);
      T current = aminval.load();
      while (!aminval.compare_exchange_weak(current, std::min(current, val)))
          ;
      return current;
  }

  template<typename T>
  NETGEN_INLINE T AtomicMax( T & maxval, T val )
  {
      std::atomic<T> & amaxval = AsAtomic(maxval);
      T current = amaxval.load();
      while (!amaxval.compare_exchange_weak(current, std::max(current, val)))
          ;
      return current;
  }

  namespace detail
  {
    template<typename T>
    struct IndexTypeHelper
    {
    private:
      template<typename T2>
      static constexpr auto check(T2* t) -> typename T2::index_type { return *t; }
      static constexpr size_t check(...);

    public:
      using type = decltype(check((T*) nullptr)); // NOLINT
    };
  
  } // namespace detail

  // Get index type of object. If object has a typedef index_type it is this type, else size_t
  template<typename T>
  using index_type = typename detail::IndexTypeHelper<T>::type;

  class MyMutex
  {
    std::atomic<bool> m;
  public:
    MyMutex() { m.store(false, std::memory_order_relaxed); }
    void lock()
    {
      bool should = false;
      while (!m.compare_exchange_weak(should, true))
        {
          should = false;
#ifdef NETGEN_ARCH_AMD64
          _mm_pause();
#endif // NETGEN_ARCH_AMD64
        }
    }
    void unlock()
    {
      m = false;
    }
  };

  class MyLock
  {
    MyMutex & mutex;
  public:
    MyLock (MyMutex & amutex) : mutex(amutex) { mutex.lock(); }
    ~MyLock () { mutex.unlock(); }
  };

  NGCORE_API int GetCompiledSIMDSize();
  NGCORE_API bool IsRangeCheckEnabled();

  NGCORE_API std::filesystem::path GetTempFilename();

} // namespace ngcore

#endif // NETGEN_CORE_UTILS_HPP

#ifndef NETGEN_CORE_UTILS_HPP
#define NETGEN_CORE_UTILS_HPP

#include <atomic>
#include <chrono>
#include <filesystem>
#include <map>
#include <ostream>
#include <optional>
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
#elif defined(__EMSCRIPTEN__)
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
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



  template <class T>
  void SaveBin (std::ostream & ost, const T & val)
  {
    const char * cp = reinterpret_cast<const char*> (&val);
    for (unsigned j = 0; j < sizeof(T); j++)
      ost.put(cp[j]);
  }
  
  
  template <class T>
  void LoadBin (std::istream & ist, T & val)
  {
    char * cp = reinterpret_cast<char*> (&val);
    for (unsigned j = 0; j < sizeof(T); j++)
      ist.get(cp[j]);
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

  
  /// min of 2 values
  template <class T>
  NETGEN_INLINE T min2 (T a, T b)
  {
    return (a < b) ? a : b;
  }
  
  /// max of 2 values
  template <class T>
  NETGEN_INLINE T max2 (T a, T b)
  {
    return (a > b) ? a : b;
}
  
  /// min of 3 values
  template <class T>
  NETGEN_INLINE T min3 (T a, T b, T c)
  {
  return (a < b) ? (a < c) ? a : c
    : (b < c) ? b : c;
  }
  
  /// max of 3 values
  template <class T>
  NETGEN_INLINE T max3 (T a, T b, T c)
  {
    ///
    return (a > b) ? ((a > c) ? a : c)
      : ((b > c) ? b : c);
  }
  
  
  /// sign of value (+1, 0, -1)
  template <class T>
  NETGEN_INLINE int sgn (T a)
  {
    return (a > 0) ? 1 : ( ( a < 0) ? -1 : 0 );
  }
  
  /// square element 
  template <class T>
  NETGEN_INLINE constexpr T sqr (const T a)
  {
    return a * a; 
  }
  
  /// element to the third power
  template <class T>
  NETGEN_INLINE T pow3 (const T a)
  {
    return a * a * a; 
  }
  
  

  NETGEN_INLINE double IfPos (double a, double b, double c) { return a>0 ? b : c; }
  NETGEN_INLINE double IfZero (double a, double b, double c) { return a==0. ? b : c; }

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

  

  
  template <int N> using IC = std::integral_constant<int,N>;  // needed for Iterate

  
  namespace detail {
    template <typename T, typename Enable = int>
    struct IsIC_trait {
      static constexpr auto check() { return false; }
    };
    
    template <typename T>
    struct IsIC_trait<T, std::enable_if_t<std::is_same_v<T, IC<T::value>> == true, int> > {
      static constexpr auto check() { return true; }  
    };
  }
  
  template <typename T>
  constexpr bool is_IC() {
    return detail::IsIC_trait<T>::check();
  }
  
  

  
  template <int NUM, typename FUNC>
  NETGEN_INLINE void Iterate (FUNC f)
  {
    if constexpr (NUM > 1) Iterate<NUM-1> (f);
  if constexpr (NUM >= 1) f(IC<NUM-1>());
  }
  
  
  template <int NUM, typename FUNC>
  NETGEN_INLINE void Switch (size_t nr, FUNC f)
  {
    if (NUM-1 == nr) f(IC<NUM-1>());
    if constexpr (NUM > 1) Switch<NUM-1> (nr, f);
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

  NGCORE_API void* GetRawSymbol( std::string func_name );

  template <typename TFunc>
  TFunc GetSymbol( std::string func_name )
  {
    return reinterpret_cast<TFunc>(GetRawSymbol(func_name));
  }

  // Class to handle/load shared libraries
  class NGCORE_API SharedLibrary
  {
    std::filesystem::path lib_name;
    std::optional<std::filesystem::path> directory_to_delete = std::nullopt;
    void *lib = nullptr;

  public:
    SharedLibrary() = default;
    SharedLibrary(const std::filesystem::path & lib_name_, std::optional<std::filesystem::path> directory_to_delete_ = std::nullopt, bool global = false );

    SharedLibrary(const SharedLibrary &) = delete;
    SharedLibrary & operator =(const SharedLibrary &) = delete;

    ~SharedLibrary();

    template <typename TFunc>
    TFunc GetSymbol( std::string func_name )
    {
      return reinterpret_cast<TFunc>(GetRawSymbol(func_name));
    }

    void Load( const std::filesystem::path & lib_name_, bool global = true);
    void Unload();
    void* GetRawSymbol( std::string func_name );
  };

} // namespace ngcore

#endif // NETGEN_CORE_UTILS_HPP

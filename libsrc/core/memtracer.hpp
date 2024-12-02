#ifndef NETGEN_CORE_MEMTRACER_HPP
#define NETGEN_CORE_MEMTRACER_HPP

#include <array>
#include <chrono>
#include <string>

#include "array.hpp"
#include "logging.hpp"
#include "paje_trace.hpp"
#include "utils.hpp"

namespace ngcore
{

  class MemoryTracer;

  namespace detail
  {
    //Type trait to check if a class implements a 'void SetMemoryTacing(int)' function
    template<typename T>
    struct has_StartMemoryTracing
    {
    private:
      template<typename T2>
      static constexpr auto check(T2*) ->
        typename std::is_same<decltype(std::declval<T2>().StartMemoryTracing()),void>::type;
      template<typename>
      static constexpr std::false_type check(...);
      using type = decltype(check<T>(nullptr)); // NOLINT
    public:
      static constexpr bool value = type::value;
    };
  } // namespace detail

  class MemoryTracer
  {
#if defined(NETGEN_TRACE_MEMORY) && !defined(__CUDA_ARCH__)
    NGCORE_API static std::vector<std::string> names;
    NGCORE_API static std::vector<int> parents;

    #if defined(NETGEN_CHECK_RANGE)
    NGCORE_API static std::atomic<size_t> total_memory;
    mutable size_t allocated_memory = 0;
    #endif // NETGEN_CHECK_RANGE

    static int CreateId(const std::string& name = "")
    {
      int id = names.size();
      names.push_back(name);
      parents.push_back(0);
      if(id==10*8*1024)
        std::cerr << "Allocated " << id << " MemoryTracer objects" << std::endl;
      return id;
    }
    mutable int id = 0;

    public:

    MemoryTracer( std::string name )
    {
      id = CreateId(name);
    }

    MemoryTracer() { }

    MemoryTracer(const MemoryTracer & tracer)
    {
      (*this) = tracer;
    }

    MemoryTracer(MemoryTracer && tracer)
    {
      (*this) = std::move(tracer);
    }

    MemoryTracer & operator=(const MemoryTracer & tracer) {
      if(tracer.id)
        id = CreateId(names[tracer.id]);
      return *this;
    }

    MemoryTracer & operator=(MemoryTracer && tracer) {
      ngcore::Swap(id, tracer.id);

      #if defined(NETGEN_CHECK_RANGE)
      ngcore::Swap(allocated_memory, tracer.allocated_memory);
      #endif // NETGEN_CHECK_RANGE

      return *this;
    }

    template <typename... TRest>
    MemoryTracer( std::string name, TRest & ... rest )
    {
      id = CreateId(name);
      Track(rest...);
    }

    #if defined(NETGEN_CHECK_RANGE)
    // check if all memory was freed when object is destroyed
    ~MemoryTracer()
    {
      NETGEN_CHECK_SAME(allocated_memory, 0);
    }
    #endif // NETGEN_CHECK_RANGE

    NETGEN_INLINE void Alloc(size_t size) const
    {
      #if defined(NETGEN_CHECK_RANGE)
      // Trace also nameless Memtracer objects if range checks are active
      if(!id && size)
        id = CreateId();
      #endif // NETGEN_CHECK_RANGE

      if(id && trace)
        trace->AllocMemory(id, size);

      #if defined(NETGEN_CHECK_RANGE)
      if(id)
      {
        allocated_memory += size;
        total_memory += size;
      }
      #endif // NETGEN_CHECK_RANGE
    }

    void Free(size_t size) const
    {
      if(id && trace)
        trace->FreeMemory(id, size);

      #if defined(NETGEN_CHECK_RANGE)
      if(id)
      {
        // check if we have at least size bytes of memory currently allocated (such that allocated_memory doesn't get negative)
        NETGEN_CHECK_RANGE(allocated_memory, static_cast<ptrdiff_t>(size), std::numeric_limits<ptrdiff_t>::max());
        allocated_memory -= size;
        total_memory -= size;
        #endif // NETGEN_CHECK_RANGE
      }
    }

    int GetId() const { return id; }

    template <typename T1, typename... TRest>
    void Track( T1 & obj, const std::string& name, TRest & ... rest ) const
    {
      Track(obj, name);
      Track(rest...);
    }

    template<typename T>
    void Track( T & obj, const std::string& name ) const
    {
      obj.GetMemoryTracer().Activate(obj, name);
      parents[obj.GetMemoryTracer().GetId()] = id;
    }

    static std::string GetName(int id)
    {
      return names[id];
    }

    std::string GetName() const
    {
      return names[id];
    }

    template<typename T>
    void Activate(T& me, const std::string& name) const
    {
      if(!id)
        {
          const_cast<MemoryTracer*>(this)->id = CreateId(name);
          if constexpr(detail::has_StartMemoryTracing<T>::value)
            me.StartMemoryTracing();
        }
      else
        SetName(name);
    }

    void SetName(const std::string& name) const
    {
      names[id] = name;
    }


    static const std::vector<std::string> & GetNames() { return names; }
    static const std::vector<int> & GetParents() { return parents; }
    static size_t GetTotalMemory()
    {
      #if defined(NETGEN_CHECK_RANGE)
      return total_memory;
      #else
      return 0;
      #endif // NETGEN_CHECK_RANGE
    }
#else // defined(NETGEN_TRACE_MEMORY) && !defined(__CUDA_ARCH__)
  public:
    MemoryTracer() {}
    MemoryTracer( std::string /* name */ ) {}
    template <typename... TRest>
    MemoryTracer( std::string /* name */, TRest & ... ) {}

    void Alloc(size_t /* size */) const {}
    void Free(size_t /* size */) const {}
    int GetId() const { return 0; }

    template <typename... TRest>
    void Track(TRest&...) const {}

    static std::string GetName(int /* id */) { return ""; }
    std::string GetName() const { return ""; }
    void SetName(std::string /* name */) const {}
    static size_t GetTotalMemory() { return 0; }
#endif // NETGEN_TRACE_MEMORY
  };
} // namespace ngcore

#endif // NETGEN_CORE_MEMTRACER_HPP

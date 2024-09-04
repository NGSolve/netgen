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

    static int CreateId(const std::string& name)
    {
      int id = names.size();
      names.push_back(name);
      parents.push_back(0);
      if(id==10*8*1024)
        std::cerr << "Allocated " << id << " MemoryTracer objects" << std::endl;
      return id;
    }
    int id;

    public:

    MemoryTracer( std::string name )
    {
      id = CreateId(name);
    }

    // not tracing
    MemoryTracer() : id(0) {}

    template <typename... TRest>
    MemoryTracer( std::string name, TRest & ... rest )
    {
      id = CreateId(name);
      Track(rest...);
    }

    NETGEN_INLINE void Alloc(size_t size) const
    {
      if(id && trace)
        trace->AllocMemory(id, size);
    }

    void Free(size_t size) const
    {
      if(id && trace)
        trace->FreeMemory(id, size);
    }

    void Swap(size_t mysize, MemoryTracer& other, size_t other_size) const
    {
      if(!trace || (id == 0 && other.id == 0))
        return;
      if(id == 0)
        return trace->ChangeMemory(other.id, mysize - other_size);
      if(other.id == 0)
        return trace->ChangeMemory(id, other_size - mysize);

      // first decrease memory, otherwise have artificial/wrong high peak memory usage
      if(mysize<other_size)
        {
          trace->ChangeMemory(other.id, mysize-other_size);
          trace->ChangeMemory(id, other_size-mysize);
        }
      else
        {
          trace->ChangeMemory(id, other_size-mysize);
          trace->ChangeMemory(other.id, mysize-other_size);
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
#else // defined(NETGEN_TRACE_MEMORY) && !defined(__CUDA_ARCH__)
  public:
    MemoryTracer() {}
    MemoryTracer( std::string /* name */ ) {}
    template <typename... TRest>
    MemoryTracer( std::string /* name */, TRest & ... ) {}

    void Alloc(size_t /* size */) const {}
    void Free(size_t /* size */) const {}
    void Swap(...) const {}
    int GetId() const { return 0; }

    template <typename... TRest>
    void Track(TRest&...) const {}

    static std::string GetName(int /* id */) { return ""; }
    std::string GetName() const { return ""; }
    void SetName(std::string /* name */) const {}
#endif // NETGEN_TRACE_MEMORY
  };
} // namespace ngcore

#endif // NETGEN_CORE_MEMTRACER_HPP

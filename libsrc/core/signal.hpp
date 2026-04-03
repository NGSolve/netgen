#ifndef NGCORE_SIGNALS_HPP
#define NGCORE_SIGNALS_HPP

#include <list>
#include <map>
#include <functional>

namespace ngcore
{
  template<typename ... ParameterTypes>
  class Signal
  {
  private:
    std::list<std::function<bool(ParameterTypes...)>> funcs;
    bool is_emitting;
  public:
    Signal() : is_emitting(true) {}

    template<typename Cls, typename FUNC>
    void Connect(Cls* self, FUNC f)
    {
      auto ptr = self->weak_from_this();
      auto func = [ptr, f](ParameterTypes... args)
        {
          if (ptr.expired())
            return false;
          f(args...);
          return true;
        };
      funcs.push_back(func);
    }

    inline void Emit(ParameterTypes ...args)
    {
      if(is_emitting)
        funcs.remove_if([&](auto& f){ return !f(args...); });
    }

    inline bool SetEmitting(bool emitting)
    {
      bool was_emitting = is_emitting;
      is_emitting = emitting;
      return was_emitting;
    }
    inline bool GetEmitting() const { return is_emitting; }
  };




  class SimpleSignal
  {
  private:
    // std::map<void*,std::function<void()>> funcs;
    std::list<std::pair<void*,std::function<void()>>> funcs;
  public:
    SimpleSignal() = default;

    template<typename FUNC>
    void Connect(void* var, FUNC f)
    {
      // funcs[var] = f;
      funcs.push_back ( { var, f } );
    }

    void Remove(void* var)
    {
      // funcs.erase(var);
      funcs.remove_if([&] (auto var_f) { return var_f.first==var; });
    }

    inline void Emit()
    {
      for (auto [key,f] : funcs)
        f();
    }
  };

  
} // namespace ngcore

#endif // NGCORE_SIGNALS_HPP

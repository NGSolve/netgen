#ifndef NGCORE_SIGNALS_HPP
#define NGCORE_SIGNALS_HPP

#include <list>
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
    void connect(Cls* self, FUNC f)
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

    inline void emit(ParameterTypes ...args)
    {
      if(is_emitting)
        funcs.remove_if([&](auto& f){ return !f(args...); });
    }

    inline void setEmitting(bool emitting) { is_emitting = emitting; }
    inline bool getEmitting() const { return is_emitting; }
  };
} // namespace ngcore

#endif // NGCORE_SIGNALS_HPP

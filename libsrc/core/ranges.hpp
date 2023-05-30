#ifndef NETGEN_CORE_RANGES_HPP
#define NETGEN_CORE_RANGES_HPP

#include <iterator>

namespace ngcore
{
  template<typename Iterator>
  class AdapterRange
  {
    Iterator _begin,_end;
  public:
    AdapterRange(Iterator abegin, Iterator aend) : _begin(abegin), _end(aend) { ; }
    Iterator begin() const { return _begin; }
    Iterator end() const { return _end; }
  };

  template<typename FUNC>
  class FilterAdapter
  {
    FUNC f;
  public:
    FilterAdapter(FUNC af) : f(af) { ; }
    FUNC GetFunction() const { return f; }
  };

  template<typename FUNC, typename Iterator>
  class FilterIterator
  {
    Iterator iter;
    Iterator end;
    FUNC f;
  public:
    FilterIterator(FUNC af, Iterator aiter, Iterator aend)
      :  iter(aiter), end(aend), f(af)
    {
      while(iter!=end && !f(*iter))
        ++iter;
    }
    inline FilterIterator& operator ++()
    {
      ++iter;
      while(iter!=end && !f(*iter))
        ++iter;
      return *this;
    }

    inline bool operator !=(FilterIterator other)
    {
      return iter != other.iter;
    }

    inline bool operator ==(FilterIterator other)
    {
      return iter == other.iter;
    }

    inline decltype(auto) operator *() const
    {
      return *iter;
    }
  };

  template<typename FUNC>
  FilterAdapter<FUNC> filter(FUNC f) { return {f}; }

  template<typename Range, typename FUNC>
  auto operator |(Range&& range, FilterAdapter<FUNC> adapter)
    -> AdapterRange<FilterIterator<FUNC,decltype(std::begin(range))>>
  {
    return {{adapter.GetFunction(),std::begin(range),std::end(range)},
        {adapter.GetFunction(), std::end(range), std::end(range)}};
  }

  template<typename FUNC, typename Iterator>
  class TransformIterator
  {
    FUNC f;
    Iterator iter;
  public:
    TransformIterator(FUNC af, Iterator aiter) : f(af), iter(aiter) { ; }

    TransformIterator& operator++() { ++iter; }
    bool operator !=(TransformIterator other) { return iter != other.iter; }
    decltype(auto) operator *() const { return f(*iter); }
  };

  template<typename FUNC>
  class TransformAdapter
  {
    FUNC f;
  public:
    TransformAdapter(FUNC af) : f(af) { ; }
    FUNC GetFunction() const { return f; }
  };

  template<typename FUNC>
  TransformAdapter<FUNC> transform(FUNC f) { return {f}; }

  template<typename Range, typename FUNC>
  auto operator |(Range&& range, TransformAdapter<FUNC> adapter)
    -> AdapterRange<TransformIterator<FUNC,decltype(std::begin(range))>>
  {
    return {{adapter.GetFunction(), std::begin(range)},
        {adapter.GetFunction(),std::end(range)}};
  }
} // namespace ngcore

#endif // NETGEN_CORE_RANGES_HPP

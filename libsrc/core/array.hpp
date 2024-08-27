#ifndef NETGEN_CORE_ARRAY_HPP
#define NETGEN_CORE_ARRAY_HPP

/**************************************************************************/
/* File:   array.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <cstring>
#include <type_traits>

#include "exception.hpp"
#include "logging.hpp"          // for logger
#include "ngcore_api.hpp"       // for NGCORE_API
#include "type_traits.hpp"      // for all_of_tmpl
#include "localheap.hpp"
#include "memtracer.hpp"
#include "utils.hpp"

namespace ngcore
{
  using std::ostream;

  template <typename ... ARGS> class Tuple 
  { 
  public:
    int Size() const { return 0; }
  };

  template <typename HEAD, typename ... TAIL>
  class Tuple<HEAD, TAIL...> : Tuple<TAIL...>
  {
    typedef Tuple<TAIL...> BASE;
    HEAD head;
  public:
    Tuple () { ; }
    Tuple (HEAD h, TAIL ... t) : Tuple<TAIL...> (t...), head(h) { ; }

    HEAD Head() const { return head; }
    Tuple<TAIL...> Tail() const { return *this; }

    int Size() const { return BASE::Size()+1; }
  };

  template <typename ... ARGS>
  ostream & operator<< (ostream & ost, Tuple<ARGS...> /* tup */)
  {
    return ost;
  }

  template <typename FIRST, typename ... ARGS>
  ostream & operator<< (ostream & ost, Tuple<FIRST, ARGS...> tup)
  {
    ost << tup.Head() << ", " << tup.Tail();
    return ost;
  }


  template <typename ... ARGS>
  Tuple<ARGS...> MakeTuple (ARGS ... args)
  {
    return Tuple<ARGS...> (args...);
  }


  template <typename AO>
  class AOWrapperIterator
  {
    const AO & ao;
    size_t ind;
  public:
    NETGEN_INLINE AOWrapperIterator (const AO &  aao, size_t ai) 
      : ao(aao), ind(ai) { ; }
    NETGEN_INLINE AOWrapperIterator operator++ (int) 
    { return AOWrapperIterator(ao, ind++); }
    NETGEN_INLINE AOWrapperIterator& operator++ ()
    { ++ind; return *this; }
    NETGEN_INLINE auto operator*() const -> decltype(ao[ind]) { return ao[ind]; }
    NETGEN_INLINE auto operator*() -> decltype(ao[ind]) { return ao[ind]; }
    NETGEN_INLINE bool operator != (AOWrapperIterator d2) { return ind != d2.ind; }
    NETGEN_INLINE bool operator == (AOWrapperIterator d2) { return ind == d2.ind; }
  };



  
  /*
    Some class which can be treated as array
   */
  template <typename T> // , typename TA = T>
  class BaseArrayObject
  {
  public:
    NETGEN_INLINE BaseArrayObject() { ; }

    NETGEN_INLINE const T & Spec() const { return static_cast<const T&> (*this); }
    NETGEN_INLINE size_t Size() const { return Spec().Size(); }
    template <typename T2>
    NETGEN_INLINE bool Contains(const T2 & el) const
    {
      for (size_t i = 0; i < Size(); i++)
        if (Spec()[i] == el)
          return true;
      return false;
    }

    static constexpr size_t ILLEGAL_POSITION = size_t(-1);
    template <typename T2>
    NETGEN_INLINE size_t Pos(const T2 & el) const
    {
      for (size_t i = 0; i < Size(); i++)
        if (Spec()[i] == el)
          return i;
      return ILLEGAL_POSITION;
    }

    template <typename T2>
    NETGEN_INLINE size_t PosSure(const T2 & el) const
    {
      for (size_t i = 0; ; i++)
        if (Spec()[i] == el)
          return i;
    }

    // NETGEN_INLINE auto & operator[] (size_t i) { return Spec()[i]; }
    NETGEN_INLINE auto operator[] (size_t i) const { return Spec()[i]; }
    // NETGEN_INLINE auto begin() const { return Spec().begin(); }
    // NETGEN_INLINE auto end() const { return Spec().end(); }
    NETGEN_INLINE auto begin () const { return AOWrapperIterator<BaseArrayObject> (*this, 0); }
    NETGEN_INLINE auto end () const { return AOWrapperIterator<BaseArrayObject> (*this, Size()); }
  };
  

  
  template <typename T>
  class AOWrapper : public BaseArrayObject<AOWrapper<T>>
  {
    T ar;
  public:
    NETGEN_INLINE AOWrapper (T aar) : ar(aar) { ; }
    NETGEN_INLINE operator T () const { return ar; }
    NETGEN_INLINE size_t Size() const { return ar.Size(); }
    NETGEN_INLINE auto operator[] (size_t i) { return ar[i]; }
    NETGEN_INLINE auto operator[] (size_t i) const { return ar[i]; }
    NETGEN_INLINE AOWrapperIterator<AOWrapper> begin () const { return AOWrapperIterator<AOWrapper> (*this, 0); }
    NETGEN_INLINE AOWrapperIterator<AOWrapper> end () const { return AOWrapperIterator<AOWrapper> (*this, Size()); }
  };

  template <typename T>
  NETGEN_INLINE AOWrapper<const T&> ArrayObject (const T & ar)
  {
    return AOWrapper<const T&> (ar);
  }

  template <typename T>
  NETGEN_INLINE AOWrapper<T> ArrayObject (T && ar)
  {
    return AOWrapper<T> (ar);
  }

  template <typename FUNC>
  auto ArrayObject (size_t s, FUNC f)
  {
    class Dummy
    {
      size_t s;
      FUNC f;
    public:
      Dummy (size_t _s, FUNC _f) : s(_s), f(_f) { ; }
      size_t Size() const { return s; }
      auto operator[] (size_t i) const { return f(i); }
    };
    return ArrayObject(Dummy(s,f));
  }

  template <typename T, typename FUNC>
  auto Substitute (const BaseArrayObject<T> & ao, FUNC f)
  {
    return ArrayObject(ao.Size(),
                       [&ao,f] (size_t i) { return f(ao[i]); });
  }
  



  /**
     nothing more but a new type for a C array.
     return value for Addr - operator of array 
  */
  template <class T>
  class CArray
  {
  protected:
    /// the data
    T * data;
  public:

    /// initialize array 
    NETGEN_INLINE CArray () { data = 0; }

    /// provide size and memory
    NETGEN_INLINE CArray (T * adata) 
      : data(adata) { ; }

    /// Access array
    NETGEN_INLINE T & operator[] (size_t i) const
    {
      return data[i]; 
    }

    NETGEN_INLINE operator T* () const { return data; }
  };


  template <typename  T>
  constexpr T IndexBASE () { return T(0); }


  class IndexFromEnd
  {
    ptrdiff_t i;
  public:
    constexpr IndexFromEnd (ptrdiff_t ai) : i(ai) { }
    IndexFromEnd operator+ (ptrdiff_t inc) const { return i+inc; }
    IndexFromEnd operator- (ptrdiff_t dec) const { return i-dec; }
    // operator ptrdiff_t () const { return i; }
    ptrdiff_t Value() const { return i; }
  };

  constexpr IndexFromEnd END(0);
  
  
  template <class T, class IndexType = size_t> class FlatArray;


  template <typename TELEM, typename IndexType>
  class ArrayIterator
  {
    FlatArray<TELEM, IndexType> ar;
    IndexType ind;
  public:
    NETGEN_INLINE ArrayIterator (FlatArray<TELEM, IndexType> aar, IndexType ai) 
      : ar(aar), ind(ai) { ; }
    NETGEN_INLINE ArrayIterator operator++ (int) 
    { return ArrayIterator(ar, ind++); }
    NETGEN_INLINE ArrayIterator operator++ ()
    { return ArrayIterator(ar, ++ind); }
    // NETGEN_INLINE const TELEM & operator*() const { return ar[ind]; }
    // NETGEN_INLINE TELEM & operator*() { return ar[ind]; }
    NETGEN_INLINE auto operator*() const -> decltype(ar[ind]) { return ar[ind]; }
    NETGEN_INLINE auto operator*() -> decltype(ar[ind]) { return ar[ind]; }
    NETGEN_INLINE bool operator != (ArrayIterator d2) { return ind != d2.ind; }
    NETGEN_INLINE bool operator == (ArrayIterator d2) { return ind == d2.ind; }
  };
  


  template <typename TSIZE>
  class ArrayRangeIterator
  {
    TSIZE ind;
  public:
    NETGEN_INLINE ArrayRangeIterator (TSIZE ai) : ind(ai) { ; }
    NETGEN_INLINE ArrayRangeIterator operator++ (int) { return ind++; }
    NETGEN_INLINE ArrayRangeIterator operator++ () { return ++ind; }
    NETGEN_INLINE TSIZE operator*() const { return ind; }
    NETGEN_INLINE TSIZE Index() { return ind; }
    NETGEN_INLINE operator TSIZE () const { return ind; }
    NETGEN_INLINE bool operator != (ArrayRangeIterator d2) { return ind != d2.ind; }
    NETGEN_INLINE bool operator == (ArrayRangeIterator d2) { return ind == d2.ind; }
  };

  /// a range of integers
  template <typename T>
  class T_Range : public BaseArrayObject <T_Range<T>>
  {
    T first, next;
  public: 
    NETGEN_INLINE T_Range () { ; }
    NETGEN_INLINE T_Range (T n) : first(0), next(n) {;}
    NETGEN_INLINE T_Range (T f, T n) : first(f), next(n) {;}
    template <typename T2>
      NETGEN_INLINE T_Range(T_Range<T2> r2) : first(r2.First()), next(r2.Next()) { ; }
    NETGEN_INLINE T First() const { return first; }
    NETGEN_INLINE T Next() const { return next; }
    NETGEN_INLINE T & First() { return first; }
    NETGEN_INLINE T & Next() { return next; }
    NETGEN_INLINE auto Size() const { return next-first; }
    NETGEN_INLINE T operator[] (size_t i) const { return first+i; }
    NETGEN_INLINE bool Contains (T i) const { return ((i >= first) && (i < next)); }
    NETGEN_INLINE T_Range Modify(int inc_beg, int inc_end) const
    { return T_Range(first+inc_beg, next+inc_end); }
    NETGEN_INLINE ArrayRangeIterator<T> begin() const { return first; }
    NETGEN_INLINE ArrayRangeIterator<T> end() const { return next; }

    NETGEN_INLINE T_Range Split (size_t nr, int tot) const
    {
      T diff = next-first;
      return T_Range (first + nr * diff / tot,
                      first + (nr+1) * diff / tot);
    }
    // NETGEN_INLINE operator IntRange () const { return IntRange(first,next); }
  };

  using IntRange = T_Range<size_t>;

  template <typename T>
  NETGEN_INLINE T_Range<T> Range (T a, T b)
  {
    return T_Range<T>(a,b);
  }

  template<typename T>
  NETGEN_INLINE auto Range (const T& ao)
    -> typename std::enable_if<has_range<T>, decltype(std::declval<T>().Range())>::type
  { return ao.Range(); }

  template <typename T>
  NETGEN_INLINE T_Range<T> Range_impl (T n, std::true_type)
  {
    return T_Range<T> (0, n);
  }

  template <typename TA>
  NETGEN_INLINE auto Range_impl (const TA & ao, std::false_type)
    -> T_Range<index_type<TA>>
  {
    return T_Range<index_type<TA>> (IndexBASE<index_type<TA>>(),
                                   IndexBASE<index_type<TA>>() + index_type<TA>(ao.Size()));
  }

  /*
    Range(obj) will create a range in using the following steps:
    
    * if obj is an integral type it will create T_Range<type(obj)>(0, obj)
    * if obj has a function Range() it will return obj.Range()
    * if obj has a typedef index_type it will return
      T_Range<index_type>(IndexBASE<index_type>(), IndexBASE<index_type>() + index_type(obj.Size()))
    * else it will return T_Range<size_t> (0, obj.Size())

   */
  template <typename T>
  auto Range(const T & x)
    -> typename std::enable_if<std::is_integral_v<T> || !has_range<T>,
                               decltype(Range_impl(x, std::is_integral<T>()))>::type {
    return Range_impl(x, std::is_integral<T>());
  }


  NETGEN_INLINE IntRange operator+ (const IntRange & range, int shift)
  {
    return IntRange (range.First()+shift, range.Next()+shift);
  }

  NETGEN_INLINE IntRange operator+ (int shift, const IntRange & range)
  {
    return IntRange (range.First()+shift, range.Next()+shift);
  }

  NETGEN_INLINE IntRange operator* (int scale, const IntRange & range)
  {
    return IntRange (scale*range.First(), scale*range.Next());
  }

  NETGEN_INLINE IntRange operator* (const IntRange & range, int scale)
  {
    return IntRange (scale*range.First(), scale*range.Next());
  }

  template <typename TI>
  inline ostream & operator<< (ostream & s, T_Range<TI> ir)
  {
    s << "[" << ir.First() << "," << ir.Next() << ")";
    return s;
  }

  template <typename ... ARGS>
  ostream & operator<< (ostream & ost, Tuple<IntRange, ARGS...> tup)
  {
    ost << tup.Head() << ", " << tup.Tail();
    return ost;
  }


  template <typename T>
  inline ostream & operator<< (ostream & ost, const BaseArrayObject<T> & array)
  {
    for (auto i : Range(array.Size()))
      ost << i << ":" << array[i] << std::endl;
    return ost;
  }


  template <typename T, typename TI, typename INDEX_ARRAY>
  class IndirectArray : public BaseArrayObject<IndirectArray<T, TI, INDEX_ARRAY> >
  {
    FlatArray<T,TI> ba;
    const INDEX_ARRAY & ia;

  public:
    NETGEN_INLINE IndirectArray (FlatArray<T,TI> aba,
                          const INDEX_ARRAY & aia)
      : ba(aba), ia(aia) { ; }
    
    NETGEN_INLINE size_t Size() const { return ia.Size(); }
    NETGEN_INLINE T & operator[] (size_t i) const { return ba[ia[i]]; }
    // NETGEN_INLINE T & operator[] (size_t i) { return ba[ia[i]]; }

    NETGEN_INLINE IndirectArray operator= (const T & val) 
    {
      for (auto i : Range(Size()))
        (*this)[i] = val;
      return IndirectArray (ba, ia);
    }

    template <typename T2>
    NETGEN_INLINE IndirectArray operator= (const BaseArrayObject<T2> & a2) 
    {
      for (auto i : Range(Size()))
	(*this)[i] = a2[i];
      return IndirectArray (ba, ia);
    }

    NETGEN_INLINE AOWrapperIterator<IndirectArray> begin() const { return { *this, 0 }; }
    NETGEN_INLINE AOWrapperIterator<IndirectArray> end() const { return { *this, Size() }; }
  };


  /**
     A simple array container.
     Array represented by size and data-pointer.
     No memory allocation and deallocation, must be provided by user.
     Helper functions for printing. 
     Optional range check by macro NETGEN_CHECK_RANGE
  */
  template <class T, class IndexType>
  class FlatArray : public BaseArrayObject<FlatArray<T,IndexType> >
  {
  protected:
    static constexpr IndexType BASE = IndexBASE<IndexType>();
    /// the size
    size_t size = 0;
    /// the data
    T * __restrict data = nullptr;
  public:
    typedef T value_type;
    typedef IndexType index_type;
    using BaseArrayObject<FlatArray>::ILLEGAL_POSITION;

    /// initialize array 
    NETGEN_INLINE FlatArray () = default;
    // { ; } // size = 0; data = 0; }

    /// copy constructor allows size-type conversion 
    NETGEN_INLINE FlatArray (const FlatArray & a2) = default;
    // : size(a2.Size()), data(a2.data) { ; } 

    /// provide size and memory
    NETGEN_INLINE FlatArray (size_t asize, T * adata) 
      : size(asize), data(adata) { ; }
    
    /// memory from local heap
    NETGEN_INLINE FlatArray(size_t asize, Allocator & lh)
      : size(asize), data(new (lh) T[asize])
    { ; }

    NETGEN_INLINE FlatArray(size_t asize, LocalHeap & lh)
      : size(asize), data (lh.Alloc<T> (asize))
    { ; }

    /// the size
    NETGEN_INLINE size_t Size() const { return size; }

    /// the data
    NETGEN_INLINE T* Data() const { return data; }

    /// Fill array with value val
    NETGEN_INLINE const FlatArray & operator= (const T & val) const
    {
      size_t hsize = size;
      T * hdata = data;
      for (size_t i = 0; i < hsize; i++)
        hdata[i] = val;
      return *this;
    }

    /// copies array
    NETGEN_INLINE const FlatArray & operator= (const FlatArray & a2) const
    {
      size_t hsize = size;
      T * hdata = data;
      for (size_t i = 0; i < hsize; i++) hdata[i] = a2.data[i];
      return *this;
    }

    template <typename T2>
    NETGEN_INLINE const FlatArray & operator= (const BaseArrayObject<T2> & a2) const
    {
      size_t hsize = size;
      T * hdata = data;
      auto p2 = a2.begin();
      for (size_t i = 0; i < hsize; i++, p2++) hdata[i] = *p2;
      return *this;
    }

    template <typename T2, std::enable_if_t<std::is_function<T2>::value>>
    NETGEN_INLINE const FlatArray & operator= (const T2 & func) const
    {
      for (size_t i = 0; i < size; i++)
        data[i] = func(i+BASE);
      return *this;
    }

//     template <typename T2>
//     const FlatArray operator= (ParallelValue<T2> val);
//     template <typename T2>
//     const FlatArray operator= (ParallelFunction<T2> val);

    /// copies pointers
    NETGEN_INLINE const FlatArray & Assign (const FlatArray & a2)
    {
      size = a2.size;
      data = a2.data;
      return *this;
    }

    /// assigns memory from local heap
    NETGEN_INLINE const FlatArray & Assign (size_t asize, LocalHeap & lh)
    {
      size = asize;
      data = lh.Alloc<T> (asize);
      return *this;
    }

    /// Access array. range check by macro NETGEN_CHECK_RANGE
    NETGEN_INLINE T & operator[] (IndexType i) const
    {
      NETGEN_CHECK_RANGE(i,BASE,size+BASE);
      return data[i-BASE]; 
    }
  
    NETGEN_INLINE T_Range<index_type> Range () const
    {
      return T_Range<index_type> (BASE, size+BASE);
    }
    
    NETGEN_INLINE const CArray<T> Addr (size_t pos) const
    {
      return CArray<T> (data+pos-BASE);
    }

    // const CArray<T> operator+ (int pos)
    // { return CArray<T> (data+pos); }
    NETGEN_INLINE T * operator+ (size_t pos) const { return data+pos; }

    /// access last element. check by macro NETGEN_CHECK_RANGE
    T & Last () const
    {
      NETGEN_CHECK_RANGE(size-1,0,size);
      return data[size-1];
    }

    /// takes sub-array starting from position pos
    NETGEN_INLINE const FlatArray<T> Part (size_t pos)
    {
      return FlatArray<T> (size-pos, data+pos);
    }

    /// takes subsize elements starting from position pos
    NETGEN_INLINE const FlatArray<T> Part (size_t pos, size_t subsize)
    {
      return FlatArray<T> (subsize, data+pos);
    }

    /// takes range starting from position start of end-start elements
    NETGEN_INLINE FlatArray<T> Range (size_t start, size_t end) const
    {
      return FlatArray<T> (end-start, data+start);
    }

    /// takes range starting from position start of end-start elements
    NETGEN_INLINE FlatArray<T> Range (size_t start, IndexFromEnd indend) const
    {
      return this->Range(start, size_t(Size()+indend.Value()));
    }
    
    /// takes range starting from position start of end-start elements
    NETGEN_INLINE FlatArray<T> Range (T_Range<size_t> range) const
    {
      return FlatArray<T> (range.Size(), data+range.First());
    }

    /// takes range starting from position start of end-start elements
    NETGEN_INLINE const FlatArray<T> operator[] (T_Range<IndexType> range) const
    {
      return FlatArray<T> (range.Size(), data+range.First());
    }
    
    template <typename TI1>
    auto operator[] (const BaseArrayObject<TI1> & ind_array) const
    {
      return IndirectArray<T, IndexType, BaseArrayObject<TI1> > (*this, ind_array);
    }

    /// first position of element elem, returns -1 if element not contained in array 
    NETGEN_INLINE size_t Pos(const T & el) const
    {
      for (size_t i = 0; i < Size(); i++)
        if (data[i] == el)
          return i;
      return ILLEGAL_POSITION;
    }

    /// does the array contain element elem ?
    NETGEN_INLINE bool Contains(const T & elem) const
    {
      return Pos(elem) != ILLEGAL_POSITION;
    }
    
    //auto begin() const { return ArrayIterator<T,IndexType> (*this, BASE); }
    // auto end() const { return ArrayIterator<T,IndexType> (*this, size+BASE); }
    NETGEN_INLINE auto begin() const { return data; }
    NETGEN_INLINE auto end() const { return data+Size(); }
  };

  template <typename T>
  FlatArray<T> View (FlatArray<T> fa) { return fa; }

  template <typename T, typename TI>
  auto Max (FlatArray<T,TI> array, typename std::remove_const<T>::type max = std::numeric_limits<T>::min()) -> T
  {
    for (auto & v : array)
      if (v > max) max = v;
    return max;
  }

  template <typename T, typename TI>
  auto Min (FlatArray<T,TI> array, typename std::remove_const<T>::type min = std::numeric_limits<T>::max()) -> T
  {
    for (auto & v : array)
      if (v < min) min = v;
    return min;
  }
  
  /// print array
  template <class T, class TIND>
  inline ostream & operator<< (ostream & s, const FlatArray<T, TIND> & a)
  {
    for (auto i : a.Range())
      s << i << ": " << a[i] << "\n";
    return s;
  }

  /// have arrays the same contents ?
  template <class T1, class T2>
  inline bool operator== (const FlatArray<T1> & a1,
                          const FlatArray<T2> & a2)
  {
    if (a1.Size () != a2.Size()) return 0;
    for (size_t i = 0; i < a1.Size(); i++)
      if (a1[i] != a2[i]) return false;
    return true;
  }
		 
  template <class T1, class T2>
  inline bool operator!= (const FlatArray<T1> & a1,
                          const FlatArray<T2> & a2)
  {
    return !(a1==a2);
  }
  

  /** 
      Dynamic array container.
   
      Array<T> is an automatically increasing array container.
      The allocated memory doubles on overflow. 
      Either the container takes care of memory allocation and deallocation,
      or the user provides one block of data.
  */
  template <class T, class IndexType = size_t>
  class Array : public FlatArray<T, IndexType>
  {
  protected:
    /// physical size of array
    size_t allocsize;
    /// that's the data we have to delete, nullptr for not owning the memory
    T * mem_to_delete;


    using FlatArray<T,IndexType>::size;
    using FlatArray<T,IndexType>::data;
    using FlatArray<T,IndexType>::BASE;

  public:
    using index_type = typename FlatArray<T, IndexType>::index_type;
    /// Generate array of logical and physical size asize
    NETGEN_INLINE explicit Array()
      : FlatArray<T,IndexType> (0, nullptr)
    {
      allocsize = 0; 
      mem_to_delete = nullptr;
    }

    NETGEN_INLINE explicit Array(size_t asize)
      : FlatArray<T,IndexType> (asize, new T[asize])
    {
      allocsize = asize; 
      mem_to_delete = data;
    }


    /// Generate array in user data
    NETGEN_INLINE Array(size_t asize, T* adata, bool ownMemory = false)
      : FlatArray<T> (asize, adata)
    {
      allocsize = asize;
      if(ownMemory)
        mem_to_delete = adata;
      else
        mem_to_delete = nullptr;
    }

    /// Generate array in user data
    template <typename ALLOCATOR>
    NETGEN_INLINE Array(size_t asize, ALLOCATOR & lh)
      : FlatArray<T> (asize, lh)
    {
      allocsize = asize; 
      mem_to_delete = nullptr;
    }

    NETGEN_INLINE Array (Array && a2) 
    {
      mt.Swap(0., a2.mt, sizeof(T) * a2.allocsize);

      size = a2.size; 
      data = a2.data;
      allocsize = a2.allocsize;
      mem_to_delete = a2.mem_to_delete;
      a2.size = 0;
      a2.allocsize = 0;
      a2.data = nullptr;
      a2.mem_to_delete = nullptr;
    }

    /// array copy 
    NETGEN_INLINE explicit Array (const Array & a2)
      : FlatArray<T,IndexType> (a2.Size(), a2.Size() ? new T[a2.Size()] : nullptr)
    {
      if constexpr (std::is_copy_assignable<T>::value)
        {
          allocsize = size;
          mem_to_delete = data;
          for (size_t i = 0; i < size; i++)
            data[i] = a2.data[i];
        }
      
// #ifdef __cpp_exceptions
#ifndef __CUDA_ARCH__
      else
        throw Exception(std::string("cannot copy-construct Array of type ") + typeid(T).name());
#endif      
    }

    
    template <typename TA>
    explicit Array (const BaseArrayObject<TA> & a2)
      : FlatArray<T,IndexType> (a2.Size(), 
                                a2.Size() ? new T[a2.Size()] : nullptr)
    {
      allocsize = size;
      mem_to_delete = data;
      /*
      for (size_t i = 0; i < size; i++)
        data[i] = a2[i];
      */
      auto p2 = a2.begin();
      for (size_t i = 0; i < size; i++, p2++)
        data[i] = *p2;
      
    }

    Array (std::initializer_list<T> list) 
      : FlatArray<T> (list.size(), 
                      list.size() ? new T[list.size()] : NULL)
    {
      allocsize = size;
      mem_to_delete = data;
      size_t cnt = 0;
      for (auto val : list)
        data[cnt++] = val;
    }

    /// array merge-copy
    explicit Array (const Array<T> & a2, const Array<T> & a3)
      : FlatArray<T> (a2.Size()+a3.Size(), 
                      a2.Size()+a3.Size() ? new T[a2.Size()+a3.Size()] : 0)
    {
      allocsize = size;
      mem_to_delete = data;
      for(size_t i = 0; i <  a2.Size(); i++)
        data[i] = a2[i];
      for (size_t i = a2.Size(), j=0; i < size; i++,j++)
        data[i] = a3[j];
    }

    /// if responsible, deletes memory
    NETGEN_INLINE ~Array()
    {
      if(mem_to_delete)
        mt.Free(sizeof(T)*allocsize);
      delete [] mem_to_delete;
    }

    // Only provide this function if T is archivable
    template<typename ARCHIVE>
    auto DoArchive(ARCHIVE& archive)
      -> typename std::enable_if_t<ARCHIVE::template is_archivable<T>, void>
    {
      if(archive.Output())
        archive << size;
      else
        {
          size_t s;
          archive & s;
          SetSize(s);
        }
      archive.Do(data, size);
    }

    /// we tell the compiler that there is no need for deleting the array ..
    NETGEN_INLINE void NothingToDelete () 
    { 
      mem_to_delete = nullptr;
    }

    /// Change logical size. If necessary, do reallocation. Keeps contents.
    NETGEN_INLINE void SetSize(size_t nsize)
    {
      if (nsize > allocsize) ReSize (nsize);
      size = nsize; 
    }

    ///
    NETGEN_INLINE void SetSize0()
    {
      size = 0; 
    }

    /// Change physical size. Keeps logical size. Keeps contents.
    NETGEN_INLINE void SetAllocSize (size_t nallocsize)
    {
      if (nallocsize > allocsize)
        ReSize (nallocsize);
    }

    /// Change physical size. Keeps logical size. Keeps contents.
    NETGEN_INLINE size_t AllocSize () const
    {
      return allocsize;
    }


    /// assigns memory from local heap
    NETGEN_INLINE const Array & Assign (size_t asize, LocalHeap & lh)
    {
      if(mem_to_delete)
        mt.Free(sizeof(T)*allocsize);
      delete [] mem_to_delete;
      size = allocsize = asize;
      data = lh.Alloc<T> (asize);
      mem_to_delete = nullptr;
      return *this;
    }

    /// Add element at end of array. reallocation if necessary.
    /// Returns index of new element.
    NETGEN_INLINE index_type Append (const T & el)
    {
      if (size == allocsize) 
        ReSize (size+1);
      data[size] = el;
      return BASE + size++;
    }

    /// Add element at end of array. reallocation not necessary.
    /// Returns index of new element.
    NETGEN_INLINE index_type AppendHaveMem (const T & el)
    {
      NETGEN_CHECK_RANGE(size, 0, allocsize);
      data[size] = el;
      return BASE + size++;
    }

    
    /// Add element at end of array. reallocation if necessary.
    /// Returns index of new element.
    NETGEN_INLINE index_type Append (T && el)
    {
      if (size == allocsize) 
        ReSize (size+1);
      data[size] = std::move(el);
      return BASE + size++;
    }

    // Add elements of initializer list to end of array. Reallocation if necessary.
    NETGEN_INLINE void Append(std::initializer_list<T> lst)
    {
      if(allocsize < size + lst.size())
        ReSize(size+lst.size());
      for(auto val : lst)
        data[size++] = val;
    }

    /// Add element at end of array. reallocation if necessary.
    NETGEN_INLINE void Insert (size_t pos, const T & el)
    {
      if (size == allocsize) 
        ReSize (size+1);
      for (size_t i = size; i > pos; i--)
        data[i] = data[i-1];
      data[pos] = el;
      size++;
    }
    
    NETGEN_INLINE Array & operator += (const T & el)
    {
      Append (el);
      return *this;
    }


    /// Append array at end of array. reallocation if necessary.
    NETGEN_INLINE void Append (FlatArray<T> source)
    {
      if(size + source.Size() >= allocsize)
        ReSize (size + source.Size() + 1);

      for(size_t i = size, j=0; j<source.Size(); i++, j++)
        data[i] = source[j];

      size += source.Size();
    }



    /// Delete element i. Move last element to position i.
    NETGEN_INLINE void DeleteElement (size_t i)
    {
      NETGEN_CHECK_RANGE(i,BASE,BASE+size);
      data[i-BASE] = std::move(data[size-1]);
      size--;
    }


    /// Delete element i. Move all remaining elements forward
    NETGEN_INLINE void RemoveElement (size_t i)
    {
      NETGEN_CHECK_RANGE(i, BASE, BASE+size);
      for(size_t j = i; j < this->size-1; j++)
	this->data[j] = this->data[j+1];
      this->size--;
    }


    /// Delete last element. 
    NETGEN_INLINE void DeleteLast ()
    {
      NETGEN_CHECK_RANGE(size-1,0,size);
      size--;
    }

    /// Deallocate memory
    NETGEN_INLINE void DeleteAll ()
    {
      if(mem_to_delete)
        mt.Free(sizeof(T)*allocsize);
      delete [] mem_to_delete;
      mem_to_delete = NULL;
      data = 0;
      size = allocsize = 0;
    }

    /// Fill array with val
    NETGEN_INLINE Array & operator= (const T & val)
    {
      FlatArray<T,IndexType>::operator= (val);
      return *this;
    }

    /// array copy
    NETGEN_INLINE Array & operator= (const Array & a2)
    {
      if constexpr (std::is_copy_assignable<T>::value)
        {
          SetSize0 ();
          SetSize (a2.Size());
          for (size_t i = 0; i < size; i++)
            data[i] = a2.data[i];
          return *this;
        }
#ifndef __CUDA_ARCH__      
      else
        throw Exception(std::string("cannot copy Array of type ") + typeid(T).name());
#endif
    }

    
    /// steal array 
    NETGEN_INLINE Array & operator= (Array && a2)
    {
      mt.Swap(sizeof(T)*allocsize, a2.mt, sizeof(T)*a2.allocsize);

      ngcore::Swap (size, a2.size);
      ngcore::Swap (data, a2.data);
      ngcore::Swap (allocsize, a2.allocsize);
      ngcore::Swap (mem_to_delete, a2.mem_to_delete);
      return *this;
    }


    /// array copy
    NETGEN_INLINE Array & operator= (const FlatArray<T> & a2)
    {
      SetSize (a2.Size());
      for (size_t i = 0; i < size; i++)
        data[i] = a2[i];
      return *this;
    }

    /*
    /// fill array with first, first+1, ... 
    Array & operator= (const IntRange & range)
    {
      SetSize (range.Size());
      for (int i = 0; i < size; i++)
        (*this)[i] = range.First()+i;
      return *this;
    }
    */
    template <typename T2>
    Array & operator= (const BaseArrayObject<T2> & a2)
    {
      size_t newsize = a2.Spec().Size();
      SetSize0 ();      
      SetSize (newsize);
      // for (size_t i = 0; i < newsize; i++)
      // (*this)[i] = a2.Spec()[i];
      size_t i = 0;
      for (auto val : a2.Spec())
        (*this)[i++] = val;
      
      return *this;
    }

    template <typename ...ARGS>
    Array & operator= (Tuple<ARGS...> tup)
    {
      SetSize (ArraySize (tup));
      StoreToArray (*this, tup);
      return *this;
    }

    Array & operator= (std::initializer_list<T> list)
    {
      *this = Array<T> (list); 
      return *this;
    }
      

//     template <typename T2>
//     Array & operator= (ParallelValue<T2> val)
//     {
//       FlatArray<T>::operator= (val);
//       return *this;
//     }
//     template <typename T2>
//     Array & operator= (ParallelFunction<T2> val)
//     {
//       FlatArray<T>::operator= (val);
//       return *this;
//     }

    
    NETGEN_INLINE void Swap (Array & b)
    {
      mt.Swap(sizeof(T) * allocsize, b.mt, sizeof(T) * b.allocsize);

      ngcore::Swap (size, b.size);
      ngcore::Swap (data, b.data);
      ngcore::Swap (allocsize, b.allocsize);
      ngcore::Swap (mem_to_delete, b.mem_to_delete);
    }

    NETGEN_INLINE void StartMemoryTracing () const
    {
      mt.Alloc(sizeof(T) * allocsize);
    }

    const MemoryTracer& GetMemoryTracer() const { return mt; }

  private:

    /// resize array, at least to size minsize. copy contents
    NETGEN_INLINE void ReSize (size_t minsize);
    MemoryTracer mt;
  };

  
  /// resize array, at least to size minsize. copy contents
  template <class T, class IndexType> 
  NETGEN_INLINE void Array<T, IndexType> :: ReSize (size_t minsize)
  {
    size_t nsize = 2 * allocsize;
    if (nsize < minsize) nsize = minsize;
    
    T * hdata = data;
    data = new T[nsize];
    mt.Alloc(sizeof(T) * nsize);

    if (hdata)
      {
        size_t mins = (nsize < size) ? nsize : size;
#if defined(__GNUG__) && __GNUC__ < 5 && !defined(__clang__)
        for (size_t i = 0; i < mins; i++) data[i] = std::move(hdata[i]);
#else
        if (std::is_trivially_copyable<T>::value)
          memcpy ((void*)data, hdata, sizeof(T)*mins);
        else
          for (size_t i = 0; i < mins; i++) data[i] = std::move(hdata[i]);
#endif
        if(mem_to_delete)
          mt.Free(sizeof(T) * allocsize);
        delete [] mem_to_delete;
      }

    mem_to_delete = data;
    allocsize = nsize;
  }

  //extern template class Array<int,int>;
  

  /**
     Array with static and dynamic memory management.
     Declares a static array which size is given by the template parameter.
     If the dynamic size fits into the static size, use static memory, 
     otherwise perform dynamic allocation
  */
  template <class T, int S> 
  class ArrayMem : public Array<T>
  {
    T mem[S];    

    using Array<T>::size;
    using Array<T>::allocsize;
    using Array<T>::data;
    using Array<T>::mem_to_delete;
    // using Array<T>::ownmem;

  public:
    /// Generate array of logical and physical size asize
    explicit ArrayMem(size_t asize = 0)    
      : Array<T> (S, mem)
    {
      size = asize;
      if (asize > S)
        {
          data = new T[asize];
          allocsize = size;
          mem_to_delete = data;
        }
    }

    /// copies from Array a2
    explicit ArrayMem(const Array<T> & a2)
      : Array<T> (S, (T*)mem)
    {
      Array<T>::operator= (a2);
    }

    /// copies from ArrayMem a2
    explicit ArrayMem(const ArrayMem & a2)
      : Array<T> (S, (T*)mem)
    {
      Array<T>::operator= (a2);
    }
  
    ArrayMem(ArrayMem && a2)
      : Array<T> (a2.Size(), (T*)mem)
    {
      if (a2.mem_to_delete)
        {
          mem_to_delete = a2.mem_to_delete;
          data = a2.data;
          allocsize = a2.allocsize;
          a2.mem_to_delete = nullptr;
          a2.data = nullptr;
          a2.size = 0;
        }
      else
        {
          allocsize = S;
          for (auto i : ngcore::Range(size))
            mem[i] = a2.mem[i];
        }
    }
    
    ArrayMem (std::initializer_list<T> list)
      : ArrayMem (list.size())
    {
      size_t cnt = 0;
      for (auto val : list)
        data[cnt++] = val;
    }
  
    template <typename T2>
    ArrayMem (const BaseArrayObject<T2> & a2)
      : ArrayMem (a2.Size())
    {
      for (size_t i : ngcore::Range(size))
        data[i] = a2[i];
    }

    
    ArrayMem & operator= (const T & val)
    {
      FlatArray<T>::operator= (val);
      return *this;
    }

    ArrayMem & operator= (ArrayMem && a2)
    {
      ngcore::Swap (mem_to_delete, a2.mem_to_delete);
      ngcore::Swap (allocsize, a2.allocsize);
      ngcore::Swap (size, a2.size);

      if (mem_to_delete==nullptr)
      {
        for (auto i : ngcore::Range(size))
          mem[i] = std::move(a2.mem[i]);
        data = mem;
      }
      else
        ngcore::Swap (data, a2.data);

      return *this;
    }


    /// array copy
    ArrayMem & operator= (const FlatArray<T> & a2)
    {
      this->SetSize (a2.Size());
      for (size_t i = 0; i < size; i++)
        (*this)[i] = a2[i];
      return *this;
    }


    template <typename T2>
    ArrayMem & operator= (const BaseArrayObject<T2> & a2)
    {
      this->SetSize (a2.Spec().Size());

      size_t i = 0;
      for (auto val : a2.Spec())
        (*this)[i++] = val;
      
      return *this;
    }

  };





  template <typename ... ARGS>
  size_t ArraySize (Tuple<ARGS...> /* tup */)  
  { return 0;}
  
  template <typename ... ARGS>
  size_t ArraySize (Tuple<int,ARGS...> tup) 
  { return 1+ArraySize(tup.Tail()); }
  
  template <typename ... ARGS>
  size_t ArraySize (Tuple<IntRange,ARGS...> tup) 
  { return tup.Head().Size()+ArraySize(tup.Tail()); }

  
  template <typename T, typename ... ARGS>
  void StoreToArray (FlatArray<T> /* a */, Tuple<ARGS...> /* tup */) { ; }
  
  template <typename T, typename ... ARGS>
  void StoreToArray (FlatArray<T> a, Tuple<int,ARGS...> tup)
  {
    a[0] = tup.Head();
    StoreToArray (a.Range(1, a.Size()), tup.Tail());
  }
  
  template <typename T, typename ... ARGS>
  void StoreToArray (FlatArray<T> a, Tuple<IntRange,ARGS...> tup)
  {
    IntRange r = tup.Head();
    a.Range(0,r.Size()) = r;
    StoreToArray (a.Range(r.Size(), a.Size()), tup.Tail());
  }

  /*
  template <typename T> template <typename ...ARGS>
  NETGEN_INLINE Array<T> & Array<T> :: operator= (Tuple<ARGS...> tup)
  {
    SetSize (ArraySize (tup));
    StoreToArray (*this, tup);
  }
  */

  /*
  /// append integers to array
  inline Array<int> & operator+= (Array<int> & array, const IntRange & range)
  {
    int oldsize = array.Size();
    int s = range.Next() - range.First();
    
    array.SetSize (oldsize+s);

    for (int i = 0; i < s; i++)
      array[oldsize+i] = range.First()+i;

    return array;
  }
  */
  

  /*
  template <typename T, typename T2>
  inline Array<T> & operator+= (Array<T> & array, const BaseArrayObject<T2> & a2)
  {
    size_t oldsize = array.Size();
    size_t s = a2.Spec().Size();
    
    array.SetSize (oldsize+s);

    for (size_t i = 0; i < s; i++)
      array[oldsize+i] = a2.Spec()[i];

    return array;
  }
  */
  
  template <typename T, typename T2>
  inline Array<T> & operator+= (Array<T> & array, const BaseArrayObject<T2> & a2)
  {
    auto oldsize = array.Size();
    auto s = a2.Spec().Size();

    array.SetSize (oldsize+s);

    for (auto val : a2.Spec())
      array[oldsize++] = val;
    
    return array;
  }

  template <typename T, typename T2>
  inline Array<T> operator+= (Array<T> && array, const BaseArrayObject<T2> & a2)
  {
    array += a2;
    return std::move(array);
  }


  /// bubble sort array
  template <class T>
  inline void BubbleSort (FlatArray<T> data)
  {
    T hv;
    for (size_t i = 0; i < data.Size(); i++)
      for (size_t j = i+1; j < data.Size(); j++)
        if (data[i] > data[j])
          {
            hv = data[i];
            data[i] = data[j];
            data[j] = hv;
          }
  }

  /// bubble sort array
  template <class T, class S>
  inline void BubbleSort (FlatArray<T> data, FlatArray<S> index)
  {
    for (size_t i = 0; i < data.Size(); i++)
      for (size_t j = i+1; j < data.Size(); j++)
	if (data[i] > data[j])
	  {
	    T hv = data[i];
	    data[i] = data[j];
	    data[j] = hv;

	    S hvs = index[i];
	    index[i] = index[j];
	    index[j] = hvs;
	  }
  }




  template <class T, typename TLESS>
  void QuickSort (FlatArray<T> data, TLESS less)
  {
    if (data.Size() <= 1) return;

    ptrdiff_t i = 0;
    ptrdiff_t j = data.Size()-1;

    T midval = data[ (i+j)/2 ];
  
    do
      {
        while (less (data[i], midval)) i++;
        while (less (midval, data[j])) j--;

        if (i <= j)
          {
	    Swap (data[i], data[j]);
            i++; j--;
          }
      }
    while (i <= j);

    QuickSort (data.Range (0, j+1), less);
    QuickSort (data.Range (i, data.Size()), less);
  }

  template <typename T>
  NETGEN_INLINE bool DefaultLess (const T & a, const T & b)
  {
    return a < b;
  }

  template <typename T>
  class DefaultLessCl
  {
  public:
    bool operator() (const T & a, const T & b) const
    {
      return a < b;
    }
  };



  template <class T>
  NETGEN_INLINE void QuickSort (FlatArray<T> data)
  {
    QuickSort (data, DefaultLessCl<T>());
  }



  template <class T, typename TLESS>
  void QuickSortI (FlatArray<T> data, FlatArray<int> index, TLESS less)
  {
    if (index.Size() <= 1) return;

    ptrdiff_t i = 0;
    ptrdiff_t j = index.Size()-1;

    int midval = index[ (i+j)/2 ];
  
    do
      {
        while (less (data[index[i]],data[midval])  ) i++;
        while (less (data[midval],  data[index[j]])) j--;

        if (i <= j)
          {
	    Swap (index[i], index[j]);
            i++; j--;
          }
      }
    while (i <= j);

    QuickSortI (data, index.Range (0, j+1), less);
    QuickSortI (data, index.Range (i, index.Size()), less);
  }


  template <class T>
  NETGEN_INLINE void QuickSortI (FlatArray<T> data, FlatArray<int> index)
  {
    QuickSortI (data, index, DefaultLessCl<T>());
  }





  template <typename T>
  NETGEN_INLINE T xxxRemoveRef (const T & x)
  {
    return x;
  }

  template <class TA1, class TA2> 
  class SumArray : public BaseArrayObject<SumArray<TA1,TA2>>
  {
    const TA1 & a1;
    const TA2 & a2;
  public:
    SumArray (const TA1 & aa1, const TA2 & aa2) : a1(aa1), a2(aa2) { ; }
    size_t Size() const { return a1.Size()+a2.Size(); }
    auto operator[] (size_t i) const -> decltype (xxxRemoveRef (a1[0])) 
    {
      return (i < a1.Size()) ? a1[i] : a2[i-a1.Size()];
    }
  };

  template <class TA1, class TA2> 
  SumArray<TA1,TA2> operator+ (const BaseArrayObject<TA1> & a1,
                               const BaseArrayObject<TA2> & a2)
  {
    return SumArray<TA1,TA2> (a1.Spec(), a2.Spec());
  }
                               

  struct HTAHelp { };
  
  // head-tail array
  template <size_t S, typename T>
  class HTArray
  {
    HTArray<S-1,T> tail;
    T head;
  public:
    constexpr HTArray () = default;
    constexpr HTArray (const HTArray &) = default;
    template <typename T2>
    constexpr HTArray (const HTArray<S,T2> & a2) : tail(a2.Tail()), head(a2.Head()) { ; }

    constexpr HTArray (T v) : tail(v), head(v) { } // all the same
    
    template <class... T2,
              std::enable_if_t<S==1+sizeof...(T2),bool> = true>
    constexpr HTArray (const T &v, T2... rest)
      : tail{HTAHelp(), v,rest...}, head(std::get<S-2>(std::tuple(rest...))) { }

    template <class... T2>
    constexpr HTArray (HTAHelp h, const T &v, T2... rest)
      : tail{h, v,rest...}, head(std::get<S-2>(std::tuple(rest...))) { }

    
    HTArray & operator= (const HTArray &) = default;

    T * Ptr () { return tail.Ptr(); }
    T & operator[] (size_t i) { return Ptr()[i]; }

    const T * Ptr () const { return tail.Ptr(); }
    const T & operator[] (size_t i) const { return Ptr()[i]; }
    template <int NR>
    T & Elem() { return (NR==S-1) ? head : tail.template Elem<NR>(); }

    auto Tail() const { return tail; }
    auto Head() const { return head; }
  };

  template <typename T>
  class HTArray<1,T>
  {
    T head;
  public:
    constexpr HTArray () = default;
    constexpr HTArray (const HTArray &) = default;
    template <typename T2>
    constexpr HTArray (const HTArray<1,T2> & a2) : head(a2.Head()) { ; }
    constexpr HTArray (T v) : head(v) { } // all the same
    template <class... T2>    
    constexpr HTArray (HTAHelp h, const T &v, T2... rest)
      : head(v) { } 

    
    HTArray & operator= (const HTArray &) = default;

    T * Ptr () { return &head; }
    T & operator[] (size_t i) { return Ptr()[i]; }

    const T * Ptr () const { return &head; }
    const T & operator[] (size_t i) const { return Ptr()[i]; }
    template <int NR>    
    T & Elem()
    {
      // assert(NR==0, "HTArray index error");
      return head;
    }

    auto Head() const { return head; }
  };

  template <typename T>
  class HTArray<0,T>
  {
    // T head; // dummy variable
  public:
    HTArray () = default;
    HTArray (const HTArray &) = default;
    template <typename T2>
    HTArray (const HTArray<0,T2> & a2) { ; }
    constexpr HTArray (T v) { } // all the same    
    HTArray & operator= (const HTArray &) = default;

    /*
    T * Ptr () { return &head; }
    T & operator[] (size_t i) { return Ptr()[i]; }

    const T * Ptr () const { return &head; }
    const T & operator[] (size_t i) const { return Ptr()[i]; }
    template <int NR>        
    T & Elem()
    {
      // assert(false, "HTArray index error");
      return head;
    }
    */
    // T * Ptr () { return (T*)(void*)&head; }
    T * Ptr () { return (T*)(void*)this; }
    T & operator[] (size_t i) { return Ptr()[i]; }
    // const T * Ptr () const { return (const T*)(const void*)&head; }
    const T * Ptr () const { return (const T*)(const void*)this; }
    const T & operator[] (size_t i) const { return Ptr()[i]; }
    template <int NR>        
    T & Elem()
    {
      throw Exception("illegal HTArray<0>::Elem<0>");
    }

  };

  template<size_t S, typename T>
  const T * operator+ (const HTArray<S,T> & ar, size_t i)
  {
    return ar.Ptr()+i;
  }
  template<size_t S, typename T>
  T * operator+ (HTArray<S,T> & ar, size_t i)
  {
    return ar.Ptr()+i;
  }










  template <typename TIA, typename TIB>
  class IteratorPair
  {
    TIA a;
    TIB b;
  public:
    IteratorPair (TIA _a, TIB _b) : a(_a), b(_b) { ; }
    
    IteratorPair & operator++() { ++a; ++b; return *this; }
    bool operator!= (const IteratorPair & it2)  { return a != it2.a; }
    
    auto operator*()
    {
      // return pair(*a,*b);
      return std::pair<decltype(*a), decltype(*b)> (*a, *b); // keep reference
    }
  };
  
  
  template <typename TA, typename TB>
  class Zip
  {
    const TA & a;
    const TB & b;
  public:
    Zip(const TA & _a, const TB & _b) : a(_a), b(_b) { ; }
    auto begin() const { return IteratorPair(a.begin(), b.begin()); }
    auto end() const { return IteratorPair(a.end(), b.end()); }
  };

  template <typename T>
  inline size_t size (const BaseArrayObject<T> & ao) { return ao.Size(); }
  
  template <typename TA>
  class Enumerate
  {
    IntRange r;
    const TA & a;
  public:
    Enumerate(const TA & _a) : r(size(_a)), a(_a) { ; }
    auto begin() const { return IteratorPair(r.begin(), a.begin()); }
    auto end() const { return IteratorPair(r.end(), a.end()); }
  };
  
  


  
}


#endif // NETGEN_CORE_ARRAY_HPP


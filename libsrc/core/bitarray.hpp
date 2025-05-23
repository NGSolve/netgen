#ifndef NETGEN_CORE_BITARRAY
#define NETGEN_CORE_BITARRAY

/**************************************************************************/
/* File:   bitarray.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <climits>
#include <cstring>
#include <ostream>

#include "array.hpp"
#include "localheap.hpp"
#include "ngcore_api.hpp"
#include "utils.hpp"

namespace ngcore
{

/**
   A compressed array of bools.

   Provides bit-operations and whole array operations.
*/
class BitArray
{
protected:
  /// number of bits
  size_t size;

  /// the data
  unsigned char * data;
  ///
  bool owns_data = true;
public:
  /// empty array
  BitArray ()
    : size(0), data(nullptr) { ; }
  /// array of asize bits
  NGCORE_API BitArray (size_t asize);
  /// array of asize bits
  NGCORE_API BitArray (size_t asize, LocalHeap & lh);
  ///
  NGCORE_API BitArray (const BitArray & ba2);
  BitArray (BitArray && ba2)
    : size(ba2.size), data(ba2.data), owns_data(ba2.owns_data)
  {
    ba2.owns_data = false;
    ba2.data = nullptr;
    mt = std::move(ba2.mt);
  }

  template <typename T>
  NETGEN_INLINE BitArray (std::initializer_list<T> list)
    : BitArray (list.size())
  {
    Clear();
    int cnt = 0;
    for (auto i = list.begin(); i < list.end(); i++, cnt++)
      if (*i) SetBit(cnt);
    StartMemoryTracing();
  }

  /// delete data
  ~BitArray ()
  {
    if (owns_data)
    {
      delete [] data;
      mt.Free(GetMemoryUsage());
    }
  }

  /// Set size, loose values
  NGCORE_API void SetSize (size_t asize);

  /// the size
  size_t Size () const { return size; }

  /// set all bits
  NGCORE_API BitArray & Set () throw();

  /// clear all bits
  NGCORE_API BitArray & Clear () throw();

  /// set bit i
  [[deprecated("Use either SetBit() or SetBitAtomic()")]]
  void Set (size_t i) { SetBitAtomic(i); }

  /// set bit i ( not thread safe )
  void SetBit (size_t i)
  {
    NETGEN_CHECK_RANGE(i, 0, size);
    data[Addr(i)] |= Mask(i);
  }

  /// set bit i ( thread safe )
  void SetBitAtomic (size_t i)
  {
    NETGEN_CHECK_RANGE(i, 0, size);
    unsigned char * p = data+Addr(i);
    unsigned char mask = Mask(i);

    AsAtomic(*p) |= mask;
  }

  /// clear bit i
  void Clear (size_t i)
  {
    NETGEN_CHECK_RANGE(i, 0, size);
    data[Addr(i)] &= ~Mask(i);
  }

  /// check bit i
  bool Test (size_t i) const
  {
    NETGEN_CHECK_RANGE(i, 0, size);
    return (data[Addr(i)] & Mask(i)) ? true : false;
  }

  /// set all bits to b
  BitArray & operator= (bool b)
  {
    if (b) Set();
    else   Clear();
    return *this;
  }

  /// check bit i
  bool operator[] (size_t i) const
  {
    NETGEN_CHECK_RANGE(i, 0, size);
    return Test(i);
  }

  NGCORE_API bool operator==(const BitArray& other) const;

  /// invert all bits
  NGCORE_API BitArray & Invert ();

  /// logical AND with ba2
  NGCORE_API BitArray & And (const BitArray & ba2);

  /// logical OR with ba2
  NGCORE_API BitArray & Or (const BitArray & ba2);

  /// copy from ba2
  NGCORE_API BitArray & operator= (const BitArray & ba2);

  NGCORE_API size_t NumSet () const;

  NGCORE_API void DoArchive(class Archive& archive);
  
  NGCORE_API auto * Data() const { return data; }

  const size_t GetMemoryUsage() const { return owns_data ? (size+CHAR_BIT-1)/CHAR_BIT : 0; }
  const MemoryTracer& GetMemoryTracer() const { return mt; }
  void StartMemoryTracing() const
  {
    mt.Alloc(GetMemoryUsage());
  }

private:
  ///
  unsigned char Mask (size_t i) const
  { return char(1) << (i % CHAR_BIT); }

  ///
  size_t Addr (size_t i) const
  { return (i / CHAR_BIT); }

  MemoryTracer mt;
};


  inline BitArray & operator|= (BitArray & me, const BitArray & you)
  {
    me.Or(you);
    return me;
  }

  inline BitArray & operator&= (BitArray & me, const BitArray & you)
  {
    me.And(you);
    return me;
  }

  inline BitArray operator| (const BitArray & a, const BitArray & b)
  {
    BitArray res = a;
    res |= b;
    return res;
  }

  inline BitArray operator& (const BitArray & a, const BitArray & b)
  {
    BitArray res = a;
    res &= b;
    return res;
  }

  inline BitArray operator~ (const BitArray & a)
  {
    BitArray res = a;
    res.Invert();
    return res;
  }

  NGCORE_API std::ostream & operator<<(std::ostream & s, const BitArray & ba);



  template <typename IndexType>
  class TBitArray : public BitArray
  {
  public:
    using BitArray::BitArray;

    void SetBit (IndexType i) { BitArray::SetBit(i-IndexBASE<IndexType>()); }
    void Clear () { BitArray::Clear(); }
    void Clear (IndexType i) { BitArray::Clear(i-IndexBASE<IndexType>()); }
    void SetBitAtomic (IndexType i) { BitArray::SetBitAtomic(i-IndexBASE<IndexType>()); }
    bool Test (IndexType i) const { return BitArray::Test(i-IndexBASE<IndexType>()); }
    
    bool operator[] (IndexType i) const { return Test(i); } 
    T_Range<IndexType> Range() const { return { IndexBASE<IndexType>(), IndexBASE<IndexType>()+Size() }; }
    NGCORE_API TBitArray & Or (const TBitArray & ba2)
    {
      BitArray::Or(ba2);
      return *this;
    }

  };

} // namespace ngcore

  
#endif // NETGEN_CORE_BITARRAY

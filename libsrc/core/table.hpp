#ifndef NETGEN_CORE_TABLE_HPP
#define NETGEN_CORE_TABLE_HPP

/**************************************************************************/
/* File:   table.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 2000                                                  */
/**************************************************************************/

#include <atomic>
#include <iostream>

#include "array.hpp"
#include "bitarray.hpp"
#include "taskmanager.hpp"
#include "ngcore_api.hpp"

namespace ngcore
{


  template <class T, class IndexType = size_t>
class FlatTable
{
protected:
  static constexpr IndexType BASE = IndexBASE<IndexType>();  
  /// number of rows
  size_t size;
  /// pointer to first in row
  size_t * index;
  /// array of data
  T * data;

public:
  FlatTable() = delete;

  NETGEN_INLINE FlatTable(size_t as, size_t * aindex, T * adata)
    : size(as), index(aindex), data(adata) { ; }

  /// Size of table
  NETGEN_INLINE size_t Size() const { return size; }

  /// Access entry
  NETGEN_INLINE const FlatArray<T> operator[] (IndexType i) const
  {
    i = i-BASE;
    return FlatArray<T> (index[i+1]-index[i], data+index[i]);
  }

  NETGEN_INLINE T * Data() const { return data; }

  NETGEN_INLINE FlatArray<T> AsArray() const
  {
    return FlatArray<T> (index[size]-index[0], data+index[0]);
  }

  NETGEN_INLINE FlatArray<size_t> IndexArray() const
  {
    return FlatArray<size_t, IndexType> (size+1, index);
  }

  /// takes range starting from position start of end-start elements
  NETGEN_INLINE FlatTable<T> Range (size_t start, size_t end) const
  {
    return FlatTable<T> (end-start, index+start-BASE, data);
  }

  /// takes range starting from position start of end-start elements
  NETGEN_INLINE FlatTable<T> Range (T_Range<size_t> range) const
  {
    return FlatTable<T> (range.Size(), index+range.First()-BASE, data);
  }

  NETGEN_INLINE T_Range<IndexType> Range () const
  {
    return T_Range<IndexType> (BASE, size+BASE);
  }
  
  class Iterator
  {
    const FlatTable & tab;
    size_t row;
  public:
    Iterator (const FlatTable & _tab, size_t _row) : tab(_tab), row(_row) { ; }
    Iterator & operator++ () { ++row; return *this; }
    FlatArray<T> operator* () const { return tab[row]; }
    bool operator!= (const Iterator & it2) { return row != it2.row; }
  };

  Iterator begin() const { return Iterator(*this, BASE); }
  Iterator end() const { return Iterator(*this, BASE+size); }
};

  NGCORE_API extern size_t * TablePrefixSum32 (FlatArray<unsigned int> entrysize);
  NGCORE_API extern size_t * TablePrefixSum64 (FlatArray<size_t> entrysize);


  NETGEN_INLINE size_t * TablePrefixSum (FlatArray<unsigned int> entrysize)
  { return TablePrefixSum32 (entrysize); }
  NETGEN_INLINE size_t * TablePrefixSum (FlatArray<int> entrysize)
  { return TablePrefixSum32 (FlatArray<unsigned> (entrysize.Size(), (unsigned int*)(int*)(entrysize.Addr(0)))); }
  NETGEN_INLINE size_t * TablePrefixSum (FlatArray<std::atomic<int>> entrysize)
  { return TablePrefixSum32 (FlatArray<unsigned> (entrysize.Size(), (unsigned int*)(std::atomic<int>*)entrysize.Addr(0))); }
  NETGEN_INLINE size_t * TablePrefixSum (FlatArray<size_t> entrysize)
  { return TablePrefixSum64 (entrysize); }


/**
    A compact Table container.
    A table contains size entries of variable size.
    The entry sizes must be known at construction.
*/
  template <class T, class IndexType = size_t>
  class Table : public FlatTable<T, IndexType>
{
protected:

  using FlatTable<T,IndexType>::size;
  using FlatTable<T,IndexType>::index;
  using FlatTable<T,IndexType>::data;

public:
  ///
  NETGEN_INLINE Table () : FlatTable<T,IndexType> (0,nullptr,nullptr) { ; }
  /// Construct table of uniform entrysize
  NETGEN_INLINE Table (size_t asize, size_t entrysize)
    : FlatTable<T,IndexType>( asize, new size_t[asize+1], new T[asize*entrysize] )
  {
    for (size_t i : IntRange(size+1))
      index[i] = i*entrysize;
  }

  /// Construct table of variable entrysize
  template <typename TI>
  NETGEN_INLINE Table (FlatArray<TI,IndexType> entrysize)
    : FlatTable<T,IndexType> (0, nullptr, nullptr)
  {
    size  = entrysize.Size();
    index = TablePrefixSum (FlatArray<TI> (entrysize.Size(), entrysize.Data()));
    size_t cnt = index[size];
    data = new T[cnt];
  }

  explicit NETGEN_INLINE Table (const Table & tab2)
    : FlatTable<T,IndexType>(0, nullptr, nullptr)
  {
    size = tab2.Size();

    index = new size_t[size+1];
    for (size_t i = 0; i <= size; i++)
      index[i] = tab2.index[i];

    size_t cnt = index[size];
    data = new T[cnt];
    for (size_t i = 0; i < cnt; i++)
      data[i] = tab2.data[i];
  }

  NETGEN_INLINE Table (Table && tab2)
    : FlatTable<T,IndexType>(0, nullptr, nullptr)
  {
    Swap (size, tab2.size);
    Swap (index, tab2.index);
    Swap (data, tab2.data);
  }

  NETGEN_INLINE Table & operator= (Table && tab2)
  {
    Swap (size, tab2.size);
    Swap (index, tab2.index);
    Swap (data, tab2.data);
    return *this;
  }



  /// Delete data
  NETGEN_INLINE ~Table ()
  {
    delete [] data;
    delete [] index;
  }

  /// Size of table
  using FlatTable<T,IndexType>::Size;

  /// number of elements in all rows
  NETGEN_INLINE size_t NElements() const { return index[size]; }

  using FlatTable<T,IndexType>::operator[];
};


/// Print table
  template <class T, typename IndexType>
  inline ostream & operator<< (ostream & s, const Table<T,IndexType> & table)
{
  for (auto i : table.Range())
    {
      s << i << ":";
      for (auto el : table[i])
	s << " " << el;
      s << "\n";
    }
  s << std::flush;
  return s;
}




  template <class T, typename IndexType=size_t>
  class TableCreator
  {
  protected:
    int mode;    // 1 .. cnt, 2 .. cnt entries, 3 .. fill table
    std::atomic<size_t> nd;
    Array<std::atomic<int>,IndexType> cnt;
    Table<T,IndexType> table;
  public:
    TableCreator()
    { nd = 0; mode = 1; }
    TableCreator (size_t acnt)
    { nd = acnt; SetMode(2); }

    Table<T,IndexType> MoveTable()
    {
      return std::move(table);
    }

    bool Done () { return mode > 3; }
    void operator++(int) { SetMode (mode+1); }

    int GetMode () const { return mode; }
    void SetMode (int amode)
    {
      mode = amode;
      if (mode == 2)
	{
	  // cnt.SetSize(nd);  // atomic has no copy
          cnt = Array<std::atomic<int>,IndexType> (nd);
          for (auto & ci : cnt) ci.store (0, std::memory_order_relaxed);
	}
      if (mode == 3)
	{
          table = Table<T,IndexType> (cnt);
          // for (auto & ci : cnt) ci = 0;
          for (auto & ci : cnt) ci.store (0, std::memory_order_relaxed);
          // cnt = 0;
	}
    }

    void SetSize (size_t _nd)
    {
      if (mode == 1)
        nd = _nd;
      else
        {
          if (nd != _nd)
            throw Exception ("cannot change size of table-creator");
        }
    }

    void Add (IndexType blocknr, const T & data)
    {
      switch (mode)
	{
	case 1:
          {
            size_t oldval = nd;
            while (blocknr+1>nd) {
              nd.compare_exchange_weak (oldval, blocknr+1);
              oldval = nd;
            }
            break;
          }
	case 2:
	  cnt[blocknr]++;
	  break;
	case 3:
          int ci = cnt[blocknr]++;
          table[blocknr][ci] = data;
	  break;
	}
    }


    void Add (IndexType blocknr, IntRange range)
    {
      switch (mode)
	{
	case 1:
          {
            size_t oldval = nd;
            while (blocknr+1>nd) {
              nd.compare_exchange_weak (oldval, blocknr+1);
              oldval = nd;
            }
            break;
          }
	case 2:
	  cnt[blocknr] += range.Size();
	  break;
	case 3:
          size_t ci = ( cnt[blocknr] += range.Size() ) - range.Size();
	  for (size_t j = 0; j < range.Size(); j++)
            table[blocknr][ci+j] = range.First()+j;
	  break;
	}
    }

    void Add (IndexType blocknr, const FlatArray<int> & dofs)
    {
      switch (mode)
	{
	case 1:
          {
            size_t oldval = nd;
            while (blocknr+1>nd) {
              nd.compare_exchange_weak (oldval, blocknr+1);
              oldval = nd;
            }
            break;
          }
	case 2:
	  cnt[blocknr] += dofs.Size();
	  break;
	case 3:
          size_t ci = ( cnt[blocknr] += dofs.Size() ) - dofs.Size();
	  for (size_t j = 0; j < dofs.Size(); j++)
            table[blocknr][ci+j] = dofs[j];
	  break;
	}
    }
  };

  class NGCORE_API FilteredTableCreator : public TableCreator<int>
  {
  protected:
    const BitArray* takedofs;
  public:
    FilteredTableCreator(const BitArray* atakedofs)
      : TableCreator<int>(), takedofs(atakedofs) { };
    FilteredTableCreator(int acnt, const BitArray* atakedofs)
      : TableCreator<int>(acnt),takedofs(atakedofs) { };
    void Add (size_t blocknr, int data);
    void Add (size_t blocknr, IntRange range);
    void Add (size_t blocknr, FlatArray<int> dofs);
  };


  /// Base class to generic DynamicTable.
  class BaseDynamicTable
  {
  protected:

    ///
    struct linestruct
    {
      ///
      int size;
      ///
      int maxsize;
      ///
      void * col;
    };

    ///
    Array<linestruct> data;
    ///
    char * oneblock;

  public:
    ///
    NGCORE_API BaseDynamicTable (int size);
    ///
    NGCORE_API BaseDynamicTable (const Array<int> & entrysizes, int elemsize);
    ///
    NGCORE_API ~BaseDynamicTable ();

    /// Changes Size of table to size, deletes data
    NGCORE_API void SetSize (int size);
    ///
    NGCORE_API void IncSize (int i, int elsize);

    NGCORE_API void DecSize (int i);
  };



  /**
      A dynamic table class.

      A DynamicTable contains entries of variable size. Entry sizes can
      be increased dynamically.
  */
  template <class T>
  class DynamicTable : public BaseDynamicTable
  {
  public:
    /// Creates table of size size
    DynamicTable (int size = 0)
      : BaseDynamicTable (size) { ; }

    /// Creates table with a priori fixed entry sizes.
    DynamicTable (const Array<int> & entrysizes)
      : BaseDynamicTable (entrysizes, sizeof(T)) { ; }

    /// Inserts element acont into row i. Does not test if already used.
    void Add (int i, const T & acont)
    {
      if (data[i].size == data[i].maxsize)
        IncSize (i, sizeof (T));
      else
        data[i].size++;
      static_cast<T*> (data[i].col) [data[i].size-1] = acont;
    }

    /// Inserts element acont into row i, iff not yet exists.
    void AddUnique (int i, const T & cont)
    {
      int es = EntrySize (i);
      int * line = const_cast<int*> (GetLine (i));
      for (int j = 0; j < es; j++)
        if (line[j] == cont)
          return;
      Add (i, cont);
    }


    /// Inserts element acont into row i. Does not test if already used.
    void AddEmpty (int i)
    {
      IncSize (i, sizeof (T));
    }

    /** Set the nr-th element in the i-th row to acont.
        Does not check for overflow. */
    void Set (int i, int nr, const T & acont)
    { static_cast<T*> (data[i].col)[nr] = acont; }


    /** Returns the nr-th element in the i-th row.
      Does not check for overflow. */
    const T & Get (int i, int nr) const
    { return static_cast<T*> (data[i].col)[nr]; }


    /** Returns pointer to the first element in row i. */
    const T * GetLine (int i) const
    { return static_cast<T*> (data[i].col); }


    /// Returns size of the table.
    int Size () const
    { return data.Size(); }

    /// Returns size of the i-th row.
    int EntrySize (int i) const
    { return data[i].size; }

    ///
    void DecEntrySize (int i)
    { DecSize(i); }

    /// Access entry i
    FlatArray<T> operator[] (int i)
    { return FlatArray<T> (data[i].size, static_cast<T*> (data[i].col)); }

    /*
    typedef const FlatArray<T> ConstFlatArray;
    /// Access entry i
    ConstFlatArray operator[] (int i) const
    { return FlatArray<T> (data[i].size, static_cast<T*> (data[i].col)); }
    */
    FlatArray<T> operator[] (int i) const
    { return FlatArray<T> (data[i].size, static_cast<T*> (data[i].col)); }
  };




  /// Print table
  template <class T>
  inline ostream & operator<< (ostream & s, const DynamicTable<T> & table)
  {
    for (int i = 0; i < table.Size(); i++)
      {
        s << i << ":";
        for (int j = 0; j < table[i].Size(); j++)
          s << " " << table[i][j];
        s << "\n";
      }
    s << std::flush;
    return s;
  }

  typedef DynamicTable<int> IntTable;

} // namespace ngcore

#endif // NETGEN_CORE_TABLE_HPP

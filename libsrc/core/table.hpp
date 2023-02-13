#ifndef NETGEN_CORE_TABLE_HPP
#define NETGEN_CORE_TABLE_HPP

/**************************************************************************/
/* File:   table.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 2000                                                  */
/**************************************************************************/

#include <atomic>
#include <iostream>
#include <optional>

#include "array.hpp"
#include "bitarray.hpp"
#include "memtracer.hpp"
#include "ngcore_api.hpp"
#include "profiler.hpp"

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
    FlatTable (const FlatTable &) = default;
    
    NETGEN_INLINE FlatTable(size_t as, size_t * aindex, T * adata)
      : size(as), index(aindex), data(adata) { ; }

    /// Size of table
    NETGEN_INLINE size_t Size() const { return size; }

    /// Access entry
    NETGEN_INLINE const FlatArray<T> operator[] (IndexType i) const
    {
      return FlatArray<T> (index[i-BASE+1]-index[i-BASE], data+index[i-BASE]);
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

    explicit NETGEN_INLINE Table (const FlatTable<T,IndexType> & tab2)
      : FlatTable<T,IndexType>(0, nullptr, nullptr)
    {
      size = tab2.Size();
      if (size == 0) return;
      
      index = new size_t[size+1];
      this->IndexArray() = tab2.IndexArray();
      // for (size_t i = 0; i <= size; i++)
      // index[i] = tab2.index[i];

      size_t cnt = index[size];
      data = new T[cnt];
      this->AsArray() = tab2.AsArray();
      /*
      for (size_t i = 0; i < cnt; i++)
        data[i] = tab2.data[i];
      */
    }
    
    explicit NETGEN_INLINE Table (const Table & tab2)
      : FlatTable<T,IndexType>(0, nullptr, nullptr)
    {
      size = tab2.Size();
      if (size == 0) return;
      
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
      tab2.mt.Free(tab2.GetMemUsage());
      Swap (size, tab2.size);
      Swap (index, tab2.index);
      Swap (data, tab2.data);
    }

    template<typename ARCHIVE>
    auto DoArchive(ARCHIVE& ar)
    {
      ar & size;
      if(size == 0)
        return;
      if(ar.Input())
        {
          index = new IndexType[size+1];
          mt.Alloc(sizeof(IndexType) * (size+1));
        }
      ar.Do(index, size+1);
      if(ar.Input())
        {
          data = new T[index[size]];
          mt.Alloc(sizeof(T) * index[size]);
        }
      ar.Do(data, index[size]);
    }

    NETGEN_INLINE Table & operator= (Table && tab2)
    {
      mt.Swap(GetMemUsage(), tab2.mt, tab2.GetMemUsage());
      Swap (size, tab2.size);
      Swap (index, tab2.index);
      Swap (data, tab2.data);
      return *this;
    }



    /// Delete data
    NETGEN_INLINE ~Table ()
    {
      mt.Free(GetMemUsage());
      delete [] data;
      delete [] index;
    }

    /// Size of table
    using FlatTable<T,IndexType>::Size;

    /// number of elements in all rows
    NETGEN_INLINE size_t NElements() const { return index[size]; }

    using FlatTable<T,IndexType>::operator[];

    NETGEN_INLINE void StartMemoryTracing (int /* mem_id */)
    {
      mt.Alloc(GetMemUsage());
    }
    const MemoryTracer& GetMemoryTracer() const { return mt; }

  private:
    size_t GetMemUsage() const { return size == 0 ? 0 : sizeof(T)*index[size] + sizeof(IndexType) * size+1; }
    MemoryTracer mt;
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

  template <typename TEntry, typename TIndex, typename TRange, typename TFunc>
  Table<TEntry, TIndex> CreateTable( const TRange & range, const TFunc & func, std::optional< size_t > cnt )
  {
      static Timer timer("CreateTable");
      RegionTimer rt(timer);
      std::unique_ptr<TableCreator<TEntry, TIndex>> pcreator;

      if(cnt)
          pcreator = std::make_unique<TableCreator<TEntry, TIndex>>(*cnt);
      else
          pcreator = std::make_unique<TableCreator<TEntry, TIndex>>();

      auto & creator = *pcreator;

      for ( ; !creator.Done(); creator++)
        ParallelForRange
          (range, [&] (auto myrange)
           {
             for (auto i : myrange)
               func(creator, i);
           }, TasksPerThread(4)
          );

    return creator.MoveTable();
  }

  template <typename TEntry, typename TIndex, typename TRange, typename TFunc>
  Table<TEntry, TIndex> CreateSortedTable( const TRange & range, const TFunc & func, std::optional< size_t > cnt )
  {
    static Timer timer("CreateSortedTable");
    RegionTimer rt(timer);
    Table<TEntry, TIndex> table = CreateTable<TEntry, TIndex>(range, func, cnt);
    ParallelForRange
      (table.Range(), [&] (auto myrange)
       {
         for (auto i : myrange)
           QuickSort(table[i]);
       }, TasksPerThread(4)
      );

    return table;
  }

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



  /**
     A dynamic table class.

     A DynamicTable contains entries of variable size. Entry sizes can
     be increased dynamically.
  */
  template <class T, class IndexType = size_t>
  class DynamicTable 
  {
  protected:
    static constexpr IndexType BASE = IndexBASE<IndexType>();

    struct linestruct
    {
      int size;
      int maxsize;
      T * col;
    };

    Array<linestruct, IndexType> data;
    T * oneblock = nullptr;
  
  public:
    /// Creates table of size size
    DynamicTable (int size = 0)
      : data(size)
    {
      for (auto & d : data)
        {
          d.maxsize = 0;
          d.size = 0;
          d.col = nullptr;
        }
      oneblock = nullptr;
    }

    /// Creates table with a priori fixed entry sizes.
    DynamicTable (const Array<int, IndexType> & entrysizes, bool setentrysize=false)
      : data(entrysizes.Size())
    {
      size_t cnt = 0;
      // size_t n = entrysizes.Size();
      
      for (auto es : entrysizes)
        cnt += es;
      oneblock = new T[cnt];
      
      cnt = 0;
      for (auto i : data.Range())
        {
          data[i].maxsize = entrysizes[i];
          if (setentrysize)
            data[i].size = entrysizes[i];
          else
            data[i].size = 0;
          data[i].col = &oneblock[cnt];
          cnt += entrysizes[i];
        }
    }

    DynamicTable (DynamicTable && tab2)
    {
      Swap (data, tab2.data);
      Swap (oneblock, tab2.oneblock);
    }

    ~DynamicTable ()
    {
      if (oneblock)
        delete [] oneblock;
      else
        for (auto & d : data)
          delete [] d.col;
    }
    
    DynamicTable & operator= (DynamicTable && tab2)
    {
      Swap (data, tab2.data);
      Swap (oneblock, tab2.oneblock);
      return *this;
    }

    /// Changes Size of table to size, deletes data
    void SetSize (int size)
    {
      for (auto & d : data)
        delete [] d.col;
    
      data.SetSize(size);
      for (auto & d : data)
        {
          d.maxsize = 0;
          d.size = 0;
          d.col = nullptr;
        }
    }      

    void ChangeSize (size_t size)
    {
      if (oneblock)
        throw Exception ("cannot change size of oneblock dynamic table");
      
      size_t oldsize = data.Size();
      if (size == oldsize) 
        return;
      
      if (size < oldsize)
        for (int i = size; i < oldsize; i++)
          delete [] data[i+BASE].col;
      
      data.SetSize(size);

      for (int i = oldsize; i < size; i++)
        {
          data[i+BASE].maxsize = 0;
          data[i+BASE].size = 0;
          data[i+BASE].col = nullptr;
        }    
    }
    

    
    ///
    void IncSize (IndexType i)
    {
      NETGEN_CHECK_RANGE(i,BASE,data.Size()+BASE);
    
      linestruct & line = data[i];
    
      if (line.size == line.maxsize)
        {
          T * p;
          if constexpr (std::is_default_constructible<T>::value)
            p = new T[(2*line.maxsize+5)];
          else
            p = reinterpret_cast<T*>(new char[(2*line.maxsize+5)*sizeof(T)]);
          for (size_t i = 0; i < line.maxsize; i++)
            p[i] = std::move(line.col[i]);
          // memcpy (p, line.col, line.maxsize * sizeof(T));
          delete [] line.col;
          line.col = p;
          line.maxsize = 2*line.maxsize+5;
        }
    
      line.size++;
    }
  
    void DecSize (IndexType i)
    {
      NETGEN_CHECK_RANGE(i,BASE,data.Size()+BASE);
      linestruct & line = data[i];
    
#ifdef NETGEN_ENABLE_CHECK_RANGE
      if (line.size == 0)
        throw Exception ("BaseDynamicTable::Dec: EntrySize < 0");
#endif
    
      line.size--;
    }
  
  
    /// Inserts element acont into row i. Does not test if already used.
    void Add (IndexType i, const T & acont)
    {
      if (data[i].size == data[i].maxsize)
        this->IncSize (i);
      else
        data[i].size++;
      data[i].col[data[i].size-1] = acont;
    }
  
    /// Inserts element acont into row i, iff not yet exists.
    void AddUnique (IndexType i, const T & cont)
    {
      int es = EntrySize (i);
      T * line = data[i].col;
      for (int j = 0; j < es; j++)
        if (line[j] == cont)
          return;
      Add (i, cont);
    }
  

    /// Inserts element acont into row i. Does not test if already used.
    void AddEmpty (IndexType i)
    {
      IncSize (i);
    }
  
    /** Set the nr-th element in the i-th row to acont.
        Does not check for overflow. */
    void Set (IndexType i, int nr, const T & acont)
    {
      data[i].col[nr] = acont;
    }


    /** Returns the nr-th element in the i-th row.
        Does not check for overflow. */
    const T & Get (IndexType i, int nr) const
    {
      return data[i].col[nr];
    }
  
  
    /** Returns pointer to the first element in row i. */
    const T * GetLine (IndexType i) const
    {
      return data[i].col;
    }
  
    /// Returns size of the table.
    size_t Size () const
    {
      return data.Size();
    }

    auto Range () const
    {
      return data.Range();
    }

    /// Returns size of the i-th row.
    int EntrySize (IndexType i) const
    {
      return data[i].size;
    }

    ///
    void DecEntrySize (IndexType i)
    {
      DecSize(i);
    }
  
    /// Access entry i
    FlatArray<T> operator[] (IndexType i)
    {
      return FlatArray<T> (data[i].size, data[i].col);
    }
  
    /*
      typedef const FlatArray<T> ConstFlatArray;
      /// Access entry i
      ConstFlatArray operator[] (int i) const
      { return FlatArray<T> (data[i].size, static_cast<T*> (data[i].col)); }
    */
    FlatArray<T> operator[] (IndexType i) const
    {
      return FlatArray<T> (data[i].size, data[i].col);
    }
  };


  /// Print table
  template <class T>
  inline ostream & operator<< (ostream & s, const DynamicTable<T> & table)
  {
    for (auto i : Range(table))
      {
        s << i << ":";
        for (int j = 0; j < table[i].Size(); j++)
          s << " " << table[i][j];
        s << "\n";
      }
    s << std::flush;
    return s;
  }


  //   Helper function to calculate coloring of a set of indices for parallel processing of independent elements/points/etc.
  //   Assigns a color to each of colors.Size() elements, such that two elements with the same color don't share a common 'dof',
  //   the mapping from element to dofs is provided by the function getDofs(int) -> iterable<int>
  //
  //   Returns the number of used colors
  template <typename Tmask>
  int ComputeColoring( FlatArray<int> colors, size_t ndofs, Tmask const & getDofs)
  {
    static Timer timer("ComputeColoring - "+Demangle(typeid(Tmask).name())); RegionTimer rt(timer);
    static_assert(sizeof(unsigned int)==4, "Adapt type of mask array");
    size_t n = colors.Size();

    Array<unsigned int> mask(ndofs);

    size_t colored_blocks = 0;

    // We are coloring with 32 colors at once and use each bit to mask conflicts
    unsigned int check = 0;
    unsigned int checkbit = 0;

    int current_color = 0;
    colors = -1;
    int maxcolor = 0;

    while(colored_blocks<n)
    {
        mask = 0;
        for (auto i : Range(n) )
        {
            if(colors[i]>-1) continue;
            check = 0;
            const auto & dofs = getDofs(i);

            // Check if adjacent dofs are already marked by current color
            for (auto dof : dofs)
                check|=mask[dof];

            // Did we find a free color?
            if(check != 0xFFFFFFFF)
            {
                checkbit = 1;
                int color = current_color;
                // find the actual color, which is free (out of 32)
                while (check & checkbit)
                {
                    color++;
                    checkbit *= 2;
                }
                colors[i] = color;
                maxcolor = color > maxcolor ? color : maxcolor;
                colored_blocks++;
                // mask all adjacent dofs with the found color
                for (auto dof : dofs)
                    mask[dof] |= checkbit;
            }
        }
        current_color+=32;
    }
    return maxcolor+1;
  }


  typedef DynamicTable<int> IntTable;

} // namespace ngcore

#endif // NETGEN_CORE_TABLE_HPP

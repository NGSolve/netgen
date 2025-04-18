#ifndef FILE_NGSTD_HASHTABLE
#define FILE_NGSTD_HASHTABLE

/**************************************************************************/
/* File:   hashtable.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <string>
#include <tuple>
#include <optional>

// #include "mpi_wrapper.hpp"
#include "ngcore_api.hpp"
#include "table.hpp"
#include "utils.hpp"

namespace ngcore
{


  template <int K> 
  class MakeTupleFromInt
  {
  public:
    template <typename I>
    auto operator()(I & i)
    { return tuple_cat(MakeTupleFromInt<K-1> ()(i), std::tie(i[K-1])); }
  };
  
  template <> 
  class MakeTupleFromInt<1>
  {
  public:
    template <typename I>
    auto operator()(I & i) { return std::tie(i[0]); }
  };
  
  

  // feature check macro for transition from INT to IVec
#define NGCORE_HAS_IVEC
  
  /// N integers
  template <int N, typename T = int>
  class IVec
  {
    /// data
    // T i[(N>0)?N:1];

    HTArray<N,T> i;
    
  public:
    ///
    constexpr NETGEN_INLINE IVec () = default;
    constexpr NETGEN_INLINE IVec (const IVec & i1) : i(i1.i) { }

    constexpr NETGEN_INLINE IVec (T ai1) : i(ai1) { }
    
    template <class... T2,
              std::enable_if_t<N==1+sizeof...(T2),bool> = true>
    constexpr IVec (const T &v, T2... rest)
      : i{v,rest...} { } 

    /*
    /// init all
    NETGEN_INLINE IVec (T ai1)
    { 
     for (int j = 0; j < N; j++) { i[j] = ai1; }
    }

    /// init i[0], i[1]
    constexpr NETGEN_INLINE IVec (T ai1, T ai2)
      : i{ai1, ai2} { ; } 

    /// init i[0], i[1], i[2]
    constexpr NETGEN_INLINE IVec (T ai1, T ai2, T ai3)
      : i{ai1, ai2, ai3} { ; } 

    /// init i[0], i[1], i[2]
    constexpr NETGEN_INLINE IVec (T ai1, T ai2, T ai3, T ai4)
      : i{ai1, ai2, ai3, ai4} { ; }
    
    /// init i[0], i[1], i[2]
    constexpr NETGEN_INLINE IVec (T ai1, T ai2, T ai3, T ai4, T ai5)
      : i{ai1, ai2, ai3, ai4, ai5} { ; }      
      
    /// init i[0], i[1], i[2]
    NETGEN_INLINE IVec (T ai1, T ai2, T ai3, T ai4, T ai5, T ai6, T ai7, T ai8, T ai9)
      : i{ai1, ai2, ai3, ai4, ai5, ai6, ai7, ai8, ai9 } { ; }            
    */
    
    template <typename ARCHIVE>
    void DoArchive(ARCHIVE& ar)
    {
      // ar.Do(i.begin(), N);
      ar.Do(i.Ptr(), N);
    }

    template <int N2, typename T2>
    NETGEN_INLINE IVec (const IVec<N2,T2> & in2)
    {
      if (N2 <= N)
        {
          for (int j = 0; j < N2; j++)
            i[j] = in2[j];
          for (int j = N2; j < N; j++)
            i[j] = 0;
        }
      else
        {
          for (int j = 0; j < N; j++)
            i[j] = in2[j];
        }
    }

    template <typename T2>
    NETGEN_INLINE IVec (const BaseArrayObject<T2> & ao)
    {
      for (int j = 0; j < N; j++)
        i[j] = ao.Spec()[j];
    }
    
    NETGEN_INLINE size_t Size() const { return N; }
    /// all ints equal ?
    NETGEN_INLINE bool operator== (const IVec & in2) const
    { 
      for (int j = 0; j < N; j++) 
	if (i[j] != in2.i[j]) return 0;
      return 1; 
    }

    /// any ints unequal ?
    NETGEN_INLINE bool operator!= (const IVec & in2) const
    {
      for (int j = 0; j < N; j++)
        if (i[j] != in2.i[j]) return 1;
      return 0;
    }

    /// sort integers
    NETGEN_INLINE IVec & Sort () & 
    {
      for (int k = 0; k < N; k++)
	for (int l = k+1; l < N; l++)
	  if (i[k] > i[l]) 
	    Swap (i[k], i[l]);
      return *this;
    }

    NETGEN_INLINE IVec Sort () &&
    {
      for (int k = 0; k < N; k++)
	for (int l = k+1; l < N; l++)
	  if (i[k] > i[l]) 
	    Swap (i[k], i[l]);
      return *this;
    }

    /// access
    NETGEN_INLINE T & operator[] (int j)
    { return i[j]; }

    /// access
    NETGEN_INLINE constexpr const T & operator[] (int j) const
    { return i[j]; }

    template <size_t J>
    constexpr T get() const { return i[J]; }
    
    operator FlatArray<T> () { return FlatArray<T> (N, &i[0]); } 

    NETGEN_INLINE IVec<N,T> & operator= (T value)
    {
      for (int j = 0; j < N; j++)
	i[j] = value;
      return *this;
    }

    template <typename T2>
    NETGEN_INLINE IVec<N,T> & operator= (IVec<N,T2> v2)
    {
      for (int j = 0; j < N; j++)
	i[j] = v2[j];
      return *this;
    }

    template <typename... Ts>
    operator std::tuple<Ts...> ()
    {
      return MakeTupleFromInt<N>()(*this);
    }

    bool Contains (T val)
    {
      for (int j = 0; j < N; j++)
        if (i[j] == val) return true;
      return false;
    }
  };

  /// sort 2 integers
  template <>
  NETGEN_INLINE IVec<2> & IVec<2>::Sort () & 
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    return *this;
  }

  template <>
  NETGEN_INLINE IVec<2> IVec<2>::Sort () &&
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    return *this;
  }

  /// sort 3 integers
  template <>
  NETGEN_INLINE IVec<3> IVec<3>::Sort () &&
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    if (i[1] > i[2]) Swap (i[1], i[2]);
    if (i[0] > i[1]) Swap (i[0], i[1]);
    return *this;
  }

  /// Print integers
  template <int N, typename T>
  inline ostream & operator<<(ostream  & s, const IVec<N,T> & i2)
  {
    for (int j = 0; j < N; j++)
      s << (int) i2[j] << " ";
    return s;
  }
  
  template <int N, typename T>
  auto begin(const IVec<N,T> & ind)
  {
    return AOWrapperIterator<IVec<N,T>> (ind, 0);
  }

  template <int N, typename T>
  auto end(const IVec<N,T> & ind)
  {
    return AOWrapperIterator<IVec<N,T>> (ind, N);    
  }





  
  template <int N, typename TI>
  NETGEN_INLINE size_t HashValue (const IVec<N,TI> & ind, size_t size)
  {
    IVec<N,size_t> lind = ind;    
    size_t sum = 0;
    for (int i = 0; i < N; i++)
      sum += lind[i];
    return sum % size;
  }

  /// hash value of 1 int
  template <typename TI>
  NETGEN_INLINE size_t HashValue (const IVec<1,TI> & ind, size_t size) 
  {
    return ind[0] % size;
  }

  /// hash value of 2 int
  template <typename TI>  
  NETGEN_INLINE size_t HashValue (const IVec<2,TI> & ind, size_t size) 
  {
    IVec<2,size_t> lind = ind;
    return (113*lind[0]+lind[1]) % size;
  }

  /// hash value of 3 int
  template <typename TI>    
  NETGEN_INLINE size_t HashValue (const IVec<3,TI> & ind, size_t size) 
  {
    IVec<3,size_t> lind = ind;
    return (113*lind[0]+59*lind[1]+lind[2]) % size;
  }

  NETGEN_INLINE size_t HashValue (size_t ind, size_t size)
  {
    return ind%size;
  }
  NETGEN_INLINE size_t HashValue (int ind, size_t size)
  {
    return size_t(ind)%size;
  }
  





  
  template <int N, typename TI>
  NETGEN_INLINE constexpr size_t HashValue2 (const IVec<N,TI> & ind, size_t mask)
  {
    IVec<N,size_t> lind = ind;    
    size_t sum = 0;
    for (int i = 0; i < N; i++)
      sum += lind[i];
    return sum & mask;
  }

  /// hash value of 1 int
  template <typename TI>
  NETGEN_INLINE constexpr size_t HashValue2 (const IVec<1,TI> & ind, size_t mask) 
  {
    return ind[0] & mask;
  }

  /// hash value of 2 int
  template <typename TI>  
  NETGEN_INLINE constexpr size_t HashValue2 (const IVec<2,TI> & ind, size_t mask) 
  {
    IVec<2,size_t> lind = ind;
    return (113*lind[0]+lind[1]) & mask;
  }

  /// hash value of 3 int
  template <typename TI>    
  NETGEN_INLINE constexpr size_t HashValue2 (const IVec<3,TI> & ind, size_t mask) 
  {
    IVec<3,size_t> lind = ind;
    return (113*lind[0]+59*lind[1]+lind[2]) & mask;
  }

  NETGEN_INLINE constexpr size_t HashValue2 (size_t ind, size_t mask)
  {
    return ind & mask;
  }
  NETGEN_INLINE constexpr size_t HashValue2 (int ind, size_t mask)
  {
    return size_t(ind) & mask;
  }
  



  
  // using ngstd::max;

  template <int D, typename T>
  NETGEN_INLINE T Max (const IVec<D,T> & i)
  {
    if (D == 0) return 0;
    T m = i[0];
    for (int j = 1; j < D; j++)
      if (i[j] > m) m = i[j];
    return m;
  }

  template <int D, typename T>
  NETGEN_INLINE T Min (const IVec<D,T> & i)
  {
    if (D == 0) return 0;
    T m = i[0];
    for (int j = 1; j < D; j++)
      if (i[j] < m) m = i[j];
    return m;
  }

  template <int D, typename T>
  NETGEN_INLINE IVec<D,T> Max (IVec<D,T> i1, IVec<D,T> i2)
  {
    IVec<D,T> tmp;
    for (int i = 0; i < D; i++)
      tmp[i] = std::max(i1[i], i2[i]);
    return tmp;
  }

  template <int D, typename T>
  NETGEN_INLINE IVec<D,T> operator+ (IVec<D,T> i1, IVec<D,T> i2)
  {
    IVec<D,T> tmp;
    for (int i = 0; i < D; i++)
      tmp[i] = i1[i]+i2[i];
    return tmp;
  }
  










  /**
     A hash-table.
     Generic identifiers are mapped to the generic type T.
     An open hashtable. The table is implemented by a DynamicTable.
     Identifiers must provide a HashValue method.
  */
  template <class T_HASH, class T>
  class HashTable
  {
    /*
    DynamicTable<T_HASH> hash;
    DynamicTable<T> cont;
    */
    DynamicTable<std::pair<T_HASH,T>> table;
  public:
    /// Constructs a hashtable of size bags.
    NETGEN_INLINE HashTable (int size)
    // : hash(size), cont(size)
      : table(size)
    { ; }
    NETGEN_INLINE ~HashTable () { ; }

    /// Sets identifier ahash to value acont
    void Set (const T_HASH & ahash, const T & acont)
    {
      int bnr = HashValue (ahash, Size());
      int pos = CheckPosition (bnr, ahash);
      if (pos != -1)
	// cont.Set (bnr, pos, acont);
        table[bnr][pos].second = acont;
      else
	{
	  // hash.Add (bnr, ahash);
	  // cont.Add (bnr, acont);
          table.Add (bnr, std::make_pair(ahash, acont));
	}        
    }

    /// get value of identifier ahash, exception if unused
    const T & Get (const T_HASH & ahash) const
    {
      int bnr = HashValue (ahash, Size());
      int pos = Position (bnr, ahash);
      // return cont.Get (bnr, pos);
      return table.Get (bnr, pos).second;
    }

    /// get value of identifier ahash, exception if unused
    const T & Get (int bnr, int pos) const
    {
      // return cont.Get (bnr, pos);
      return table.Get (bnr, pos).second;
    }

    /// is identifier used ?
    bool Used (const T_HASH & ahash) const
    {
      // return (CheckPosition (HashValue (ahash, hash.Size()), ahash) != -1);
      return (CheckPosition (HashValue (ahash, table.Size()), ahash) != -1);
    }

    /// is identifier used ?
    bool Used (const T_HASH & ahash, int & bnr, int & pos) const
    {
      // bnr = HashValue (ahash, hash.Size());
      bnr = HashValue (ahash, Size());
      pos = CheckPosition (bnr, ahash);
      return (pos != -1);
    }


    /// number of hash entries
    size_t Size () const
    {
      // return hash.Size();
      return table.Size();
    }

    /// size of hash entry
    size_t EntrySize (int bnr) const
    {
      // return hash[bnr].Size();
      return table[bnr].Size();
    }

    /// get identifier and value of entry bnr, position colnr
    void GetData (int bnr, int colnr, T_HASH & ahash, T & acont) const
    {
      // ahash = hash[bnr][colnr];
      // acont = cont[bnr][colnr];
      ahash = table[bnr][colnr].first;
      acont = table[bnr][colnr].second;
    }

    /// set identifier and value of entry bnr, position colnr
    void SetData (int bnr, int colnr, const T_HASH & ahash, const T & acont)
    {
      // hash[bnr][colnr] = ahash;
      // cont[bnr][colnr] = acont;
      table[bnr][colnr] = std::make_pair(ahash, acont);
    }    

    /// returns position of index. returns -1 on unused
    int CheckPosition (int bnr, const T_HASH & ind) const
    {
      /*
      for (int i = 0; i < hash[bnr].Size(); i++)
	if (hash[bnr][i] == ind)
	  return i;
      */
      for (int i = 0; i < table[bnr].Size(); i++)
	if (table[bnr][i].first == ind)
	  return i;
      return -1;
    }

    /// returns position of index. exception on unused
    int Position (int bnr, const T_HASH & ind) const
    {
      for (int i = 0; i < table[bnr].Size(); i++)
	if (table[bnr][i].first == ind)
	  return i;
      throw Exception ("Ask for unused hash-value");
    }

    T & operator[] (T_HASH ahash)
    {
      int bnr, pos;
      if (Used (ahash, bnr, pos))
        return table[bnr][pos].second;
      else
        {
	  // hash.Add (bnr, ahash);
	  // cont.Add (bnr, T(0));
          table.Add (bnr, std::make_pair(ahash, T(0)));
          // return cont[bnr][cont[bnr].Size()-1];
          return table[bnr][table[bnr].Size()-1].second;
        }
    }

    const T & operator[] (T_HASH ahash) const
    {
      return Get(ahash);
    }

    class Iterator
    {
      const HashTable & ht;
      int bnr;
      int pos;
    public:
      Iterator (const HashTable & aht, int abnr, int apos)
        : ht(aht), bnr(abnr), pos(apos) { ; }
      std::pair<T_HASH,T> operator* () const
      {
        T_HASH hash; 
        T data;
        ht.GetData (bnr, pos, hash, data);
        return std::pair<T_HASH,T> (hash, data);
      }

      Iterator & operator++() 
      {
        pos++;
        if (pos == ht.EntrySize(bnr))
          {
            pos = 0;
            bnr++;
            for ( ; bnr < ht.Size(); bnr++)
              if (ht.EntrySize(bnr) != 0) break;
          }
        return *this;
      }
      
      bool operator!= (const Iterator & it2) { return bnr != it2.bnr || pos != it2.pos; }
    };

    Iterator begin () const 
    {
      int i = 0;
      for ( ; i < Size(); i++)
        if (EntrySize(i) != 0) break;
      return Iterator(*this, i,0); 
    }
    Iterator end () const { return Iterator(*this, Size(),0); }
  };



  inline size_t RoundUp2 (size_t i)
  {
    size_t res = 1;
    while (res < i) res *= 2; // hope it will never be too large 
    return res; 
  }

  template <typename T>
  constexpr inline T InvalidHash() { return T(-1); }

  template <typename T_HASH>
  struct CHT_trait
  {
    constexpr static inline T_HASH Invalid() { return InvalidHash<T_HASH>(); }
    constexpr static inline size_t HashValue (const T_HASH & hash, size_t mask) { return HashValue2(hash, mask); }
  };

  template <typename T1, typename T2>
  struct CHT_trait<std::tuple<T1,T2>>
  {
    constexpr static inline std::tuple<T1,T2> Invalid() { return { CHT_trait<T1>::Invalid(), CHT_trait<T2>::Invalid() } ; }
    constexpr static inline size_t HashValue (const std::tuple<T1,T2> & hash, size_t mask)
    {
      return (CHT_trait<T1>::HashValue(std::get<0>(hash), mask) + CHT_trait<T2>::HashValue(std::get<1>(hash),mask)) & mask;
    }
  };

  

  /**
     A closed hash-table.
     All information is stored in one fixed array.
     The array should be allocated with the double size of the expected number of entries.
  */
  template <class T_HASH, class T>
  class ClosedHashTable
  {
  protected:
    ///
    size_t size;
    size_t mask;
    ///
    size_t used = 0;
    ///
    Array<T_HASH> hash;
    ///
    Array<T> cont;
    ///
    // T_HASH invalid = -1;
    // static constexpr T_HASH invalid = InvalidHash<T_HASH>();
    static constexpr T_HASH invalid = CHT_trait<T_HASH>::Invalid();
  public:
    ///
    ClosedHashTable (size_t asize = 128)
      : size(RoundUp2(asize)), hash(size), cont(size)
    {
      mask = size-1;
      // hash = T_HASH(invalid);
      // hash = InvalidHash<T_HASH>();
      hash = CHT_trait<T_HASH>::Invalid();
    }

    ClosedHashTable (ClosedHashTable && ht2) = default;

    /// allocate on local heap
    ClosedHashTable (size_t asize, LocalHeap & lh)
      : size(RoundUp2(asize)), mask(size-1), hash(size, lh), cont(size, lh)
    {
      // hash = T_HASH(invalid);
      hash = InvalidHash<T_HASH>();
    }

    ClosedHashTable & operator= (ClosedHashTable && ht2) = default;

    /// 
    size_t Size() const
    {
      return size;
    }

    /// is position used
    bool UsedPos (size_t pos) const
    {
      return ! (hash[pos] == invalid); 
    }

    /// number of used elements
    size_t UsedElements () const
    {
      return used;
    }

    size_t Position (const T_HASH ind) const
    {
      // size_t i = HashValue2(ind, mask);
      size_t i = CHT_trait<T_HASH>::HashValue(ind, mask);
      while (true)
	{
	  if (hash[i] == ind) return i;
	  if (hash[i] == invalid) return size_t(-1);
          i = (i+1) & mask;          
	}
    }

    void DoubleSize()
    {
      ClosedHashTable tmp(2*Size());
      for (auto both : *this)
        tmp[both.first] = both.second;
      *this = std::move(tmp);
    }
    
    // returns true if new position is created
    bool PositionCreate (const T_HASH ind, size_t & apos)
    {
      if (UsedElements()*2 > Size()) DoubleSize();
      
      // size_t i = HashValue2 (ind, mask);
      size_t i = CHT_trait<T_HASH>::HashValue (ind, mask);

      while (true)
	{
	  if (hash[i] == invalid)
	    { 
	      hash[i] = ind; 
	      apos = i;
              used++;
	      return true;
	    }
	  if (hash[i] == ind) 
	    { 
	      apos = i; 
	      return false; 
	    }
          i = (i+1) & mask;
	}
    }


    ///
    void Set (const T_HASH & ahash, const T & acont)
    {
      size_t pos;
      PositionCreate (ahash, pos);
      hash[pos] = ahash;
      cont[pos] = acont;
    }

    ///
    const T & Get (const T_HASH & ahash) const
    {
      size_t pos = Position (ahash);
      if (pos == size_t(-1))
        throw Exception (std::string("illegal key: ") + ToString(ahash) );
      return cont[pos];
    }

    ///
    bool Used (const T_HASH & ahash) const
    {
      return (Position (ahash) != size_t(-1));
    }

    inline std::optional<T> GetIfUsed (const T_HASH & ahash) const
    {
      size_t pos = Position (ahash);
      if (pos != size_t(-1))
        return cont[pos];
      else
        return std::nullopt;
    }
    

    void SetData (size_t pos, const T_HASH & ahash, const T & acont)
    {
      hash[pos] = ahash;
      cont[pos] = acont;
    }

    void GetData (size_t pos, T_HASH & ahash, T & acont) const
    {
      ahash = hash[pos];
      acont = cont[pos];
    }
  
    void SetData (size_t pos, const T & acont)
    {
      cont[pos] = acont;
    }

    void GetData (size_t pos, T & acont) const
    {
      acont = cont[pos];
    }

    T GetData (size_t pos) const
    {
      return cont[pos];
    }

    std::pair<T_HASH,T> GetBoth (size_t pos) const
    {
      return std::pair<T_HASH,T> (hash[pos], cont[pos]);
    }

    const T & operator[] (T_HASH key) const { return Get(key); }
    T & operator[] (T_HASH key)
    {
      size_t pos;
      PositionCreate(key, pos);
      return cont[pos];
    }
    
    void SetSize (size_t asize)
    {
      size = asize;
      hash.Alloc(size);
      cont.Alloc(size);

      // for (size_t i = 0; i < size; i++)
      // hash[i] = invalid;
      hash = T_HASH(invalid);
    }

    void Delete (T_HASH key)
    {
      size_t pos = Position(key);
      if (pos == size_t(-1)) return;
      hash[pos] = invalid; used--;
      
      while (1)
        {
          size_t nextpos = pos+1;
          if (nextpos == size) nextpos = 0;
          if (hash[nextpos] == invalid) break;
          
          auto key = hash[nextpos];
          auto val = cont[nextpos];
          hash[pos] = invalid; used--;
          
          Set (key, val);
          pos = nextpos;
        }
    }

    void DeleteData()
    {
      hash = T_HASH(invalid);
      used = 0;
    }

    template <typename ARCHIVE>
    void DoArchive (ARCHIVE& ar)
    {
      ar & hash & cont;
      ar & size & mask & used;
    }    

    struct EndIterator { };
    
    class Iterator
    {
      const ClosedHashTable & tab;
      size_t nr;
    public:
      Iterator (const ClosedHashTable & _tab, size_t _nr)
        : tab(_tab), nr(_nr)
      {
        while (nr < tab.Size() && !tab.UsedPos(nr)) nr++;
      }
      Iterator & operator++()
      {
        nr++;
        while (nr < tab.Size() && !tab.UsedPos(nr)) nr++;
        return *this;
      }

      bool operator!= (EndIterator it2) { return nr != tab.Size(); }
      
      auto operator* () const { return tab.GetBoth(nr); }
    };

    Iterator begin() const { return Iterator(*this, 0); }
    EndIterator end() const { return EndIterator(); }
  };

  template <class T_HASH, class T>  
  ostream & operator<< (ostream & ost,
                        const ClosedHashTable<T_HASH,T> & tab)
  {
    /*
    for (size_t i = 0; i < tab.Size(); i++)
      if (tab.UsedPos(i))
        {
          T_HASH key;
          T val;
          tab.GetData (i, key, val);
          ost << key << ": " << val << ", ";
        }
    */
    for (auto [key,val] : tab)
      ost << key << ": " << val << ", ";      
    return ost;
  }

  template <typename TI>
  NETGEN_INLINE size_t HashValue (const IVec<3,TI> ind)
  {
    IVec<3,size_t> lind = ind;
    return 113*lind[0]+59*lind[1]+lind[2];
  }

  template <typename TI>  
  NETGEN_INLINE size_t HashValue (const IVec<2,TI> ind)
  {
    IVec<2,size_t> lind = ind;
    return 113*lind[0]+lind[1];
  }

  template <typename TI>  
  NETGEN_INLINE size_t HashValue (const IVec<1,TI> ind)
  {
    return ind[0];
  }


  template <typename TKEY, typename T>
  class ParallelHashTable
  {
    class ClosedHT
    {
      Array<TKEY> keys;
      Array<T> values;
      size_t used;
      
    public:
      ClosedHT(size_t asize = 256) : keys(asize), values(asize), used(0)
      {
        keys = TKEY(-1);
      }

      size_t Size () const { return keys.Size(); }
      size_t Used () const { return used; }

      ClosedHT & operator= (ClosedHT&&) = default;

      void Resize()
      {
        ClosedHT tmp(keys.Size()*2);
        for (size_t i = 0; i < keys.Size(); i++)
          if (keys[i] != TKEY(-1))
            {
              TKEY hkey = keys[i];
              T hval = values[i];
              size_t hhash = HashValue(hkey);
              size_t hhash2 = hhash / 256;
              tmp.DoSave(hkey, [hval] (T & v) { v = hval; }, hhash2);
            }
        (*this) = std::move(tmp);
      }
      
      template <typename TFUNC>
      auto Do (TKEY key, TFUNC func, size_t hash)
      {
        if (used > keys.Size()/2)
          Resize();
        return DoSave (key, func, hash);
      }
      
      template <typename TFUNC>
      auto DoSave (TKEY key, TFUNC func, size_t hash)
      {
        size_t pos = hash & (keys.Size()-1);
        while (1)
          {
            if (keys[pos] == key)
              break;
            if (keys[pos] == TKEY(-1))
              {
                keys[pos] = key;
                values[pos] = T(0);
                used++;
                break;
              }
            pos++;
            if (pos == keys.Size()) pos = 0;
          }
        return func(values[pos]);
      }
      
      T Get (TKEY key, size_t hash)
      {
        size_t pos = hash & (keys.Size()-1);
        while (1)
          {
            if (keys[pos] == key)
              return values[pos];
            if (keys[pos] == TKEY(-1))
              throw Exception ("ParallelHashTable::Get of unused key");
            pos++;
            if (pos == keys.Size()) pos = 0;
          }
      }

      size_t GetCosts (TKEY key, size_t hash)
      {
        size_t pos = hash & (keys.Size()-1);
        size_t costs = 1;
        while (1)
          {
            if (keys[pos] == key)
              return costs;
            if (keys[pos] == TKEY(-1))
              throw Exception ("ParallelHashTable::Get of unused key");
            costs++;
            pos++;
            if (pos == keys.Size()) pos = 0;
          }
      }


      template <typename TFUNC>
      void Iterate (TFUNC func) const
      {
        for (size_t i = 0; i < keys.Size(); i++)
          if (keys[i] != TKEY(-1))
            func(keys[i], values[i]);
      }
        
      void Print (ostream & ost) const
      {
        for (size_t i = 0; i < keys.Size(); i++)
          if (keys[i] != TKEY(-1))
            ost << keys[i] << ": " << values[i] << ", ";
      }
    };

    Array<ClosedHT> hts;
    class alignas(64) MyMutex64 : public MyMutex { };
    
    Array<MyMutex64> locks;

  public:
    ParallelHashTable() : hts(256), locks(256) { ; }
    size_t NumBuckets() const { return hts.Size(); }
    auto & Bucket(size_t nr) { return hts[nr]; }
    size_t BucketSize(size_t nr) const { return hts[nr].Size(); }
    size_t Used (size_t nr) const { return hts[nr].Used(); } 
    size_t Used() const
    {
      size_t used = 0;
      for (auto & ht : hts)
        used += ht.Used();
      return used;
    }  
    template <typename TFUNC>
    auto Do (TKEY key, TFUNC func)
    {
      size_t hash = HashValue(key);
      size_t hash1 = hash % 256;
      size_t hash2 = hash / 256;
      
      // locks[hash1].lock();
      // hts[hash1].Do (key, func, hash2);
      // locks[hash1].unlock();
      MyLock lock(locks[hash1]);
      return hts[hash1].Do (key, func, hash2);
    }
    
    T Get (TKEY key)
    {
      size_t hash = HashValue(key);
      size_t hash1 = hash % 256;
      size_t hash2 = hash / 256;
      
      return hts[hash1].Get (key, hash2);
    }

    auto GetCosts (TKEY key)
    {
      size_t hash = HashValue(key);
      size_t hash1 = hash % 256;
      size_t hash2 = hash / 256;
      
      return hts[hash1].GetCosts (key, hash2);
    }

    
    template <typename TFUNC>
    void Iterate(TFUNC func) const
    {
      for (auto & bucket : hts)
        bucket.Iterate(func);
    }

    template <typename TFUNC>
    void Iterate(size_t nr, TFUNC func) const
    {
      hts[nr].Iterate(func);
    }


    template <typename FUNC>
    void IterateParallel (FUNC func)
    {
      Array<size_t> base(NumBuckets());
      size_t sum = 0;
      for (size_t i = 0; i < NumBuckets(); i++)
        {
          base[i] = sum;
          sum += Used(i); 
        }
      ParallelFor(NumBuckets(),
                  [&] (size_t nr)
                  {
                    size_t cnt = base[nr];
                    Iterate(nr,
                            [&cnt, func] (TKEY key, T val)
                            {
                              func(cnt, key, val);
                              cnt++;
                            });
                  });
    }
    

    

    void Print (ostream & ost) const
    {
      for (size_t i : Range(hts))
        if (hts[i].Used() > 0)
          {
            ost << i << ": ";
            hts[i].Print(ost);
          }
    }
  };

  template <typename TKEY, typename T>
  inline ostream & operator<< (ostream & ost, const ParallelHashTable<TKEY,T> & ht)
  {
    ht.Print(ost);
    return ost;
  }









  template <class T, class IndexType>
  class CompressedTable
  {
    Table<T, size_t> table;
    ClosedHashTable<IndexType, size_t> idmap;
    
  public:
    CompressedTable (Table<T, size_t> && atable, ClosedHashTable<IndexType, size_t> && aidmap)
      : table(std::move(atable)), idmap(std::move(aidmap)) { }

    FlatArray<T> operator[](IndexType id) const
    {
      if (auto nr = idmap.GetIfUsed(id))
        return table[*nr];
      else
        return { 0, nullptr };
    }
    auto & GetTable() { return table; }
  };


  template <class T, typename IndexType>
  class CompressedTableCreator
  {
  protected:
    int mode;    // 1 .. cnt, 2 .. cnt entries, 3 .. fill table
    size_t nd;   // number of entries;
    ClosedHashTable<IndexType, size_t> idmap;
    Array<int,size_t> cnt;
    Table<T,size_t> table;
  public:
    CompressedTableCreator()
    { nd = 0; mode = 1; }

    CompressedTable<T,IndexType> MoveTable()
    {
      return { std::move(table), std::move(idmap) };
    }

    bool Done () { return mode > 3; }
    void operator++(int) { SetMode (mode+1); }

    int GetMode () const { return mode; }
    void SetMode (int amode)
    {
      mode = amode;
      if (mode == 2)
	{
          cnt.SetSize(nd);  
          cnt = 0;
	}
      if (mode == 3)
	{
          table = Table<T,size_t> (cnt);
          cnt = 0;
	}
    }

    void Add (IndexType blocknr, const T & data)
    {
      switch (mode)
	{
	case 1:
          {
            if (!idmap.Used (blocknr))
              idmap[blocknr] = nd++;
            break;
          }
	case 2:
	  cnt[idmap.Get(blocknr)]++;
	  break;
	case 3:
          size_t cblock = idmap.Get(blocknr);
          int ci = cnt[cblock]++;
          table[cblock][ci] = data;
	  break;
	}
    }
  };


  









  
} // namespace ngcore


/*
#ifdef PARALLEL
namespace ngcore {
  template<int S, typename T>
  class MPI_typetrait<ngcore::IVec<S, T> >
  {
  public:
    /// gets the MPI datatype
    static MPI_Datatype MPIType () 
    { 
      static MPI_Datatype MPI_T = 0;
      if (!MPI_T)
	{
	  MPI_Type_contiguous ( S, MPI_typetrait<T>::MPIType(), &MPI_T);
	  MPI_Type_commit ( &MPI_T );
	}
      return MPI_T;
    }
  };
}
#endif
*/

namespace ngcore
{
  template<typename T> struct MPI_typetrait;
  
  template<int S, typename T>
  struct MPI_typetrait<IVec<S, T> > {
    static auto MPIType () {
      return MPI_typetrait<std::array<T,S>>::MPIType();
    }
  };
}



namespace std
{
  // structured binding support
  template <auto N, typename T>
  struct tuple_size<ngcore::IVec<N,T>> : std::integral_constant<std::size_t, N> {};
  template<size_t N, auto M, typename T> struct tuple_element<N,ngcore::IVec<M,T>> { using type = T; };
}

#endif

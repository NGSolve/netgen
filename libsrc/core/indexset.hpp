#ifndef NETGEN_CORE_INDEXSET_HPP
#define NETGEN_CORE_INDEXSET_HPP

#include "array.hpp"
#include "bitarray.hpp"

namespace ngcore
{
  /// set of small non-negative integers: O(1) insert and lookup, O(members) clear and iteration
  class IndexSet
  {
    Array<int> set;
    BitArray flags;
  public:
    IndexSet (size_t maxind = 0) { SetMaxIndex (maxind); }

    /// grow index range to at least maxind
    void SetMaxIndex (size_t maxind)
    {
      if (maxind > flags.Size())
        {
          flags.SetSize (2*maxind);
          flags.Clear();
        }
    }

    bool Contains (int ind) const { return flags.Test (ind); }

    void Add (int ind)
    {
      if (!flags.Test (ind))
        {
          set.Append (ind);
          flags.SetBit (ind);
        }
    }

    void Clear ()
    {
      for (int ind : set)
        flags.Clear (ind);
      set.SetSize (0);
    }

    FlatArray<int> GetArray () const { return set; }
    size_t Size () const { return set.Size(); }
    auto begin () const { return set.begin(); }
    auto end () const { return set.end(); }
  };
}

#endif // NETGEN_CORE_INDEXSET_HPP

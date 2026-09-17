#ifndef NETGEN_CORE_DYNSTRIDEARRAY_HPP
#define NETGEN_CORE_DYNSTRIDEARRAY_HPP

/*
  Array of records with a run-time width:
    slot i = [ Header ][ Tail_0 x width ][ Tail_1 x width ] ...
  All tails share one width, slots are `stride` bytes apart.
  Header may be void (pure uniform-width table).
  operator[] and iteration hand out DynStrideView / const views.
*/

#include <cstring>
#include <new>
#include <type_traits>

#include "array.hpp"
#include "memtracer.hpp"

namespace ngcore
{
  namespace detail
  {
    template <class T>
    constexpr size_t SizeOrZero () { if constexpr (std::is_void_v<T>) return 0; else return sizeof(T); }
    template <class T>
    constexpr size_t AlignOrOne () { if constexpr (std::is_void_v<T>) return 1; else return alignof(T); }
    constexpr size_t RoundUpTo (size_t n, size_t a) { return (n + a - 1) / a * a; }
  }

  template <class... Ts> struct TailList { };

  template <class Header, class... Tails>
  struct DynStrideLayout
  {
    static constexpr size_t NT = sizeof...(Tails);
    static constexpr size_t tail_size[NT] = { sizeof(Tails)... };
    static constexpr size_t tail_align[NT] = { alignof(Tails)... };
    static constexpr size_t max_align =
      std::max({ detail::AlignOrOne<Header>(), alignof(Tails)... });

    size_t width = 0;
    size_t offset[NT] = {};
    size_t stride = 0;

    constexpr DynStrideLayout () : DynStrideLayout(0) { }
    constexpr explicit DynStrideLayout (size_t awidth) : width(awidth)
    {
      size_t pos = detail::SizeOrZero<Header>();
      for (size_t k = 0; k < NT; k++)
        {
          pos = detail::RoundUpTo (pos, tail_align[k]);
          offset[k] = pos;
          pos += width * tail_size[k];
        }
      stride = detail::RoundUpTo (pos, max_align);
    }
  };


  template <class Header, bool IS_CONST, class... Tails>
  class DynStrideView
  {
  public:
    using Layout = DynStrideLayout<Header, Tails...>;
    using CharPtr = std::conditional_t<IS_CONST, const char*, char*>;
    template <class T> using Ptr = std::conditional_t<IS_CONST, const T, T>;
    template <size_t K> using Tail = std::tuple_element_t<K, std::tuple<Tails...>>;

  private:
    CharPtr p;
    const Layout * layout;

  public:
    NETGEN_INLINE DynStrideView (CharPtr ap, const Layout * alayout) : p(ap), layout(alayout) { }
    NETGEN_INLINE DynStrideView (const DynStrideView<Header, false, Tails...> & v2)
      : p(v2.Data()), layout(&v2.GetLayout()) { }

    NETGEN_INLINE size_t Width () const { return layout->width; }
    NETGEN_INLINE CharPtr Data () const { return p; }
    NETGEN_INLINE const Layout & GetLayout () const { return *layout; }

    template <class H = Header, class = std::enable_if_t<!std::is_void_v<H>>>
    NETGEN_INLINE Ptr<H> & Head () const { return *reinterpret_cast<Ptr<H>*> (p); }

    template <size_t K = 0>
    NETGEN_INLINE Ptr<Tail<K>> * TailPtr () const
    { return reinterpret_cast<Ptr<Tail<K>>*> (p + layout->offset[K]); }

    // first n entries of tail K, default full width
    template <size_t K = 0>
    NETGEN_INLINE FlatArray<Ptr<Tail<K>>> TailArray (size_t n) const
    {
      NETGEN_CHECK_RANGE(n, 0, layout->width+1);
      return FlatArray<Ptr<Tail<K>>> (n, TailPtr<K>());
    }
    template <size_t K = 0>
    NETGEN_INLINE FlatArray<Ptr<Tail<K>>> TailArray () const
    { return FlatArray<Ptr<Tail<K>>> (layout->width, TailPtr<K>()); }

    // single-tail record without header: the view acts as a FlatArray of its tail
    template <size_t N = sizeof...(Tails), class = std::enable_if_t<N == 1 && std::is_void_v<Header>>>
    NETGEN_INLINE size_t Size () const { return layout->width; }
    template <size_t N = sizeof...(Tails), class = std::enable_if_t<N == 1 && std::is_void_v<Header>>>
    NETGEN_INLINE decltype(auto) operator[] (size_t k) const
    {
      NETGEN_CHECK_RANGE(k, 0, layout->width);
      return (TailPtr<0>()[k]);
    }
    template <size_t N = sizeof...(Tails), class = std::enable_if_t<N == 1 && std::is_void_v<Header>>>
    NETGEN_INLINE auto begin () const { return TailPtr<0>(); }
    template <size_t N = sizeof...(Tails), class = std::enable_if_t<N == 1 && std::is_void_v<Header>>>
    NETGEN_INLINE auto end () const { return TailPtr<0>() + layout->width; }

    // copy header and min(width) tail entries from another view
    template <bool C2, class V = DynStrideView, class = std::enable_if_t<!IS_CONST, V>>
    void Assign (const DynStrideView<Header, C2, Tails...> & v2) const
    {
      if constexpr (!std::is_void_v<Header>)
        std::memcpy (p, v2.Data(), sizeof(Header));
      size_t w = std::min (Width(), v2.Width());
      CopyTails (v2, w, std::make_index_sequence<sizeof...(Tails)>());
    }

  private:
    template <bool C2, size_t... K>
    void CopyTails (const DynStrideView<Header, C2, Tails...> & v2, size_t w, std::index_sequence<K...>) const
    {
      ((std::memcpy ((void*)TailPtr<K>(), v2.template TailPtr<K>(), w*sizeof(Tail<K>))), ...);
    }
  };


  template <class Header, class TAILS, class IndexType = size_t>
  class FlatDynStrideArray;
  template <class Header, class TAILS, class IndexType = size_t>
  class DynStrideArray;


  /// non-owning: slots [0,size) at data, layout owned elsewhere
  template <class Header, class IndexType, class... Tails>
  class FlatDynStrideArray<Header, TailList<Tails...>, IndexType>
  {
  public:
    using Layout = DynStrideLayout<Header, Tails...>;
    using View = DynStrideView<Header, false, Tails...>;
    using ConstView = DynStrideView<Header, true, Tails...>;
    using RawView = View;
    using RawConstView = ConstView;
    static constexpr size_t NT = sizeof...(Tails);
    using index_type = IndexType;
    static constexpr IndexType BASE = IndexBASE<IndexType>();

  protected:
    char * data;
    size_t size;
    const Layout * layout;

  public:
    NETGEN_INLINE FlatDynStrideArray (size_t asize, char * adata, const Layout * alayout)
      : data(adata), size(asize), layout(alayout) { }

    NETGEN_INLINE size_t Size () const { return size; }
    NETGEN_INLINE size_t Width () const { return layout->width; }
    NETGEN_INLINE size_t Stride () const { return layout->stride; }
    NETGEN_INLINE const Layout & GetLayout () const { return *layout; }
    NETGEN_INLINE char * Data () { return data; }
    NETGEN_INLINE const char * Data () const { return data; }
    NETGEN_INLINE T_Range<IndexType> Range () const { return T_Range<IndexType> (BASE, BASE+size); }

    NETGEN_INLINE View operator[] (IndexType i)
    {
      NETGEN_CHECK_RANGE(i, BASE, BASE+size);
      return View (data + (i-BASE)*layout->stride, layout);
    }
    NETGEN_INLINE ConstView operator[] (IndexType i) const
    {
      NETGEN_CHECK_RANGE(i, BASE, BASE+size);
      return ConstView (data + (i-BASE)*layout->stride, layout);
    }

    NETGEN_INLINE View First () { return (*this)[BASE]; }
    NETGEN_INLINE ConstView First () const { return (*this)[BASE]; }
    NETGEN_INLINE View Last () { return (*this)[BASE+size-1]; }
    NETGEN_INLINE ConstView Last () const { return (*this)[BASE+size-1]; }

    // 0-based sub-array
    NETGEN_INLINE FlatDynStrideArray<Header, TailList<Tails...>, size_t> Range (T_Range<size_t> r) const
    {
      NETGEN_CHECK_RANGE(r.Next(), 0, size+1);
      return { r.Size(), data + r.First()*layout->stride, layout };
    }
    template <typename TI = IndexType, typename = std::enable_if_t<!std::is_same_v<TI,size_t>>>
    NETGEN_INLINE FlatDynStrideArray<Header, TailList<Tails...>, size_t> Range (T_Range<IndexType> r) const
    {
      return Range (T_Range<size_t> (r.First()-BASE, r.Next()-BASE));
    }
    NETGEN_INLINE auto Range (size_t start, size_t end) const
    {
      return Range (T_Range<size_t> (start, end));
    }
    NETGEN_INLINE auto Range (size_t start, IndexFromEnd indend) const
    {
      return Range (start, size_t(size+indend.Value()));
    }
    template <typename TI = IndexType, typename = std::enable_if_t<!std::is_same_v<TI,size_t>>>
    NETGEN_INLINE auto Range (IndexType from, IndexType next) const
    {
      return Range (T_Range<IndexType> (from, next));
    }
    template <typename TI = IndexType, typename = std::enable_if_t<!std::is_same_v<TI,size_t>>>
    NETGEN_INLINE auto Range (IndexType from, IndexFromEnd indend) const
    {
      return Range (from, IndexType(BASE+size+indend.Value()));
    }

    template <class RAW, class V>
    class Iterator
    {
      typename RAW::CharPtr p;
      const Layout * layout;
    public:
      NETGEN_INLINE Iterator (typename RAW::CharPtr ap, const Layout * alayout) : p(ap), layout(alayout) { }
      NETGEN_INLINE Iterator & operator++ () { p += layout->stride; return *this; }
      NETGEN_INLINE V operator* () const { return V (p, layout); }
      NETGEN_INLINE bool operator!= (const Iterator & it2) const { return p != it2.p; }
      NETGEN_INLINE bool operator== (const Iterator & it2) const { return p == it2.p; }
    };

    NETGEN_INLINE auto begin () { return Iterator<RawView, View> (data, layout); }
    NETGEN_INLINE auto end () { return Iterator<RawView, View> (data + size*layout->stride, layout); }
    NETGEN_INLINE auto begin () const { return Iterator<RawConstView, ConstView> (data, layout); }
    NETGEN_INLINE auto end () const { return Iterator<RawConstView, ConstView> (data + size*layout->stride, layout); }
  };


  template <class Header, class IndexType, class... Tails>
  class DynStrideArray<Header, TailList<Tails...>, IndexType>
    : public FlatDynStrideArray<Header, TailList<Tails...>, IndexType>
  {
    static_assert (std::is_void_v<Header> || std::is_trivially_copyable_v<Header>,
                   "DynStrideArray: header must be trivially copyable");
    static_assert ((std::is_trivially_copyable_v<Tails> && ...),
                   "DynStrideArray: tails must be trivially copyable");
    static_assert (sizeof...(Tails) > 0, "DynStrideArray: needs at least one tail");

    using FLAT = FlatDynStrideArray<Header, TailList<Tails...>, IndexType>;
  public:
    using typename FLAT::Layout;
    using typename FLAT::RawView;
    using typename FLAT::RawConstView;
    using typename FLAT::View;
    using typename FLAT::ConstView;
    using FLAT::NT;
    using FLAT::BASE;
    using FLAT::Size;
    using FLAT::Width;
    using FLAT::operator[];

  private:
    using FLAT::data;
    using FLAT::size;
    size_t allocsize = 0;   // in slots
    Layout ownlayout;
    MemoryTracer mt;

  public:
    DynStrideArray () : FLAT(0, nullptr, &ownlayout) { }
    explicit DynStrideArray (size_t awidth) : FLAT(0, nullptr, &ownlayout), ownlayout(awidth) { }
    DynStrideArray (size_t asize, size_t awidth) : FLAT(0, nullptr, &ownlayout), ownlayout(awidth)
    {
      Allocate (asize);
      size = asize;
      Init (0, size);
    }

    DynStrideArray (const DynStrideArray & a2) : FLAT(0, nullptr, &ownlayout), ownlayout(a2.ownlayout)
    {
      Allocate (a2.size);
      size = a2.size;
      if (size) std::memcpy (data, a2.data, size*ownlayout.stride);
    }

    DynStrideArray (DynStrideArray && a2)
      : FLAT(a2.size, a2.data, &ownlayout), allocsize(a2.allocsize), ownlayout(a2.ownlayout)
    {
      a2.data = nullptr; a2.size = 0; a2.allocsize = 0;
    }

    ~DynStrideArray () { Free(); }

    DynStrideArray & operator= (const DynStrideArray & a2)
    {
      if (&a2 == this) return *this;
      Free();
      ownlayout = a2.ownlayout;
      Allocate (a2.size);
      size = a2.size;
      if (size) std::memcpy (data, a2.data, size*ownlayout.stride);
      return *this;
    }

    DynStrideArray & operator= (DynStrideArray && a2)
    {
      Swap (a2);
      return *this;
    }

    void Swap (DynStrideArray & a2)
    {
      std::swap (data, a2.data);
      std::swap (size, a2.size);
      std::swap (allocsize, a2.allocsize);
      std::swap (ownlayout, a2.ownlayout);
    }

    NETGEN_INLINE size_t AllocSize () const { return allocsize; }

    // keeps contents, new slots default-initialized
    void SetSize (size_t nsize)
    {
      if (nsize > allocsize) ReSize (nsize);
      if (nsize > size) Init (size, nsize);
      size = nsize;
    }
    NETGEN_INLINE void SetSize0 () { size = 0; }

    void SetAllocSize (size_t nallocsize)
    {
      if (nallocsize > allocsize) ReSize (nallocsize);
    }

    // re-stride to a new width, keeps min(old,new) tail entries per slot
    void SetWidth (size_t nwidth)
    {
      if (nwidth == ownlayout.width) return;
      Layout nlayout (nwidth);
      size_t nalloc = std::max (allocsize, size);
      char * ndata = nalloc ? AllocBytes (nalloc*nlayout.stride) : nullptr;
      for (size_t i = 0; i < size; i++)
        RawView (ndata + i*nlayout.stride, &nlayout)
          .Assign (RawConstView (data + i*ownlayout.stride, &ownlayout));
      if (nwidth > ownlayout.width)
        for (size_t i = 0; i < size; i++)
          InitTails (ndata + i*nlayout.stride, nlayout, ownlayout.width, nwidth,
                     std::make_index_sequence<NT>());
      Free();
      data = ndata;
      allocsize = nalloc;
      ownlayout = nlayout;
      if (allocsize) mt.Alloc (allocsize*ownlayout.stride);
    }

    // append default-initialized slot, returns its index
    IndexType Append ()
    {
      if (size == allocsize) ReSize (size+1);
      Init (size, size+1);
      return BASE + size++;
    }

    // append copy of a raw view, growing the width if needed
    template <bool C2>
    IndexType Append (const DynStrideView<Header, C2, Tails...> & v2)
    {
      if (v2.Width() > ownlayout.width) SetWidth (v2.Width());
      if (size == allocsize) ReSize (size+1);
      RawView (data + size*ownlayout.stride, &ownlayout).Assign (v2);
      if (v2.Width() < ownlayout.width)
        InitTails (data + size*ownlayout.stride, ownlayout, v2.Width(), ownlayout.width,
                   std::make_index_sequence<NT>());
      return BASE + size++;
    }

    // delete slot i, move last slot into it (Array semantics)
    void DeleteElement (IndexType ind)
    {
      NETGEN_CHECK_RANGE(ind, BASE, BASE+size);
      size_t i = ind-BASE;
      if (i+1 < size)
        std::memcpy (data + i*ownlayout.stride, data + (size-1)*ownlayout.stride, ownlayout.stride);
      size--;
    }

    // delete slot i, move all following slots forward
    void RemoveElement (IndexType ind)
    {
      NETGEN_CHECK_RANGE(ind, BASE, BASE+size);
      size_t i = ind-BASE;
      if (i+1 < size)
        std::memmove (data + i*ownlayout.stride, data + (i+1)*ownlayout.stride, (size-i-1)*ownlayout.stride);
      size--;
    }

    template <typename FUNC>
    void RemoveElementIf (FUNC func)
    {
      size_t keep = 0;
      for (size_t j = 0; j < size; j++)
        if (!func (View (data + j*ownlayout.stride, &ownlayout)))
          {
            if (keep != j)
              std::memcpy (data + keep*ownlayout.stride, data + j*ownlayout.stride, ownlayout.stride);
            keep++;
          }
      size = keep;
    }

    NETGEN_INLINE void DeleteLast ()
    {
      NETGEN_CHECK_RANGE(size-1, 0, size);
      size--;
    }

    void DeleteAll ()
    {
      Free();
      size = 0;
    }

    const MemoryTracer & GetMemoryTracer () const { return mt; }
    void StartMemoryTracing () const { mt.Alloc (allocsize*ownlayout.stride); }

  private:
    static char * AllocBytes (size_t nbytes)
    {
      return static_cast<char*> (::operator new (nbytes, std::align_val_t(Layout::max_align)));
    }
    static void FreeBytes (char * p)
    {
      ::operator delete (p, std::align_val_t(Layout::max_align));
    }

    void Allocate (size_t nalloc)
    {
      allocsize = nalloc;
      data = nalloc ? AllocBytes (nalloc*ownlayout.stride) : nullptr;
      if (nalloc) mt.Alloc (nalloc*ownlayout.stride);
    }

    void Free ()
    {
      if (data)
        {
          mt.Free (allocsize*ownlayout.stride);
          FreeBytes (data);
        }
      data = nullptr;
      allocsize = 0;
    }

    void ReSize (size_t minsize)
    {
      size_t nsize = std::max (2*allocsize, minsize);
      char * ndata = AllocBytes (nsize*ownlayout.stride);
      if (data)
        {
          std::memcpy (ndata, data, size*ownlayout.stride);
          Free();
        }
      data = ndata;
      allocsize = nsize;
      mt.Alloc (allocsize*ownlayout.stride);
    }

    // default-construct header and tails of slots [first, next)
    void Init (size_t first, size_t next)
    {
      for (size_t i = first; i < next; i++)
        {
          char * p = data + i*ownlayout.stride;
          if constexpr (!std::is_void_v<Header>)
            new (p) Header();
          InitTails (p, ownlayout, 0, ownlayout.width, std::make_index_sequence<NT>());
        }
    }

    template <size_t... K>
    static void InitTails (char * p, const Layout & l, size_t from, size_t to, std::index_sequence<K...>)
    {
      ((InitTail<K> (p + l.offset[K], from, to)), ...);
    }

    template <size_t K>
    static void InitTail (char * p, size_t from, size_t to)
    {
      using T = std::tuple_element_t<K, std::tuple<Tails...>>;
      for (size_t j = from; j < to; j++)
        new (p + j*sizeof(T)) T();
    }
  };

} // namespace ngcore

#endif // NETGEN_CORE_DYNSTRIDEARRAY_HPP

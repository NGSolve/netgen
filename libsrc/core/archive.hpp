#ifndef NETGEN_CORE_ARCHIVE_HPP
#define NETGEN_CORE_ARCHIVE_HPP

#include <algorithm>
#include <any>
#include <array>                // for array
#include <complex>              // for complex
#include <cstring>              // for size_t, strlen
#include <filesystem>           // for path
#include <fstream>              // for ifstream, ofstream
#include <functional>           // for function
#include <map>                  // for map
#include <memory>               // for shared_ptr
#include <optional>             // for optional
#include <string>               // for string
#include <type_traits>          // for declval, enable_if_t, false_type, is_co...
#include <cstddef>              // for std::byte
#include <set>                  // for set
#include <typeinfo>             // for type_info
#include <utility>              // for move, swap, pair
#include <vector>               // for vector

#include "exception.hpp"        // for UnreachableCodeException, Exception
#include "ngcore_api.hpp"       // for NGCORE_API
#include "type_traits.hpp"      // for all_of_tmpl
#include "utils.hpp"            // for Demangle, unlikely
#include "version.hpp"          // for VersionInfo

#ifdef NETGEN_PYTHON
namespace pybind11
{
  class object;
}
#endif // NETGEN_PYTHON

namespace ngcore
{
  template <typename T>
  struct Shallow {
    T val;
    Shallow() = default;
    Shallow(T aval) : val(aval) { ; }
    operator T&() { return val; }
  };

  // Helper to detect shared_from_this
  template <typename T>
  class has_shared_from_this2
    {
    private:
      // typedef T* T_ptr;
      template <typename C> static std::true_type test(decltype(((C*)nullptr)->shared_from_this()));
      template <typename C> static std::false_type test(...);
      
    public:
      // If the test returns true_type, then T has shared_from_this
      static constexpr bool value = decltype(test<T>(0))::value;
  };
  
  

  
  template <typename T, typename = void>
  class has_shallow_archive : public std::false_type {};
  
  template <typename T>
  class has_shallow_archive<T, std::void_t<decltype(T::shallow_archive)>>
    : public std::is_same<decltype(T::shallow_archive), std::true_type> {};
  

  
#ifdef NETGEN_PYTHON
  pybind11::object CastAnyToPy(const std::any& a);
#endif // NETGEN_PYTHON

  class NGCORE_API Archive;
  namespace detail
  {
    template <class T, class Tuple, size_t... Is>
    T* construct_from_tuple(Tuple&& tuple, std::index_sequence<Is...> ) {
      // return new T{std::get<Is>(std::forward<Tuple>(tuple))...};
      return new T{std::get<Is>(std::move(tuple))...};
    }

    template <class T, class Tuple>
    T* construct_from_tuple(Tuple&& tuple) {
      return construct_from_tuple<T>(std::forward<Tuple>(tuple),
                                     std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>{}
                                     );
    }

    // create new pointer of type T if it is default constructible, else throw
    template<typename T, typename... TArgs>
    T* constructIfPossible(std::tuple<TArgs...> args)
    {
      if constexpr(std::is_constructible_v<T, TArgs...>)
        return construct_from_tuple<T>(args);
          throw Exception(std::string(Demangle(typeid(T).name())) +
                          " is not constructible!");
    }

    template <typename T> T *constructIfPossible()
    {
      if constexpr(std::is_constructible_v<T>)
        return new T();
      throw Exception(std::string(Demangle(typeid(T).name())) +
                      " is not default constructible!");
    }

    //Type trait to check if a class implements a 'void DoArchive(Archive&)' function
    template<typename T>
    struct has_DoArchive
    {
    private:
      template<typename T2>
      static constexpr auto check(T2*) ->
        typename std::is_same<decltype(std::declval<T2>().DoArchive(std::declval<Archive&>())),void>::type;
      template<typename>
      static constexpr std::false_type check(...);
      using type = decltype(check<T>(nullptr)); // NOLINT
    public:
      NGCORE_API static constexpr bool value = type::value;
    };

    // Check if class is archivable
    template<typename T>
    struct is_Archivable_struct
    {
    private:
      template<typename T2>
      static constexpr auto check(T2*) ->
        typename std::is_same<decltype(std::declval<Archive>() & std::declval<T2&>()),Archive&>::type;
      template<typename>
      static constexpr std::false_type check(...);
      using type = decltype(check<T>(nullptr)); // NOLINT
    public:
      NGCORE_API static constexpr bool value = type::value;
    };

    template <typename T>
    struct has_GetCArgs
    {
      template <typename C> static std::true_type check( decltype( sizeof(&C::GetCArgs )) ) { return std::true_type(); }
      template <typename> static std::false_type check(...) { return std::false_type(); }
      typedef decltype( check<T>(sizeof(char)) ) type;
      static constexpr type value = type();
    };
    template<typename T>
    constexpr bool has_GetCArgs_v = has_GetCArgs<T>::value;

    template<typename T,
    typename std::enable_if<!has_GetCArgs_v<T>>::type* = nullptr>
    std::tuple<> GetCArgs(T&val) { return {}; }

    template<typename T,
    typename std::enable_if<has_GetCArgs_v<T>>::type* = nullptr>
    auto GetCArgs(T&val) {
      return val.GetCArgs();
    }

    template<typename T>
    using TCargs = decltype(GetCArgs<T>(*static_cast<T*>(nullptr)));


    struct ClassArchiveInfo
    {
      // create new object of this type and return a void* pointer that is points to the location
      // of the (base)class given by type_info
      // std::function<void*(const std::type_info&)> creator;
      void* (*creator)(const std::type_info&, Archive&);
      // This caster takes a void* pointer to the type stored in this info and casts it to a
      // void* pointer pointing to the (base)class type_info
      // std::function<void*(const std::type_info&, void*)> upcaster;
      void* (*upcaster) (const std::type_info&, void*);
      // This caster takes a void* pointer to the (base)class type_info and returns void* pointing
      // to the type stored in this info
      // std::function<void*(const std::type_info&, void*)> downcaster;
      void* (*downcaster)(const std::type_info&, void*);

      // Archive constructor arguments
      // std::function<void(Archive&, void*)> cargs_archiver;
      void (*cargs_archiver)(Archive&, void*);

#ifdef NETGEN_PYTHON
      // std::function<pybind11::object(const std::any&)> anyToPyCaster;
      pybind11::object (*anyToPyCaster)(const std::any&);
#endif // NETGEN_PYTHON
    };
  } // namespace detail

  template<typename T>
  constexpr bool is_archivable = detail::is_Archivable_struct<T>::value;

  
  template <typename T, typename ... Trest>
  constexpr size_t TotSize () 
  {
    if constexpr (sizeof...(Trest) == 0)
                   return sizeof(T);
    else
      return sizeof(T) + TotSize<Trest...> ();
  }
  
  
  // Base Archive class
  class NGCORE_API Archive
  {
    const bool is_output;
    // how many different shared_ptr/pointer have been (un)archived
    int shared_ptr_count{0}, ptr_count{0};
    // maps for archived shared pointers and pointers
    std::map<void*, int> shared_ptr2nr{}, ptr2nr{};
    // vectors for storing the unarchived (shared) pointers
    std::vector<std::shared_ptr<void>> nr2shared_ptr{};
    std::vector<void*> nr2ptr{};
  protected:
    bool shallow_to_python = false;
    std::map<std::string, VersionInfo> version_map = GetLibraryVersions();
  public:
    template<typename T>
      static constexpr bool is_archivable = detail::is_Archivable_struct<T>::value;

    Archive() = delete;
    Archive(const Archive&) = delete;
    Archive(Archive&&) = delete;
    Archive (bool ais_output) : is_output(ais_output) { ; }

    virtual ~Archive() { ; }

    // If the object is pickled, all shallow archived objects will be pickled as a list,
    // instead of written as a binary archive. This allows pickle to serialize every object only
    // once and put them together correctly afterwards. Therefore all objects that may live in
    // Python should be archived using this Shallow function. If Shallow is called from C++ code
    // it archives the object normally.
#ifdef NETGEN_PYTHON
    template<typename T>
    Archive& Shallow(T& val); // implemented in python_ngcore.hpp
#else // NETGEN_PYTHON
    template<typename T>
    Archive& Shallow(T& val)
    {
      static_assert(detail::is_any_pointer<T>, "ShallowArchive must be given pointer type!");
        *this & val;
      return *this;
    }
#endif // NETGEN_PYTHON

#ifdef NETGEN_PYTHON
    virtual void ShallowOutPython(const pybind11::object& /*unused*/)
    { throw UnreachableCodeException{}; }
    virtual void ShallowInPython(pybind11::object &)
    { throw UnreachableCodeException{}; }
#endif // NETGEN_PYTHON

    Archive& operator=(const Archive&) = delete;
    Archive& operator=(Archive&&) = delete;

    bool Output () const { return is_output; }
    bool Input () const { return !is_output; }
    const VersionInfo& GetVersion(const std::string& library)
    { return version_map[library]; }

    // only used for PyArchive
    virtual void NeedsVersion(const std::string& /*unused*/, const std::string& /*unused*/) {}

    // Pure virtual functions that have to be implemented by In-/OutArchive
    virtual Archive & operator & (std::byte & d) = 0;    
    virtual Archive & operator & (float & d) = 0;
    virtual Archive & operator & (double & d) = 0;
    virtual Archive & operator & (int & i) = 0;
    virtual Archive & operator & (long & i) = 0;
    virtual Archive & operator & (size_t & i) = 0;
    virtual Archive & operator & (short & i) = 0;
    virtual Archive & operator & (unsigned char & i) = 0;
    virtual Archive & operator & (bool & b) = 0;
    virtual Archive & operator & (std::string & str) = 0;
    virtual Archive & operator & (char *& str) = 0;

    Archive & operator & (VersionInfo & version)
    {
        if(Output())
            (*this) << version.to_string();
        else
        {
            std::string s;
            (*this) & s;
            version = VersionInfo(s);
        }
        return *this;
    }

    // Archive std classes ================================================
    template<typename T>
    Archive& operator & (std::complex<T>& c)
    {
      if(Output())
          (*this) << c.real() << c.imag();
      else
        {
          T tmp;
          (*this) & tmp;
          c.real(tmp);
          (*this) & tmp;
          c.imag(tmp);
        }
      return (*this);
    }
    template<typename T>
    Archive& operator & (std::vector<T>& v)
    {
      size_t size;
      if(Output())
          size = v.size();
      (*this) & size;
      if(Input())
        v.resize(size);
      Do(&v[0], size);
      return (*this);
    }
 
    // archive implementation for enums
    template<typename T>
    auto operator & (T& val) -> std::enable_if_t<std::is_enum<T>::value, Archive&>
    {
      int enumval;
      if(Output())
        enumval = int(val);
      *this & enumval;
      if(Input())
        val = T(enumval);
      return *this;
    }

    // vector<bool> has special implementation (like a bitarray) therefore
    // it needs a special overload (this could probably be more efficient, but we
    // don't use it that often anyway)
    Archive& operator& (std::vector<bool>& v)
    {
      size_t size;
      if(Output())
        size = v.size();
      (*this) & size;
      if(Input())
        {
          v.resize(size);
          bool b;
          for(size_t i=0; i<size; i++)
            {
              (*this) & b;
              v[i] = b;
            }
        }
      else
        {
          for(bool b : v)
            (*this) & b;
        }
      return *this;
    }
    template<typename T1, typename T2>
    Archive& operator& (std::map<T1, T2>& map)
    {
      if(Output())
        {
          (*this) << size_t(map.size());
          for(auto& pair : map)
              (*this) << pair.first << pair.second;
        }
      else
        {
          size_t size = 0;
          (*this) & size;
          T1 key; T2 val;
          for(size_t i = 0; i < size; i++)
            {
              T1 key; T2 val;
              (*this) & key & val;
              map[key] = val;
            }
        }
      return (*this);
    }
    template<typename T>
    Archive& operator& (std::optional<T>& opt)
    {
        bool has_value = opt.has_value();
        (*this) & has_value;
        if(has_value)
          {
            if(Output())
                (*this) << *opt;
            else
            {
                T value;
                (*this) & value;
                opt = value;
            }
          }
      return (*this);
    }
    template <typename T>
    Archive& operator&(std::set<T> &s)
    {
      auto size = s.size();
      (*this) & size;
      if(Output())
        for(const auto & val : s)
          (*this) << val;
      else
      {
          for(size_t i=0; i<size; i++)
          {
              T val;
              (*this) & val;
              s.insert(val);
          }
      }
      return *this;
    }

    // Archive arrays =====================================================
    // this functions can be overloaded in Archive implementations for more efficiency
    template <typename T, typename = std::enable_if_t<is_archivable<T>>>
    Archive & Do (T * data, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & data[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (std::byte * d, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & d[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (double * d, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & d[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (int * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (long * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (size_t * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (short * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (unsigned char * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; }; // NOLINT

    virtual Archive & Do (bool * b, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & b[j]; }; return *this; }; // NOLINT

    // Archive a class implementing a (void DoArchive(Archive&)) method =======
    template<typename T, typename=std::enable_if_t<detail::has_DoArchive<T>::value>>
    Archive& operator & (T& val)
    {
      val.DoArchive(*this); return *this;
    }



    
    // pack elements to binary
    template <typename ... Types>
      Archive & DoPacked (Types & ... args)
    {
      if (true) // (isbinary)
        {
          constexpr size_t totsize = TotSize<Types...>(); // (args...);
          std::byte mem[totsize];
          if (is_output)
            {
              CopyToBin (&mem[0], args...);
              Do(&mem[0], totsize);
            }
          else
            {
              Do(&mem[0], totsize);
              CopyFromBin (&mem[0], args...);
            }
        }
      // else
      // cout << "DoPacked of non-binary called --> individual pickling" << endl;
      return *this;
    }
    
    
    template <typename T, typename ... Trest>
      constexpr void CopyToBin (std::byte * ptr, T & first, Trest & ...rest) const
    {
      memcpy (ptr, &first, sizeof(first));
      CopyToBin(ptr+sizeof(first), rest...);
    }
    constexpr void CopyToBin (std::byte * ptr) const { }
    
    template <typename T, typename ... Trest>
      constexpr void CopyFromBin (std::byte * ptr, T & first, Trest & ...rest) const
    {
      memcpy (&first, ptr, sizeof(first));
      CopyFromBin(ptr+sizeof(first), rest...);
    }
    constexpr void CopyFromBin (std::byte * ptr) const { }


      

    template <typename T>
    Archive& operator & (ngcore::Shallow<T>& shallow)
    {
      this->Shallow(shallow.val);
      return *this;
    }
      

    // Archive shared_ptrs =================================================
    template <typename T>
    Archive& operator & (std::shared_ptr<T>& ptr)
    {
      if constexpr(has_shallow_archive<T>::value)
        if (shallow_to_python)
          {
            Shallow (ptr);
            return *this;
          }
          
      if(Output())
        {
          // save -2 for nullptr
          if(!ptr)
            return (*this) << -2;

          void* reg_ptr = ptr.get();
          bool neededDowncast = false;
          // Downcasting is only possible for our registered classes
          if(typeid(T) != typeid(*ptr))
            {
              if(!IsRegistered(Demangle(typeid(*ptr).name())))
                  throw Exception(std::string("Archive error: Polymorphic type ")
                                  + Demangle(typeid(*ptr).name())
                                  + " not registered for archive");
              reg_ptr = GetArchiveRegister(Demangle(typeid(*ptr).name())).downcaster(typeid(T), ptr.get());
              // if there was a true downcast we have to store more information
              if(reg_ptr != static_cast<void*>(ptr.get()))
                neededDowncast = true;
            }
          auto pos = shared_ptr2nr.find(reg_ptr);
          // if not found store -1 and the pointer
          if(pos == shared_ptr2nr.end())
            {
              auto p = ptr.get();
              (*this) << -1;
              (*this) & neededDowncast & p;
              // if we did downcast we store the true type as well
              if(neededDowncast)
                (*this) << Demangle(typeid(*ptr).name());
              shared_ptr2nr[reg_ptr] = shared_ptr_count++;
              return *this;
            }
          // if found store the position and if it has to be downcasted and how
          (*this) << pos->second << neededDowncast;
          if(neededDowncast)
            (*this) << Demangle(typeid(*ptr).name());
        }
      else // Input
        {
          int nr;
          (*this) & nr;
          // -2 restores a nullptr
          if(nr == -2)
            {
              ptr = nullptr;
              return *this;
            }
          // -1 restores a new shared ptr by restoring the inner pointer and creating a shared_ptr to it
          if (nr == -1)
            {
              T* p = nullptr;
              bool neededDowncast;
              (*this) & neededDowncast & p;
              ptr = std::shared_ptr<T>(p);
              // if we did downcast we need to store a shared_ptr<void> to the true object
              if(neededDowncast)
                {
                  std::string name;
                  (*this) & name;
                  auto info = GetArchiveRegister(name);
                  // for this we use an aliasing constructor to create a shared pointer sharing lifetime
                  // with our shared ptr, but pointing to the true object
                  nr2shared_ptr.push_back(std::shared_ptr<void>(std::static_pointer_cast<void>(ptr),
                                                                info.downcaster(typeid(T),
                                                                                ptr.get())));
                }
              else
                  nr2shared_ptr.push_back(ptr);
            }
          else
            {
              auto other = nr2shared_ptr[nr];
              bool neededDowncast;
              (*this) & neededDowncast;
              if(neededDowncast)
                {
                  // if there was a downcast we can expect the class to be registered (since archiving
                  // wouldn't have worked else)
                  std::string name;
                  (*this) & name;
                  auto info = GetArchiveRegister(name);
                  // same trick as above, create a shared ptr sharing lifetime with
                  // the shared_ptr<void> in the register, but pointing to our object
                  ptr = std::static_pointer_cast<T>(std::shared_ptr<void>(other,
                                                                          info.upcaster(typeid(T),
                                                                               other.get())));
                }
              else
                {
                  ptr = std::static_pointer_cast<T>(other);
                }
            }
        }
      return *this;
    }

    // Archive pointers =======================================================
    template <typename T>
    Archive & operator& (T *& p)
    {
      if (Output())
        {
          // if the pointer is null store -2
          if (!p)
              return (*this) << -2;
          auto reg_ptr = static_cast<void*>(p);
          if(typeid(T) != typeid(*p))
            {
              if(!IsRegistered(Demangle(typeid(*p).name())))
                throw Exception(std::string("Archive error: Polymorphic type ")
                                + Demangle(typeid(*p).name())
                                + " not registered for archive");
              reg_ptr = GetArchiveRegister(Demangle(typeid(*p).name())).downcaster(typeid(T), static_cast<void*>(p));
            }
          auto pos = ptr2nr.find(reg_ptr);
          // if the pointer is not found in the map create a new entry
          if (pos == ptr2nr.end())
            {
              ptr2nr[reg_ptr] = ptr_count++;
              if(typeid(*p) == typeid(T))
                if (std::is_constructible<T>::value)
                  return (*this) << -1 & (*p);
                else
                  {
                    if (IsRegistered(Demangle(typeid(*p).name())))
                    {
                      (*this) << -3 << Demangle(typeid(*p).name());
                      GetArchiveRegister(Demangle(typeid(*p).name())).
                        cargs_archiver(*this, p);
                      return (*this) & (*p);
                    }
                    else
                      throw Exception(std::string("Archive error: Class ") +
                                      Demangle(typeid(*p).name()) + " does not provide a default constructor!");
                  }
              else
                {
                  // if a pointer to a base class is archived, the class hierarchy must be registered
                  // to avoid compile time issues we allow this behaviour only for "our" classes that
                  // implement a void DoArchive(Archive&) member function
                  // To recreate the object we need to store the true type of it
                  if(!IsRegistered(Demangle(typeid(*p).name())))
                    throw Exception(std::string("Archive error: Polymorphic type ")
                                    + Demangle(typeid(*p).name())
                                    + " not registered for archive");
                  (*this) << -3 << Demangle(typeid(*p).name());
                  GetArchiveRegister(Demangle(typeid(*p).name())).
                    cargs_archiver(*this, p);
                  return (*this) & (*p);
                }
            }
          else
            {
              (*this) & pos->second;
              bool downcasted = !(reg_ptr == static_cast<void*>(p) );
              // store if the class has been downcasted and the name
              (*this) << downcasted << Demangle(typeid(*p).name());
            }
        }
      else
        {
          int nr;
          (*this) & nr;
          if (nr == -2) // restore a nullptr
              p = nullptr;
          else if (nr == -1) // create a new pointer of standard type (no virtual or multiple inheritance,...)
            {
              p = detail::constructIfPossible<T>();
              nr2ptr.push_back(p);
              (*this) & *p;
            }
          else if(nr == -3) // restore one of our registered classes that can have multiple inheritance,...
            {
              // As stated above, we want this special behaviour only for our classes that implement DoArchive
              std::string name;
              (*this) & name;
              auto info = GetArchiveRegister(name);
              // the creator creates a new object of type name, and returns a void* pointing
              // to T (which may have an offset)
              p = static_cast<T*>(info.creator(typeid(T), *this));
              // we store the downcasted pointer (to be able to find it again from
              // another class in a multiple inheritance tree)
              nr2ptr.push_back(info.downcaster(typeid(T),p));
              (*this) & *p;
            }
          else
            {
              bool downcasted;
              std::string name;
              (*this) & downcasted & name;
              if(downcasted)
                {
                  // if the class has been downcasted we can assume it is in the register
                  auto info = GetArchiveRegister(name);
                  p = static_cast<T*>(info.upcaster(typeid(T), nr2ptr[nr]));
                }
              else
                p = static_cast<T*>(nr2ptr[nr]);
            }
        }
      return *this;
    }

    Archive& operator&(std::tuple<>&) { return *this; }

    template <typename... T>
    Archive& operator&(std::tuple<T...> &t)
    {
      // call operator& for each element of the tuple
      std::apply([this](auto&... arg) { std::make_tuple(((*this) & arg).IsParallel()...);}, t);
      return *this;
    }

    // const ptr
    template<typename T>
    Archive& operator &(const T*& t)
    {
      return (*this) & const_cast<T*&>(t); // NOLINT
    }

    // Write a read only variable
    template <typename T>
    Archive & operator << (const T & t)
    {
      T ht(t);
      (*this) & ht;
      return *this;
    }

    virtual void FlushBuffer() {}

    bool parallel = false;
    bool IsParallel() const { return parallel; }
    void SetParallel (bool _parallel) { parallel = _parallel; }
    
  private:
  template<typename T, typename Bases>
    friend class RegisterClassForArchive;

#ifdef NETGEN_PYTHON
    friend pybind11::object CastAnyToPy(const std::any&);
#endif // NETGEN_PYTHON

    // Returns ClassArchiveInfo of Demangled typeid
    static const detail::ClassArchiveInfo& GetArchiveRegister(const std::string& classname);
    // Set ClassArchiveInfo for Demangled typeid, this is done by creating an instance of
    // RegisterClassForArchive<type, bases...>
    static void SetArchiveRegister(const std::string& classname, const detail::ClassArchiveInfo& info);
    static bool IsRegistered(const std::string& classname);

    // Helper class for up-/downcasting
    template<typename T, typename ... Bases>
    struct Caster{};

    template<typename T>
    struct Caster<T, std::tuple<>>
    {
      static void* tryUpcast (const std::type_info& /*unused*/, T* /*unused*/)
      {
        throw Exception("Upcast not successful, some classes are not registered properly for archiving!");
      }
      static void* tryDowncast (const std::type_info& /*unused*/, void* /*unused*/)
      {
        throw Exception("Downcast not successful, some classes are not registered properly for archiving!");
      }
    };

    template<typename T, typename B1>
    struct Caster<T,B1>
    {
      static void* tryUpcast(const std::type_info& ti, T* p)
      {
        try {
          return GetArchiveRegister(Demangle(typeid(B1).name()))
            .upcaster(ti, static_cast<void *>(dynamic_cast<B1 *>(p)));
        } catch (const Exception &) {
        throw Exception("Upcast not successful, some classes are not "
                                "registered properly for archiving!");
        }
      }

      static void* tryDowncast(const std::type_info& ti, void* p)
      {
        if(typeid(B1) == ti)
          return dynamic_cast<T*>(static_cast<B1*>(p));
        try
          {
            return dynamic_cast<T*>(static_cast<B1*>(GetArchiveRegister(Demangle(typeid(B1).name())).
                                                     downcaster(ti, p)));
        } catch (const Exception &) {
            throw Exception("Downcast not successful, some classes are not "
                            "registered properly for archiving!");
        }
      }
    };

    template<typename T, typename B1, typename ... Brest>
    struct Caster<T,std::tuple<B1, Brest...>>
    {
      static void* tryUpcast(const std::type_info& ti, T* p)
      {
        try
          { return GetArchiveRegister(Demangle(typeid(B1).name())).
              upcaster(ti, static_cast<void*>(dynamic_cast<B1*>(p))); }
        catch(const Exception&)
          { return Caster<T, std::tuple<Brest...>>::tryUpcast(ti, p); }
      }

      static void* tryDowncast(const std::type_info& ti, void* p)
      {
        if(typeid(B1) == ti)
          return dynamic_cast<T*>(static_cast<B1*>(p));
        try
          {
            return dynamic_cast<T*>(static_cast<B1*>(GetArchiveRegister(Demangle(typeid(B1).name())).
                                                     downcaster(ti, p)));
          }
        catch(const Exception&)
          {
            return Caster<T, std::tuple<Brest...>>::tryDowncast(ti, p);
          }
      }
    };
  };

  // BinaryOutArchive ======================================================================
  class NGCORE_API BinaryOutArchive : public Archive
  {
    static constexpr size_t BUFFERSIZE = 1024;
    std::array<char,BUFFERSIZE> buffer{};
    size_t ptr = 0;
  protected:
    std::shared_ptr<std::ostream> stream;
  public:
    BinaryOutArchive() = delete;
    BinaryOutArchive(const BinaryOutArchive&) = delete;
    BinaryOutArchive(BinaryOutArchive&&) = delete;
    BinaryOutArchive(std::shared_ptr<std::ostream>&& astream)
      : Archive(true), stream(std::move(astream))
    { }
    BinaryOutArchive(const std::filesystem::path& filename)
      : BinaryOutArchive(std::make_shared<std::ofstream>(filename)) {}
    ~BinaryOutArchive () override { FlushBuffer(); }

    BinaryOutArchive& operator=(const BinaryOutArchive&) = delete;
    BinaryOutArchive& operator=(BinaryOutArchive&&) = delete;

    using Archive::operator&;
    Archive & operator & (std::byte & d) override
    { return Write(d); }
    Archive & operator & (float & f) override
    { return Write(f); }
    Archive & operator & (double & d) override
    { return Write(d); }
    Archive & operator & (int & i) override
    { return Write(i); }
    Archive & operator & (short & i) override
    { return Write(i); }
    Archive & operator & (long & i) override
    {
      // for platform independence
      if constexpr (sizeof(long) == 8)
        return Write(i);
      else
        return Write(static_cast<int64_t>(i));
    }
    Archive & operator & (size_t & i) override
    {
      // for platform independence
      if constexpr (sizeof(size_t) == 8)
        return Write(i);
      else
        return Write(static_cast<uint64_t>(i));
    }
    Archive & operator & (unsigned char & i) override
    { return Write(i); }
    Archive & operator & (bool & b) override
    { return Write(b); }

    Archive & operator & (std::string & str) override
    {
      int len = str.length();
      (*this) & len;
      FlushBuffer();
      if(len)
        stream->write (&str[0], len);
      return *this;
    }
    Archive & operator & (char *& str) override
    {
      long len = str ? static_cast<long>(strlen (str)) : -1;
      (*this) & len;
      FlushBuffer();
      if(len > 0)
        stream->write (&str[0], len); // NOLINT
      return *this;
    }
    void FlushBuffer() override
    {
      if (ptr > 0)
        {
          stream->write(&buffer[0], ptr);
          ptr = 0;
        }
    }
    Archive & Do (std::byte * d, size_t n) override
    {
      FlushBuffer();
      stream->write(reinterpret_cast<char*>(d), n*sizeof(std::byte)); return *this;
    } 

  private:
    template <typename T>
    Archive & Write (T x)
    {
      static_assert(sizeof(T) < BUFFERSIZE, "Cannot write large types with this function!");
      if (unlikely(ptr > BUFFERSIZE-sizeof(T)))
        {
          stream->write(&buffer[0], ptr);
          ptr = 0;
        }
      memcpy(&buffer[ptr], &x, sizeof(T));
      ptr += sizeof(T);
      return *this;
    }
  };

  // BinaryInArchive ======================================================================
  class NGCORE_API BinaryInArchive : public Archive
  {
  protected:
    std::shared_ptr<std::istream> stream;
  public:
    BinaryInArchive (std::shared_ptr<std::istream>&& astream)
      : Archive(false), stream(std::move(astream))
    { }
    BinaryInArchive (const std::filesystem::path& filename)
      : BinaryInArchive(std::make_shared<std::ifstream>(filename)) { ; }

    using Archive::operator&;
    Archive & operator & (std::byte & d) override
    { Read(d); return *this; }
    Archive & operator & (float & f) override
    { Read(f); return *this; }
    Archive & operator & (double & d) override
    { Read(d); return *this; }
    Archive & operator & (int & i) override
    { Read(i); return *this; }
    Archive & operator & (short & i) override
    { Read(i); return *this; }
    Archive & operator & (long & i) override
    {
      // for platform independence
      if constexpr (sizeof(long) == 8)
        Read(i);
      else
      {
        int64_t tmp = 0;
        Read(tmp);
        i = tmp;
      }
      return *this;
    }
    Archive & operator & (size_t & i) override
    {
      // for platform independence
      if constexpr (sizeof(long) == 8)
        Read(i);
      else
      {
        uint64_t tmp = 0;
        Read(tmp);
        i = tmp;
      }
      return *this;
    }
    Archive & operator & (unsigned char & i) override
    { Read(i); return *this; }
    Archive & operator & (bool & b) override
    { Read(b); return *this; }
    Archive & operator & (std::string & str) override
    {
      int len;
      (*this) & len;
      str.resize(len);
      if(len)
        stream->read(&str[0], len); // NOLINT
      return *this;
    }
    Archive & operator & (char *& str) override
    {
      long len;
      (*this) & len;
      if(len == -1)
        str = nullptr;
      else
        {
          str = new char[len+1]; // NOLINT
          stream->read(&str[0], len); // NOLINT
          str[len] = '\0'; // NOLINT
        }
      return *this;
    }

    Archive & Do (std::byte * d, size_t n) override
    { stream->read(reinterpret_cast<char*>(d), n*sizeof(std::byte)); return *this; } // NOLINT
    Archive & Do (double * d, size_t n) override
    { stream->read(reinterpret_cast<char*>(d), n*sizeof(double)); return *this; } // NOLINT
    Archive & Do (int * i, size_t n) override
    { stream->read(reinterpret_cast<char*>(i), n*sizeof(int)); return *this; } // NOLINT
    Archive & Do (size_t * i, size_t n) override
    {
      // for platform independence
      if constexpr (sizeof(long) == 8)
        stream->read(reinterpret_cast<char*>(i), n*sizeof(size_t)); // NOLINT
      else
        for(size_t j = 0; j < n; j++)
          (*this) & i[j];
      return *this;
    }

  private:
    template<typename T>
    inline void Read(T& val)
    { stream->read(reinterpret_cast<char*>(&val), sizeof(T)); } // NOLINT
  };

  // TextOutArchive ======================================================================
  class NGCORE_API TextOutArchive : public Archive
  {
  protected:
    std::shared_ptr<std::ostream> stream;
  public:
    TextOutArchive (std::shared_ptr<std::ostream>&& astream)
      : Archive(true), stream(std::move(astream))
    { }
    TextOutArchive (const std::filesystem::path& filename) :
      TextOutArchive(std::make_shared<std::ofstream>(filename)) { }

    using Archive::operator&;
    Archive & operator & (std::byte & d) override
    { *stream << int(d) << ' '; return *this; }
    Archive & operator & (float & f) override
    { *stream << f << '\n'; return *this; }
    Archive & operator & (double & d) override
    { *stream << d << '\n'; return *this; }
    Archive & operator & (int & i) override
    { *stream << i << '\n'; return *this; }
    Archive & operator & (short & i) override
    { *stream << i << '\n'; return *this; }
    Archive & operator & (long & i) override
    { *stream << i << '\n'; return *this; }
    Archive & operator & (size_t & i) override
    { *stream << i << '\n'; return *this; }
    Archive & operator & (unsigned char & i) override
    { *stream << int(i) << '\n'; return *this; }
    Archive & operator & (bool & b) override
    { *stream << (b ? 't' : 'f') << '\n'; return *this; }
    Archive & operator & (std::string & str) override
    {
      int len = str.length();
      *stream << len << '\n';
      if(len)
        {
          stream->write(&str[0], len); // NOLINT
          *stream << '\n';
        }
      return *this;
    }
    Archive & operator & (char *& str) override
    {
      long len = str ? static_cast<long>(strlen (str)) : -1;
      *this & len;
      if(len > 0)
        {
          stream->write (&str[0], len); // NOLINT
          *stream << '\n';
        }
      return *this;
    }
  };

  // TextInArchive ======================================================================
  class NGCORE_API TextInArchive : public Archive
  {
  protected:
    std::shared_ptr<std::istream> stream;
  public:
    TextInArchive (std::shared_ptr<std::istream>&& astream) :
      Archive(false), stream(std::move(astream))
    { }
    TextInArchive (const std::filesystem::path& filename)
      : TextInArchive(std::make_shared<std::ifstream>(filename)) {}

    using Archive::operator&;
    Archive & operator & (std::byte & d) override
    { int tmp; *stream >> tmp; d = std::byte(tmp); return *this; }
    Archive & operator & (float & f) override
    { *stream >> f; return *this; }
    Archive & operator & (double & d) override
    { *stream >> d; return *this; }
    Archive & operator & (int & i) override
    { *stream >> i; return *this; }
    Archive & operator & (short & i) override
    { *stream >> i; return *this; }
    Archive & operator & (long & i) override
    { *stream >> i; return *this; }
    Archive & operator & (size_t & i) override
    { *stream >> i; return *this; }
    Archive & operator & (unsigned char & i) override
    { int _i; *stream >> _i; i = _i; return *this; }
    Archive & operator & (bool & b) override
    { char c; *stream >> c; b = (c=='t'); return *this; }
    Archive & operator & (std::string & str) override
    {
      // Ignore \r (carriage return) characters when reading strings
      // this is necessary for instance when a file was written on Windows and is read on Unix

      int len;
      *stream >> len;
      char ch;
      stream->get(ch); // read newline character
      if(ch == '\r') // windows line endings -> read \n as well
        stream->get(ch);
      str.resize(len);
      if(len)
        stream->get(&str[0], len+1, '\0');

      // remove all \r characters from the string, check if size changed
      // if so, read the remaining characters
      str.erase(std::remove(str.begin(), str.end(), '\r'), str.cend());
      size_t chars_to_read = len-str.size();
      while (chars_to_read>0)
      {
        auto old_size = str.size();
        str.resize(len);

        stream->get(&str[old_size], chars_to_read+1, '\0');
        str.erase(std::remove(str.begin()+old_size, str.end(), '\r'), str.cend());
        chars_to_read = len - str.size();
      }
      return *this;
    }
    Archive & operator & (char *& str) override
    {
      long len;
      (*this) & len;
      char ch;
      if(len == -1)
        {
          str = nullptr;
          return (*this);
        }
      str = new char[len+1]; // NOLINT
      if(len)
        {
          stream->get(ch); // \n
          if(ch == '\r') // windows line endings, read \n as well
            stream->get(ch);
          stream->get(&str[0], len+1, '\0'); // NOLINT
        }
      str[len] = '\0'; // NOLINT
      return *this;
    }
  };

  // HashArchive =================================================================
  // This class enables to easily create hashes for archivable objects by xoring
  // threw its data

  class NGCORE_API HashArchive : public Archive
  {
    size_t hash_value = 0;
    char* h;
    int offset = 0;
  public:
    HashArchive() : Archive(true)
      { h = (char*)&hash_value; }

    using Archive::operator&;
    Archive & operator & (std::byte & d) override { return ApplyHash(d); }    
    Archive & operator & (float & f) override { return ApplyHash(f); }
    Archive & operator & (double & d) override { return ApplyHash(d); }
    Archive & operator & (int & i) override { return ApplyHash(i); }
    Archive & operator & (short & i) override { return ApplyHash(i); }
    Archive & operator & (long & i) override { return ApplyHash(i); }
    Archive & operator & (size_t & i) override { return ApplyHash(i); }
    Archive & operator & (unsigned char & i) override { return ApplyHash(i); }
    Archive & operator & (bool & b) override { return ApplyHash(b); }
    Archive & operator & (std::string & str) override
    { for(auto c : str) ApplyHash(c);  return *this; }
    Archive & operator & (char *& str) override
    { char* s = str; while(*s != '\0') ApplyHash(*(s++)); return *this; }

    // HashArchive can be used in const context
    template<typename T>
      Archive & operator& (const T& val) const
    { return (*this) & const_cast<T&>(val); }

    size_t GetHash() const { return hash_value; }

  private:
    template<typename T>
      Archive& ApplyHash(T val)
    {
      size_t n = sizeof(T);
      char* pval = (char*)&val;
      for(size_t i = 0; i < n; i++)
        {
          h[offset++] ^= pval[i];
          offset %= 8;
        }
      return *this;
    }
  };

} // namespace ngcore

#endif // NETGEN_CORE_ARCHIVE_HPP

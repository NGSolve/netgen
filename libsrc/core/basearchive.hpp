#ifndef NG_BASEARCHIVE_HPP
#define NG_BASEARCHIVE_HPP

namespace ngcore
{
  class Archive;
  std::string demangle(const char* typeinfo);

  // create new pointer of type T if it is default constructible, else throw
  template<typename T>
  T* constructIfPossible_impl(...)
  { throw std::runtime_error(std::string(demangle(typeid(T).name())) + " is not default constructible!"); }

  template<typename T, typename= typename std::enable_if<std::is_constructible<T>::value>::type>
  T* constructIfPossible_impl(int) { return new T; }

  template<typename T>
  T* constructIfPossible() { return constructIfPossible_impl<T>(int{}); }

  // Type trait to check if a class implements a 'void DoArchive(Archive&)' function
  template<typename T>
  struct has_DoArchive
  {
  private:
    template<typename T2>
    static constexpr auto check(T2*) ->
      typename std::is_same<decltype(std::declval<T2>().DoArchive(std::declval<Archive&>())),void>::type;
    template<typename>
    static constexpr std::false_type check(...);
    typedef decltype(check<T>(0)) type;
  public:
    static constexpr bool value = type::value;
  };

  // Info stored by registering a class using the RegisterClassForArchive struct in the map
  // stored in GetArchiveRegister
  struct ClassArchiveInfo
  {
    // create new object of this type and return a void* pointer that is points to the location
    // of the (base)class given by type_info
    std::function<void*(const std::type_info&)> creator;
    // This caster takes a void* pointer to the type stored in this info and casts it to a
    // void* pointer pointing to the (base)class type_info
    std::function<void*(const std::type_info&, void*)> upcaster;
    // This caster takes a void* pointer to the (base)class type_info and returns void* pointing
    // to the type stored in this info
    std::function<void*(const std::type_info&, void*)> downcaster;
  };

  // Returns a map of from the mangled typeids to the ClassArchiveInfo
  std::map<std::string, ClassArchiveInfo>& GetArchiveRegister();

  // Helper class for up-/downcasting
  template<typename T, typename ... Bases>
  struct Caster
  {
    static void* tryUpcast(const std::type_info& ti, T* p);
    static void* tryDowncast(const std::type_info& ti, void* p);
  };

  // Base Archive class
  class Archive
  {
    bool is_output;
    // how many different shared_ptr/pointer have been (un)archived
    int shared_ptr_count, ptr_count;
    // maps for archived shared pointers and pointers
    std::map<void*, int> shared_ptr2nr, ptr2nr;
    // vectors for storing the unarchived (shared) pointers
    std::vector<std::shared_ptr<void>> nr2shared_ptr;
    std::vector<void*> nr2ptr;
  public:
    Archive (bool ais_output) : is_output(ais_output), shared_ptr_count(0), ptr_count(0) { ; }
    virtual ~Archive() { ; }

    bool Output () { return is_output; }
    bool Input () { return !is_output; }

    // Pure virtual functions that have to be implemented by In-/OutArchive
    virtual Archive & operator & (double & d) = 0;
    virtual Archive & operator & (int & i) = 0;
    virtual Archive & operator & (long & i) = 0;
    virtual Archive & operator & (size_t & i) = 0;
    virtual Archive & operator & (short & i) = 0;
    virtual Archive & operator & (unsigned char & i) = 0;
    virtual Archive & operator & (bool & b) = 0;
    virtual Archive & operator & (std::string & str) = 0;
    virtual Archive & operator & (char *& str) = 0;


    // Archive std classes ================================================
    template<typename T>
    Archive& operator & (std::complex<T>& c)
    {
      if(is_output)
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
      if(is_output)
          size = v.size();
      (*this) & size;
      if(!is_output)
        v.resize(size);
      Do(&v[0], size);
      return (*this);
    }
    // Archive arrays =====================================================
    // this functions can be overloaded in Archive implementations for more efficiency
    template <typename T>
    Archive & Do (T * data, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & data[j]; }; return *this; };

    virtual Archive & Do (double * d, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & d[j]; }; return *this; };

    virtual Archive & Do (int * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (long * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (size_t * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (short * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (unsigned char * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (bool * b, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & b[j]; }; return *this; };

    // Archive a class implementing a (void DoArchive(Archive&)) method =======
    template<typename T, typename=std::enable_if_t<has_DoArchive<T>::value>>
    Archive& operator & (T& val)
    {
      val.DoArchive(*this); return *this;
    }

    // Archive shared_ptrs =================================================
    template <typename T>
    Archive& operator & (std::shared_ptr<T>& ptr)
    {
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
              if(GetArchiveRegister().count(demangle(typeid(*ptr).name())) == 0)
                  throw std::runtime_error(std::string("Archive error: Polymorphic type ")
                                           + demangle(typeid(*ptr).name())
                                           + " not registered for archive");
              reg_ptr = GetArchiveRegister()[demangle(typeid(*ptr).name())].downcaster(typeid(T), ptr.get());
              // if there was a true downcast we have to store more information
              if(reg_ptr != (void*) ptr.get())
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
                (*this) << demangle(typeid(*ptr).name());
              shared_ptr2nr[reg_ptr] = shared_ptr_count++;
              return *this;
            }
          // if found store the position and if it has to be downcasted and how
          (*this) << pos->second << neededDowncast;
          if(neededDowncast)
            (*this) << demangle(typeid(*ptr).name());
          return (*this);
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
          else if (nr == -1)
            {
              T* p;
              bool neededDowncast;
              (*this) & neededDowncast & p;
              ptr = std::shared_ptr<T>(p);
              // if we did downcast we need to store a shared_ptr<void> to the true object
              if(neededDowncast)
                {
                  std::string name;
                  (*this) & name;
                  auto info = GetArchiveRegister()[name];
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
                  auto info = GetArchiveRegister()[name];
                  // same trick as above, create a shared ptr sharing lifetime with
                  // the shared_ptr<void> in the register, but pointing to our object
                  ptr = std::static_pointer_cast<T>(std::shared_ptr<void>(other,
                                                                          info.upcaster(typeid(T),
                                                                               other.get())));
                }
              else
                ptr = std::static_pointer_cast<T>(other);
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
          void* reg_ptr = (void*)p;
          if(typeid(T) != typeid(*p))
            {
              if(GetArchiveRegister().count(demangle(typeid(*p).name())) == 0)
                throw std::runtime_error(std::string("Archive error: Polymorphic type ")
                                         + demangle(typeid(*p).name())
                                         + " not registered for archive");
              else
                reg_ptr = GetArchiveRegister()[demangle(typeid(*p).name())].downcaster(typeid(T), (void*) p);
            }
          auto pos = ptr2nr.find(reg_ptr);
          // if the pointer is not found in the map create a new entry
          if (pos == ptr2nr.end())
            {
              ptr2nr[reg_ptr] = ptr_count++;
              if(typeid(*p) == typeid(T))
                if (std::is_constructible<T>::value)
                               {
                                 return (*this) << -1 & (*p);
                               }
                else
                  throw std::runtime_error(std::string("Archive error: Class ") +
                                           demangle(typeid(*p).name()) + " does not provide a default constructor!");
              else
                {
                  // if a pointer to a base class is archived, the class hierarchy must be registered
                  // to avoid compile time issues we allow this behaviour only for "our" classes that
                  // implement a void DoArchive(Archive&) member function
                  // To recreate the object we need to store the true type of it
                  if(GetArchiveRegister().count(demangle(typeid(*p).name())) == 0)
                    throw std::runtime_error(std::string("Archive error: Polymorphic type ")
                                             + demangle(typeid(*p).name())
                                             + " not registered for archive");
                  else
                    return (*this) << -3 << demangle(typeid(*p).name()) & (*p);
                }
            }
          else
            {
              (*this) & pos->second;
              bool downcasted = !(reg_ptr == (void*) p);
              // store if the class has been downcasted and the name
              (*this) << downcasted << demangle(typeid(*p).name());
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
              p = constructIfPossible<T>();
              nr2ptr.push_back(p);
              (*this) & *p;
            }
          else if(nr == -3) // restore one of our registered classes that can have multiple inheritance,...
            {
              // As stated above, we want this special behaviour only for our classes that implement DoArchive
              std::string name;
              (*this) & name;
              auto info = GetArchiveRegister()[name];
              // the creator creates a new object of type name, and returns a void* pointing
              // to T (which may have an offset)
              p = (T*) info.creator(typeid(T));
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
                    auto info = GetArchiveRegister()[name];
                    p = (T*) info.upcaster(typeid(T), nr2ptr[nr]);
                }
              else
                p = (T*) nr2ptr[nr];
            }
        }
      return *this;
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
  };

  template<typename T, typename ... Bases>
  class RegisterClassForArchive
  {
  public:
    RegisterClassForArchive()
    {
      ClassArchiveInfo info;
      info.creator = [this,&info](const std::type_info& ti) -> void*
                     { return typeid(T) == ti ? constructIfPossible<T>()
                         : Caster<T, Bases...>::tryUpcast(ti, constructIfPossible<T>()); };
      info.upcaster = [this](const std::type_info& ti, void* p) -> void*
                    { return typeid(T) == ti ? p : Caster<T, Bases...>::tryUpcast(ti, (T*) p); };
      info.downcaster = [this](const std::type_info& ti, void* p) -> void*
                        { return typeid(T) == ti ? p : Caster<T, Bases...>::tryDowncast(ti, p); };
      GetArchiveRegister()[std::string(demangle(typeid(T).name()))] = info;
    }
  };

  template<typename T>
  struct Caster<T>
  {
    static void* tryUpcast (const std::type_info& ti, T* p)
    {
      throw std::runtime_error("Upcast not successful, some classes are not registered properly for archiving!");
    }
    static void* tryDowncast (const std::type_info& ti, void* p)
    {
      throw std::runtime_error("Downcast not successful, some classes are not registered properly for archiving!");
    }
  };

  template<typename T, typename B1, typename ... Brest>
  struct Caster<T,B1,Brest...>
  {
    static void* tryUpcast(const std::type_info& ti, T* p)
    {
      try
        { return GetArchiveRegister()[demangle(typeid(B1).name())].upcaster(ti, (void*) (dynamic_cast<B1*>(p))); }
      catch(std::exception)
        { return Caster<T, Brest...>::tryUpcast(ti, p); }
    }

    static void* tryDowncast(const std::type_info& ti, void* p)
    {
      if(typeid(B1) == ti)
        return dynamic_cast<T*>((B1*) p);
      try
        { return GetArchiveRegister()[demangle(typeid(B1).name())].downcaster(ti, (void*) ((B1*)p)); }
      catch(std::exception)
        { return Caster<T, Brest...>::tryDowncast(ti, p); }
    }
  };
}

#endif // NG_BASEARCHIVE_HPP

#ifndef NG_BASEARCHIVE_HPP
#define NG_BASEARCHIVE_HPP

namespace ngcore
{
  class Archive;

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
          if constexpr(has_DoArchive<T>::value)
                        {
                          if(GetArchiveRegister().count(std::string(typeid(*ptr).name())) == 0)
                            throw std::runtime_error(std::string("Archive error: Polymorphic type ")
                                                     + typeid(*ptr).name()
                                                     + " not registered for archive");
                          else
                            reg_ptr = GetArchiveRegister()[typeid(*ptr).name()].downcaster(typeid(T), ptr.get());
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
              if(neededDowncast)
                (*this) << std::string(typeid(*ptr).name());
              shared_ptr2nr[reg_ptr] = shared_ptr_count++;
              return *this;
            }
          // if found store the position and if it has to be downcasted and how
          (*this) << pos->second << neededDowncast;
          if(neededDowncast)
              (*this) << std::string(typeid(*ptr).name());
          return (*this);
        }
      else // Input
        {
          int nr;
          (*this) & nr;
          if(nr == -2)
            {
              ptr = nullptr;
              return *this;
            }
          else if (nr == -1)
            {
              T* p;
              bool neededDowncast;
              (*this) & neededDowncast & p;
              ptr = std::shared_ptr<T>(p);
              if(neededDowncast)
                {
                  std::string name;
                  (*this) & name;
                  auto info = GetArchiveRegister()[name];
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
                  if constexpr(has_DoArchive<T>::value)
                                {
                                  std::string name;
                                  (*this) & name;
                                  auto info = GetArchiveRegister()[name];
                                  ptr = std::static_pointer_cast<T>(std::shared_ptr<void>(other,
                                                                                          info.upcaster(typeid(T),
                                                                                               other.get())));
                                }
                  else
                    throw std::runtime_error("Shouldn't get here...");
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
            {
              int m2 = -2;
              (*this) & m2;
              return *this;
            }
          void* reg_ptr = (void*)p;
          if constexpr(has_DoArchive<T>::value)
                        {
                          if(GetArchiveRegister().count(std::string(typeid(*p).name())) == 0)
                            throw std::runtime_error(std::string("Archive error: Polimorphic type ")
                                                     + typeid(*p).name()
                                                     + " not registered for archive");
                          else
                            reg_ptr = GetArchiveRegister()[typeid(*p).name()].downcaster(typeid(T), p);
                        }
          auto pos = ptr2nr.find(reg_ptr);
          // if the pointer is not found in the map create a new entry
          if (pos == ptr2nr.end())
            {
              ptr2nr[reg_ptr] = ptr_count++;
              if(typeid(*p) == typeid(T))
                if constexpr (std::is_constructible<T>::value)
                               {
                                 return (*this) << -1 & (*p);
                               }
                else
                  throw std::runtime_error(std::string("Archive error: Class ") +
                                           typeid(*p).name() + " does not provide a default constructor!");
              else
                {
                  // We want this special behaviour only for our classes that implement DoArchive
                  if constexpr(has_DoArchive<T>::value)
                                {
                                  if(GetArchiveRegister().count(std::string(typeid(*p).name())) == 0)
                                    throw std::runtime_error(std::string("Archive error: Polimorphic type ")
                                                             + typeid(*p).name()
                                                             + " not registered for archive");
                                  else
                                    return (*this) << -3 << std::string(typeid(*p).name()) & (*p);
                                }
                  else
                    throw std::runtime_error(std::string("Archive error: Class ")
                                             + typeid(*p).name()
                                             + " is polymorphic but not registered for archiving");
                }
            }
          else
            {
              (*this) & pos->second;
              (*this) << std::string(typeid(*p).name());
            }
        }
      else
        {
          int nr;
          (*this) & nr;
          if (nr == -2)
            {
              p = nullptr;
            }
          else if (nr == -1)
            {
              if constexpr (std::is_constructible<T>::value)
                             {
                               p = new T;
                               nr2ptr.push_back(p);
                               (*this) & *p;
                             }
              else
                throw std::runtime_error("Class isn't registered properly");

            }
          else if(nr == -3)
            {
              // We want this special behaviour only for our classes that implement DoArchive
              if constexpr(has_DoArchive<T>::value)
                            {
                              std::string name;
                              (*this) & name;
                              auto info = GetArchiveRegister()[name];
                              p = (T*) info.creator(typeid(T));
                              nr2ptr.push_back(info.downcaster(typeid(T),p));
                              (*this) & *p;
                            }
              else
                throw std::runtime_error("Class isn't registered properly");
            }
          else
            {
              std::string name;
              (*this) & name;
              if constexpr(has_DoArchive<T>::value)
                            {
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
      static_assert(std::is_constructible_v<T>, "Class registered for archive must be default constructible");
      ClassArchiveInfo info;
      info.creator = [this,&info](const std::type_info& ti) -> void*
                     { return typeid(T) == ti ? new T : Caster<T, Bases...>::tryUpcast(ti, new T); };
      info.upcaster = [this](const std::type_info& ti, void* p) -> void*
                    { return typeid(T) == ti ? p : Caster<T, Bases...>::tryUpcast(ti, (T*) p); };
      info.downcaster = [this](const std::type_info& ti, void* p) -> void*
                        { return typeid(T) == ti ? p : Caster<T, Bases...>::tryDowncast(ti, p); };
      GetArchiveRegister()[std::string(typeid(T).name())] = info;
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
        { return GetArchiveRegister()[typeid(B1).name()].upcaster(ti, (void*) (dynamic_cast<B1*>(p))); }
      catch(std::exception)
        { return Caster<T, Brest...>::tryUpcast(ti, p); }
    }

    static void* tryDowncast(const std::type_info& ti, void* p)
    {
      if(typeid(B1) == ti)
        return dynamic_cast<T*>((B1*) p);
      try
        { return GetArchiveRegister()[typeid(B1).name()].downcaster(ti, (void*) ((B1*)p)); }
      catch(std::exception)
        { return Caster<T, Brest...>::tryDowncast(ti, p); }
    }
  };
}

#endif // NG_BASEARCHIVE_HPP

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

  std::map<std::string, std::function<void*()>>& GetArchiveRegister();

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
        v.reserve(size);
      Do(&v[0], size);
      return (*this);
    }
    // Archive arrays =====================================================
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
          auto pos = shared_ptr2nr.find((void*) ptr.get());
          // if not found store -1 and the pointer
          if(pos == shared_ptr2nr.end())
            {
              shared_ptr2nr[(void*) ptr.get()] = shared_ptr_count++;
              auto p = ptr.get();
              return (*this) << -1 & p;
            }
          // if found store the position
          return (*this) << pos->second;
        }
      else // Input
        {
          int nr;
          (*this) & nr;
          if(nr == -2)
            ptr = nullptr;
          else if (nr == -1)
            {
              T* p;
              (*this) & p;
              ptr = std::shared_ptr<T>(p);
              nr2shared_ptr.push_back(ptr);
            }
          else
            ptr = std::reinterpret_pointer_cast<T>(nr2shared_ptr[nr]);
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
          auto pos = ptr2nr.find( (void*) p);
          // if the pointer is not found in the map create a new entry
          if (pos == ptr2nr.end())
            {
              ptr2nr[(void*) p] = ptr_count++;
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
                                  if(GetArchiveRegister().count(typeid(*p).name()) == 0)
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
            }
        }
      else
        {
          int nr;
          (*this) & nr;
          // cout << "in, got nr " << nr << endl;
          if (nr == -2)
            {
              p = nullptr;
            }
          else if (nr == -1)
            {
              if constexpr (std::is_constructible<T>::value)
                             {
                               p = new T;
                               // cout << "create new ptr, p = " << p << endl;
                               (*this) & *p;
                               nr2ptr.push_back(p);
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
                              p = reinterpret_cast<T*>(GetArchiveRegister()[name]());
                              nr2ptr.push_back(p);
                            }
              else
                throw std::runtime_error("Class isn't registered properly");
            }
          else
            {
              p = (T*)nr2ptr[nr];
              // cout << "reuse ptr " << nr << ": " << p << endl;
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
  };

  template<typename T>
  void RegisterClassForArchive()
  {
    static_assert(std::is_constructible_v<T>, "Class registered for archive must be default constructible");
    GetArchiveRegister()[std::string(typeid(T).name())] = []() -> void* { return new T; };
  }
}

#endif // NG_BASEARCHIVE_HPP

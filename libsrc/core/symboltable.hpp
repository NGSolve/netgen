#ifndef NETGEN_CORE_SYMBOLTABLE_HPP
#define NETGEN_CORE_SYMBOLTABLE_HPP

#include <ostream>
#include <string>
#include <vector>

#include "exception.hpp"
#include "ngcore_api.hpp"

namespace ngcore
{
  /**
      A symbol table.

      The symboltable provides a mapping from string identifiers
      to the generic type T. The strings are copied.
      Complexity by name access is linear, by index is constant.
  */
  template <class T>
  class SymbolTable
  {
    std::vector<std::string> names;
    std::vector<T> data;
  public:
    using value_type = T;
    using reference = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;

    /// Creates a symboltable
    SymbolTable () = default;
    SymbolTable (const SymbolTable<T> &) = default;
    SymbolTable (SymbolTable<T> &&) noexcept = default;

    ~SymbolTable() = default;

    SymbolTable& operator=(const SymbolTable<T>&) = default;
    SymbolTable& operator=(SymbolTable<T>&&) = default;

    template<typename ARCHIVE>
    auto DoArchive(ARCHIVE& ar)
      -> typename std::enable_if_t<ARCHIVE::template is_archivable<T>, void>
    {
      ar & names & data;
    }

    /// INDEX of symbol name, throws exception if unused
    size_t Index (std::string_view name) const
    {
      for (size_t i = 0; i < names.size(); i++)
        if (names[i] == name) return i;
      throw RangeException("SymbolTable", name);
    }

    /// Index of symbol name, returns -1 if unused
    int CheckIndex (std::string_view name) const
    {
      for (int i = 0; i < names.size(); i++)
        if (names[i] == name) return i;
      return -1;
    }

    /// number of identifiers
    size_t Size() const
    {
      return data.size();
    }

    /// Returns reference to element. exception for unused identifier
    reference operator[] (std::string_view name)
    {
      return data[Index (name)];
    }

    const_reference operator[] (std::string_view name) const
    {
      return data[Index (name)];
    }

    /// Returns reference to i-th element, range check only in debug build
    reference operator[] (size_t i)
    {
      NETGEN_CHECK_RANGE(i, 0, data.size());
      return data[i];
    }

    /// Returns const reference to i-th element, range check only in debug build
    const_reference operator[] (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, data.size());
      return data[i];
    }

    /// Returns name of i-th element, range check only in debug build
    const std::string & GetName (size_t i) const
    {
      NETGEN_CHECK_RANGE(i, 0, names.size());
      return names[i];
    }

    /// Associates el to the string name, overrides if name is used
    void Set (std::string_view name, const T & el)
    {
      int i = CheckIndex (name);
      if (i >= 0)
        data[i] = el;
      else
        {
          data.push_back(el);
          names.push_back(std::string(name));
        }
    }


    
    /*
    bool Used (const std::string & name) const
    {
      return CheckIndex(name) >= 0;
    }
    */
    
    bool Used (std::string_view name) const
    {
      return CheckIndex(name) >= 0;
    }
    
    /// Deletes symboltable
    inline void DeleteAll ()
    {
      names.clear();
      data.clear();
    }

    // Adds all elements from other symboltable
    SymbolTable<T>& Update(const SymbolTable<T>& tbl2)
    {
      for (size_t i = 0; i < tbl2.Size(); i++)
        Set (tbl2.GetName(i), tbl2[i]);
      return *this;
    }
  };

  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const SymbolTable<T> & st)
  {
    for (int i = 0; i < st.Size(); i++)
      ost << st.GetName(i) << " : " << st[i] << std::endl;
    return ost;
  }
} // namespace ngcore

#endif // NETGEN_CORE_SYMBOLTABLE_HPP

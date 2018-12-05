#ifndef NG_ARCHIVE_HPP
#define NG_ARCHIVE_HPP

namespace ngcore
{
  // BinaryOutArchive ======================================================================
  class BinaryOutArchive : public Archive
  {
    size_t ptr = 0;
    enum { BUFFERSIZE = 1024 };
    char buffer[BUFFERSIZE];
    std::shared_ptr<std::ostream> fout;
  public:
    BinaryOutArchive(std::shared_ptr<std::ostream> afout) : Archive(true), fout(afout)
    {
      (*this) & GetLibraryVersions();
    }
    BinaryOutArchive(std::string filename)
      : BinaryOutArchive(std::make_shared<std::ofstream>(filename)) {}
    virtual ~BinaryOutArchive () { FlushBuffer(); }

    const VersionInfo& getVersion(const std::string& library)
    { return GetLibraryVersions()[library]; }

    using Archive::operator&;
    virtual Archive & operator & (double & d)
    { return Write(d); }
    virtual Archive & operator & (int & i)
    { return Write(i); }
    virtual Archive & operator & (short & i)
    { return Write(i); }
    virtual Archive & operator & (long & i)
    { return Write(i); }
    virtual Archive & operator & (size_t & i)
    { return Write(i); }
    virtual Archive & operator & (unsigned char & i)
    { return Write(i); }
    virtual Archive & operator & (bool & b)
    { return Write(b); }
    virtual Archive & operator & (std::string & str)
    {
      if (ptr > 0) FlushBuffer();
      int len = str.length();
      fout->write (reinterpret_cast<char*>(&len), sizeof(int));
      fout->write (&str[0], len);
      return *this;
    }
    virtual Archive & operator & (char *& str)
    {
      if (ptr > 0) FlushBuffer();
      int len = strlen (str);
      fout->write (reinterpret_cast<char*>(&len), sizeof(int));
      fout->write (&str[0], len);
      return *this;
    }
    void FlushBuffer()
    {
      if (ptr > 0)
        {
          fout->write(&buffer[0], ptr);
          ptr = 0;
        }
    }

  private:
    template <typename T>
    Archive & Write (T x)
    {
      if (unlikely(ptr > BUFFERSIZE-sizeof(T)))
        {
          fout->write(&buffer[0], ptr);
          * (T*) (&buffer[0]) = x;
          ptr = sizeof(T);
          return *this;
        }
      * (T*) (&buffer[ptr]) = x;
      ptr += sizeof(T);
      return *this;
    }
  };

  // BinaryInArchive ======================================================================
  class BinaryInArchive : public Archive
  {
    std::map<std::string, VersionInfo> vinfo;
    std::shared_ptr<std::istream> fin;
  public:
    BinaryInArchive (std::shared_ptr<std::istream> afin) : Archive(false), fin(afin)
    {
      (*this) & vinfo;
    }
    BinaryInArchive (std::string filename)
      : BinaryInArchive(std::make_shared<std::ifstream>(filename)) { ; }

    const VersionInfo& getVersion(const std::string& library)
    { return vinfo[library]; }

    using Archive::operator&;
    virtual Archive & operator & (double & d)
    { Read(d); return *this; }
    virtual Archive & operator & (int & i)
    { Read(i); return *this; }
    virtual Archive & operator & (short & i)
    { Read(i); return *this; }
    virtual Archive & operator & (long & i)
    { Read(i); return *this; }
    virtual Archive & operator & (size_t & i)
    { Read(i); return *this; }
    virtual Archive & operator & (unsigned char & i)
    { Read(i); return *this; }
    virtual Archive & operator & (bool & b)
    { Read(b); return *this; }
    virtual Archive & operator & (std::string & str)
    {
      int len;
      Read(len);
      str.resize(len);
      fin->read(&str[0], len);
      return *this;
    }
    virtual Archive & operator & (char *& str)
    {
      int len;
      Read(len);
      str = new char[len+1];
      fin->read(&str[0], len);
      str[len] = '\0';
      return *this;
    }

    virtual Archive & Do (double * d, size_t n)
    { fin->read(reinterpret_cast<char*>(d), n*sizeof(double)); return *this; }
    virtual Archive & Do (int * i, size_t n)
    { fin->read(reinterpret_cast<char*>(i), n*sizeof(int)); return *this; }
    virtual Archive & Do (size_t * i, size_t n)
    { fin->read(reinterpret_cast<char*>(i), n*sizeof(size_t)); return *this; }

  private:
    template<typename T>
    inline void Read(T& val)
    { fin->read(reinterpret_cast<char*>(&val), sizeof(T)); }
  };

  // TextOutArchive ======================================================================
  class TextOutArchive : public Archive
  {
    std::shared_ptr<std::ostream> fout;
  public:
    TextOutArchive (std::shared_ptr<std::ostream> afout) : Archive(true), fout(afout)
    {
      (*this) & GetLibraryVersions();
    }
    TextOutArchive (std::string filename) :
      TextOutArchive(std::make_shared<std::ofstream>(filename.c_str())) { }

    const VersionInfo& getVersion(const std::string& library)
    { return GetLibraryVersions()[library]; }

    using Archive::operator&;
    virtual Archive & operator & (double & d)
    { *fout << d << '\n'; return *this; }
    virtual Archive & operator & (int & i)
    { *fout << i << '\n'; return *this; }
    virtual Archive & operator & (short & i)
    { *fout << i << '\n'; return *this; }
    virtual Archive & operator & (long & i)
    { *fout << i << '\n'; return *this; }
    virtual Archive & operator & (size_t & i)
    { *fout << i << '\n'; return *this; }
    virtual Archive & operator & (unsigned char & i)
    { *fout << int(i) << '\n'; return *this; }
    virtual Archive & operator & (bool & b)
    { *fout << (b ? 't' : 'f') << '\n'; return *this; }
    virtual Archive & operator & (std::string & str)
    {
      int len = str.length();
      *fout << len << '\n';
      if(len)
        {
          fout->write(&str[0], len);
          *fout << '\n';
        }
      return *this;
    }
    virtual Archive & operator & (char *& str)
    {
      int len = strlen (str);
      *fout << len << '\n';
      if(len)
        {
          fout->write (&str[0], len);
          *fout << '\n';
        }
      return *this;
    }
  };

  // TextInArchive ======================================================================
  class TextInArchive : public Archive
  {
    std::map<std::string, VersionInfo> vinfo;
    std::shared_ptr<std::istream> fin;
  public:
    TextInArchive (std::shared_ptr<std::istream> afin) : Archive(false), fin(afin)
    {
      (*this) & vinfo;
    }
    TextInArchive (std::string filename)
      : TextInArchive(std::make_shared<std::ifstream>(filename)) {}

    const VersionInfo& getVersion(const std::string& library)
    { return vinfo[library]; }

    using Archive::operator&;
    virtual Archive & operator & (double & d)
    { *fin >> d; return *this; }
    virtual Archive & operator & (int & i)
    { *fin >> i; return *this; }
    virtual Archive & operator & (short & i)
    { *fin >> i; return *this; }
    virtual Archive & operator & (long & i)
    { *fin >> i; return *this; }
    virtual Archive & operator & (size_t & i)
    { *fin >> i; return *this; }
    virtual Archive & operator & (unsigned char & i)
    { int _i; *fin >> _i; i = _i; return *this; }
    virtual Archive & operator & (bool & b)
    { char c; *fin >> c; b = (c=='t'); return *this; }
    virtual Archive & operator & (std::string & str)
    {
      int len;
      *fin >> len;
      char ch;
      fin->get(ch); // '\n'
      str.resize(len);
      if(len)
        fin->get(&str[0], len+1, '\0');
      return *this;
    }
    virtual Archive & operator & (char *& str)
    {
      int len;
      *fin >> len;
      char ch;
      fin->get(ch); // '\n'
      str = new char[len+1];
      if(len)
        fin->get(&str[0], len, '\0');
      str[len] = 0;
      return *this;
    }
  };
}
#endif // NG_ARCHIVE_HPP

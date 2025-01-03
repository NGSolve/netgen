#ifndef WRITEUSER
#define WRITEUSER

/**************************************************************************/
/* File:    writeuser.hh                                                  */
/* Authors: many                                                          */
/* Date:    10. Dec. 97                                                   */
/**************************************************************************/

#include <filesystem>
#include <functional>
#include <optional>

#include <meshing.hpp>

namespace netgen {

using namespace std::filesystem;

typedef std::function<void (const Mesh & mesh, const filesystem::path & filename)> FWrite;
typedef std::function<void (Mesh & mesh, const filesystem::path & filename)> FRead;
typedef std::function<bool (const filesystem::path & filename)> FTest;

struct UserFormatRegister {
  struct UserFormatEntry {
    string format;
    Array<string> extensions;
    optional<FRead> read;
    optional<FWrite> write;
  };
  static void Register(UserFormatEntry && entry);

  static const bool HaveFormat(string format);
  DLL_HEADER static const UserFormatEntry & Get(string format);

  template<typename TFunc>
  static void IterateFormats(TFunc func, bool need_read=false, bool need_write=false) {
    Array<string> import_formats;
    for(const auto & e: getFormats())
    if((!need_read || e.second.read) && (!need_write || e.second.write))
      import_formats.Append(e.second.format);
    QuickSort(import_formats);
    for(auto format : import_formats)
      func(Get(format));
  }

  DLL_HEADER static std::map<std::string, UserFormatEntry>& getFormats();
};

struct RegisterUserFormat {
  RegisterUserFormat(string format, Array<string> extensions, optional<FRead> read, optional<FWrite> write, FTest ftest = [](const filesystem::path & ){return true;})
  {
    UserFormatRegister::Register({format, std::move(extensions), std::move(read), std::move(write)});
  }
};

DLL_HEADER void ReadFile(Mesh & mesh, const filesystem::path & filename);
DLL_HEADER void ReadUserFormat(Mesh & mesh, const filesystem::path & filename, const string & format = "");


extern bool DLL_HEADER WriteUserFormat (const string & format,
                                        const Mesh & mesh,
                                        const filesystem::path & filename);

extern void DLL_HEADER RegisterUserFormats (NgArray<const char*> & names,
                                 NgArray<const char*> & extensions);

}

#endif


#include "ngcore_api.hpp"
#include "utils.hpp"
#include "logging.hpp"
#include "simd_generic.hpp"

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN
#else // WIN32
#include <cxxabi.h>
#include <dlfcn.h>
#endif //WIN32
//
#include <array>
#include <filesystem>
#include <iostream>
#include <regex>
#include <string>
#include <thread>

#include "ngstream.hpp"


namespace ngcore
{
    namespace detail
    {
        // see https://github.com/RobotLocomotion/drake/blob/master/common/nice_type_name.cc
        static const auto demangle_regexes =
            std::array<std::pair<std::regex, std::string>, 8>{
                // Remove unwanted keywords and following space. (\b is word boundary.)
                std::make_pair(std::regex("\\b(class|struct|enum|union) "), ""),
                // Tidy up anonymous namespace.
                {std::regex("[`(]anonymous namespace[')]"), "(anonymous)"},
                // Replace Microsoft __int64 with long long.
                {std::regex("\\b__int64\\b"), "long long"},
                // Temporarily replace spaces we want to keep with "!". (\w is
                // alphanumeric or underscore.)
                {std::regex("(\\w) (\\w)"), "$1!$2"},
                {std::regex(" "), ""},  // Delete unwanted spaces.
                // Some compilers throw in extra namespaces like "__1" or "__cxx11".
                // Delete them.
                {std::regex("\\b__[[:alnum:]_]+::"), ""},
                {std::regex("!"), " "},  // Restore wanted spaces.

                // Recognize std::string's full name and abbreviate.
                {std::regex("\\bstd::basic_string<char,std::char_traits<char>,"
                        "std::allocator<char>>"), "std::string"}
            };
        std::string CleanupDemangledName( std::string s )
        {
            for(const auto & [r, sub] : demangle_regexes)
                s = std::regex_replace (s,r,sub);
#ifdef EMSCRIPTEN
            // for some reason regex_replace is not working at all
            std::string temp = s;
            s = "";
            for(auto c : temp)
              if(c!=' ')
                s+=c;
#endif // EMSCRIPTEN

            return s;
        }
    } // namespace detail

  // parallel netgen
  int id = 0, ntasks = 1;

#ifdef WIN32
  // windows does demangling in typeid(T).name()
  NGCORE_API std::string Demangle(const char* typeinfo) {
      std::string name = typeinfo;
      return detail::CleanupDemangledName(name);
  }
#else
  NGCORE_API std::string Demangle(const char* typeinfo)
  {
    int status=0;
    try
      {
        char *s = abi::__cxa_demangle(typeinfo, nullptr, nullptr, &status);
        std::string result;
        if (s == nullptr)
          result = typeinfo;
        else
          {
            result = s;
            free(s);
          }
        result = detail::CleanupDemangledName(result);
        return result;
      }
    catch( const std::exception & e )
      {
        GetLogger("utils")->warn("{}:{} cannot demangle {}, status: {}, error:{}", __FILE__, __LINE__, typeinfo, status, e.what());
      }
    std::string name = typeinfo;
    return detail::CleanupDemangledName(name);
  }
#endif

  double seconds_per_tick = [] () noexcept
  {
      auto tick_start = GetTimeCounter();
      double tstart = WallTime();
      double tend = WallTime()+0.001;

      // wait for 1ms and compare wall time with time counter
      while(WallTime()<tend);

      auto tick_end = GetTimeCounter();
      tend = WallTime();

      return (tend-tstart)/static_cast<double>(tick_end-tick_start);
  }();

  const std::chrono::time_point<TClock> wall_time_start = TClock::now();

  int printmessage_importance = getenv("NG_MESSAGE_LEVEL") ? atoi(getenv("NG_MESSAGE_LEVEL")) : 0;
  bool NGSOStream :: glob_active = true;

  NGCORE_API int GetCompiledSIMDSize()
  {
      return GetDefaultSIMDSize();
  }

  NGCORE_API bool IsRangeCheckEnabled()
  {
#ifdef NETGEN_ENABLE_CHECK_RANGE
      return true;
#else
      return false;
#endif
  }

  NGCORE_API std::filesystem::path GetTempFilename()
  {
      static int counter = 0;
      auto path = std::filesystem::temp_directory_path();
      path += ".temp_netgen_file_"+ToString(counter++)+"_"+ToString(GetTimeCounter());
      return path;
  }


  SharedLibrary :: SharedLibrary(const std::filesystem::path & lib_name_, std::optional<std::filesystem::path> directory_to_delete_, bool global )
      : lib_name(lib_name_),directory_to_delete(directory_to_delete_)
  {
    Load(lib_name, global);
  }

  SharedLibrary :: ~SharedLibrary()
  {
    Unload();
    if(directory_to_delete)
      for([[maybe_unused]] auto i : Range(5))
      {
        // on Windows, a (detached?) child process of the compiler/linker might still block the directory
        // wait for it to finish (up to a second)
        try
        {
          std::filesystem::remove_all(*directory_to_delete);
          directory_to_delete = std::nullopt;
          break;
        }
        catch(const std::exception &e)
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
      }
    if(directory_to_delete)
      std::cerr << "Could not delete " << directory_to_delete->string() << std::endl;
  }

  void SharedLibrary :: Load( const std::filesystem::path & lib_name_, bool global )
  {
    Unload();
    lib_name = lib_name_;
#ifdef WIN32
    lib = LoadLibrary(lib_name.wstring().c_str());
    if (!lib) throw std::runtime_error(std::string("Could not load library ") + lib_name.string());
#else // WIN32
    auto flags = RTLD_NOW;
    if (global) flags |= RTLD_GLOBAL;
    lib = dlopen(lib_name.c_str(), flags);
    if(lib == nullptr) throw std::runtime_error(dlerror());
#endif // WIN32
  }

  void SharedLibrary :: Unload() {
    if(lib)
    {
#ifdef WIN32
      FreeLibrary((HMODULE)lib);
#else // WIN32
      int rc = dlclose(lib);
      if(rc != 0) std::cerr << "Failed to close library " << lib_name << std::endl;
#endif // WIN32
    }
  }

  void* SharedLibrary :: GetRawSymbol( std::string func_name )
  {
#ifdef WIN32
    void* func = GetProcAddress((HMODULE)lib, func_name.c_str());
    if(func == nullptr)
      throw std::runtime_error(std::string("Could not find function ") + func_name + " in library " + lib_name.string());
#else // WIN32
    void* func = dlsym(lib, func_name.c_str());
    if(func == nullptr)
        throw std::runtime_error(dlerror());
#endif // WIN32

    return func;
  }

  void* GetRawSymbol( std::string func_name )
  {
    void * func = nullptr;
#ifdef WIN32
    throw std::runtime_error("GetRawSymbol not implemented on WIN32");
#else // WIN32
    func = dlsym(RTLD_DEFAULT, func_name.c_str());
    if(func == nullptr)
        throw std::runtime_error(dlerror());
#endif // WIN32
    return func;
  }


} // namespace ngcore


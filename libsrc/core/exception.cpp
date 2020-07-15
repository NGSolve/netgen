#include "exception.hpp"
#include "utils.hpp"

namespace ngcore
{
  Exception :: Exception(const std::string& s)
    : m_what(s) {}
  
  Exception :: Exception(const char* s)
    : m_what(s) {}


  void ThrowException(const std::string & s)
  {
    throw Exception (s);
  }
  
  void ThrowException(const char * s)
  {
    throw Exception (s);
  }
} // namespace ngcore


// ********* STUFF FOR GETBACKTRACE ***************************
#ifdef __GNUC__

#include <execinfo.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <dlfcn.h>
#include <array>
#include <memory>
#include <cxxabi.h>
#include <signal.h>
#include <vector>

namespace ngcore
{
  namespace detail
  {
    static int exec(std::string cmd, std::string & out) {
      std::array<char, 128> buffer;
      FILE *pipe = popen(cmd.c_str(), "r");

      if (!pipe)
        throw std::runtime_error("popen() failed!");

      out = "";
      while (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
        out += buffer.data();

      int error_code = pclose(pipe);
      return error_code;
    }

#ifdef __APPLE__
    // Split output line from backtrace_symbols to recover function name and offset
    // then use `nm` command line tool to get the address of the function
    // then use `add42line` command line tool to map function address + offset to line in source code
    static std::string TranslateBacktrace( std::string s, std::string libname )
    {
      // example line
      // 1   libngcore.dylib                     0x000000010ddb298c _ZL21ngcore_signal_handleri + 316
      constexpr char reset_shell[] = "\033[0m";
      constexpr char green[] = "\033[32m";
      constexpr char yellow[] = "\033[33m";

      std::istringstream in(s);

      std::string libname1, funcname, addr, plus_sign;
      size_t i,offset;

      in >> i >> libname1 >> addr >> funcname >> plus_sign >> std::hex >> offset;

      std::stringstream out;

      if(!funcname.empty() && !libname.empty())
      {
        std::string nm_command = "nm " + libname + " | grep \"" + funcname + "$\" | cut -f 1 -d ' '";
        std::string output;
        auto exit_code = exec(nm_command, output);
        auto fptr = std::strtoul(output.c_str(), 0, 16);
        if(fptr == 0)
            return out.str()+'\n';
        std::stringstream offset_s;
        offset_s << "0x" << std::hex << fptr+offset - 5;
        std::string addr2line_command = std::string("atos -o ") + libname + " --fullPath " + offset_s.str();
        exit_code = exec(addr2line_command, output);
        if(exit_code==0)
          out << " at " << green << output << reset_shell;
        else
          out << '\n';
      }
      else
        out << s << '\n';

      return out.str();
    }
#else // __APPLE__

    // Split output line from backtrace_symbols to recover function name and offset
    // then use `nm` command line tool to get the address of the function
    // then use `addr2line` command line tool to map function address + offset to line in source code
    static std::string TranslateBacktrace( std::string s, std::string /*dummy*/ )
    {
      // example line:
      // /home/mhochsteger/install/ngs_clang/bin/../lib/libngcore.so(_ZN6ngcore11TaskManager4LoopEi+0x1e0) [0x7f2991fe1030]
      constexpr char reset_shell[] = "\033[0m";
      constexpr char green[] = "\033[32m";
      constexpr char yellow[] = "\033[33m";

      auto brace_open_pos = s.find('(');
      auto brace_close_pos = s.find(')', brace_open_pos);
      auto plus_pos = s.find('+', brace_open_pos);
      auto bracket_open_pos = s.find('[');
      auto bracket_close_pos = s.find(']');

      auto libname = s.substr(0, brace_open_pos);
      auto funcname = s.substr(brace_open_pos+1, plus_pos - brace_open_pos - 1);
      auto offset = std::strtoul(s.substr(plus_pos+1, brace_close_pos - plus_pos - 1).c_str(), 0, 16);
      auto position = std::strtoul(s.substr(bracket_open_pos+1, bracket_close_pos - bracket_open_pos - 1).c_str(), 0, 16);
      std::stringstream out;

      if(!funcname.empty())
      {
        std::vector<char> buffer(10240);
        int status;
        size_t size = buffer.size();
        abi::__cxa_demangle(funcname.c_str(), &buffer[0], &size, &status);
        out << "in " << yellow << &buffer[0] << reset_shell << '\n';

        std::string nm_command = "nm " + libname + " | grep " + funcname + " | cut -f 1 -d ' '";
        std::string output;
        auto exit_code = exec(nm_command, output);
        auto fptr = std::strtoul(output.c_str(), 0, 16);

        std::stringstream offset_s;
        offset_s << "0x" << std::hex << fptr+offset - 5;
        std::string addr2line_command = std::string("addr2line -i -p -e ") + libname + " " + offset_s.str();
        exit_code = exec(addr2line_command, output);
        if(exit_code==0)
          {
            std::stringstream soutput(output);
            std::string s;
            while(soutput)
              {
                if(getline(soutput, s))
                    out << "\t   at " << green << s << reset_shell << '\n';
              }
          }
        else
          out << '\n';
      }
      else
        out << s << '\n';

      return out.str();
    }
#endif // __APPLE__

  } // namespace detail

  std::string GetBackTrace()
  {
    if(!getenv("NG_BACKTRACE"))
        return "";
    std::cerr << "Collecting backtrace..." << std::endl;
    std::stringstream result;
    void *bt[100];
    int bt_size;
    char **bt_syms;
    int i;

    bt_size = backtrace(bt, 100);
    bt_syms = backtrace_symbols(bt, bt_size);
    Dl_info info;
    for (i = 1; i < bt_size-1; i++)
      {
        dladdr(bt[i], &info);
        size_t len = strlen(bt_syms[i]);
        result << '#'<< i << '\t' << detail::TranslateBacktrace( bt_syms[i], info.dli_fname );
      }
    free(bt_syms);
    return result.str();
  }

} // namespace ngcore

static void ngcore_signal_handler(int sig)
{
  static bool first_call = true;
  if(!first_call)
      exit(1); // avoid endless recursions if signals are caused by this handler
  first_call = false;

  switch(sig)
    {
      case SIGABRT:
          std::cerr << "Caught SIGABRT: usually caused by abort() or assert()" << std::endl;
          break;
      case SIGILL:
          std::cerr << "Caught SIGILL: illegal instruction" << std::endl;
          break;
      case SIGSEGV:
          std::cerr << "Caught SIGSEGV: segmentation fault" << std::endl;
          break;
    }

  std::cerr << ngcore::GetBackTrace() << std::endl;
  exit(1);
}

// register signal handler when library is loaded
static bool dummy = []()
{
    if(getenv("NG_BACKTRACE"))
    {
      signal(SIGABRT, ngcore_signal_handler);
      signal(SIGILL, ngcore_signal_handler);
      signal(SIGSEGV, ngcore_signal_handler);
    }
    return true;
}();

#else // __GNUC__

namespace ngcore
{
  std::string GetBackTrace()
  {
    return std::string();
  }
} // namespace ngcore

#endif // __GNUC__

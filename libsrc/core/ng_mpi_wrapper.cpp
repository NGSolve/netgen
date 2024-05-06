#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "ng_mpi.hpp"
#include "ngstream.hpp"
#include "utils.hpp"

using std::cerr;
using std::cout;
using std::endl;

namespace ngcore {

static std::unique_ptr<SharedLibrary> mpi_lib, ng_mpi_lib;

void InitMPI(std::filesystem::path mpi_lib_path,
             std::filesystem::path ng_libs_dir) {
  if (ng_mpi_lib) return;
  cout << IM(3) << "InitMPI" << endl;

  typedef void (*get_version_handle)(char *, int *);
  typedef int (*init_handle)(int *, char ***);
  typedef int (*mpi_initialized_handle)(int *);
  typedef void (*ng_init_handle)();

  init_handle mpi_init;
  mpi_initialized_handle mpi_initialized;
  get_version_handle get_version;

  mpi_lib = std::make_unique<SharedLibrary>(mpi_lib_path, std::nullopt, true);

  try {
    mpi_init = GetSymbol<init_handle>("MPI_Init");
    mpi_initialized = GetSymbol<mpi_initialized_handle>("MPI_Initialized");
    get_version = GetSymbol<get_version_handle>("MPI_Get_library_version");
  } catch (std::runtime_error &e) {
    cerr << "Could not load MPI symbols: " << e.what() << endl;
    throw e;
  }

  int flag = 0;
  mpi_initialized(&flag);
  if (!flag) {
    typedef const char *pchar;
    int argc = 1;
    pchar args[] = {"netgen", nullptr};
    pchar *argv = &args[0];
    cout << IM(5) << "Calling MPI_Init" << endl;
    mpi_init(&argc, (char ***)argv);
  }

  char version_string[65536];
  int resultlen = 0;
  get_version(version_string, &resultlen);
  mpi_library_version = version_string;
  cout << IM(7) << "MPI version: " << version_string << endl;

  std::string libname = "";
  if (mpi_library_version.substr(0, 8) == "Open MPI") {
    cout << IM(5) << "Have Open MPI" << endl;
    libname = std::string("libng_openmpi") + NETGEN_SHARED_LIBRARY_SUFFIX;
  } else if (mpi_library_version.substr(0, 5) == "MPICH") {
    cout << IM(5) << "Have MPICH" << endl;
    libname = std::string("libng_mpich") + NETGEN_SHARED_LIBRARY_SUFFIX;
  } else
    cerr << "Unknown MPI version, skipping init: " << version_string<< endl;

  if (libname.size()) {
    ng_mpi_lib = std::make_unique<SharedLibrary>(libname);
    auto ng_init = ng_mpi_lib->GetSymbol<ng_init_handle>("ng_init_mpi");
    ng_init();
  }
}

static std::runtime_error no_mpi() {
  return std::runtime_error("MPI not enabled");
}

std::string mpi_library_version = "";

#if defined(NG_PYTHON) && defined(NG_MPI4PY)
decltype(NG_MPI_CommFromMPI4Py) NG_MPI_CommFromMPI4Py =
    [](py::handle, NG_MPI_Comm &) -> bool { throw no_mpi(); };
decltype(NG_MPI_CommToMPI4Py) NG_MPI_CommToMPI4Py =
    [](NG_MPI_Comm) -> py::handle { throw no_mpi(); };
#endif

#include "ng_mpi_generated_dummy_init.hpp"

}  // namespace ngcore

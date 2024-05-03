#include <iostream>
#include <stdexcept>
#include <filesystem>

#include "ng_mpi.hpp"
#include "ngstream.hpp"
#include "utils.hpp"

using std::cout;
using std::endl;

namespace ngcore {

static std::unique_ptr<SharedLibrary> mpi_lib, ng_mpi_lib;

void InitMPI(std::filesystem::path mpi_lib_path, std::filesystem::path ng_libs_dir) {
  if (ng_mpi_lib) return;
  cout << IM(3) << "InitMPI" << endl;

  typedef void (*get_version_handle)(char *, int *);
  typedef int (*init_handle)(int *, char ***);
  typedef int (*mpi_initialized_handle)(int *);
  typedef void (*ng_init_handle)();

  init_handle mpi_init;
  mpi_initialized_handle mpi_initialized;
  get_version_handle get_version;
  try {
    mpi_init = GetSymbol<init_handle>("MPI_Init");
    cout << IM(3) << "MPI already loaded " << mpi_init << endl;
    mpi_initialized = GetSymbol<mpi_initialized_handle>("MPI_Initialized");
    get_version = GetSymbol<get_version_handle>("MPI_Get_library_version");
  } catch (std::runtime_error &e) {
    cout << IM(3) << "MPI not loaded" << endl;
    mpi_lib = std::make_unique<SharedLibrary>(mpi_lib_path,
                                              std::nullopt, true);
    mpi_init = mpi_lib->GetSymbol<init_handle>("MPI_Init");
    mpi_initialized =
        mpi_lib->GetSymbol<mpi_initialized_handle>("MPI_Initialized");
    get_version =
        mpi_lib->GetSymbol<get_version_handle>("MPI_Get_library_version");
  }

  int flag = 0;
  mpi_initialized(&flag);
  if (!flag) {
    cout << IM(3) << "Calling MPI_Init" << endl;
    mpi_init(nullptr, nullptr);
  }

  char version_string[65536];
  int resultlen = 0;
  get_version(version_string, &resultlen);
  mpi_library_version = version_string;
  cout << IM(3) << "MPI version: " << version_string << endl;

  std::string libname = "";
  if (mpi_library_version.substr(0, 8) == "Open MPI") {
    cout << IM(3) << "Have Open MPI" << endl;
    libname = std::string("libng_openmpi") + NETGEN_SHARED_LIBRARY_SUFFIX;
  } else if (mpi_library_version.substr(0, 5) == "MPICH") {
    cout << IM(3) << "Have MPICH" << endl;
    libname = std::string("libng_mpich.so") + NETGEN_SHARED_LIBRARY_SUFFIX;
  } else
    cout << IM(3) << "Unknown MPI" << endl;

  if (libname.size()) {
    cout << IM(3) << "loading " << libname << endl;
    cout << IM(3) << "NG_MPI_INT before " << NG_MPI_INT.value << endl;
    ng_mpi_lib = std::make_unique<SharedLibrary>(libname);
    auto ng_init = ng_mpi_lib->GetSymbol<ng_init_handle>("ng_init_mpi");
    cout << IM(3) << "have ng_init " << ng_init << endl;
    ng_init();
    cout << IM(3) << "NG_MPI_INT after  " << NG_MPI_INT.value << endl;

    int size, rank;
    NG_MPI_Comm_size(NG_MPI_COMM_WORLD, &size);
    NG_MPI_Comm_rank(NG_MPI_COMM_WORLD, &rank);
    cout << IM(3) << "Hello from " << rank << " of " << size << endl;
  }
}

static std::runtime_error no_mpi() {
  return std::runtime_error("MPI not enabled");
}

std::string mpi_library_version = "";

#include "ng_mpi_generated_dummy_init.hpp"

}  // namespace ngcore

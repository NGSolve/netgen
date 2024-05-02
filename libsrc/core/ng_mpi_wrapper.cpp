#include <mpi.h>

#include <iostream>
#include <stdexcept>

#include "ng_mpi.hpp"
#include "utils.hpp"

using std::cout;
using std::endl;

namespace ngcore {

static std::unique_ptr<SharedLibrary> mpi_lib, ng_mpi_lib;

void InitMPI() {
  if (ng_mpi_lib) return;
  cout << "InitMPI" << endl;

  typedef void (*get_version_handle)(char *, int *);
  typedef int (*init_handle)(int *, char ***);
  typedef int (*mpi_initialized_handle)(int *);
  typedef void (*ng_init_handle)();

  init_handle mpi_init;
  mpi_initialized_handle mpi_initialized;
  get_version_handle get_version;
  try {
    mpi_init = GetSymbol<init_handle>("MPI_Init");
    cout << "MPI already loaded " << mpi_init << endl;
    mpi_initialized = GetSymbol<mpi_initialized_handle>("MPI_Initialized");
    get_version = GetSymbol<get_version_handle>("MPI_Get_library_version");
  } catch (std::runtime_error &e) {
    cout << "MPI not loaded" << endl;
    mpi_lib = std::make_unique<SharedLibrary>("libmpi.so", std::nullopt, true);
    mpi_init = mpi_lib->GetSymbol<init_handle>("MPI_Init");
    mpi_initialized =
        mpi_lib->GetSymbol<mpi_initialized_handle>("MPI_Initialized");
    get_version =
        mpi_lib->GetSymbol<get_version_handle>("MPI_Get_library_version");
  }

  int flag = 0;
  mpi_initialized(&flag);
  if (!flag) {
    cout << "Calling MPI_Init" << endl;
    mpi_init(nullptr, nullptr);
  }

  char version_string[65536];
  int resultlen = 0;
  get_version(version_string, &resultlen);
  mpi_library_version = version_string;
  cout << "MPI version: " << version_string << endl;

  std::string libname = "";
  if (mpi_library_version.substr(0, 8) == "Open MPI") {
    cout << "Have Open MPI" << endl;
    libname = "/opt/netgen/lib/libng_openmpi.so";
  } else if (mpi_library_version.substr(0, 5) == "MPICH") {
    cout << "Have MPICH" << endl;
    libname = "/opt/netgen/lib/libng_mpich.so";
  } else
    cout << "Unknown MPI" << endl;

  if (libname.size()) {
    cout << "loading " << libname << endl;
    cout << "NG_MPI_INT before " << NG_MPI_INT.value << endl;
    ng_mpi_lib = std::make_unique<SharedLibrary>(libname);
    auto ng_init = ng_mpi_lib->GetSymbol<ng_init_handle>("ng_init_mpi");
    cout << "have ng_init " << ng_init << endl;
    ng_init();
    cout << "NG_MPI_INT after  " << NG_MPI_INT.value << endl;
  }
}

// NG_MPI_Comm NG_MPI_COMM_WORLD = 0;

// NG_MPI_Datatype NG_MPI_INT = 0;
// NG_MPI_Datatype NG_MPI_SHORT = 0;
// NG_MPI_Datatype NG_MPI_CHAR = 0;
// NG_MPI_Datatype NG_MPI_UINT64_T = 0;
// NG_MPI_Datatype NG_MPI_DOUBLE = 0;
// NG_MPI_Datatype NG_MPI_C_BOOL = 0;

static std::runtime_error no_mpi() {
  return std::runtime_error("MPI not enabled");
}

std::string mpi_library_version = "";

#include "ng_mpi_generated_dummy_init.hpp"

}  // namespace ngcore

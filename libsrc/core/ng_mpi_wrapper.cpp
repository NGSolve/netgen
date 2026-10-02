#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "ng_mpi.hpp"
#include "array.hpp"
#include "ngstream.hpp"
#ifdef NG_PYTHON
#include "python_ngcore.hpp"
#endif // NG_PYTHON
#include "utils.hpp"

using std::cerr;
using std::cout;
using std::endl;

namespace ngcore {

static std::unique_ptr<SharedLibrary> mpi_lib, ng_mpi_lib;
static bool need_mpi_finalize = false;

struct MPIFinalizer {
  ~MPIFinalizer() {
    if (need_mpi_finalize) {
      cout << IM(5) << "Calling MPI_Finalize" << endl;
      NG_MPI_Finalize();
    }
  }
} mpi_finalizer;

bool MPI_Loaded() { return ng_mpi_lib != nullptr; }

void InitMPI(std::optional<std::filesystem::path> mpi_lib_path) {
  if (ng_mpi_lib) return;

  cout << IM(3) << "InitMPI" << endl;

  std::string vendor = "";
  std::string mpi4py_lib_file = "";

  if (mpi_lib_path) {
    // Dynamic load of given shared MPI library
    // Then call MPI_Init, read the library version and set the vender name
    try {
      typedef int (*init_handle)(int *, char ***);
      typedef int (*mpi_initialized_handle)(int *);
      mpi_lib =
          std::make_unique<SharedLibrary>(*mpi_lib_path, std::nullopt, true);
      auto mpi_init = mpi_lib->GetSymbol<init_handle>("MPI_Init");
      auto mpi_initialized =
          mpi_lib->GetSymbol<mpi_initialized_handle>("MPI_Initialized");

      int flag = 0;
      mpi_initialized(&flag);
      if (!flag) {
        int argc = 1;
        char name[] = "netgen";
        char *args[] = {name, nullptr};
        char **argv = args;
        cout << IM(5) << "Calling MPI_Init" << endl;
        mpi_init(&argc, &argv);   // was passing argv instead of &argv
        need_mpi_finalize = true;
      }

      char c_version_string[65536];
      c_version_string[0] = '\0';
      int result_len = 0;
      typedef void (*get_version_handle)(char *, int *);
      auto get_version =
          mpi_lib->GetSymbol<get_version_handle>("MPI_Get_library_version");
      get_version(c_version_string, &result_len);
      vendor = c_version_string;   // library version string, searched for vendor names below
    } catch (std::runtime_error &e) {
      cerr << "Could not load MPI: " << e.what() << endl;
      throw e;
    }
  } else {
#ifdef NG_PYTHON
    // Use mpi4py to init MPI library and get the vendor name
    auto mpi4py = py::module::import("mpi4py.MPI");
    vendor = mpi4py.attr("get_vendor")()[py::int_(0)].cast<std::string>();

#ifndef WIN32
    // Load mpi4py library (it exports all MPI symbols) to have all MPI symbols
    // available before the ng_mpi wrapper is loaded This is not necessary on
    // windows as the matching mpi dll is linked to the ng_mpi wrapper directly
    mpi4py_lib_file = mpi4py.attr("__file__").cast<std::string>();
    mpi_lib =
        std::make_unique<SharedLibrary>(mpi4py_lib_file, std::nullopt, true);
    try {
      char c_version_string[65536];
      c_version_string[0] = '\0';
      int result_len = 0;
      typedef void (*get_version_handle)(char *, int *);
      mpi_lib->GetSymbol<get_version_handle>("MPI_Get_library_version")(c_version_string, &result_len);
      vendor += std::string(" / ") + c_version_string;
    } catch (std::runtime_error &) {
    }
#endif  // WIN32
#endif // NG_PYTHON
  }

  // The wrapper variant is chosen by ABI, not by vendor name: every MPI library
  // is either Open MPI-ABI (predefined handles are exported symbols) or
  // MPICH-ABI (predefined handles are integer constants: MPICH, Cray MPICH,
  // MVAPICH, Intel MPI, ...). Vendor names only select the dedicated variants.
  auto mentions = [&](const char *name) {
    return vendor.find(name) != std::string::npos;
  };
  Array<std::string> candidates;
  if (mentions("Microsoft MPI"))
    candidates.Append("ng_microsoft_mpi");
  else if (mentions("Intel(R) MPI") || mentions("Intel MPI"))
    candidates = {"ng_intel_mpi", "ng_mpich"};
  else if (mentions("Open MPI") || mentions("Spectrum MPI"))
    candidates.Append("ng_openmpi");
  else if (mentions("MPICH") || mentions("MVAPICH"))
    candidates.Append("ng_mpich");
  else {
    bool openmpi_abi = false;
    if (mpi_lib)
      try {
        mpi_lib->GetSymbol<void *>("ompi_mpi_comm_world");
        openmpi_abi = true;
      } catch (std::runtime_error &) {
      }
    candidates.Append(openmpi_abi ? "ng_openmpi" : "ng_mpich");
  }

  // Load the ng_mpi wrapper and call ng_init_mpi to set all function pointers
  typedef void (*ng_init_handle)();
  std::string errors;
  for (auto name : candidates) {
    std::string ng_lib_name = name + NETGEN_SHARED_LIBRARY_SUFFIX;
    try {
      ng_mpi_lib = std::make_unique<SharedLibrary>(ng_lib_name);
      break;
    } catch (std::runtime_error &e) {
      errors += std::string("\n  ") + ng_lib_name + ": " + e.what();
    }
  }
  if (!ng_mpi_lib)
    throw std::runtime_error("Could not load the MPI wrapper library for \"" + vendor +
                             "\" (is Netgen built with USE_MPI=ON?):" + errors);
  ng_mpi_lib->GetSymbol<ng_init_handle>("ng_init_mpi")();
  std::cout << IM(3) << "MPI wrapper loaded: " << candidates[0] << " for " << vendor << endl;
}

static std::runtime_error no_mpi() {
  return std::runtime_error("MPI not enabled");
}

#ifdef NG_PYTHON
decltype(NG_MPI_CommFromMPI4Py) NG_MPI_CommFromMPI4Py =
    [](py::handle py_obj, NG_MPI_Comm &ng_comm) -> bool {
  // If this gets called, it means that we want to convert an mpi4py
  // communicator to a Netgen MPI communicator, but the Netgen MPI wrapper
  // runtime was not yet initialized.

  // store the current address of this function
  auto old_converter_address = NG_MPI_CommFromMPI4Py;

  // initialize the MPI wrapper runtime, this sets all the function pointers
  InitMPI();

  // if the initialization was successful, the function pointer should have
  // changed
  // -> call the actual conversion function
  if (NG_MPI_CommFromMPI4Py != old_converter_address)
    return NG_MPI_CommFromMPI4Py(py_obj, ng_comm);

  // otherwise, something strange happened
  throw no_mpi();
};
decltype(NG_MPI_CommToMPI4Py) NG_MPI_CommToMPI4Py =
    [](NG_MPI_Comm) -> py::handle { throw no_mpi(); };
#endif  // NG_PYTHON

#include "ng_mpi_generated_dummy_init.hpp"

}  // namespace ngcore

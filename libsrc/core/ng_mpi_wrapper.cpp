#ifdef PARALLEL

#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "ng_mpi.hpp"
#include "ngstream.hpp"
#include "python_ngcore.hpp"
#include "utils.hpp"

using std::cerr;
using std::cout;
using std::endl;

#ifndef NG_MPI_WRAPPER
#define MPI4PY_LIMITED_API 1
#define MPI4PY_LIMITED_API_SKIP_MESSAGE 1
#define MPI4PY_LIMITED_API_SKIP_SESSION 1
#include "mpi4py_pycapi.h"  // mpi4py < 4.0.0
#endif // NG_MPI_WRAPPER

namespace ngcore {

#ifdef NG_MPI_WRAPPER
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
        typedef const char *pchar;
        int argc = 1;
        pchar args[] = {"netgen", nullptr};
        pchar *argv = &args[0];
        cout << IM(5) << "Calling MPI_Init" << endl;
        mpi_init(&argc, (char ***)argv);
        need_mpi_finalize = true;
      }

      char c_version_string[65536];
      c_version_string[0] = '\0';
      int result_len = 0;
      typedef void (*get_version_handle)(char *, int *);
      auto get_version =
          mpi_lib->GetSymbol<get_version_handle>("MPI_Get_library_version");
      get_version(c_version_string, &result_len);
      std::string version = c_version_string;

      if (version.substr(0, 8) == "Open MPI")
        vendor = "Open MPI";
      else if (version.substr(0, 5) == "MPICH")
        vendor = "MPICH";
      else if (version.substr(0, 13) == "Microsoft MPI")
        vendor = "Microsoft MPI";
      else if (version.substr(0, 12) == "Intel(R) MPI")
        vendor = "Intel MPI";
      else
        throw std::runtime_error(
            std::string("Unknown MPI version: " + version));
    } catch (std::runtime_error &e) {
      cerr << "Could not load MPI: " << e.what() << endl;
      throw e;
    }
  } else {
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
#endif  // WIN32
  }

  std::string ng_lib_name = "";
  if (vendor == "Open MPI")
    ng_lib_name = "ng_openmpi";
  else if (vendor == "MPICH")
    ng_lib_name = "ng_mpich";
  else if (vendor == "Microsoft MPI")
    ng_lib_name = "ng_microsoft_mpi";
  else if (vendor == "Intel MPI")
    ng_lib_name = "ng_intel_mpi";
  else
    throw std::runtime_error("Unknown MPI vendor: " + vendor);

  ng_lib_name += NETGEN_SHARED_LIBRARY_SUFFIX;

  // Load the ng_mpi wrapper and call ng_init_mpi to set all function pointers
  typedef void (*ng_init_handle)();
  ng_mpi_lib = std::make_unique<SharedLibrary>(ng_lib_name);
  ng_mpi_lib->GetSymbol<ng_init_handle>("ng_init_mpi")();
  std::cout << IM(3) << "MPI wrapper loaded, vendor: " << vendor << endl;
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
#else  // NG_MPI_WRAPPER

static bool imported_mpi4py = false;
#ifdef NG_PYTHON
decltype(NG_MPI_CommFromMPI4Py) NG_MPI_CommFromMPI4Py =
    [](py::handle src, NG_MPI_Comm &dst) -> bool {
  if (!imported_mpi4py) {
    import_mpi4py__MPI();
    imported_mpi4py = true;
  }
  PyObject *py_src = src.ptr();
  auto type = Py_TYPE(py_src);
  if (PyObject_TypeCheck(py_src, &PyMPIComm_Type)) {
    dst = *PyMPIComm_Get(py_src);
    return !PyErr_Occurred();
  }
  return false;
};

decltype(NG_MPI_CommToMPI4Py) NG_MPI_CommToMPI4Py =
    [](NG_MPI_Comm src) -> py::handle {
  if (!imported_mpi4py) {
    import_mpi4py__MPI();
    imported_mpi4py = true;
  }
  return py::handle(PyMPIComm_New(src));
};

#endif  // NG_PYTHON

bool MPI_Loaded() { return true; }
void InitMPI(std::optional<std::filesystem::path>) {}

#endif  // NG_MPI_WRAPPER

}  // namespace ngcore

#endif  // PARALLEL

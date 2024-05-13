#define OMPI_SKIP_MPICXX

#include "ng_mpi.hpp"

#include <mpi.h>

#include <type_traits>

#include "ngcore_api.hpp"
#include "pybind11/pytypes.h"

#if defined(NG_PYTHON) && defined(NG_MPI4PY)
#include <mpi4py.h>

#include "python_ngcore.hpp"

namespace py = pybind11;
#endif

#ifdef MSMPI_VER
int MPI_Comm_create_group(MPI_Comm arg0, MPI_Group arg1, int arg2,
                          MPI_Comm* arg3) {
  throw std::runtime_error(
      "MPI_Comm_create_group not supported on Microsoft MPI");
}
static MPI_Datatype MPI_CXX_DOUBLE_COMPLEX;
#endif  // MSMPI_VER

namespace ngcore {

static_assert(sizeof(MPI_Status) <= sizeof(NG_MPI_Status), "Size mismatch");
static_assert(alignof(MPI_Status) <= alignof(NG_MPI_Status), "Size mismatch");

int mpi2ng(int value) { return value; }
void* mpi2ng(void* ptr) { return ptr; }

NG_MPI_Status* mpi2ng(MPI_Status* status) {
  return reinterpret_cast<NG_MPI_Status*>(status);
}

#if !defined(MPICH) && !defined(MSMPI_VER)
NG_MPI_Comm mpi2ng(MPI_Comm comm) { return reinterpret_cast<uintptr_t>(comm); }
#endif

template <size_t size, size_t stride>
void gather_strided_array(size_t count, char* data) {
  static_assert(size <= stride, "Size must be less than or equal to stride");
  if constexpr (size < stride) {
    char* dst = data;
    char* src = data;
    for (auto i : Range(count)) {
      memcpy(dst, src, size);
      dst += size;
      src += stride;
    }
  }
}

template <typename T>
T cast_ng2mpi(uintptr_t obj) {
  if constexpr (std::is_pointer_v<T>)
    return reinterpret_cast<T>(obj);
  else
    return static_cast<T>(obj);
}

template <typename T>
T cast_ng2mpi(uintptr_t* ptr) {
  if constexpr (std::is_pointer_v<T>)
    return reinterpret_cast<T>(ptr);
  else
    return static_cast<T>(ptr);
}

template <typename T, typename TSrc>
T* cast_ng2mpi(TSrc* ptr, int count) {
  gather_strided_array<sizeof(T), sizeof(TSrc)>(count,
                                                reinterpret_cast<char*>(ptr));
  return reinterpret_cast<T*>(ptr);
}

MPI_Comm ng2mpi(NG_MPI_Comm comm) {
  static_assert(sizeof(MPI_Comm) <= sizeof(comm.value), "Size mismatch");
  static_assert(alignof(MPI_Comm) <= alignof(NG_MPI_Comm), "Size mismatch");
  return cast_ng2mpi<MPI_Comm>(comm.value);
}

MPI_Group ng2mpi(NG_MPI_Group group) {
  static_assert(sizeof(MPI_Group) <= sizeof(group.value), "Size mismatch");
  static_assert(alignof(MPI_Group) <= alignof(NG_MPI_Group), "Size mismatch");
  return cast_ng2mpi<MPI_Group>(group.value);
}

MPI_Comm* ng2mpi(NG_MPI_Comm* comm) {
  return cast_ng2mpi<MPI_Comm*>(&comm->value);
}
MPI_Group* ng2mpi(NG_MPI_Group* group) {
  return cast_ng2mpi<MPI_Group*>(&group->value);
}
MPI_Datatype* ng2mpi(NG_MPI_Datatype* type) {
  return cast_ng2mpi<MPI_Datatype*>(&type->value);
}
MPI_Datatype* ng2mpi(NG_MPI_Datatype* type, int count) {
  return cast_ng2mpi<MPI_Datatype>(&type->value, count);
}
MPI_Request* ng2mpi(NG_MPI_Request* request) {
  return cast_ng2mpi<MPI_Request*>(&request->value);
}
MPI_Request* ng2mpi(NG_MPI_Request* request, int count) {
  return cast_ng2mpi<MPI_Request>(&request->value, count);
}
MPI_Status* ng2mpi(NG_MPI_Status* status) {
  return reinterpret_cast<MPI_Status*>(status);
}
MPI_Aint* ng2mpi(NG_MPI_Aint* aint) {
  return reinterpret_cast<MPI_Aint*>(aint);
}
MPI_Aint* ng2mpi(NG_MPI_Aint* aint, int count) {
  return cast_ng2mpi<MPI_Aint>(aint, count);
}

MPI_Datatype ng2mpi(NG_MPI_Datatype type) {
  static_assert(sizeof(MPI_Datatype) <= sizeof(type.value), "Size mismatch");
  return cast_ng2mpi<MPI_Datatype>(type.value);
}

MPI_Request ng2mpi(NG_MPI_Request request) {
  static_assert(sizeof(MPI_Request) <= sizeof(request.value), "Size mismatch");
  return cast_ng2mpi<MPI_Request>(request.value);
}

MPI_Op ng2mpi(NG_MPI_Op op) {
  static_assert(sizeof(MPI_Op) <= sizeof(op.value), "Size mismatch");
  return cast_ng2mpi<MPI_Op>(op.value);
}

MPI_Aint ng2mpi(NG_MPI_Aint aint) {
  static_assert(sizeof(MPI_Aint) <= sizeof(aint.value), "Size mismatch");
  return cast_ng2mpi<MPI_Aint>(aint.value);
}

void* ng2mpi(void* ptr) { return ptr; }
char* ng2mpi(char* ptr) { return ptr; }
char*** ng2mpi(char*** ptr) { return ptr; }
int* ng2mpi(int* ptr) { return ptr; }
int ng2mpi(int value) { return value; }

}  // namespace ngcore

using namespace ngcore;

extern "C" {
NGCORE_API_EXPORT void ng_init_mpi();
}

static bool imported_mpi4py = false;

void ng_init_mpi() {
#if defined(NG_PYTHON) && defined(NG_MPI4PY)
  NG_MPI_CommFromMPI4Py = [](py::handle src, NG_MPI_Comm& dst) -> bool {
    if (!imported_mpi4py) {
      import_mpi4py();
      imported_mpi4py = true;
    }
    PyObject* py_src = src.ptr();
    auto type = Py_TYPE(py_src);
    if (PyObject_TypeCheck(py_src, &PyMPIComm_Type)) {
      dst = mpi2ng(*PyMPIComm_Get(py_src));
      return !PyErr_Occurred();
    }
    return false;
  };
  NG_MPI_CommToMPI4Py = [](NG_MPI_Comm src) -> py::handle {
    if (!imported_mpi4py) {
      import_mpi4py();
      imported_mpi4py = true;
    }
    return py::handle(PyMPIComm_New(ng2mpi(src)));
  };
#endif

#include "ng_mpi_generated_init.hpp"
}

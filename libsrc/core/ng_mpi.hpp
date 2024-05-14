#ifndef NG_MPI_HPP_INCLUDED
#define NG_MPI_HPP_INCLUDED

#ifdef PARALLEL

#include <cstdint>
#include <filesystem>
#include <optional>

#include "ngcore_api.hpp"

#if defined(NG_PYTHON) && defined(NG_MPI4PY)
#include <pybind11/pybind11.h>

namespace py = pybind11;
#endif

#ifndef NG_MPI_WRAPPER
#include <mpi.h>
#if defined(NG_PYTHON) && defined(NG_MPI4PY)
#include <mpi4py.h>
#endif
#endif  // NG_MPI_WRAPPER

namespace ngcore {

NGCORE_API bool MPI_Loaded();
NGCORE_API void InitMPI(
    std::optional<std::filesystem::path> mpi_lib_path = std::nullopt);

#ifdef NG_MPI_WRAPPER
inline void not_implemented() { throw std::runtime_error("Not implemented"); }

struct NG_MPI_Status {
  uintptr_t data[4];
};

struct NG_MPI_Comm {
  uintptr_t value;
  NG_MPI_Comm() { value = 0; }
  NG_MPI_Comm(uintptr_t value_) : value(value_) {}
  NG_MPI_Comm(const NG_MPI_Comm &comm) : value(comm.value) {}

  void operator=(int value_) { value = value_; }
  void operator=(uintptr_t value_) { value = value_; }
  bool operator==(const NG_MPI_Comm &comm) const { return value == comm.value; }
  bool operator!=(const NG_MPI_Comm &comm) const { return value != comm.value; }
};

struct NG_MPI_Datatype {
  uintptr_t value = 0;
  NG_MPI_Datatype() = default;
  NG_MPI_Datatype(uintptr_t value_) : value(value_) {}
  operator bool() const { return value != 0; }
  void operator=(NG_MPI_Datatype type) { value = type.value; }
  void operator=(uintptr_t value_) { value = value_; }
  void operator=(void *value_) { value = reinterpret_cast<uintptr_t>(value_); }
};

struct NG_MPI_Request {
  uintptr_t value = 0;
  NG_MPI_Request() = default;
  NG_MPI_Request(uintptr_t value_) : value(value_) {}
  void operator=(uintptr_t value_) { value = value_; }
};

struct NG_MPI_Op {
  uintptr_t value;
  NG_MPI_Op(uintptr_t value_) : value(value_) {}
  void operator=(uintptr_t value_) { value = value_; }
  void operator=(void *value_) { value = reinterpret_cast<uintptr_t>(value_); }
};

struct NG_MPI_Group {
  uintptr_t value = 0;
  NG_MPI_Group(uintptr_t value_) : value(value_) {}
  NG_MPI_Group() = default;
};

struct NG_MPI_Aint {
  intptr_t value = 0;
  NG_MPI_Aint(intptr_t value_) : value(value_) {}
  NG_MPI_Aint() = default;
};

#else
using NG_MPI_Status = MPI_Status;
using NG_MPI_Comm  = MPI_Comm;
using NG_MPI_Datatype = MPI_Datatype;
using NG_MPI_Request = MPI_Request;
using NG_MPI_Op = MPI_Op;
using NG_MPI_Group = MPI_Group;
using NG_MPI_Aint = MPI_Aint;
#endif

#include "ng_mpi_generated_declarations.hpp"

#if defined(NG_PYTHON) && defined(NG_MPI4PY)
NGCORE_API extern bool (*NG_MPI_CommFromMPI4Py)(py::handle, NG_MPI_Comm &);
NGCORE_API extern py::handle (*NG_MPI_CommToMPI4Py)(NG_MPI_Comm);
#endif

}  // namespace ngcore

#endif  // PARALLEL
#endif  // NG_MPI_HPP_INCLUDED

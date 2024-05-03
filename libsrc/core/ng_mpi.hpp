#ifndef NG_MPI_HPP_INCLUDED
#define NG_MPI_HPP_INCLUDED

#include <cstdint>
#include <filesystem>
#include <string>

#include "ngcore_api.hpp"

namespace ngcore {

NGCORE_API void InitMPI(
    std::filesystem::path mpi_lib_path = std::string("libmpi") +
                                         NETGEN_SHARED_LIBRARY_SUFFIX,
    std::filesystem::path ng_libs_dir = "");
NGCORE_API extern std::string mpi_library_version;

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
  uintptr_t value = 0;
  NG_MPI_Aint(uintptr_t value_) : value(value_) {}
  NG_MPI_Aint() = default;
};

#include "ng_mpi_generated_declarations.hpp"

}  // namespace ngcore

#endif  // NG_MPI_HPP_INCLUDED

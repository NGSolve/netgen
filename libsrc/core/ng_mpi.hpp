#ifndef NG_MPI_HPP_INCLUDED
#define NG_MPI_HPP_INCLUDED

#include <cstdint>
#include <stdexcept>
#include <string>

#include "ngcore_api.hpp"

namespace ngcore {

NGCORE_API void InitMPI();
NGCORE_API extern std::string mpi_library_version;

inline void not_implemented() { throw std::runtime_error("Not implemented"); }

struct NG_MPI_Status {
  uint64_t data[4];
};

struct NG_MPI_Comm {
  uintptr_t value;
  NG_MPI_Comm() { value = 0;}
  NG_MPI_Comm(uintptr_t v) : value(v) {}
};

struct NG_MPI_Datatype {
  uintptr_t value;
  NG_MPI_Datatype(uintptr_t v) : value(v) {}
  operator bool() const { return value != 0; }
};

struct NG_MPI_Request {
  uintptr_t value = 0;
  NG_MPI_Request(uintptr_t v) : value(v) {}
  NG_MPI_Request() = default;
};

struct NG_MPI_Op {
  uintptr_t value;
  NG_MPI_Op(uintptr_t v) : value(v) {}
};

struct NG_MPI_Group {
  uintptr_t value = 0;
  NG_MPI_Group(uintptr_t v) : value(v) {}
  NG_MPI_Group() = default;
};


#include "ng_mpi_generated_declarations.hpp"

}  // namespace ngcore

#endif  // NG_MPI_HPP_INCLUDED

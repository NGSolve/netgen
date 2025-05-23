#ifndef NG_MPI_NATIVE_HPP
#define NG_MPI_NATIVE_HPP

#include <mpi.h>

#include "mpi_wrapper.hpp"
#include "ng_mpi.hpp"

namespace ngcore {

MPI_Comm NG_MPI_Native(NG_MPI_Comm comm) {
  return reinterpret_cast<MPI_Comm>(comm.value);
}

MPI_Comm NG_MPI_Native(NgMPI_Comm comm) {
  return reinterpret_cast<MPI_Comm>(static_cast<NG_MPI_Comm>(comm).value);
}

}  // namespace ngcore

#endif  // NG_MPI_NATIVE_HPP

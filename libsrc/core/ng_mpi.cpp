#include "ng_mpi.hpp"

#include <mpi.h>

#include <type_traits>

#include "ngcore_api.hpp"

namespace ngcore {

template <typename T>
uintptr_t mpi2ng(T t) {
  if constexpr (std::is_pointer_v<T>)
    return reinterpret_cast<uintptr_t>(t);
  else
    return static_cast<uintptr_t>(t);
}

template <typename T>
T cast_ng2mpi(uintptr_t t) {
  if constexpr (std::is_pointer_v<T>)
    return reinterpret_cast<T>(t);
  else
    return static_cast<T>(t);
}

template <typename T>
T cast_ng2mpi(uintptr_t* t) {
  if constexpr (std::is_pointer_v<T>)
    return reinterpret_cast<T>(t);
  else
    return static_cast<T>(t);
}

MPI_Comm ng2mpi(NG_MPI_Comm c) {
  static_assert(sizeof(MPI_Comm) <= sizeof(c.value), "Size mismatch");
  return cast_ng2mpi<MPI_Comm>(c.value);
}

MPI_Comm* ng2mpi(NG_MPI_Comm* c) { return cast_ng2mpi<MPI_Comm*>(&c->value); }

MPI_Datatype ng2mpi(NG_MPI_Datatype c) {
  static_assert(sizeof(MPI_Datatype) <= sizeof(c.value), "Size mismatch");
  return cast_ng2mpi<MPI_Datatype>(c.value);
}

MPI_Request ng2mpi(NG_MPI_Request c) {
  static_assert(sizeof(MPI_Request) <= sizeof(c.value), "Size mismatch");
  return cast_ng2mpi<MPI_Request>(c.value);
}

int* ng2mpi(int* c) { return c; }

}  // namespace ngcore

using namespace ngcore;

NGCORE_API_EXPORT extern "C" void ng_init_mpi();

void ng_init_mpi() {
#include "ng_mpi_generated_init.hpp"
}

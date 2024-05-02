#include "ng_mpi.hpp"

#include <mpi.h>

#include <type_traits>

#include "ngcore_api.hpp"

namespace ngcore {

// template <typename T>
// uintptr_t mpi2ng(T t) {
//   if constexpr (std::is_pointer_v<T>)
//     return reinterpret_cast<uintptr_t>(t);
//   else
//     return static_cast<uintptr_t>(t);
// }

static_assert(sizeof(MPI_Status) <= sizeof(NG_MPI_Status), "Size mismatch");
static_assert(alignof(MPI_Status) <= alignof(NG_MPI_Status), "Size mismatch");

int mpi2ng(int v) { return v; }
void* mpi2ng(void*p) { return p; }

// TODO: When we are dealing with arrays of multiple MPI_Status, we need to copy them together in continuous memory
NG_MPI_Status* mpi2ng(MPI_Status*p) { return reinterpret_cast<NG_MPI_Status*>(p); }

NG_MPI_Comm mpi2ng(MPI_Comm c) { return reinterpret_cast<uintptr_t>(c); }

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
  static_assert(alignof(MPI_Comm) <= alignof(NG_MPI_Comm), "Size mismatch");
  return cast_ng2mpi<MPI_Comm>(c.value);
}

MPI_Group ng2mpi(NG_MPI_Group c) {
  static_assert(sizeof(MPI_Group) <= sizeof(c.value), "Size mismatch");
  static_assert(alignof(MPI_Group) <= alignof(NG_MPI_Group), "Size mismatch");
  return cast_ng2mpi<MPI_Group>(c.value);
}

MPI_Comm* ng2mpi(NG_MPI_Comm* c) { return cast_ng2mpi<MPI_Comm*>(&c->value); }
MPI_Group* ng2mpi(NG_MPI_Group* c) { return cast_ng2mpi<MPI_Group*>(&c->value); }
MPI_Datatype* ng2mpi(NG_MPI_Datatype* c) { return cast_ng2mpi<MPI_Datatype*>(&c->value); }
MPI_Request* ng2mpi(NG_MPI_Request* c) { return cast_ng2mpi<MPI_Request*>(&c->value); }
MPI_Status* ng2mpi(NG_MPI_Status* c) { return reinterpret_cast<MPI_Status*>(c); }

MPI_Datatype ng2mpi(NG_MPI_Datatype c) {
  static_assert(sizeof(MPI_Datatype) <= sizeof(c.value), "Size mismatch");
  return cast_ng2mpi<MPI_Datatype>(c.value);
}

MPI_Request ng2mpi(NG_MPI_Request c) {
  static_assert(sizeof(MPI_Request) <= sizeof(c.value), "Size mismatch");
  return cast_ng2mpi<MPI_Request>(c.value);
}

void* ng2mpi(void* c) { return c; }
char* ng2mpi(char* c) { return c; }
char*** ng2mpi(char*** c) { return c; }
int* ng2mpi(int* c) { return c; }
int ng2mpi(int c) { return c; }

}  // namespace ngcore

using namespace ngcore;

NGCORE_API_EXPORT extern "C" void ng_init_mpi();

void ng_init_mpi() {
#include "ng_mpi_generated_init.hpp"
}

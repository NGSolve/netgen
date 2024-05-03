#include "ng_mpi.hpp"

#include <mpi.h>

#include <type_traits>

#include "ngcore_api.hpp"

namespace ngcore {

static_assert(sizeof(MPI_Status) <= sizeof(NG_MPI_Status), "Size mismatch");
static_assert(alignof(MPI_Status) <= alignof(NG_MPI_Status), "Size mismatch");

int mpi2ng(int value) { return value; }
void* mpi2ng(void* ptr) { return ptr; }

// TODO: When we are dealing with arrays of multiple MPI_Statuses, we need to
// copy them together in continuous memory
NG_MPI_Status* mpi2ng(MPI_Status* status) {
  return reinterpret_cast<NG_MPI_Status*>(status);
}

#ifndef MPICH
NG_MPI_Comm mpi2ng(MPI_Comm comm) { return reinterpret_cast<uintptr_t>(comm); }
#endif

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
MPI_Request* ng2mpi(NG_MPI_Request* request) {
  return cast_ng2mpi<MPI_Request*>(&request->value);
}
MPI_Status* ng2mpi(NG_MPI_Status* status) {
  return reinterpret_cast<MPI_Status*>(status);
}
MPI_Aint* ng2mpi(NG_MPI_Aint* aint) {
  return reinterpret_cast<MPI_Aint*>(aint);
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

NGCORE_API_EXPORT extern "C" void ng_init_mpi();

void ng_init_mpi() {
#include "ng_mpi_generated_init.hpp"
}

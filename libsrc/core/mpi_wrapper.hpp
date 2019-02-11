#ifndef NGCORE_MPIWRAPPER_HPP
#define NGCORE_MPIWRAPPER_HPP

#ifdef PARALLEL
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#endif


namespace ngcore
{

#ifdef PARALLEL

  template <class T> struct MPI_typetrait  { };
  
  template <> struct MPI_typetrait<int> {
    static MPI_Datatype MPIType () { return MPI_INT; } };

  template <> struct MPI_typetrait<short> {
    static MPI_Datatype MPIType () { return MPI_SHORT; } };

  template <> struct MPI_typetrait<char> {
    static MPI_Datatype MPIType () { return MPI_CHAR; } };  

  template <> struct MPI_typetrait<size_t> {
    static MPI_Datatype MPIType () { return MPI_UNIT64_T; } };

  template <> struct MPI_typetrait<double> {
    static MPI_Datatype MPIType () { return MPI_DOUBLE; } };

  template <> struct MPI_typetrait<bool> {
    static MPI_Datatype MPIType () { return MPI_C_BOOL; } };

  
  template <class T, class T2 = decltype(MPI_typetrait<T>::MPIType())>
  inline MPI_Datatype GetMPIType () {
    return MPI_typetrait<T>::MPIType();
  }


  class NgMPI_Comm
  {
    MPI_Comm comm;
    int * refcount;
  public:
    NgMPI_Comm (MPI_Comm _comm, bool owns = false)
      : comm(_comm)
    {
      if (!owns)
        refcount = nullptr;
      else
        refcount = new int{1};
    }
    
    NgMPI_Comm (const NgMPI_Comm & c)
      : comm(c.comm), refcount(c.refcount)
    {
      if (refcount) (*refcount)++;
    }

    NgMPI_Comm (NgMPI_Comm && c)
      : comm(c.comm), refcount(c.refcount)
    {
      c.refcount = nullptr;
    }
    
    ~NgMPI_Comm()
    {
      if (refcount)
        if (--(*refcount) == 0)
          MPI_Comm_free(&comm);
    }
    
    operator MPI_Comm() const { return comm; }

    auto Rank() const { int r; MPI_Comm_rank(comm, &r); return r; }
    auto Size() const { int s; MPI_Comm_size(comm, &s); return s; }    
  };

  
#else
  
  class NgMPI_Comm
  {
    
  public:
    NgMPI_Comm (int _comm, bool owns = false)
    { ; }

    size_t Rank() const { return 0; }
    size_t Size() const { return 1; }
  };  
  
#endif






  

  
}

#endif


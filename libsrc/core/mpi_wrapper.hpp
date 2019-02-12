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

  template <> struct MPI_typetrait<unsigned char> {
    static MPI_Datatype MPIType () { return MPI_CHAR; } };  

  template <> struct MPI_typetrait<size_t> {
    static MPI_Datatype MPIType () { return MPI_UINT64_T; } };

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
    int rank, size;
  public:
    NgMPI_Comm (MPI_Comm _comm, bool owns = false)
      : comm(_comm)
    {
      if (!owns)
        refcount = nullptr;
      else
        refcount = new int{1};
      
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &size);
    }
    
    NgMPI_Comm (const NgMPI_Comm & c)
      : comm(c.comm), refcount(c.refcount), rank(c.rank), size(c.size)
    {
      if (refcount) (*refcount)++;
    }

    NgMPI_Comm (NgMPI_Comm && c)
      : comm(c.comm), refcount(c.refcount), rank(c.rank), size(c.size)
    {
      c.refcount = nullptr;
    }
    
    ~NgMPI_Comm()
    {
      if (refcount)
        if (--(*refcount) == 0)
          MPI_Comm_free(&comm);
    }

    NgMPI_Comm & operator= (const NgMPI_Comm & c)
    {
      if (refcount)
        if (--(*refcount) == 0)
          MPI_Comm_free(&comm);

      refcount = c.refcount;
      if (refcount) (*refcount)++;      
      size = c.size;
      rank = c.rank;
      return *this;
    }
    
    
    
    operator MPI_Comm() const { return comm; }

    int Rank() const { return rank; } // int r; MPI_Comm_rank(comm, &r); return r; }
    int Size() const { return size; } // int s; MPI_Comm_size(comm, &s); return s; }


    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    void Send( T & val, int dest, int tag) {
      MPI_Send (&val, 1, GetMPIType<T>(), dest, tag, comm);
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    void MyMPI_Recv (T & val, int src, int tag) {
      MPI_Recv (&val, 1, GetMPIType<T>(), src, tag, comm, MPI_STATUS_IGNORE);
    }

    
  };

  
#else
  class MPI_Comm {
    int nr;
  public:
    MPI_Comm (int _nr = 0) : nr(_nr) { ; }
    operator int() const { return nr; }
    bool operator== (MPI_Comm c2) const { return nr == c2.nr; }
  };
  static MPI_Comm MPI_COMM_WORLD = 12345, MPI_COMM_NULL = 10000;
  
  class NgMPI_Comm
  {
    
  public:
    NgMPI_Comm (MPI_Comm _comm, bool owns = false)
    { ; }

    size_t Rank() const { return 0; }
    size_t Size() const { return 1; }

    operator MPI_Comm() const { return MPI_Comm(); }

    template<typename T>
    void Send( T & val, int dest, int tag) { ; }
    
    template<typename T>
    void MyMPI_Recv (T & val, int src, int tag) { ; }
  };  
  
#endif






  

  
}

#endif


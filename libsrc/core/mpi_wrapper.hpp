#ifndef NGCORE_MPIWRAPPER_HPP
#define NGCORE_MPIWRAPPER_HPP

#ifdef PARALLEL
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#endif

#include "array.hpp"
#include "exception.hpp"

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
  protected:
    MPI_Comm comm;
    bool valid_comm;
    int * refcount;
    int rank, size;
  public:
    NgMPI_Comm ()
      : valid_comm(false), refcount(nullptr), rank(0), size(1)
    { ; }

    NgMPI_Comm (MPI_Comm _comm, bool owns = false)
      : comm(_comm), valid_comm(true)
    {
      if (!owns)
        refcount = nullptr;
      else
        refcount = new int{1};
      
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &size);
    }
    
    NgMPI_Comm (const NgMPI_Comm & c)
      : comm(c.comm), valid_comm(c.valid_comm), refcount(c.refcount),
        rank(c.rank), size(c.size)
    {
      if (refcount) (*refcount)++;
    }

    NgMPI_Comm (NgMPI_Comm && c)
      : comm(c.comm), valid_comm(c.valid_comm), refcount(c.refcount),
        rank(c.rank), size(c.size)
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
      comm = c.comm;
      valid_comm = c.valid_comm;
      size = c.size;
      rank = c.rank;
      return *this;
    }
    
    class InvalidCommException : public Exception {
    public:
      InvalidCommException() : Exception("Do not have a valid communicator") { ; }
    };
    
    operator MPI_Comm() const {
      if (!valid_comm) throw InvalidCommException();
      return comm;
    }

    int Rank() const { return rank; }
    int Size() const { return size; }
    void Barrier() const {
      if (size > 1) MPI_Barrier (comm);
    }
    

    /** --- blocking P2P --- **/

    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    void Send (T & val, int dest, int tag) const {
      MPI_Send (&val, 1, GetMPIType<T>(), dest, tag, comm);
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    void Send(FlatArray<T> s, int dest, int tag) const {
      MPI_Send (s.Data(), s.Size(), GetMPIType<T>(), dest, tag, comm);
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    void Recv (T & val, int src, int tag) const {
      MPI_Recv (&val, 1, GetMPIType<T>(), src, tag, comm, MPI_STATUS_IGNORE);
    }

    template <typename T, typename T2 = decltype(GetMPIType<T>())>
    void Recv (FlatArray <T> s, int src, int tag) const {
      MPI_Recv (s.Data(), s.Size(), GetMPIType<T> (), src, tag, comm, MPI_STATUS_IGNORE);
    }
    
    template <typename T, typename T2 = decltype(GetMPIType<T>())>
    void Recv (Array <T> & s, int src, int tag) const
    {
      MPI_Status status;
      int len;
      const MPI_Datatype MPI_T  = GetMPIType<T> ();
      MPI_Probe (src, tag, comm, &status);
      MPI_Get_count (&status, MPI_T, &len);
      s.SetSize (len);
      MPI_Recv (s.Data(), len, MPI_T, src, tag, comm, MPI_STATUS_IGNORE);
    }

    /** --- non-blocking P2P --- **/
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    MPI_Request ISend (T & val, int dest, int tag) const
    {
      MPI_Request request;
      MPI_Isend (&val, 1, GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }

    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    MPI_Request ISend (const FlatArray<T> & s, int dest, int tag) const
    {
      MPI_Request request;
      MPI_Isend (s.Data(), s.Size(), GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }

    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    MPI_Request IRecv (T & val, int dest, int tag) const
    {
      MPI_Request request;
      MPI_Irecv (&val, 1, GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    MPI_Request IRecv (const FlatArray<T> & s, int src, int tag) const
    { 
      MPI_Request request;
      MPI_Irecv (s.Data(), s.Size(), GetMPIType<T>(), src, tag, comm, &request);
      return request;
    }

    
    /** --- collectives --- **/

    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    T Reduce (T d, const MPI_Op & op, int root = 0) const
    {
      if (size == 1) return d;
      
      T global_d;
      MPI_Reduce (&d, &global_d, 1, GetMPIType<T>(), op, root, comm);
      return global_d;
    }
    
    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    T AllReduce (T d, const MPI_Op & op) const
    {
      if (size == 1) return d;
      
      T global_d;
      MPI_Allreduce ( &d, &global_d, 1, GetMPIType<T>(), op, comm);
      return global_d;
    }

    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    void Bcast (T & s, int root = 0) const {
      if (size == 1) return ;
      MPI_Bcast (&s, 1, GetMPIType<T>(), root, comm);
    }
    
    void Bcast (std::string & s, int root = 0) const 
    {
      if (size == 1) return;
      int len = s.length();
      Bcast (len, root);
      if (rank != 0) s.resize (len);
      MPI_Bcast (&s[0], len, MPI_CHAR, root, comm);
    }

    NgMPI_Comm SubCommunicator (FlatArray<int> procs) const
    {
      MPI_Comm subcomm;
      MPI_Group gcomm, gsubcomm;
      MPI_Comm_group(comm, &gcomm);
      MPI_Group_incl(gcomm, procs.Size(), procs.Data(), &gsubcomm);
      MPI_Comm_create_group(comm, gsubcomm, 4242, &subcomm);
      return NgMPI_Comm(subcomm, true);
    }

  }; // class NgMPI_Comm

  NETGEN_INLINE void MyMPI_WaitAll (FlatArray<MPI_Request> requests)
  {
    if (!requests.Size()) return;
    MPI_Waitall (requests.Size(), requests.Data(), MPI_STATUSES_IGNORE);
  }
  
  NETGEN_INLINE int MyMPI_WaitAny (FlatArray<MPI_Request> requests)
  {
    int nr;
    MPI_Waitany (requests.Size(), requests.Data(), &nr, MPI_STATUS_IGNORE);
    return nr;
  }

#else // PARALLEL
  class MPI_Comm {
    int nr;
  public:
    MPI_Comm (int _nr = 0) : nr(_nr) { ; }
    operator int() const { return nr; }
    bool operator== (MPI_Comm c2) const { return nr == c2.nr; }
  };
  static MPI_Comm MPI_COMM_WORLD = 12345, MPI_COMM_NULL = 10000;

  typedef int MPI_Op;
  typedef int MPI_Request;
  
  enum { MPI_SUM = 0, MPI_MIN = 1, MPI_MAX = 2 };
  
  class NgMPI_Comm
  {
    
  public:
    NgMPI_Comm () { ; } 
    NgMPI_Comm (MPI_Comm _comm, bool owns = false) { ; }

    size_t Rank() const { return 0; }
    size_t Size() const { return 1; }
    void Barrier() const { ; } 
    operator MPI_Comm() const { return MPI_Comm(); }

    template<typename T>
    void Send( T & val, int dest, int tag) const { ; }
    
    template<typename T>
    void Send(FlatArray<T> s, int dest, int tag) const { ; }

    template<typename T>
    void Recv (T & val, int src, int tag) const { ; }

    template <typename T>
    void Recv (FlatArray <T> s, int src, int tag) const { ; }

    template <typename T>
    void Recv (Array <T> & s, int src, int tag) const { ; }

    template<typename T>
    MPI_Request ISend (T & val, int dest, int tag) const { return 0; } 
    
    template<typename T>
    MPI_Request ISend (const FlatArray<T> & s, int dest, int tag) const { return 0; }

    template<typename T>
    MPI_Request IRecv (T & val, int dest, int tag) const { return 0; } 
    
    template<typename T>
    MPI_Request IRecv (const FlatArray<T> & s, int src, int tag) const { return 0; }

    template <typename T>
    T Reduce (T d, const MPI_Op & op, int root = 0) const { return d; }
    
    template <typename T>
    T AllReduce (T d, const MPI_Op & op) const { return d; }

    template <typename T>
    void Bcast (T & s, int root = 0) const { ; } 

    NgMPI_Comm SubCommunicator (FlatArray<int> procs) const
    { return *this; }
  };  

  NETGEN_INLINE void MyMPI_WaitAll (FlatArray<MPI_Request> requests) { ; }

#endif // PARALLEL

} // namespace ngcore

#endif // NGCORE_MPIWRAPPER_HPP


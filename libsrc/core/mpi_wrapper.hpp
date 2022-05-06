#ifndef NGCORE_MPIWRAPPER_HPP
#define NGCORE_MPIWRAPPER_HPP

#ifdef PARALLEL
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#endif

#include "array.hpp"
#include "table.hpp"
#include "exception.hpp"
#include "profiler.hpp"
#include "ngstream.hpp"

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

  template <> struct MPI_typetrait<signed char> {
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
  
  template <class T>
  inline MPI_Datatype GetMPIType (T &) {
    return GetMPIType<T>();
  }


  inline void MyMPI_WaitAll (FlatArray<MPI_Request> requests)
  {
    static Timer t("MPI - WaitAll"); RegionTimer reg(t);    
    if (!requests.Size()) return;
    MPI_Waitall (requests.Size(), requests.Data(), MPI_STATUSES_IGNORE);
  }
  
  inline int MyMPI_WaitAny (FlatArray<MPI_Request> requests)
  {
    int nr;
    MPI_Waitany (requests.Size(), requests.Data(), &nr, MPI_STATUS_IGNORE);
    return nr;
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
      int flag;
      MPI_Initialized (&flag);
      if (!flag)
        {
          valid_comm = false;
          refcount = nullptr;
          rank = 0;
          size = 1;
          return;
        }

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

    bool ValidCommunicator() const
    {
      return valid_comm;
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
      static Timer t("MPI - Barrier"); RegionTimer reg(t);
      if (size > 1) MPI_Barrier (comm);
    }
    

    /** --- blocking P2P --- **/

    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    void Send (T & val, int dest, int tag) const {
      MPI_Send (&val, 1, GetMPIType<T>(), dest, tag, comm);
    }

    void Send (const std::string & s, int dest, int tag) const {
      MPI_Send( const_cast<char*> (&s[0]), s.length(), MPI_CHAR, dest, tag, comm);
    }
    
    template<typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Send(FlatArray<T,TI> s, int dest, int tag) const {
      MPI_Send (s.Data(), s.Size(), GetMPIType<T>(), dest, tag, comm);
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    void Recv (T & val, int src, int tag) const {
      MPI_Recv (&val, 1, GetMPIType<T>(), src, tag, comm, MPI_STATUS_IGNORE);
    }

    void Recv (std::string & s, int src, int tag) const {    
      MPI_Status status;
      int len;
      MPI_Probe (src, tag, comm, &status);
      MPI_Get_count (&status, MPI_CHAR, &len);
      // s.assign (len, ' ');
      s.resize (len);
      MPI_Recv( &s[0], len, MPI_CHAR, src, tag, comm, MPI_STATUS_IGNORE);
    }
    

    template <typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Recv (FlatArray <T,TI> s, int src, int tag) const {
      MPI_Recv (s.Data(), s.Size(), GetMPIType<T> (), src, tag, comm, MPI_STATUS_IGNORE);
    }
    
    template <typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Recv (Array <T,TI> & s, int src, int tag) const
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
    MPI_Request ISend (FlatArray<T> s, int dest, int tag) const
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
    MPI_Request IRecv (FlatArray<T> s, int src, int tag) const
    { 
      MPI_Request request;
      MPI_Irecv (s.Data(), s.Size(), GetMPIType<T>(), src, tag, comm, &request);
      return request;
    }

    
    /** --- collectives --- **/

    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    T Reduce (T d, const MPI_Op & op, int root = 0) const
    {
      static Timer t("MPI - Reduce"); RegionTimer reg(t);          
      if (size == 1) return d;
      
      T global_d;
      MPI_Reduce (&d, &global_d, 1, GetMPIType<T>(), op, root, comm);
      return global_d;
    }
    
    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    T AllReduce (T d, const MPI_Op & op) const
    {
      static Timer t("MPI - AllReduce"); RegionTimer reg(t);
      if (size == 1) return d;
      
      T global_d;
      MPI_Allreduce ( &d, &global_d, 1, GetMPIType<T>(), op, comm);
      return global_d;
    }

    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    void AllReduce (FlatArray<T> d, const MPI_Op & op) const
    {
      static Timer t("MPI - AllReduce Array"); RegionTimer reg(t);
      if (size == 1) return;
      
      MPI_Allreduce (MPI_IN_PLACE, d.Data(), d.Size(), GetMPIType<T>(), op, comm);
    }
    
    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    void Bcast (T & s, int root = 0) const {
      if (size == 1) return;
      static Timer t("MPI - Bcast"); RegionTimer reg(t);
      MPI_Bcast (&s, 1, GetMPIType<T>(), root, comm);
    }

    
    template <class T>
    void Bcast (Array<T> & d, int root = 0)
    {
      if (size == 1) return;
      
      int ds = d.Size();
      Bcast (ds, root);
      if (Rank() != root) d.SetSize (ds);
      if (ds != 0)
        MPI_Bcast (d.Data(), ds, GetMPIType<T>(), root, comm);
    }

    
    void Bcast (std::string & s, int root = 0) const 
    {
      if (size == 1) return;
      int len = s.length();
      Bcast (len, root);
      if (rank != 0) s.resize (len);
      MPI_Bcast (&s[0], len, MPI_CHAR, root, comm);
    }

    template <typename T>
    void AllToAll (FlatArray<T> send, FlatArray<T> recv) const
    {
      MPI_Alltoall (send.Data(), 1, GetMPIType<T>(), 
                    recv.Data(), 1, GetMPIType<T>(), comm);
    }


    template <typename T>
    void ScatterRoot (FlatArray<T> send) const
    {
      if (size == 1) return;
      MPI_Scatter (send.Data(), 1, GetMPIType<T>(),
                   MPI_IN_PLACE, -1, GetMPIType<T>(), 0, comm);
    }
    
    template <typename T>
    void Scatter (T & recv) const
    {
      if (size == 1) return;      
      MPI_Scatter (NULL, 0, GetMPIType<T>(),
                   &recv, 1, GetMPIType<T>(), 0, comm);
    }

    template <typename T>
    void GatherRoot (FlatArray<T> recv) const
    {
      recv[0] = T(0);
      if (size == 1) return;      
      MPI_Gather (MPI_IN_PLACE, 1, GetMPIType<T>(), 
                  recv.Data(), 1, GetMPIType<T>(), 0, comm);
    }

    template <typename T>
    void Gather (T send) const
    {
      if (size == 1) return;            
      MPI_Gather (&send, 1, GetMPIType<T>(), 
                  NULL, 1, GetMPIType<T>(), 0, comm);
    }

    
    template <typename T>
    void AllGather (T val, FlatArray<T> recv) const
    {
      if (size == 1)
        {
          recv[0] = val;
          return;
        }
      MPI_Allgather (&val, 1, GetMPIType<T>(),
                     recv.Data(), 1, GetMPIType<T>(), 
                     comm);
    }
    


    template <typename T>
    void ExchangeTable (DynamicTable<T> & send_data, 
                        DynamicTable<T> & recv_data, int tag)
    {
      Array<int> send_sizes(size);
      Array<int> recv_sizes(size);
      
      for (int i = 0; i < size; i++)
        send_sizes[i] = send_data[i].Size();
      
      AllToAll (send_sizes, recv_sizes);
    
      recv_data = DynamicTable<T> (recv_sizes, true);
      
      Array<MPI_Request> requests;
      for (int dest = 0; dest < size; dest++)
        if (dest != rank && send_data[dest].Size())
          requests.Append (ISend (FlatArray<T>(send_data[dest]), dest, tag));
      
      for (int dest = 0; dest < size; dest++)
        if (dest != rank && recv_data[dest].Size())
          requests.Append (IRecv (FlatArray<T>(recv_data[dest]), dest, tag));

      MyMPI_WaitAll (requests);
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





  
  class MyMPI
  {
    bool initialized_by_me;
  public:
    MyMPI(int argc, char ** argv) 
    {
      int is_init = -1;
      MPI_Initialized(&is_init);
      if (!is_init)
        {
          MPI_Init (&argc, &argv);
          initialized_by_me = true;
        }
      else
        initialized_by_me = false;
      
      NgMPI_Comm comm(MPI_COMM_WORLD);
      NGSOStream::SetGlobalActive (comm.Rank() == 0);
      
      if (comm.Size() > 1)
        TaskManager::SetNumThreads (1);
    }
    
    ~MyMPI()
    {
      if (initialized_by_me)
        MPI_Finalize ();
    }
  };
  


  

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
  typedef int MPI_Datatype;  
  typedef int MPI_Request;
  
  enum { MPI_SUM = 0, MPI_MIN = 1, MPI_MAX = 2, MPI_LOR = 4711 };
  
  class NgMPI_Comm
  {
    
  public:
    NgMPI_Comm () { ; } 
    NgMPI_Comm (MPI_Comm _comm, bool owns = false) { ; }

    size_t Rank() const { return 0; }
    size_t Size() const { return 1; }
    bool ValidCommunicator() const { return false; }
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
    MPI_Request ISend (FlatArray<T> s, int dest, int tag) const { return 0; }

    template<typename T>
    MPI_Request IRecv (T & val, int dest, int tag) const { return 0; } 
    
    template<typename T>
    MPI_Request IRecv (FlatArray<T> s, int src, int tag) const { return 0; }

    template <typename T>
    T Reduce (T d, const MPI_Op & op, int root = 0) const { return d; }
    
    template <typename T>
    T AllReduce (T d, const MPI_Op & op) const { return d; }

    template <typename T>
    void AllReduce (FlatArray<T> d, const MPI_Op & op) const { ; }
    
    template <typename T>
    void Bcast (T & s, int root = 0) const { ; } 

    template <class T>
    void Bcast (Array<T> & d, int root = 0) { ; } 

    template <typename T>
    void AllGather (T val, FlatArray<T> recv) const
    {
      recv[0] = val;
    }

    template <typename T>
    void ExchangeTable (DynamicTable<T> & send_data, 
                        DynamicTable<T> & recv_data, int tag) { ; }

    
    NgMPI_Comm SubCommunicator (FlatArray<int> procs) const
    { return *this; }
  };  

  inline void MyMPI_WaitAll (FlatArray<MPI_Request> requests) { ; }
  inline int MyMPI_WaitAny (FlatArray<MPI_Request> requests) { return 0; }

  class MyMPI
  {
  public:
    MyMPI(int argc, char ** argv) { ; }
  };

  
#endif // PARALLEL

} // namespace ngcore

#endif // NGCORE_MPIWRAPPER_HPP


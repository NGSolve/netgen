#ifndef NGCORE_MPIWRAPPER_HPP
#define NGCORE_MPIWRAPPER_HPP

#include <array>

#include <complex>

#include "array.hpp"
#include "table.hpp"
#include "exception.hpp"
#include "profiler.hpp"
#include "ngstream.hpp"
#include "ng_mpi.hpp"

namespace ngcore
{

#ifdef PARALLEL

  template <class T> struct MPI_typetrait  { };
  
  template <> struct MPI_typetrait<int> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_INT; } };

  template <> struct MPI_typetrait<short> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_SHORT; } };

  template <> struct MPI_typetrait<char> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_CHAR; } };

  template <> struct MPI_typetrait<signed char> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_CHAR; } };
  
  template <> struct MPI_typetrait<unsigned char> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_CHAR; } };

  template <> struct MPI_typetrait<size_t> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_UINT64_T; } };

  template <> struct MPI_typetrait<double> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_DOUBLE; } };

  template <> struct MPI_typetrait<std::complex<double>> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_CXX_DOUBLE_COMPLEX; } };

  template <> struct MPI_typetrait<bool> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_C_BOOL; } };


  template<typename T, size_t S>
  struct MPI_typetrait<std::array<T,S>>
  {
    static NG_MPI_Datatype MPIType ()
    { 
      static NG_MPI_Datatype NG_MPI_T = 0;
      if (!NG_MPI_T)
	{
	  NG_MPI_Type_contiguous ( S, MPI_typetrait<T>::MPIType(), &NG_MPI_T);
	  NG_MPI_Type_commit ( &NG_MPI_T );
	}
      return NG_MPI_T;
    }
  };
  
  template <class T, class T2 = decltype(MPI_typetrait<T>::MPIType())>
  inline NG_MPI_Datatype GetMPIType () {
    return MPI_typetrait<T>::MPIType();
  }

  template <class T>
  inline NG_MPI_Datatype GetMPIType (T &) {
    return GetMPIType<T>();
  }

  class NgMPI_Request
  {
    NG_MPI_Request request;
  public:
    NgMPI_Request (NG_MPI_Request requ) : request{requ} { }
    NgMPI_Request (const NgMPI_Request&) = delete;    
    NgMPI_Request (NgMPI_Request&&) = default;
    ~NgMPI_Request () { NG_MPI_Wait (&request, NG_MPI_STATUS_IGNORE); }
    void Wait() {  NG_MPI_Wait (&request, NG_MPI_STATUS_IGNORE); }
    operator NG_MPI_Request() &&
    {
      auto tmp = request;
      request = NG_MPI_REQUEST_NULL;
      return tmp;
    }
  };

  class NgMPI_Requests
  {
    Array<NG_MPI_Request> requests;
  public:
    NgMPI_Requests() = default;
    ~NgMPI_Requests() { WaitAll(); }

    void Reset() { requests.SetSize0(); }
    
    NgMPI_Requests & operator+= (NgMPI_Request && r)
    {
      requests += NG_MPI_Request(std::move(r));
      return *this;
    }

    NgMPI_Requests & operator+= (NG_MPI_Request r)
    {
      requests += r;
      return *this;
    }
    
    void WaitAll()
    {
      static Timer t("NgMPI - WaitAll"); RegionTimer reg(t);    
      if (!requests.Size()) return;
      NG_MPI_Waitall (requests.Size(), requests.Data(), NG_MPI_STATUSES_IGNORE);
    }

    int WaitAny ()
    {
      int nr;
      NG_MPI_Waitany (requests.Size(), requests.Data(), &nr, NG_MPI_STATUS_IGNORE);
      return nr;
    }
  };
  
  [[deprecated("use requests.WaitAll instread")]]
  inline void MyMPI_WaitAll (FlatArray<NG_MPI_Request> requests)
  {
    static Timer t("MPI - WaitAll"); RegionTimer reg(t);    
    if (!requests.Size()) return;
    NG_MPI_Waitall (requests.Size(), requests.Data(), NG_MPI_STATUSES_IGNORE);
  }

  [[deprecated("use requests.WaitAny instread")]]  
  inline int MyMPI_WaitAny (FlatArray<NG_MPI_Request> requests)
  {
    int nr;
    NG_MPI_Waitany (requests.Size(), requests.Data(), &nr, NG_MPI_STATUS_IGNORE);
    return nr;
  }

  

  class NgMPI_Comm
  {
  protected:
    NG_MPI_Comm comm;
    bool valid_comm;
    int * refcount;
    int rank, size;
  public:
    NgMPI_Comm ()
      : valid_comm(false), refcount(nullptr), rank(0), size(1)
    { ; }

    NgMPI_Comm (NG_MPI_Comm _comm, bool owns = false)
      : comm(_comm), valid_comm(true)
    {
      int flag;
      NG_MPI_Initialized (&flag);
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
      
      NG_MPI_Comm_rank(comm, &rank);
      NG_MPI_Comm_size(comm, &size);
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
          NG_MPI_Comm_free(&comm);
    }

    bool ValidCommunicator() const
    {
      return valid_comm;
    }
    
    NgMPI_Comm & operator= (const NgMPI_Comm & c)
    {
      if (refcount)
        if (--(*refcount) == 0)
          NG_MPI_Comm_free(&comm);

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
    
    operator NG_MPI_Comm() const {
      if (!valid_comm) throw InvalidCommException();
      return comm;
    }

    int Rank() const { return rank; }
    int Size() const { return size; }
    void Barrier() const {
      static Timer t("MPI - Barrier"); RegionTimer reg(t);
      if (size > 1) NG_MPI_Barrier (comm);
    }
    

    /** --- blocking P2P --- **/

    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    void Send (T & val, int dest, int tag) const {
      NG_MPI_Send (&val, 1, GetMPIType<T>(), dest, tag, comm);
    }

    void Send (const std::string & s, int dest, int tag) const {
      NG_MPI_Send( const_cast<char*> (&s[0]), s.length(), NG_MPI_CHAR, dest, tag, comm);
    }
    
    template<typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Send(FlatArray<T,TI> s, int dest, int tag) const {
      NG_MPI_Send (s.Data(), s.Size(), GetMPIType<T>(), dest, tag, comm);
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    void Recv (T & val, int src, int tag) const {
      NG_MPI_Recv (&val, 1, GetMPIType<T>(), src, tag, comm, NG_MPI_STATUS_IGNORE);
    }

    void Recv (std::string & s, int src, int tag) const {    
      NG_MPI_Status status;
      int len;
      NG_MPI_Probe (src, tag, comm, &status);
      NG_MPI_Get_count (&status, NG_MPI_CHAR, &len);
      // s.assign (len, ' ');
      s.resize (len);
      NG_MPI_Recv( &s[0], len, NG_MPI_CHAR, src, tag, comm, NG_MPI_STATUS_IGNORE);
    }
    

    template <typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Recv (FlatArray <T,TI> s, int src, int tag) const {
      NG_MPI_Recv (s.Data(), s.Size(), GetMPIType<T> (), src, tag, comm, NG_MPI_STATUS_IGNORE);
    }
    
    template <typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Recv (Array <T,TI> & s, int src, int tag) const
    {
      NG_MPI_Status status;
      int len;
      const NG_MPI_Datatype NG_MPI_T  = GetMPIType<T> ();
      NG_MPI_Probe (src, tag, comm, &status);
      NG_MPI_Get_count (&status, NG_MPI_T, &len);
      s.SetSize (len);
      NG_MPI_Recv (s.Data(), len, NG_MPI_T, src, tag, comm, NG_MPI_STATUS_IGNORE);
    }

    /** --- non-blocking P2P --- **/

    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    [[nodiscard]] NG_MPI_Request ISend (T & val, int dest, int tag) const
    {
      NG_MPI_Request request;
      NG_MPI_Isend (&val, 1, GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    [[nodiscard]] NG_MPI_Request ISend (FlatArray<T> s, int dest, int tag) const
    {
      NG_MPI_Request request;
      NG_MPI_Isend (s.Data(), s.Size(), GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())> 
    [[nodiscard]] NG_MPI_Request IRecv (T & val, int dest, int tag) const
    {
      NG_MPI_Request request;
      NG_MPI_Irecv (&val, 1, GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }
    
    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    [[nodiscard]] NG_MPI_Request IRecv (FlatArray<T> s, int src, int tag) const
    { 
      NG_MPI_Request request;
      NG_MPI_Irecv (s.Data(), s.Size(), GetMPIType<T>(), src, tag, comm, &request);
      return request;
    }

    
    /** --- collectives --- **/

    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    T Reduce (T d, const NG_MPI_Op & op, int root = 0) const
    {
      static Timer t("MPI - Reduce"); RegionTimer reg(t);          
      if (size == 1) return d;
      
      T global_d;
      NG_MPI_Reduce (&d, &global_d, 1, GetMPIType<T>(), op, root, comm);
      return global_d;
    }
    
    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    T AllReduce (T d, const NG_MPI_Op & op) const
    {
      static Timer t("MPI - AllReduce"); RegionTimer reg(t);
      if (size == 1) return d;
      
      T global_d;
      NG_MPI_Allreduce ( &d, &global_d, 1, GetMPIType<T>(), op, comm);
      return global_d;
    }

    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    void AllReduce (FlatArray<T> d, const NG_MPI_Op & op) const
    {
      static Timer t("MPI - AllReduce Array"); RegionTimer reg(t);
      if (size == 1) return;
      
      NG_MPI_Allreduce (NG_MPI_IN_PLACE, d.Data(), d.Size(), GetMPIType<T>(), op, comm);
    }
    
    template <typename T, typename T2 = decltype(GetMPIType<T>())> 
    void Bcast (T & s, int root = 0) const {
      if (size == 1) return;
      static Timer t("MPI - Bcast"); RegionTimer reg(t);
      NG_MPI_Bcast (&s, 1, GetMPIType<T>(), root, comm);
    }


    template <class T, size_t S>
    void Bcast (std::array<T,S> & d, int root = 0) const
    {
      if (size == 1) return;
      if (S != 0)
        NG_MPI_Bcast (&d[0], S, GetMPIType<T>(), root, comm);
    }
    
    
    template <class T>
    void Bcast (Array<T> & d, int root = 0) const
    {
      if (size == 1) return;
      
      int ds = d.Size();
      Bcast (ds, root);
      if (Rank() != root) d.SetSize (ds);
      if (ds != 0)
        NG_MPI_Bcast (d.Data(), ds, GetMPIType<T>(), root, comm);
    }

    
    void Bcast (std::string & s, int root = 0) const 
    {
      if (size == 1) return;
      int len = s.length();
      Bcast (len, root);
      if (rank != 0) s.resize (len);
      NG_MPI_Bcast (&s[0], len, NG_MPI_CHAR, root, comm);
    }


    
    template <class T, size_t S>
    [[nodiscard]] NgMPI_Request IBcast (std::array<T,S> & d, int root = 0) const
    {
      NG_MPI_Request request;      
      NG_MPI_Ibcast (&d[0], S, GetMPIType<T>(), root, comm, &request);
      return request;
    }

    template <class T>
     [[nodiscard]] NgMPI_Request IBcast (FlatArray<T> d, int root = 0) const
    {
      NG_MPI_Request request;      
      int ds = d.Size();
      NG_MPI_Ibcast (d.Data(), ds, GetMPIType<T>(), root, comm, &request);
      return request;
    }

    
    template <typename T>
    void AllToAll (FlatArray<T> send, FlatArray<T> recv) const
    {
      NG_MPI_Alltoall (send.Data(), 1, GetMPIType<T>(),
                       recv.Data(), 1, GetMPIType<T>(), comm);
    }


    template <typename T>
    void ScatterRoot (FlatArray<T> send) const
    {
      if (size == 1) return;
      NG_MPI_Scatter (send.Data(), 1, GetMPIType<T>(),
                      NG_MPI_IN_PLACE, -1, GetMPIType<T>(), 0, comm);
    }
    
    template <typename T>
    void Scatter (T & recv) const
    {
      if (size == 1) return;      
      NG_MPI_Scatter (NULL, 0, GetMPIType<T>(),
                      &recv, 1, GetMPIType<T>(), 0, comm);
    }

    template <typename T>
    void GatherRoot (FlatArray<T> recv) const
    {
      recv[0] = T(0);
      if (size == 1) return;      
      NG_MPI_Gather (NG_MPI_IN_PLACE, 1, GetMPIType<T>(),
                     recv.Data(), 1, GetMPIType<T>(), 0, comm);
    }

    template <typename T>
    void Gather (T send) const
    {
      if (size == 1) return;            
      NG_MPI_Gather (&send, 1, GetMPIType<T>(),
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
      NG_MPI_Allgather (&val, 1, GetMPIType<T>(),
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
      
      NgMPI_Requests requests;
      for (int dest = 0; dest < size; dest++)
        if (dest != rank && send_data[dest].Size())
          requests += ISend (FlatArray<T>(send_data[dest]), dest, tag);
      
      for (int dest = 0; dest < size; dest++)
        if (dest != rank && recv_data[dest].Size())
          requests += IRecv (FlatArray<T>(recv_data[dest]), dest, tag);

      requests.WaitAll();
    }
    



    
    NgMPI_Comm SubCommunicator (FlatArray<int> procs) const
    {
      NG_MPI_Comm subcomm;
      NG_MPI_Group gcomm, gsubcomm;
      NG_MPI_Comm_group(comm, &gcomm);
      NG_MPI_Group_incl(gcomm, procs.Size(), procs.Data(), &gsubcomm);
      NG_MPI_Comm_create_group(comm, gsubcomm, 4242, &subcomm);
      return NgMPI_Comm(subcomm, true);
    }

  }; // class NgMPI_Comm

#else // PARALLEL
  class NG_MPI_Comm {
    int nr;
  public:
    NG_MPI_Comm (int _nr = 0) : nr(_nr) { ; }
    operator int() const { return nr; }
    bool operator== (NG_MPI_Comm c2) const { return nr == c2.nr; }
  };
  static NG_MPI_Comm NG_MPI_COMM_WORLD = 12345, NG_MPI_COMM_NULL = 10000;

  typedef int NG_MPI_Op;
  typedef int NG_MPI_Datatype;  
  typedef int NG_MPI_Request;
  
  enum { NG_MPI_SUM = 0, NG_MPI_MIN = 1, NG_MPI_MAX = 2, NG_MPI_LOR = 4711 };

  inline void NG_MPI_Type_contiguous ( int, NG_MPI_Datatype, NG_MPI_Datatype*) { ; } 
  inline void NG_MPI_Type_commit ( NG_MPI_Datatype * ) { ; }

  template <class T> struct MPI_typetrait  {
    static NG_MPI_Datatype MPIType () { return -1; }    
  };
  template <class T, class T2=void>
  inline NG_MPI_Datatype GetMPIType () { return -1; }

  class NgMPI_Request {
  public:
    NgMPI_Request() = default;
    NgMPI_Request(NgMPI_Request &&) { ; }    
    NgMPI_Request(NG_MPI_Request &&) { ; }
  };
  class NgMPI_Requests
  {
  public:
    NgMPI_Requests & operator+= (NgMPI_Request &&) { return *this; }
    NgMPI_Requests & operator+= (NG_MPI_Request r) { return *this; }
    void Reset() { ; }
    void WaitAll() { ; }
    int WaitAny() { return 0; }
  };
  
  class NgMPI_Comm
  {
    
  public:
    NgMPI_Comm () { ; } 
    NgMPI_Comm (NG_MPI_Comm _comm, bool owns = false) { ; }

    size_t Rank() const { return 0; }
    size_t Size() const { return 1; }
    bool ValidCommunicator() const { return false; }
    void Barrier() const { ; } 
    operator NG_MPI_Comm() const { return NG_MPI_Comm(); }

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
    NG_MPI_Request ISend (T & val, int dest, int tag) const { return 0; } 
    
    template<typename T>
    NG_MPI_Request ISend (FlatArray<T> s, int dest, int tag) const { return 0; }

    template<typename T>
    NG_MPI_Request IRecv (T & val, int dest, int tag) const { return 0; } 
    
    template<typename T>
    NG_MPI_Request IRecv (FlatArray<T> s, int src, int tag) const { return 0; }

    template <typename T>
    T Reduce (T d, const NG_MPI_Op & op, int root = 0) const { return d; }
    
    template <typename T>
    T AllReduce (T d, const NG_MPI_Op & op) const { return d; }

    template <typename T>
    void AllReduce (FlatArray<T> d, const NG_MPI_Op & op) const { ; }
    
    template <typename T>
    void Bcast (T & s, int root = 0) const { ; } 

    template <class T, size_t S>
    void Bcast (std::array<T,S> & d, int root = 0) const {}
    
    template <class T>
    void Bcast (Array<T> & d, int root = 0) const { ; } 

    template <class T, size_t S>
    NG_MPI_Request IBcast (std::array<T,S> & d, int root = 0) const { return 0; }

    template <class T>
    NG_MPI_Request IBcast (FlatArray<T> d, int root = 0) const { return 0; } 
    
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

  inline void MyMPI_WaitAll (FlatArray<NG_MPI_Request> requests) { ; }
  inline int MyMPI_WaitAny (FlatArray<NG_MPI_Request> requests) { return 0; }

#endif // PARALLEL

} // namespace ngcore

#endif // NGCORE_MPIWRAPPER_HPP


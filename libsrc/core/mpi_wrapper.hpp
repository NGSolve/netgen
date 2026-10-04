#ifndef NGCORE_MPIWRAPPER_HPP
#define NGCORE_MPIWRAPPER_HPP

#include <array>

#include <complex>
#include <limits>
#include <memory>

#include "array.hpp"
#include "table.hpp"
#include "exception.hpp"
#include "profiler.hpp"
#include "ngstream.hpp"
#include "ng_mpi.hpp"

namespace ngcore
{

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

  template <> struct MPI_typetrait<std::byte> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_BYTE; } };

  template <> struct MPI_typetrait<size_t> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_UINT64_T; } };

  template <> struct MPI_typetrait<double> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_DOUBLE; } };

  template <> struct MPI_typetrait<float> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_FLOAT; } };

  template <> struct MPI_typetrait<std::complex<double>> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_CXX_DOUBLE_COMPLEX; } };

  template <> struct MPI_typetrait<bool> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_C_BOOL; } };

  template <> struct MPI_typetrait<long> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_LONG; } };

  template <> struct MPI_typetrait<long long> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_LONG_LONG; } };

  template <> struct MPI_typetrait<unsigned> {
    static NG_MPI_Datatype MPIType () { return NG_MPI_UNSIGNED; } };

  template <class T> struct MPI_typetrait<const T> : MPI_typetrait<T> { };


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

  /// message size in the int count MPI expects
  inline int MPI_Count (size_t n)
  {
    if (n > size_t(std::numeric_limits<int>::max()))
      throw Exception("MPI message with " + ToString(n) + " entries exceeds the int count limit");
    return int(n);
  }

  /// owns a request, waits for it on destruction
  class NgMPI_Request
  {
    NG_MPI_Request request;
  public:
    NgMPI_Request () : request(NG_MPI_REQUEST_NULL) { }
    NgMPI_Request (NG_MPI_Request requ) : request{requ} { }
    NgMPI_Request (const NgMPI_Request&) = delete;
    NgMPI_Request (NgMPI_Request && r) : request(r.request) { r.request = NG_MPI_REQUEST_NULL; }
    NgMPI_Request & operator= (NgMPI_Request && r)
    {
      Wait();
      request = r.request;
      r.request = NG_MPI_REQUEST_NULL;
      return *this;
    }
    ~NgMPI_Request () { Wait(); }

    void Wait()
    {
      if (request.value != NG_MPI_REQUEST_NULL.value)
        NG_MPI_Wait (&request, NG_MPI_STATUS_IGNORE);
      request = NG_MPI_REQUEST_NULL;
    }
    /// hands the request over, the caller has to wait for it
    operator NG_MPI_Request() &&
    {
      auto tmp = request;
      request = NG_MPI_REQUEST_NULL;
      return tmp;
    }
  };

  /// a set of requests, waits for all of them on destruction
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
      requests.SetSize0();
    }

    int WaitAny ()
    {
      int nr;
      NG_MPI_Waitany (requests.Size(), requests.Data(), &nr, NG_MPI_STATUS_IGNORE);
      return nr;
    }
  };


  /*
    A communicator. Without a loaded and initialized MPI it is the serial
    communicator (rank 0 of 1) and all operations are local. Copies share
    the MPI communicator, which is freed with the last copy if owned.
    Counts are int as in MPI, larger messages throw.
  */
  class NgMPI_Comm
  {
  protected:
    NG_MPI_Comm comm;
    bool valid_comm;
    std::shared_ptr<NG_MPI_Comm> owner;   // frees the communicator with the last copy
    int rank, size;

    static constexpr int tag_subcomm = 4242;
  public:
    NgMPI_Comm ()
      : comm(), valid_comm(false), rank(0), size(1)
    { ; }

    NgMPI_Comm (NG_MPI_Comm _comm, bool owns = false)
      : comm(_comm), valid_comm(true), rank(0), size(1)
    {
      int flag = 0;
      if (MPI_Loaded())
        NG_MPI_Initialized (&flag);
      if (!flag)
        {
          valid_comm = false;
          return;
        }
      if (owns)
        owner = std::shared_ptr<NG_MPI_Comm> (new NG_MPI_Comm(_comm),
                                              [] (NG_MPI_Comm * c) { NG_MPI_Comm_free (c); delete c; });
      NG_MPI_Comm_rank (comm, &rank);
      NG_MPI_Comm_size (comm, &size);
    }

    NgMPI_Comm (const NgMPI_Comm &) = default;
    NgMPI_Comm (NgMPI_Comm &&) = default;
    NgMPI_Comm & operator= (const NgMPI_Comm &) = default;
    NgMPI_Comm & operator= (NgMPI_Comm &&) = default;
    ~NgMPI_Comm() = default;

    bool ValidCommunicator() const { return valid_comm; }

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
    void Send (const T & val, int dest, int tag) const {
      NG_MPI_Send (const_cast<T*>(&val), 1, GetMPIType<T>(), dest, tag, comm);
    }

    void Send (const std::string & s, int dest, int tag) const {
      NG_MPI_Send (const_cast<char*>(s.data()), MPI_Count(s.length()), NG_MPI_CHAR, dest, tag, comm);
    }

    template<typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Send (FlatArray<T,TI> s, int dest, int tag) const {
      NG_MPI_Send (const_cast<std::remove_const_t<T>*>(s.Data()), MPI_Count(s.Size()), GetMPIType<T>(), dest, tag, comm);
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
      s.resize (len);
      NG_MPI_Recv (s.data(), len, NG_MPI_CHAR, src, tag, comm, NG_MPI_STATUS_IGNORE);
    }

    template <typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Recv (FlatArray<T,TI> s, int src, int tag) const {
      NG_MPI_Recv (s.Data(), MPI_Count(s.Size()), GetMPIType<T>(), src, tag, comm, NG_MPI_STATUS_IGNORE);
    }

    /// receives an array of any size
    template <typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void Recv (Array<T,TI> & s, int src, int tag) const
    {
      NG_MPI_Status status;
      int len;
      const NG_MPI_Datatype NG_MPI_T = GetMPIType<T>();
      NG_MPI_Probe (src, tag, comm, &status);
      NG_MPI_Get_count (&status, NG_MPI_T, &len);
      s.SetSize (len);
      NG_MPI_Recv (s.Data(), len, NG_MPI_T, src, tag, comm, NG_MPI_STATUS_IGNORE);
    }

    /** --- non-blocking P2P, the returned request waits on destruction --- **/

    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    [[nodiscard]] NgMPI_Request ISend (const T & val, int dest, int tag) const
    {
      NG_MPI_Request request;
      NG_MPI_Isend (const_cast<T*>(&val), 1, GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }

    template<typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    [[nodiscard]] NgMPI_Request ISend (FlatArray<T,TI> s, int dest, int tag) const
    {
      NG_MPI_Request request;
      NG_MPI_Isend (const_cast<std::remove_const_t<T>*>(s.Data()), MPI_Count(s.Size()), GetMPIType<T>(), dest, tag, comm, &request);
      return request;
    }

    template<typename T, typename T2 = decltype(GetMPIType<T>())>
    [[nodiscard]] NgMPI_Request IRecv (T & val, int src, int tag) const
    {
      NG_MPI_Request request;
      NG_MPI_Irecv (&val, 1, GetMPIType<T>(), src, tag, comm, &request);
      return request;
    }

    template<typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    [[nodiscard]] NgMPI_Request IRecv (FlatArray<T,TI> s, int src, int tag) const
    {
      NG_MPI_Request request;
      NG_MPI_Irecv (s.Data(), MPI_Count(s.Size()), GetMPIType<T>(), src, tag, comm, &request);
      return request;
    }


    /** --- collectives --- **/

    /// result on root, the input on the other ranks
    template <typename T, typename T2 = decltype(GetMPIType<T>())>
    T Reduce (T d, const NG_MPI_Op & op, int root = 0) const
    {
      static Timer t("MPI - Reduce"); RegionTimer reg(t);
      if (size == 1) return d;
      T global_d = d;
      NG_MPI_Reduce (&d, &global_d, 1, GetMPIType<T>(), op, root, comm);
      return global_d;
    }

    template <typename T, typename T2 = decltype(GetMPIType<T>())>
    T AllReduce (T d, const NG_MPI_Op & op) const
    {
      static Timer t("MPI - AllReduce"); RegionTimer reg(t);
      if (size == 1) return d;
      T global_d = d;
      NG_MPI_Allreduce (&d, &global_d, 1, GetMPIType<T>(), op, comm);
      return global_d;
    }

    template <typename T, typename TI, typename T2 = decltype(GetMPIType<T>())>
    void AllReduce (FlatArray<T,TI> d, const NG_MPI_Op & op) const
    {
      static Timer t("MPI - AllReduce Array"); RegionTimer reg(t);
      if (size == 1) return;
      NG_MPI_Allreduce (NG_MPI_IN_PLACE, d.Data(), MPI_Count(d.Size()), GetMPIType<T>(), op, comm);
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

    /// the array is resized on the other ranks
    template <class T, class TI>
    void Bcast (Array<T,TI> & d, int root = 0) const
    {
      if (size == 1) return;
      int ds = MPI_Count(d.Size());
      Bcast (ds, root);
      if (rank != root) d.SetSize (ds);
      if (ds != 0)
        NG_MPI_Bcast (d.Data(), ds, GetMPIType<T>(), root, comm);
    }

    void Bcast (std::string & s, int root = 0) const
    {
      if (size == 1) return;
      int len = MPI_Count(s.length());
      Bcast (len, root);
      if (rank != root) s.resize (len);
      if (len != 0)
        NG_MPI_Bcast (s.data(), len, NG_MPI_CHAR, root, comm);
    }

    template <class T, size_t S>
    [[nodiscard]] NgMPI_Request IBcast (std::array<T,S> & d, int root = 0) const
    {
      NG_MPI_Request request;
      NG_MPI_Ibcast (&d[0], S, GetMPIType<T>(), root, comm, &request);
      return request;
    }

    template <class T, class TI>
    [[nodiscard]] NgMPI_Request IBcast (FlatArray<T,TI> d, int root = 0) const
    {
      NG_MPI_Request request;
      NG_MPI_Ibcast (d.Data(), MPI_Count(d.Size()), GetMPIType<T>(), root, comm, &request);
      return request;
    }

    /// one entry to and from every rank
    template <typename T>
    void AllToAll (FlatArray<T> send, FlatArray<T> recv) const
    {
      if (size == 1) { recv[0] = send[0]; return; }
      NG_MPI_Alltoall (send.Data(), 1, GetMPIType<T>(),
                       recv.Data(), 1, GetMPIType<T>(), comm);
    }

    /// one entry per rank from root; send is used on root only
    template <typename T>
    void Scatter (FlatArray<T> send, T & recv, int root = 0) const
    {
      if (size == 1) { recv = send[0]; return; }
      NG_MPI_Scatter (rank == root ? send.Data() : nullptr, 1, GetMPIType<T>(),
                      &recv, 1, GetMPIType<T>(), root, comm);
    }

    /// one entry per rank to root; recv is filled on root only
    template <typename T>
    void Gather (T send, FlatArray<T> recv, int root = 0) const
    {
      if (size == 1) { recv[0] = send; return; }
      NG_MPI_Gather (&send, 1, GetMPIType<T>(),
                     rank == root ? recv.Data() : nullptr, 1, GetMPIType<T>(), root, comm);
    }

    template <typename T>
    void AllGather (T val, FlatArray<T> recv) const
    {
      if (size == 1) { recv[0] = val; return; }
      NG_MPI_Allgather (&val, 1, GetMPIType<T>(),
                        recv.Data(), 1, GetMPIType<T>(), comm);
    }

    /// arrays of any size to root; recv has one row per rank, filled on root only
    template <typename T, typename TI>
    void GatherTable (FlatArray<T,TI> send, Table<T> & recv, int root = 0) const
    {
      static Timer t("MPI - GatherTable"); RegionTimer reg(t);
      int n = MPI_Count(send.Size());
      Array<int> sizes(size);
      Gather (n, sizes, root);
      if (rank != root)
        {
          NG_MPI_Gatherv (const_cast<std::remove_const_t<T>*>(send.Data()), n, GetMPIType<T>(),
                          nullptr, nullptr, nullptr, GetMPIType<T>(), root, comm);
          return;
        }
      recv = Table<T> (sizes);
      Array<int> displ(size);
      for (int i = 0, cnt = 0; i < size; i++)
        { displ[i] = cnt; cnt += sizes[i]; }
      if (size == 1)
        { recv[0] = send; return; }
      NG_MPI_Gatherv (const_cast<std::remove_const_t<T>*>(send.Data()), n, GetMPIType<T>(),
                      recv.AsArray().Data(), sizes.Data(), displ.Data(), GetMPIType<T>(), root, comm);
    }

    /// arrays of any size to every rank; recv has one row per rank
    template <typename T, typename TI>
    void AllGatherTable (FlatArray<T,TI> send, Table<T> & recv) const
    {
      static Timer t("MPI - AllGatherTable"); RegionTimer reg(t);
      int n = MPI_Count(send.Size());
      Array<int> sizes(size);
      AllGather (n, sizes);
      recv = Table<T> (sizes);
      if (size == 1)
        { recv[0] = send; return; }
      Array<int> displ(size);
      for (int i = 0, cnt = 0; i < size; i++)
        { displ[i] = cnt; cnt += sizes[i]; }
      NG_MPI_Allgatherv (const_cast<std::remove_const_t<T>*>(send.Data()), n, GetMPIType<T>(),
                         recv.AsArray().Data(), sizes.Data(), displ.Data(), GetMPIType<T>(), comm);
    }

    /// row i of send_data goes to rank i, row i of recv_data comes from rank i (own row is not copied)
    template <typename T>
    void ExchangeTable (DynamicTable<T> & send_data,
                        DynamicTable<T> & recv_data, int tag) const
    {
      Array<int> send_sizes(size);
      Array<int> recv_sizes(size);

      for (int i = 0; i < size; i++)
        send_sizes[i] = MPI_Count(send_data[i].Size());

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

    /// communicator of the given ranks, collective over them
    NgMPI_Comm SubCommunicator (FlatArray<int> procs) const
    {
      NG_MPI_Comm subcomm;
      NG_MPI_Group gcomm, gsubcomm;
      NG_MPI_Comm_group (comm, &gcomm);
      NG_MPI_Group_incl (gcomm, MPI_Count(procs.Size()), procs.Data(), &gsubcomm);
      NG_MPI_Comm_create_group (comm, gsubcomm, tag_subcomm, &subcomm);
      NG_MPI_Group_free (&gsubcomm);
      NG_MPI_Group_free (&gcomm);
      return NgMPI_Comm (subcomm, true);
    }


    /** --- deprecated role-split collectives --- **/

    template <typename T>
    [[deprecated("use Scatter(send, recv, root) on all ranks")]]
    void ScatterRoot (FlatArray<T> send) const
    {
      if (size == 1) return;
      NG_MPI_Scatter (send.Data(), 1, GetMPIType<T>(),
                      NG_MPI_IN_PLACE, -1, GetMPIType<T>(), 0, comm);
    }

    template <typename T>
    [[deprecated("use Scatter(send, recv, root) on all ranks")]]
    void Scatter (T & recv) const
    {
      if (size == 1) return;
      NG_MPI_Scatter (NULL, 0, GetMPIType<T>(),
                      &recv, 1, GetMPIType<T>(), 0, comm);
    }

    template <typename T>
    [[deprecated("use Gather(send, recv, root) on all ranks")]]
    void GatherRoot (FlatArray<T> recv) const
    {
      recv[0] = T(0);
      if (size == 1) return;
      NG_MPI_Gather (NG_MPI_IN_PLACE, 1, GetMPIType<T>(),
                     recv.Data(), 1, GetMPIType<T>(), 0, comm);
    }

    template <typename T>
    [[deprecated("use Gather(send, recv, root) on all ranks")]]
    void Gather (T send) const
    {
      if (size == 1) return;
      NG_MPI_Gather (&send, 1, GetMPIType<T>(),
                     NULL, 1, GetMPIType<T>(), 0, comm);
    }

  }; // class NgMPI_Comm


} // namespace ngcore

#endif // NGCORE_MPIWRAPPER_HPP

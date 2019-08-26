#ifndef FILE_PARALLEL
#define FILE_PARALLEL



#ifdef VTRACE
#include "vt_user.h"
#else
  #define VT_USER_START(n)
  #define VT_USER_END(n)
  #define VT_TRACER(n)
#endif


namespace netgen
{
  // using ngcore::id;
  // using ngcore::ntasks;

#ifndef PARALLEL
  /** without MPI, we need a dummy typedef **/
  // typedef int MPI_Comm;
#endif

  /** This is the "standard" communicator that will be used for netgen-objects. **/
  // extern DLL_HEADER NgMPI_Comm ng_comm;

#ifdef OLD
#ifdef PARALLEL
  inline int MyMPI_GetNTasks (MPI_Comm comm /* = ng_comm */)
  {
    int ntasks;
    MPI_Comm_size(comm, &ntasks);
    return ntasks;
  }
  inline int MyMPI_GetId (MPI_Comm comm /* = ng_comm */)
  {
    int id;
    MPI_Comm_rank(comm, &id);
    return id;
  }
#else
  // enum { MPI_COMM_WORLD = 12345, MPI_COMM_NULL = 0};
  inline int MyMPI_GetNTasks (MPI_Comm comm /* = ng_comm */) { return 1; }
  inline int MyMPI_GetId (MPI_Comm comm /* = ng_comm */) { return 0; }
#endif
#endif
  
  /*
#ifdef PARALLEL
  // For python wrapping of communicators
  struct PyMPI_Comm {
    MPI_Comm comm;
    bool owns_comm;
    PyMPI_Comm (MPI_Comm _comm, bool _owns_comm = false) : comm(_comm), owns_comm(_owns_comm) { }
    PyMPI_Comm (const PyMPI_Comm & c) = delete;
    ~PyMPI_Comm () {
      if (owns_comm)
	MPI_Comm_free(&comm);
    }
    inline int Rank() const { return MyMPI_GetId(comm); }
    inline int Size() const { return MyMPI_GetNTasks(comm); }
  };
#else
  // dummy without MPI
  struct PyMPI_Comm {
    MPI_Comm comm = 0;
    PyMPI_Comm (MPI_Comm _comm, bool _owns_comm = false) { }
    ~PyMPI_Comm () { }
    inline int Rank() const { return 0; }
    inline int Size() const { return 1; }
  };
#endif
  */
  
#ifdef PARALLEL
  template <class T>
  inline MPI_Datatype MyGetMPIType ( ) 
  { cerr << "ERROR in GetMPIType() -- no type found" << endl;return 0; }
  template <>
  inline MPI_Datatype MyGetMPIType<int> ( )
  { return MPI_INT; }
  template <>
  inline MPI_Datatype MyGetMPIType<double> ( ) 
  { return MPI_DOUBLE; }
  template <>
  inline MPI_Datatype MyGetMPIType<char> ( ) 
  { return MPI_CHAR; }
  template<>
  inline MPI_Datatype MyGetMPIType<size_t> ( ) 
  { return MPI_UINT64_T; }
#else
  typedef int MPI_Datatype;
  template <class T> inline MPI_Datatype MyGetMPIType ( ) { return 0; }
#endif

#ifdef PARALLEL
  enum { MPI_TAG_CMD = 110 };
  enum { MPI_TAG_MESH = 210 };
  enum { MPI_TAG_VIS = 310 };

  inline void MyMPI_Send (int i, int dest, int tag, MPI_Comm comm /* = ng_comm */)
  {
    int hi = i;
    MPI_Send( &hi, 1, MPI_INT, dest, tag, comm);
  }

  inline void MyMPI_Recv (int & i, int src, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Status status;
    MPI_Recv( &i, 1, MPI_INT, src, tag, comm, &status);
  }



  inline void MyMPI_Send (const string & s, int dest, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Send( const_cast<char*> (s.c_str()), s.length(), MPI_CHAR, dest, tag, comm);
  }

  inline void MyMPI_Recv (string & s, int src, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_CHAR, &len);
    s.assign (len, ' ');
    MPI_Recv( &s[0], len, MPI_CHAR, src, tag, comm, &status);
  }

 


  template <class T, int BASE>
  inline void MyMPI_Send (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Send( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, comm);
  }

  template <class T, int BASE>
  inline void MyMPI_Recv ( NgFlatArray<T, BASE> s, int src, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Status status;
    MPI_Recv( &s.First(), s.Size(), MyGetMPIType<T>(), src, tag, comm, &status);
  }

  template <class T, int BASE>
  inline void MyMPI_Recv ( NgArray <T, BASE> & s, int src, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, comm, &status);
    MPI_Get_count (&status, MyGetMPIType<T>(), &len);

    s.SetSize (len);
    MPI_Recv( &s.First(), len, MyGetMPIType<T>(), src, tag, comm, &status);
  }

  template <class T, int BASE>
  inline int MyMPI_Recv ( NgArray <T, BASE> & s, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Status status;
    int len;
    MPI_Probe (MPI_ANY_SOURCE, tag, comm, &status);

    int src = status.MPI_SOURCE;

    MPI_Get_count (&status, MyGetMPIType<T>(), &len);

    s.SetSize (len);
    MPI_Recv( &s.First(), len, MyGetMPIType<T>(), src, tag, comm, &status);

    return src;
  }


  /*
  template <class T, int BASE>
  inline void MyMPI_ISend (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Request & request)
  {
    MPI_Isend( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, MPI_COMM_WORLD, & request);
  }


  template <class T, int BASE>
  inline void MyMPI_IRecv (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Request & request)
  {
    MPI_Irecv( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, MPI_COMM_WORLD, & request);
  }
  */

  template <class T, int BASE>
  inline MPI_Request MyMPI_ISend (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Request request;
    MPI_Isend( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, comm, &request);
    return request;
  }


  template <class T, int BASE>
  inline MPI_Request MyMPI_IRecv (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Request request;
    MPI_Irecv( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, comm, &request);
    return request;
  }

  /*
  template <class T, int BASE>
  inline void MyMPI_ISend (NgFlatArray<T, BASE> s, int dest, int tag)
  {
    MPI_Request request;
    MPI_Isend( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, MPI_COMM_WORLD, &request);
    MPI_Request_free (&request);
  }


  template <class T, int BASE>
  inline void MyMPI_IRecv (NgFlatArray<T, BASE> s, int dest, int tag)
  {
    MPI_Request request;
    MPI_Irecv( &s.First(), s.Size(), MyGetMPIType<T>(), dest, tag, MPI_COMM_WORLD, &request);
    MPI_Request_free (&request);
  }
  */



  /*
    send a table entry to each of the processes in the group ...
    receive-table entries will be set
   */

  /*
  template <typename T>
  inline void MyMPI_ExchangeTable (TABLE<T> & send_data, 
				   TABLE<T> & recv_data, int tag,
				   MPI_Comm comm = MPI_COMM_WORLD)
  {
    int ntasks, rank;
    MPI_Comm_size(comm, &ntasks);
    MPI_Comm_rank(comm, &rank);

    NgArray<MPI_Request> requests;
    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank)
	requests.Append (MyMPI_ISend (send_data[dest], dest, tag, comm));

    for (int i = 0; i < ntasks-1; i++)
      {
	MPI_Status status;
	MPI_Probe (MPI_ANY_SOURCE, tag, comm, &status);
	int size, src = status.MPI_SOURCE;
	MPI_Get_count (&status, MPI_INT, &size);
	recv_data.SetEntrySize (src, size, sizeof(T));
	requests.Append (MyMPI_IRecv (recv_data[src], src, tag, comm));
      }
    MPI_Barrier (comm);
    MPI_Waitall (requests.Size(), &requests[0], MPI_STATUS_IGNORE);
  }
  */

  template <typename T>
  inline void MyMPI_ExchangeTable (TABLE<T> & send_data, 
				   TABLE<T> & recv_data, int tag,
				   const NgMPI_Comm & comm /* = ng_comm */)
  {
    /*
    int rank = MyMPI_GetId(comm);
    int ntasks = MyMPI_GetNTasks(comm);
    */
    int rank = comm.Rank();
    int ntasks = comm.Size();
    
    NgArray<int> send_sizes(ntasks);
    NgArray<int> recv_sizes(ntasks);
    for (int i = 0; i < ntasks; i++)
      send_sizes[i] = send_data[i].Size();
    
    MPI_Alltoall (&send_sizes[0], 1, MPI_INT, 
		  &recv_sizes[0], 1, MPI_INT, comm);

      // in-place is buggy !
//    MPI_Alltoall (MPI_IN_PLACE, 1, MPI_INT, 
//		  &recv_sizes[0], 1, MPI_INT, comm);


    for (int i = 0; i < ntasks; i++)
      recv_data.SetEntrySize (i, recv_sizes[i], sizeof(T));
    
    NgArray<MPI_Request> requests;
    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank && send_data[dest].Size())
	requests.Append (MyMPI_ISend (send_data[dest], dest, tag, comm));

    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank && recv_data[dest].Size())
	requests.Append (MyMPI_IRecv (recv_data[dest], dest, tag, comm));

    // MPI_Barrier (comm);
    MPI_Waitall (requests.Size(), &requests[0], MPI_STATUS_IGNORE);
  }







  extern void MyMPI_SendCmd (const char * cmd);
  extern string MyMPI_RecvCmd ();




  template <class T>
  inline void MyMPI_Bcast (T & s, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Bcast (&s, 1, MyGetMPIType<T>(), 0, comm);
  }

  template <class T>
  inline void MyMPI_Bcast (NgArray<T, 0> & s, NgMPI_Comm comm /* = ng_comm */)
  {
    int size = s.Size();
    MyMPI_Bcast (size, comm);
    // if (MyMPI_GetId(comm) != 0) s.SetSize (size);
    if (comm.Rank() != 0) s.SetSize (size);
    MPI_Bcast (&s[0], size, MyGetMPIType<T>(), 0, comm);
  }

  template <class T>
  inline void MyMPI_Bcast (NgArray<T, 0> & s, int root, MPI_Comm comm /* = ng_comm */)
  {
    int id;
    MPI_Comm_rank(comm, &id);

    int size = s.Size();
    MPI_Bcast (&size, 1, MPI_INT, root, comm);
    if (id != root) s.SetSize (size);
    if ( !size ) return;
    MPI_Bcast (&s[0], size, MyGetMPIType<T>(), root, comm);
  }

  template <class T, class T2>
  inline void MyMPI_Allgather (const T & send, NgFlatArray<T2> recv, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Allgather( const_cast<T*> (&send), 1, MyGetMPIType<T>(), &recv[0], 1, MyGetMPIType<T2>(), comm);
  }

  template <class T, class T2>
  inline void MyMPI_Alltoall (NgFlatArray<T> send, NgFlatArray<T2> recv, MPI_Comm comm /* = ng_comm */)
  {
    MPI_Alltoall( &send[0], 1, MyGetMPIType<T>(), &recv[0], 1, MyGetMPIType<T2>(), comm);
  }

//   template <class T, class T2>
//   inline void MyMPI_Alltoall_Block (NgFlatArray<T> send, NgFlatArray<T2> recv, int blocklen, MPI_Comm comm = ng_comm)
//   {
//     MPI_Alltoall( &send[0], blocklen, MyGetMPIType<T>(), &recv[0], blocklen, MyGetMPIType<T2>(), comm);
//   }



  /*
  inline void MyMPI_Send (  int *& s, int len,  int dest, int tag)
  {
    int hlen = len;
    MPI_Send( &hlen, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    MPI_Send( s, len, MPI_INT, dest, tag, MPI_COMM_WORLD);
  }


  inline void MyMPI_Recv ( int *& s, int & len, int src, int tag)
  {
    MPI_Status status;
    MPI_Recv( &len, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
    if ( s ) 
      delete [] s;
    s = new int [len];
    MPI_Recv( s, len, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
  }



  inline void MyMPI_Send ( double * s, int len,  int dest, int tag)
  {
     MPI_Send( &len, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
     MPI_Send( s, len, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  }


  inline void MyMPI_Recv ( double *& s, int & len, int src, int tag)
  {
    MPI_Status status;
    MPI_Recv( &len, 1, MPI_INT, src, tag, MPI_COMM_WORLD, &status);
    if ( s )
      delete [] s;
    s = new double [len];
    MPI_Recv( s, len, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);
  }
  */

#endif // PARALLEL

}

#endif

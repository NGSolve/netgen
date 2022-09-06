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

#ifdef OLD
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
#endif

  
  enum { MPI_TAG_CMD = 110 };
  enum { MPI_TAG_MESH = 210 };
  enum { MPI_TAG_VIS = 310 };

#ifdef PARALLEL

  [[deprecated("mympi_send int, use comm.Send instead")]]            
  inline void MyMPI_Send (int i, int dest, int tag, MPI_Comm comm)
  {
    int hi = i;
    MPI_Send( &hi, 1, MPI_INT, dest, tag, comm);
  }
  
  [[deprecated("mympi_revc int, use comm.Recv instead")]]            
  inline void MyMPI_Recv (int & i, int src, int tag, MPI_Comm comm)
  {
    MPI_Status status;
    MPI_Recv( &i, 1, MPI_INT, src, tag, comm, &status);
  }

  [[deprecated("mympi_send string, use comm.Send instead")]]              
  inline void MyMPI_Send (const string & s, int dest, int tag, MPI_Comm comm)
  {
    MPI_Send( const_cast<char*> (s.c_str()), s.length(), MPI_CHAR, dest, tag, comm);
  }

  [[deprecated("mympi_revc string, use comm.Recv instead")]]              
  inline void MyMPI_Recv (string & s, int src, int tag, MPI_Comm comm)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count (&status, MPI_CHAR, &len);
    s.assign (len, ' ');
    MPI_Recv( &s[0], len, MPI_CHAR, src, tag, comm, &status);
  }

 

  template <class T, int BASE>
  [[deprecated("mympi_send ngflatarray, use comm.send instead")]]              
  inline void MyMPI_Send (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Comm comm)
  {
    MPI_Send( &s.First(), s.Size(), GetMPIType<T>(), dest, tag, comm);
  }

  template <class T, int BASE>
  [[deprecated("mympi_recv ngflatarray, use comm.Recv instead")]]                
  inline void MyMPI_Recv ( NgFlatArray<T, BASE> s, int src, int tag, MPI_Comm comm)
  {
    MPI_Status status;
    MPI_Recv( &s.First(), s.Size(), GetMPIType<T>(), src, tag, comm, &status);
  }

  template <class T, int BASE>
  [[deprecated("use ngcore - Array instead")]]
  inline void MyMPI_Recv ( NgArray <T, BASE> & s, int src, int tag, MPI_Comm comm)
  {
    MPI_Status status;
    int len;
    MPI_Probe (src, tag, comm, &status);
    MPI_Get_count (&status, GetMPIType<T>(), &len);

    s.SetSize (len);
    MPI_Recv( &s.First(), len, GetMPIType<T>(), src, tag, comm, &status);
  }

  template <class T, int BASE>
  [[deprecated("use ngcore - Array instead")]]
  inline int MyMPI_Recv ( NgArray <T, BASE> & s, int tag, MPI_Comm comm)
  {
    MPI_Status status;
    int len;
    MPI_Probe (MPI_ANY_SOURCE, tag, comm, &status);

    int src = status.MPI_SOURCE;

    MPI_Get_count (&status, GetMPIType<T>(), &len);

    s.SetSize (len);
    MPI_Recv( &s.First(), len, GetMPIType<T>(), src, tag, comm, &status);

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
  [[deprecated("mympi_isend ngflatarray, use comm.send instead")]]
  [[deprecated("use ngcore - Array instead")]]
  inline MPI_Request MyMPI_ISend (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Comm comm)
  {
    MPI_Request request;
    MPI_Isend( &s.First(), s.Size(), GetMPIType<T>(), dest, tag, comm, &request);
    return request;
  }

  template <class T, int BASE>
  [[deprecated("mympi_irecv ngflatarray, use comm.recv instead")]]                
  inline MPI_Request MyMPI_IRecv (NgFlatArray<T, BASE> s, int dest, int tag, MPI_Comm comm)
  {
    MPI_Request request;
    MPI_Irecv( &s.First(), s.Size(), GetMPIType<T>(), dest, tag, comm, &request);
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

  template <typename T>
  [[deprecated("do we need that ? ")]]       
  inline void MyMPI_ExchangeTable (TABLE<T> & send_data, 
				   TABLE<T> & recv_data, int tag,
				   const NgMPI_Comm & comm)
  {
    int rank = comm.Rank();
    int ntasks = comm.Size();
    
    Array<int> send_sizes(ntasks);
    Array<int> recv_sizes(ntasks);
    for (int i = 0; i < ntasks; i++)
      send_sizes[i] = send_data[i].Size();

    comm.AllToAll (send_sizes, recv_sizes);
    
    for (int i = 0; i < ntasks; i++)
      recv_data.SetEntrySize (i, recv_sizes[i], sizeof(T));
    
    Array<MPI_Request> requests;
    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank && send_data[dest].Size())
        requests.Append (comm.ISend (FlatArray<T>(send_data[dest]), dest, tag));

    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank && recv_data[dest].Size())
        requests.Append (comm.IRecv (FlatArray<T>(recv_data[dest]), dest, tag));

    MyMPI_WaitAll (requests);
  }


  template <typename T>
  [[deprecated("do we need that ? ")]]       
  inline void MyMPI_ExchangeTable (DynamicTable<T> & send_data, 
				   DynamicTable<T> & recv_data, int tag,
				   const NgMPI_Comm & comm)
  {
    int rank = comm.Rank();
    int ntasks = comm.Size();
    
    Array<int> send_sizes(ntasks);
    Array<int> recv_sizes(ntasks);
    for (int i = 0; i < ntasks; i++)
      send_sizes[i] = send_data[i].Size();

    comm.AllToAll (send_sizes, recv_sizes);
    
    // for (int i = 0; i < ntasks; i++)
    // recv_data.SetEntrySize (i, recv_sizes[i], sizeof(T));
    recv_data = DynamicTable<T> (recv_sizes, true);
    
    Array<MPI_Request> requests;
    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank && send_data[dest].Size())
        requests.Append (comm.ISend (FlatArray<T>(send_data[dest]), dest, tag));

    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank && recv_data[dest].Size())
        requests.Append (comm.IRecv (FlatArray<T>(recv_data[dest]), dest, tag));

    MyMPI_WaitAll (requests);
  }



  
  [[deprecated("do we still send commands?")]]                      
  DLL_HEADER void MyMPI_SendCmd (const char * cmd);
  [[deprecated("do we still send commands?")]]                        
  extern string MyMPI_RecvCmd ();


  template <class T>
  [[deprecated("use comm.BCast instead")]]                      
  inline void MyMPI_Bcast (T & s, MPI_Comm comm)
  {
    MPI_Bcast (&s, 1, GetMPIType<T>(), 0, comm);
  }

  template <class T>
  [[deprecated("use comm.BCast instead")]]                        
  inline void MyMPI_Bcast (NgArray<T, 0> & s, NgMPI_Comm comm)
  {
    int size = s.Size();
    // MyMPI_Bcast (size, comm);
    comm.Bcast(size);
    // if (MyMPI_GetId(comm) != 0) s.SetSize (size);
    if (comm.Rank() != 0) s.SetSize (size);
    MPI_Bcast (&s[0], size, GetMPIType<T>(), 0, comm);
  }

  template <class T>
  [[deprecated("use comm.BCast instead")]]                        
  inline void MyMPI_Bcast (NgArray<T, 0> & s, int root, MPI_Comm comm)
  {
    int id;
    MPI_Comm_rank(comm, &id);

    int size = s.Size();
    MPI_Bcast (&size, 1, MPI_INT, root, comm);
    if (id != root) s.SetSize (size);
    if ( !size ) return;
    MPI_Bcast (&s[0], size, GetMPIType<T>(), root, comm);
  }

  template <class T, class T2>
  [[deprecated("mympi_allgather deprecated, use comm.allgather")]]                    
  inline void MyMPI_Allgather (const T & send, NgFlatArray<T2> recv, MPI_Comm comm)
  {
    MPI_Allgather( const_cast<T*> (&send), 1, GetMPIType<T>(), &recv[0], 1, GetMPIType<T2>(), comm);
  }

  template <class T, class T2>
  [[deprecated("mympi_alltoall deprecated, use comm.alltoall")]]                  
  inline void MyMPI_Alltoall (NgFlatArray<T> send, NgFlatArray<T2> recv, MPI_Comm comm)
  {
    MPI_Alltoall( &send[0], 1, GetMPIType<T>(), &recv[0], 1, GetMPIType<T2>(), comm);
  }


#else
  template <typename T>
  [[deprecated("do we need that ? ")]]       
  inline void MyMPI_ExchangeTable (TABLE<T> & send_data, 
				   TABLE<T> & recv_data, int tag,
				   const NgMPI_Comm & comm)
  {
    ;
  }

  template <typename T>
  [[deprecated("do we need that ? ")]]       
  inline void MyMPI_ExchangeTable (DynamicTable<T> & send_data, 
				   DynamicTable<T> & recv_data, int tag,
				   const NgMPI_Comm & comm)
  { ; } 
#endif // PARALLEL

}

#endif

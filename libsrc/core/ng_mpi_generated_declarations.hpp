#ifdef NG_MPI_WRAPPER
NGCORE_API extern double (*NG_MPI_Wtime)();
NGCORE_API extern int (*NG_MPI_Allgather)(void*, int, NG_MPI_Datatype, void*, int, NG_MPI_Datatype, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Allreduce)(void*, void*, int, NG_MPI_Datatype, NG_MPI_Op, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Alltoall)(void*, int, NG_MPI_Datatype, void*, int, NG_MPI_Datatype, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Barrier)(NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Bcast)(void*, int, NG_MPI_Datatype, int, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Ibcast)(void*, int, NG_MPI_Datatype, int, NG_MPI_Comm, NG_MPI_Request*);
NGCORE_API extern int (*NG_MPI_Comm_c2f)(NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Comm_create)(NG_MPI_Comm, NG_MPI_Group, NG_MPI_Comm*);
NGCORE_API extern int (*NG_MPI_Comm_create_group)(NG_MPI_Comm, NG_MPI_Group, int, NG_MPI_Comm*);
NGCORE_API extern int (*NG_MPI_Comm_free)(NG_MPI_Comm*);
NGCORE_API extern int (*NG_MPI_Comm_group)(NG_MPI_Comm, NG_MPI_Group*);
NGCORE_API extern int (*NG_MPI_Comm_rank)(NG_MPI_Comm, int*);
NGCORE_API extern int (*NG_MPI_Comm_size)(NG_MPI_Comm, int*);
NGCORE_API extern int (*NG_MPI_Finalize)();
NGCORE_API extern int (*NG_MPI_Gather)(void*, int, NG_MPI_Datatype, void*, int, NG_MPI_Datatype, int, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Gatherv)(void*, int, NG_MPI_Datatype, void*, int*, int*, NG_MPI_Datatype, int, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Get_count)(NG_MPI_Status*, NG_MPI_Datatype, int*);
NGCORE_API extern int (*NG_MPI_Get_processor_name)(char*, int*);
NGCORE_API extern int (*NG_MPI_Group_incl)(NG_MPI_Group, int, int*, NG_MPI_Group*);
NGCORE_API extern int (*NG_MPI_Init)(int*, char***);
NGCORE_API extern int (*NG_MPI_Init_thread)(int*, char***, int, int*);
NGCORE_API extern int (*NG_MPI_Initialized)(int*);
NGCORE_API extern int (*NG_MPI_Iprobe)(int, int, NG_MPI_Comm, int*, NG_MPI_Status*);
NGCORE_API extern int (*NG_MPI_Irecv)(void*, int, NG_MPI_Datatype, int, int, NG_MPI_Comm, NG_MPI_Request*);
NGCORE_API extern int (*NG_MPI_Isend)(void*, int, NG_MPI_Datatype, int, int, NG_MPI_Comm, NG_MPI_Request*);
NGCORE_API extern int (*NG_MPI_Probe)(int, int, NG_MPI_Comm, NG_MPI_Status*);
NGCORE_API extern int (*NG_MPI_Query_thread)(int*);
NGCORE_API extern int (*NG_MPI_Recv)(void*, int, NG_MPI_Datatype, int, int, NG_MPI_Comm, NG_MPI_Status*);
NGCORE_API extern int (*NG_MPI_Recv_init)(void*, int, NG_MPI_Datatype, int, int, NG_MPI_Comm, NG_MPI_Request*);
NGCORE_API extern int (*NG_MPI_Reduce)(void*, void*, int, NG_MPI_Datatype, NG_MPI_Op, int, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Reduce_local)(void*, void*, int, NG_MPI_Datatype, NG_MPI_Op);
NGCORE_API extern int (*NG_MPI_Request_free)(NG_MPI_Request*);
NGCORE_API extern int (*NG_MPI_Scatter)(void*, int, NG_MPI_Datatype, void*, int, NG_MPI_Datatype, int, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Send)(void*, int, NG_MPI_Datatype, int, int, NG_MPI_Comm);
NGCORE_API extern int (*NG_MPI_Send_init)(void*, int, NG_MPI_Datatype, int, int, NG_MPI_Comm, NG_MPI_Request*);
NGCORE_API extern int (*NG_MPI_Startall)(int, NG_MPI_Request*);
NGCORE_API extern int (*NG_MPI_Type_commit)(NG_MPI_Datatype*);
NGCORE_API extern int (*NG_MPI_Type_contiguous)(int, NG_MPI_Datatype, NG_MPI_Datatype*);
NGCORE_API extern int (*NG_MPI_Type_create_resized)(NG_MPI_Datatype, NG_MPI_Aint, NG_MPI_Aint, NG_MPI_Datatype*);
NGCORE_API extern int (*NG_MPI_Type_create_struct)(int, int*, NG_MPI_Aint*, NG_MPI_Datatype*, NG_MPI_Datatype*);
NGCORE_API extern int (*NG_MPI_Type_free)(NG_MPI_Datatype*);
NGCORE_API extern int (*NG_MPI_Type_get_extent)(NG_MPI_Datatype, NG_MPI_Aint*, NG_MPI_Aint*);
NGCORE_API extern int (*NG_MPI_Type_indexed)(int, int*, int*, NG_MPI_Datatype, NG_MPI_Datatype*);
NGCORE_API extern int (*NG_MPI_Type_size)(NG_MPI_Datatype, int*);
NGCORE_API extern int (*NG_MPI_Wait)(NG_MPI_Request*, NG_MPI_Status*);
NGCORE_API extern int (*NG_MPI_Waitall)(int, NG_MPI_Request*, NG_MPI_Status*);
NGCORE_API extern int (*NG_MPI_Waitany)(int, NG_MPI_Request*, int*, NG_MPI_Status*);
NGCORE_API extern NG_MPI_Comm NG_MPI_COMM_NULL;
NGCORE_API extern NG_MPI_Comm NG_MPI_COMM_WORLD;
NGCORE_API extern NG_MPI_Datatype NG_MPI_CHAR;
NGCORE_API extern NG_MPI_Datatype NG_MPI_CXX_DOUBLE_COMPLEX;
NGCORE_API extern NG_MPI_Datatype NG_MPI_C_BOOL;
NGCORE_API extern NG_MPI_Datatype NG_MPI_DATATYPE_NULL;
NGCORE_API extern NG_MPI_Datatype NG_MPI_DOUBLE;
NGCORE_API extern NG_MPI_Datatype NG_MPI_FLOAT;
NGCORE_API extern NG_MPI_Datatype NG_MPI_INT;
NGCORE_API extern NG_MPI_Datatype NG_MPI_SHORT;
NGCORE_API extern NG_MPI_Datatype NG_MPI_UINT64_T;
NGCORE_API extern NG_MPI_Op NG_MPI_LOR;
NGCORE_API extern NG_MPI_Op NG_MPI_MAX;
NGCORE_API extern NG_MPI_Op NG_MPI_MIN;
NGCORE_API extern NG_MPI_Op NG_MPI_SUM;
NGCORE_API extern NG_MPI_Request NG_MPI_REQUEST_NULL;
NGCORE_API extern NG_MPI_Status* NG_MPI_STATUSES_IGNORE;
NGCORE_API extern NG_MPI_Status* NG_MPI_STATUS_IGNORE;
NGCORE_API extern int NG_MPI_ANY_SOURCE;
NGCORE_API extern int NG_MPI_ANY_TAG;
NGCORE_API extern int NG_MPI_MAX_PROCESSOR_NAME;
NGCORE_API extern int NG_MPI_PROC_NULL;
NGCORE_API extern int NG_MPI_ROOT;
NGCORE_API extern int NG_MPI_SUBVERSION;
NGCORE_API extern int NG_MPI_THREAD_MULTIPLE;
NGCORE_API extern int NG_MPI_THREAD_SINGLE;
NGCORE_API extern int NG_MPI_VERSION;
NGCORE_API extern void* NG_MPI_IN_PLACE;
#else  // NG_MPI_WRAPPER
#define NG_MPI_Wtime MPI_Wtime
#define NG_MPI_Allgather MPI_Allgather
#define NG_MPI_Allreduce MPI_Allreduce
#define NG_MPI_Alltoall MPI_Alltoall
#define NG_MPI_Barrier MPI_Barrier
#define NG_MPI_Bcast MPI_Bcast
#define NG_MPI_Ibcast MPI_Ibcast
#define NG_MPI_Comm_c2f MPI_Comm_c2f
#define NG_MPI_Comm_create MPI_Comm_create
#define NG_MPI_Comm_create_group MPI_Comm_create_group
#define NG_MPI_Comm_free MPI_Comm_free
#define NG_MPI_Comm_group MPI_Comm_group
#define NG_MPI_Comm_rank MPI_Comm_rank
#define NG_MPI_Comm_size MPI_Comm_size
#define NG_MPI_Finalize MPI_Finalize
#define NG_MPI_Gather MPI_Gather
#define NG_MPI_Gatherv MPI_Gatherv
#define NG_MPI_Get_count MPI_Get_count
#define NG_MPI_Get_processor_name MPI_Get_processor_name
#define NG_MPI_Group_incl MPI_Group_incl
#define NG_MPI_Init MPI_Init
#define NG_MPI_Init_thread MPI_Init_thread
#define NG_MPI_Initialized MPI_Initialized
#define NG_MPI_Iprobe MPI_Iprobe
#define NG_MPI_Irecv MPI_Irecv
#define NG_MPI_Isend MPI_Isend
#define NG_MPI_Probe MPI_Probe
#define NG_MPI_Query_thread MPI_Query_thread
#define NG_MPI_Recv MPI_Recv
#define NG_MPI_Recv_init MPI_Recv_init
#define NG_MPI_Reduce MPI_Reduce
#define NG_MPI_Reduce_local MPI_Reduce_local
#define NG_MPI_Request_free MPI_Request_free
#define NG_MPI_Scatter MPI_Scatter
#define NG_MPI_Send MPI_Send
#define NG_MPI_Send_init MPI_Send_init
#define NG_MPI_Startall MPI_Startall
#define NG_MPI_Type_commit MPI_Type_commit
#define NG_MPI_Type_contiguous MPI_Type_contiguous
#define NG_MPI_Type_create_resized MPI_Type_create_resized
#define NG_MPI_Type_create_struct MPI_Type_create_struct
#define NG_MPI_Type_free MPI_Type_free
#define NG_MPI_Type_get_extent MPI_Type_get_extent
#define NG_MPI_Type_indexed MPI_Type_indexed
#define NG_MPI_Type_size MPI_Type_size
#define NG_MPI_Wait MPI_Wait
#define NG_MPI_Waitall MPI_Waitall
#define NG_MPI_Waitany MPI_Waitany
#define NG_MPI_COMM_NULL MPI_COMM_NULL
#define NG_MPI_COMM_WORLD MPI_COMM_WORLD
#define NG_MPI_CHAR MPI_CHAR
#define NG_MPI_CXX_DOUBLE_COMPLEX MPI_CXX_DOUBLE_COMPLEX
#define NG_MPI_C_BOOL MPI_C_BOOL
#define NG_MPI_DATATYPE_NULL MPI_DATATYPE_NULL
#define NG_MPI_DOUBLE MPI_DOUBLE
#define NG_MPI_FLOAT MPI_FLOAT
#define NG_MPI_INT MPI_INT
#define NG_MPI_SHORT MPI_SHORT
#define NG_MPI_UINT64_T MPI_UINT64_T
#define NG_MPI_LOR MPI_LOR
#define NG_MPI_MAX MPI_MAX
#define NG_MPI_MIN MPI_MIN
#define NG_MPI_SUM MPI_SUM
#define NG_MPI_REQUEST_NULL MPI_REQUEST_NULL
#define NG_MPI_STATUSES_IGNORE MPI_STATUSES_IGNORE
#define NG_MPI_STATUS_IGNORE MPI_STATUS_IGNORE
#define NG_MPI_ANY_SOURCE MPI_ANY_SOURCE
#define NG_MPI_ANY_TAG MPI_ANY_TAG
#define NG_MPI_MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#define NG_MPI_PROC_NULL MPI_PROC_NULL
#define NG_MPI_ROOT MPI_ROOT
#define NG_MPI_SUBVERSION MPI_SUBVERSION
#define NG_MPI_THREAD_MULTIPLE MPI_THREAD_MULTIPLE
#define NG_MPI_THREAD_SINGLE MPI_THREAD_SINGLE
#define NG_MPI_VERSION MPI_VERSION
#define NG_MPI_IN_PLACE MPI_IN_PLACE
#endif // NG_MPI_WRAPPER

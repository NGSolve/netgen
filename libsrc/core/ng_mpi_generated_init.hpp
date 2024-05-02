NG_MPI_Init = [](int* arg0, char*** arg1)->int { return MPI_Init( ng2mpi(arg0),  ng2mpi(arg1)); };
NG_MPI_Initialized = [](int* arg0)->int { return MPI_Initialized( ng2mpi(arg0)); };
NG_MPI_Finalize = []()->int { return MPI_Finalize(); };
NG_MPI_Barrier = [](NG_MPI_Comm arg0)->int { return MPI_Barrier( ng2mpi(arg0)); };
NG_MPI_Comm_free = [](NG_MPI_Comm* arg0)->int { return MPI_Comm_free( ng2mpi(arg0)); };
NG_MPI_Comm_rank = [](NG_MPI_Comm arg0, int* arg1)->int { return MPI_Comm_rank( ng2mpi(arg0),  ng2mpi(arg1)); };
NG_MPI_Comm_size = [](NG_MPI_Comm arg0, int* arg1)->int { return MPI_Comm_size( ng2mpi(arg0),  ng2mpi(arg1)); };
NG_MPI_Type_commit = [](NG_MPI_Datatype* arg0)->int { return MPI_Type_commit( ng2mpi(arg0)); };
NG_MPI_Waitall = [](int arg0, NG_MPI_Request* arg1, NG_MPI_Status* arg2)->int { return MPI_Waitall( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2)); };
NG_MPI_Waitany = [](int arg0, NG_MPI_Request* arg1, int* arg2, NG_MPI_Status* arg3)->int { return MPI_Waitany( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2),  ng2mpi(arg3)); };
NG_MPI_Wait = [](NG_MPI_Request* arg0, NG_MPI_Status* arg1)->int { return MPI_Wait( ng2mpi(arg0),  ng2mpi(arg1)); };
NG_MPI_Send = [](void* arg0, int arg1, NG_MPI_Datatype arg2, int arg3, int arg4, NG_MPI_Comm arg5)->int { return MPI_Send( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2),  ng2mpi(arg3),  ng2mpi(arg4),  ng2mpi(arg5)); };
NG_MPI_Probe = [](int arg0, int arg1, NG_MPI_Comm arg2, NG_MPI_Status* arg3)->int { return MPI_Probe( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2),  ng2mpi(arg3)); };
NG_MPI_Recv = [](void* arg0, int arg1, NG_MPI_Datatype arg2, int arg3, int arg4, NG_MPI_Comm arg5, NG_MPI_Status* arg6)->int { return MPI_Recv( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2),  ng2mpi(arg3),  ng2mpi(arg4),  ng2mpi(arg5),  ng2mpi(arg6)); };
NG_MPI_Get_count = [](NG_MPI_Status* arg0, NG_MPI_Datatype arg1, int* arg2)->int { return MPI_Get_count( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2)); };
NG_MPI_Bcast = [](void* arg0, int arg1, NG_MPI_Datatype arg2, int arg3, NG_MPI_Comm arg4)->int { return MPI_Bcast( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2),  ng2mpi(arg3),  ng2mpi(arg4)); };
NG_MPI_Comm_group = [](NG_MPI_Comm arg0, NG_MPI_Group* arg1)->int { return MPI_Comm_group( ng2mpi(arg0),  ng2mpi(arg1)); };
NG_MPI_Group_incl = [](NG_MPI_Group arg0, int arg1, int* arg2, NG_MPI_Group* arg3)->int { return MPI_Group_incl( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2),  ng2mpi(arg3)); };
NG_MPI_Comm_create_group = [](NG_MPI_Comm arg0, NG_MPI_Group arg1, int arg2, NG_MPI_Comm* arg3)->int { return MPI_Comm_create_group( ng2mpi(arg0),  ng2mpi(arg1),  ng2mpi(arg2),  ng2mpi(arg3)); };
NG_MPI_Get_processor_name = [](char* arg0, int* arg1)->int { return MPI_Get_processor_name( ng2mpi(arg0),  ng2mpi(arg1)); };
NG_MPI_COMM_WORLD = mpi2ng(MPI_COMM_WORLD);
NG_MPI_STATUS_IGNORE = mpi2ng(MPI_STATUS_IGNORE);
NG_MPI_STATUSES_IGNORE = mpi2ng(MPI_STATUSES_IGNORE);
NG_MPI_INT = mpi2ng(MPI_INT);
NG_MPI_SHORT = mpi2ng(MPI_SHORT);
NG_MPI_CHAR = mpi2ng(MPI_CHAR);
NG_MPI_UINT64_T = mpi2ng(MPI_UINT64_T);
NG_MPI_DOUBLE = mpi2ng(MPI_DOUBLE);
NG_MPI_C_BOOL = mpi2ng(MPI_C_BOOL);
NG_MPI_IN_PLACE = mpi2ng(MPI_IN_PLACE);
NG_MPI_MAX_PROCESSOR_NAME = mpi2ng(MPI_MAX_PROCESSOR_NAME);

functions = [
        ("double", "MPI_Wtime"),
        ("int", "MPI_Allgather", "void*", "int", "MPI_Datatype", "void*", "int", "MPI_Datatype", "MPI_Comm"),
        ("int", "MPI_Allreduce", "void*", "void*", "int", "MPI_Datatype", "MPI_Op", "MPI_Comm"),
        ("int", "MPI_Alltoall", "void*", "int", "MPI_Datatype", "void*", "int", "MPI_Datatype", "MPI_Comm"),
        ("int", "MPI_Barrier", "MPI_Comm"),
        ("int", "MPI_Bcast", "void*", "int", "MPI_Datatype", "int", "MPI_Comm"),
        ("int", "MPI_Comm_create_group", "MPI_Comm", "MPI_Group", "int", "MPI_Comm*"),
        ("int", "MPI_Comm_free", "MPI_Comm*"),
        ("int", "MPI_Comm_group", "MPI_Comm", "MPI_Group*"),
        ("int", "MPI_Comm_rank", "MPI_Comm", "int*"),
        ("int", "MPI_Comm_size", "MPI_Comm", "int*"),
        ("int", "MPI_Finalize"),
        ("int", "MPI_Gather", "void*", "int", "MPI_Datatype", "void*", "int", "MPI_Datatype", "int", "MPI_Comm"),
        ("int", "MPI_Get_count", "MPI_Status*", "MPI_Datatype", "int*"),
        ("int", "MPI_Get_processor_name", "char*", "int*"),
        ("int", "MPI_Group_incl", "MPI_Group", "int", "int*", "MPI_Group*"),
        ("int", "MPI_Init", "int*", "char***"),
        ("int", "MPI_Init_thread", "int*", "char***", "int", "int*"),
        ("int", "MPI_Initialized", "int*"),
        ("int", "MPI_Iprobe", "int", "int", "MPI_Comm", "int*", "MPI_Status*"),
        ("int", "MPI_Irecv", "void*", "int", "MPI_Datatype", "int", "int", "MPI_Comm", "MPI_Request*"),
        ("int", "MPI_Isend", "void*", "int", "MPI_Datatype", "int", "int", "MPI_Comm", "MPI_Request*"),
        ("int", "MPI_Probe", "int", "int", "MPI_Comm", "MPI_Status*"),
        ("int", "MPI_Query_thread", "int*"),
        ("int", "MPI_Recv", "void*", "int", "MPI_Datatype", "int", "int", "MPI_Comm", "MPI_Status*"),
        ("int", "MPI_Reduce", "void*", "void*", "int", "MPI_Datatype", "MPI_Op", "int", "MPI_Comm"),
        ("int", "MPI_Reduce_local", "void*", "void*", "int", "MPI_Datatype", "MPI_Op"),
        ("int", "MPI_Request_free", "MPI_Request*"),
        ("int", "MPI_Scatter", "void*", "int", "MPI_Datatype", "void*", "int", "MPI_Datatype", "int", "MPI_Comm"),
        ("int", "MPI_Send", "void*", "int", "MPI_Datatype", "int", "int", "MPI_Comm"),
        ("int", "MPI_Type_commit", "MPI_Datatype*"),
        ("int", "MPI_Type_contiguous", "int", "MPI_Datatype", "MPI_Datatype*"),
        ("int", "MPI_Type_create_resized", "MPI_Datatype", "MPI_Aint", "MPI_Aint", "MPI_Datatype*"),
        ("int", "MPI_Type_create_struct", "int", "int*:0", "MPI_Aint*:0", "MPI_Datatype*:0", "MPI_Datatype*"),
        ("int", "MPI_Type_free", "MPI_Datatype*"),
        ("int", "MPI_Type_get_extent", "MPI_Datatype", "MPI_Aint*", "MPI_Aint*"),
        ("int", "MPI_Type_indexed", "int", "int*:0", "int*:0", "MPI_Datatype", "MPI_Datatype*"),
        ("int", "MPI_Type_size", "MPI_Datatype*", "int*"),
        ("int", "MPI_Wait", "MPI_Request*", "MPI_Status*"),
        ("int", "MPI_Waitall", "int", "MPI_Request*:0", "MPI_Status*"),
        ("int", "MPI_Waitany", "int", "MPI_Request*:0", "int*", "MPI_Status*"),
        ]

constants = [
        ("MPI_Comm", "MPI_COMM_WORLD"),
        ("MPI_Datatype", "MPI_CHAR"),
        ("MPI_Datatype", "MPI_CXX_DOUBLE_COMPLEX"),
        ("MPI_Datatype", "MPI_C_BOOL"),
        ("MPI_Datatype", "MPI_DATATYPE_NULL"),
        ("MPI_Datatype", "MPI_DOUBLE"),
        ("MPI_Datatype", "MPI_INT"),
        ("MPI_Datatype", "MPI_SHORT"),
        ("MPI_Datatype", "MPI_UINT64_T"),
        ("MPI_Op", "MPI_LOR"),
        ("MPI_Op", "MPI_MAX"),
        ("MPI_Op", "MPI_MIN"),
        ("MPI_Op", "MPI_SUM"),
        ("MPI_Status*", "MPI_STATUSES_IGNORE"),
        ("MPI_Status*", "MPI_STATUS_IGNORE"),
        ("int", "MPI_ANY_SOURCE"),
        ("int", "MPI_ANY_TAG"),
        ("int", "MPI_MAX_PROCESSOR_NAME"),
        ("int", "MPI_PROC_NULL"),
        ("int", "MPI_ROOT"),
        ("int", "MPI_SUBVERSION"),
        ("int", "MPI_THREAD_MULTIPLE"),
        ("int", "MPI_THREAD_SINGLE"),
        ("int", "MPI_VERSION"),
        ("void*", "MPI_IN_PLACE"),
]

def get_args(f, counts=False):
    args = []
    for arg in f[2:]:
        has_count = ':' in arg
        if has_count:
            s, count = arg.split(':')
            count = int(count)
        else:
            s = arg
            count = None
        if s.startswith("MPI_"):
            s = "NG_" + s
        if counts:
            args.append((s, count))
        else:
            args.append(s)
    return args

def generate_declarations():
    code = ""
    for f in functions:
        ret = f[0]
        name = f[1]
        args = ", ".join(get_args(f))
        code += f"NGCORE_API extern {ret} (*NG_{name})({args});\n"

    for typ, name in constants:
        if typ.startswith("MPI_"):
            typ = "NG_" + typ
        code += f"NGCORE_API extern {typ} NG_{name};\n"

    with open("ng_mpi_generated_declarations.hpp", "w") as f:
        f.write(code)

def generate_dummy_init():
    code = ""
    for f in functions:
        ret = f[0]
        name = f[1]
        args = ", ".join(get_args(f))
        code += f"decltype(NG_{name}) NG_{name} = []({args})->{ret} {{ throw no_mpi(); }};\n"

    for typ, name in constants:
        if typ.startswith("MPI_"):
            typ = "NG_" + typ
        code += f"{typ} NG_{name} = 0;\n"

    with open("ng_mpi_generated_dummy_init.hpp", "w") as f:
        f.write(code)

def generate_init():
    code = ""
    for f in functions:
        ret = f[0]
        name = f[1]
        args = get_args(f, counts=True)
        in_args  =''
        call_args = ''
        for i in range(len(args)):
            arg, count = args[i]
            if i > 0:
                in_args += ', '
                call_args += ', '
            in_args += arg + f" arg{i}"
            if not arg.startswith("NG_"):
                # plain type (like int, int *, etc.), just pass the argument along
                call_args += f" arg{i}"
            elif count is None:
                # MPI type (by value or pointer), but just one object, no arrays
                call_args += f" ng2mpi(arg{i})"
            else:
                # arrays of MPI types, we need to copy them due to incompatible size
                call_args += f" ng2mpi(arg{i}, arg{count})"
        code += f"NG_{name} = []({in_args})->{ret} {{ return {name}({call_args}); }};\n"

    for _, name in constants:
        code += f"NG_{name} = mpi2ng({name});\n"

    with open("ng_mpi_generated_init.hpp", "w") as f:
        f.write(code)

if __name__ == "__main__":
    generate_declarations()
    generate_dummy_init()
    generate_init()

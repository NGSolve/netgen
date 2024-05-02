functions = [
        ("int", "MPI_Init", "int*", "char***"),
        ("int", "MPI_Initialized", "int*"),
        ("int", "MPI_Finalize"),
        ("int", "MPI_Barrier", "MPI_Comm"),
        ("int", "MPI_Comm_free", "MPI_Comm*"),
        ("int", "MPI_Comm_rank", "MPI_Comm", "int*"),
        ("int", "MPI_Comm_size", "MPI_Comm", "int*"),
        ("int", "MPI_Type_commit", "MPI_Datatype*"),
        ("int", "MPI_Waitall", "int", "MPI_Request*", "MPI_Status*"),
        ("int", "MPI_Waitany", "int", "MPI_Request*", "int*", "MPI_Status*"),
        ("int", "MPI_Wait", "MPI_Request*", "MPI_Status*"),
        ("int", "MPI_Send", "void*", "int", "MPI_Datatype", "int", "int", "MPI_Comm"),
        ("int", "MPI_Probe", "int", "int", "MPI_Comm", "MPI_Status*"),
        ("int", "MPI_Recv", "void*", "int", "MPI_Datatype", "int", "int", "MPI_Comm", "MPI_Status*"),
        ("int", "MPI_Get_count", "MPI_Status*", "MPI_Datatype", "int*"),
        ("int", "MPI_Bcast", "void*", "int", "MPI_Datatype", "int", "MPI_Comm"),
        ("int", "MPI_Comm_group", "MPI_Comm", "MPI_Group*"),
        ("int", "MPI_Group_incl", "MPI_Group", "int", "int*", "MPI_Group*"),
        ("int", "MPI_Comm_create_group", "MPI_Comm", "MPI_Group", "int", "MPI_Comm*"),
        ("int", "MPI_Get_processor_name", "char*", "int*"),
        ]

constants = [
        ("MPI_Comm", "MPI_COMM_WORLD"),
        ("MPI_Status*", "MPI_STATUS_IGNORE"),
        ("MPI_Status*", "MPI_STATUSES_IGNORE"),
        ("MPI_Datatype", "MPI_INT"),
        ("MPI_Datatype", "MPI_SHORT"),
        ("MPI_Datatype", "MPI_CHAR"),
        ("MPI_Datatype", "MPI_UINT64_T"),
        ("MPI_Datatype", "MPI_DOUBLE"),
        ("MPI_Datatype", "MPI_C_BOOL"),
        ("void*", "MPI_IN_PLACE"),
        ("int", "MPI_MAX_PROCESSOR_NAME"),
]

def get_args(f):
    return ["NG_"+a if a.startswith("MPI_") else a for a in f[2:]]

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
        args = get_args(f)
        in_args  =''
        call_args = ''
        for i, a in enumerate(args):
            if i > 0:
                in_args += ', '
                call_args += ', '
            in_args += a + f" arg{i}"
            call_args += f" ng2mpi(arg{i})"
        code += f"NG_{name} = []({in_args})->{ret} {{ return {name}({call_args}); }};\n"

    for _, name in constants:
        code += f"NG_{name} = mpi2ng({name});\n"

    with open("ng_mpi_generated_init.hpp", "w") as f:
        f.write(code)

if __name__ == "__main__":
    generate_declarations()
    generate_dummy_init()
    generate_init()

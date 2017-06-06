add_definitions(-DNGINTERFACE_EXPORTS)
add_library(gen OBJECT
        array.cpp bitarray.cpp dynamicmem.cpp flags.cpp 
        hashtabl.cpp mystring.cpp ngexception.cpp optmem.cpp parthreads.cpp 
        profiler.cpp seti.cpp sort.cpp spbita2d.cpp symbolta.cpp table.cpp 
        mpi_interface.cpp gzstream.cpp
        )

set_target_properties( gen PROPERTIES POSITION_INDEPENDENT_CODE ON )

install( FILES ngexception.hpp DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel )

install(FILES
  archive_base.hpp array.hpp autodiff.hpp autoptr.hpp bitarray.hpp
  dynamicmem.hpp flags.hpp hashtabl.hpp mpi_interface.hpp myadt.hpp
  ngsimd.hpp mystring.hpp netgenout.hpp ngexception.hpp ngpython.hpp
  optmem.hpp parthreads.hpp profiler.hpp seti.hpp sort.hpp
  spbita2d.hpp stack.hpp symbolta.hpp table.hpp template.hpp
  gzstream.h
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/general COMPONENT netgen_devel
)

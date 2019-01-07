add_definitions(-DNGINTERFACE_EXPORTS)
add_library(gen OBJECT
        array.cpp bitarray.cpp dynamicmem.cpp flags.cpp 
        hashtabl.cpp mystring.cpp optmem.cpp parthreads.cpp
        seti.cpp sort.cpp spbita2d.cpp table.cpp
        mpi_interface.cpp gzstream.cpp
        )

set_target_properties( gen PROPERTIES POSITION_INDEPENDENT_CODE ON )

install(FILES
  array.hpp autodiff.hpp autoptr.hpp bitarray.hpp
  dynamicmem.hpp flags.hpp hashtabl.hpp mpi_interface.hpp myadt.hpp
  ngsimd.hpp mystring.hpp netgenout.hpp ngpython.hpp
  optmem.hpp parthreads.hpp seti.hpp sort.hpp
  spbita2d.hpp stack.hpp table.hpp template.hpp
  gzstream.h
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/general COMPONENT netgen_devel
)

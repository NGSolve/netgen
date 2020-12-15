add_definitions(-DNGINTERFACE_EXPORTS)
add_library(gen INTERFACE)
set(sdir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(gen INTERFACE
  ${sdir}/ngarray.cpp ${sdir}/ngbitarray.cpp ${sdir}/dynamicmem.cpp
        ${sdir}/hashtabl.cpp ${sdir}/mystring.cpp ${sdir}/optmem.cpp ${sdir}/parthreads.cpp
        ${sdir}/seti.cpp ${sdir}/sort.cpp ${sdir}/spbita2d.cpp ${sdir}/table.cpp
        ${sdir}/mpi_interface.cpp ${sdir}/gzstream.cpp
        )

install(FILES
  ngarray.hpp autodiff.hpp autoptr.hpp ngbitarray.hpp
  dynamicmem.hpp hashtabl.hpp mpi_interface.hpp myadt.hpp
  mystring.hpp netgenout.hpp ngpython.hpp
  optmem.hpp parthreads.hpp seti.hpp sort.hpp
  spbita2d.hpp stack.hpp table.hpp template.hpp
  gzstream.h
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/general COMPONENT netgen_devel
)

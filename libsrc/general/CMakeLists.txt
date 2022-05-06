target_sources(nglib PRIVATE
    dynamicmem.cpp
    gzstream.cpp
    hashtabl.cpp
    mystring.cpp
    ngarray.cpp
    ngbitarray.cpp
    optmem.cpp
    parthreads.cpp
    seti.cpp
    sort.cpp
    spbita2d.cpp
    table.cpp
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

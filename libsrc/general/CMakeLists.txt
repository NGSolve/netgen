target_sources(nglib PRIVATE
    gzstream.cpp
    hashtabl.cpp
    mystring.cpp
    optmem.cpp
    parthreads.cpp
    seti.cpp
    sort.cpp
    table.cpp
)

install(FILES
  autodiff.hpp
  hashtabl.hpp myadt.hpp
  mystring.hpp netgenout.hpp ngpython.hpp
  optmem.hpp parthreads.hpp seti.hpp sort.hpp
  table.hpp template.hpp
  gzstream.h
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/general COMPONENT netgen_devel
)

target_sources(nglib PRIVATE
    gzstream.cpp
    optmem.cpp
    seti.cpp
    sort.cpp
    table.cpp
)

install(FILES
  myadt.hpp
  netgenout.hpp ngpython.hpp
  optmem.hpp seti.hpp sort.hpp
  table.hpp template.hpp
  gzstream.h
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/general COMPONENT netgen_devel
)

target_sources(nglib PRIVATE
    gzstream.cpp
    optmem.cpp
    table.cpp
)

install(FILES
  myadt.hpp
  netgenout.hpp ngpython.hpp
  optmem.hpp
  table.hpp template.hpp
  gzstream.h
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/general COMPONENT netgen_devel
)

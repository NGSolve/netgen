target_sources(nglib PRIVATE
    gzstream.cpp
    table.cpp
)

install(FILES
  myadt.hpp
  ngpython.hpp
  table.hpp template.hpp
  gzstream.h
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/general COMPONENT netgen_devel
)

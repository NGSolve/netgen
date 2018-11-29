add_definitions(-DNGINTERFACE_EXPORTS)
add_library(ngcore OBJECT basearchive.cpp)

set_target_properties(ngcore PROPERTIES POSITION_INDEPENDENT_CODE ON )

install(FILES ngcore.hpp archive.hpp basearchive.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/core COMPONENT netgen_devel
)

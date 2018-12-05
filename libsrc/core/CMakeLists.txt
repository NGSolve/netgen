add_definitions(-DNGINTERFACE_EXPORTS)
add_library(ngcore basearchive.cpp)

set_target_properties(ngcore PROPERTIES POSITION_INDEPENDENT_CODE ON )

install(TARGETS ngcore DESTINATION ${NG_INSTALL_DIR} COMPONENT netgen)

install(FILES ngcore.hpp archive.hpp basearchive.hpp version.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/core COMPONENT netgen_devel
)

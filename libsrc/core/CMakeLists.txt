add_library(ngcore basearchive.cpp)

set_target_properties(ngcore PROPERTIES POSITION_INDEPENDENT_CODE ON )
target_compile_definitions(ngcore PRIVATE -DNGCORE_EXPORTS)

install(TARGETS ngcore DESTINATION ${NG_INSTALL_DIR} COMPONENT netgen)

install(FILES ngcore.hpp archive.hpp basearchive.hpp version.hpp type_traits.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/core COMPONENT netgen_devel
)

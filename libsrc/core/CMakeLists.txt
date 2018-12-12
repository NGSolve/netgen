add_library(ngcore SHARED archive.cpp version.cpp)

target_compile_definitions(ngcore PRIVATE -DNGCORE_EXPORTS)

install(TARGETS ngcore DESTINATION ${NG_INSTALL_DIR} COMPONENT netgen)

install(FILES ngcore.hpp archive.hpp type_traits.hpp version.hpp ngcore_api.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/core COMPONENT netgen_devel)

if(ENABLE_CPP_CORE_GUIDELINES_CHECK)
  set_target_properties(ngcore PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif(ENABLE_CPP_CORE_GUIDELINES_CHECK)

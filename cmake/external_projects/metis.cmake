set(METIS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_metis)
set(METIS_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/metis)

ExternalProject_Add(project_metis
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
  URL "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"
  URL_MD5 5465e67079419a69e0116de24fce58fe
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  CMAKE_ARGS
         -DGKLIB_PATH=${METIS_SRC_DIR}/GKlib
         -DCMAKE_INSTALL_PREFIX=${METIS_DIR}
  UPDATE_COMMAND "" # Disable update
  BUILD_IN_SOURCE 1
  )

set_vars( NETGEN_CMAKE_ARGS METIS_DIR )

list(APPEND NETGEN_DEPENDENCIES project_metis)

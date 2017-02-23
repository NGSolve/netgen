set(METIS_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/parmetis/metis)
set(PARMETIS_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/parametis)
ExternalProject_Add(parmetis
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
  URL "http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz"
  URL_MD5 f69c479586bf6bb7aff6a9bc0c739628
  CMAKE_CACHE_ARGS -DMETIS_PATH:string=${METIS_DIR} -DGKLIB_PATH:string=${METIS_DIR}/GKlib
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  UPDATE_COMMAND "" # Disable update
  BUILD_IN_SOURCE 1
  )
set_vars( NETGEN_CMAKE_ARGS METIS_DIR )

list(APPEND NETGEN_DEPENDENCIES parmetis)

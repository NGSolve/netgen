set(METIS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_metis)
set(METIS_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/metis)

ExternalProject_Add(project_metis
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
  URL https://bitbucket.org/petsc/pkg-metis/get/v5.1.0-p12.tar.gz
  URL_MD5 6cd66f75f88dfa2cf043de011f85d8bc
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  CMAKE_ARGS
         -DCMAKE_POLICY_VERSION_MINIMUM=3.5
         -DGKLIB_PATH=${METIS_SRC_DIR}/GKlib
         -DCMAKE_INSTALL_PREFIX=${METIS_DIR}
	 -DCMAKE_POSITION_INDEPENDENT_CODE=ON
	 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
	 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  UPDATE_COMMAND "" # Disable update
  BUILD_IN_SOURCE 1
  )

set_vars( NETGEN_CMAKE_ARGS METIS_DIR )

list(APPEND NETGEN_DEPENDENCIES project_metis)

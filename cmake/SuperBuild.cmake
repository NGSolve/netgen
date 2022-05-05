include (ExternalProject)

option( BUILD_ZLIB "Build and link static version of zlib (useful for pip binaries)" OFF )
option( BUILD_OCC "Build and link static version of occ (useful for pip binaries)" OFF )
set_property (DIRECTORY PROPERTY EP_PREFIX dependencies)

set (NETGEN_DEPENDENCIES)
set (LAPACK_DEPENDENCIES)
set (NETGEN_CMAKE_ARGS "" CACHE INTERNAL "")
set (SUBPROJECT_CMAKE_ARGS "" CACHE INTERNAL "")

set (SUBPROJECT_ARGS
    LIST_SEPARATOR |
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
)

# only show output on failure in ci-builds
if(DEFINED ENV{CI})
    set (SUBPROJECT_ARGS
        LOG_DOWNLOAD ON
        LOG_BUILD ON
        LOG_INSTALL ON
        LOG_CONFIGURE ON
    )
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.14.0")
        set (SUBPROJECT_ARGS
            ${SUBPROJECT_ARGS}
            LOG_OUTPUT_ON_FAILURE ON
            LOG_MERGED_STDOUTERR ON
        )
    endif()
endif()


set (NETGEN_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} )

macro(set_vars VAR_OUT)
  foreach(varname ${ARGN})
    if(NOT "${${varname}}" STREQUAL "")
      string(REPLACE ";" "|" varvalue "${${varname}}" )
      set(${VAR_OUT} "${${VAR_OUT}};-D${varname}=${varvalue}" CACHE INTERNAL "")
    endif()
  endforeach()
endmacro()
#######################################################################

set_vars(SUBPROJECT_CMAKE_ARGS CMAKE_OSX_DEPLOYMENT_TARGET)
set_vars(SUBPROJECT_CMAKE_ARGS CMAKE_OSX_SYSROOT)
set_vars(SUBPROJECT_CMAKE_ARGS CMAKE_C_COMPILER)
set_vars(SUBPROJECT_CMAKE_ARGS CMAKE_CXX_COMPILER)
set_vars(SUBPROJECT_CMAKE_ARGS CMAKE_BUILD_TYPE)

set(SUBPROJECT_CMAKE_ARGS "${SUBPROJECT_CMAKE_ARGS};-DCMAKE_POSITION_INDEPENDENT_CODE=ON" CACHE INTERNAL "")

if(USE_CCACHE)
  find_program(CCACHE_FOUND NAMES ccache ccache.bat)
  if(CCACHE_FOUND)
      set(SUBPROJECT_CMAKE_ARGS "${SUBPROJECT_CMAKE_ARGS};-DCMAKE_CXX_COMPILER_LAUNCHER=${CCACHE_FOUND}" CACHE INTERNAL "")
  endif()
endif()

#######################################################################
set (DEPS_DOWNLOAD_URL "https://github.com/NGSolve/ngsolve_dependencies/releases/download/v1.0.0" CACHE STRING INTERNAL)
set (OCC_DOWNLOAD_URL_WIN "${DEPS_DOWNLOAD_URL}/occ75_win64.zip" CACHE STRING INTERNAL)
set (TCLTK_DOWNLOAD_URL_WIN "${DEPS_DOWNLOAD_URL}/tcltk_win64.zip" CACHE STRING INTERNAL)
set (ZLIB_DOWNLOAD_URL_WIN "${DEPS_DOWNLOAD_URL}/zlib_win64.zip" CACHE STRING INTERNAL)
set (CGNS_DOWNLOAD_URL_WIN "${DEPS_DOWNLOAD_URL}/cgns_win64.zip" CACHE STRING INTERNAL)
set (CGNS_DOWNLOAD_URL_MAC "${DEPS_DOWNLOAD_URL}/cgns_mac.zip" CACHE STRING INTERNAL)


if(UNIX)
  message("Checking for write permissions in install directory...")
  execute_process(COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX})
  execute_process(COMMAND test -w ${CMAKE_INSTALL_PREFIX} RESULT_VARIABLE res)
  if(res)
    message(WARNING "No write access at install directory, please set correct permissions")
  endif()
endif(UNIX)

if(USE_OCC)
if(BUILD_OCC)
  set(OCC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/occ)

  ExternalProject_Add(project_occ
    URL https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V7_6_1.zip
    URL_MD5 e891d85cad61c5cc7ccba3d0110f0c8c
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    ${SUBPROJECT_ARGS}
    CMAKE_ARGS
         -DCMAKE_INSTALL_PREFIX=${OCC_DIR}
         -DCMAKE_PREFIX_PATH=${OCC_DIR}
         -DBUILD_LIBRARY_TYPE:STRING=Static
         -DBUILD_MODULE_FoundationClasses:BOOL=ON
         -DBUILD_MODULE_ModelingData:BOOL=ON
         -DBUILD_MODULE_ModelingAlgorithms:BOOL=ON
         -DBUILD_MODULE_DataExchange:BOOL=ON
         -DBUILD_MODULE_Visualization:BOOL=OFF
         -DBUILD_MODULE_ApplicationFramework:BOOL=OFF
         -DBUILD_MODULE_Draw:BOOL=OFF
         -DUSE_FREETYPE:BOOL=OFF
         -DUSE_OPENGL:BOOL=OFF
         -DUSE_XLIB:BOOL=OFF
         -DBUILD_DOC_Overview:BOOL=OFF
         ${SUBPROJECT_CMAKE_ARGS}
    UPDATE_COMMAND ""
    )

  list(APPEND NETGEN_DEPENDENCIES project_occ)
  set(OpenCascade_ROOT ${OCC_DIR})
else(BUILD_OCC)
    if(WIN32 AND NOT OCC_INCLUDE_DIR AND NOT OpenCASCADE_DIR)
        # we can download prebuilt occ binaries for windows
        ExternalProject_Add(win_download_occ
          ${SUBPROJECT_ARGS}
          URL ${OCC_DOWNLOAD_URL_WIN}
          UPDATE_COMMAND "" # Disable update
          BUILD_IN_SOURCE 1
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}
          )
        list(APPEND NETGEN_DEPENDENCIES win_download_occ)
    else()
        find_package(OpenCascade NAMES OpenCasCade OpenCASCADE opencascade REQUIRED)
    endif()
endif(BUILD_OCC)
endif(USE_OCC)

if(BUILD_ZLIB)
  set(ZLIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/zlib)
  ExternalProject_Add(project_zlib
    ${SUBPROJECT_ARGS}
    URL https://github.com/madler/zlib/archive/refs/tags/v1.2.11.zip
    URL_MD5 9d6a627693163bbbf3f26403a3a0b0b1
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    CMAKE_ARGS
         -DCMAKE_INSTALL_PREFIX=${ZLIB_DIR}
         ${SUBPROJECT_CMAKE_ARGS}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    )

  list(APPEND NETGEN_DEPENDENCIES project_zlib)
  list(APPEND NETGEN_CMAKE_PREFIX_PATH ${ZLIB_DIR})
  if(WIN32)
    # force linking the static library
    set(ZLIB_INCLUDE_DIRS ${ZLIB_DIR}/include)
    set(ZLIB_LIBRARIES ${ZLIB_DIR}/lib/zlibstatic.lib)
  endif(WIN32)
else()
    include(cmake/external_projects/zlib.cmake)
endif()

#######################################################################
if (USE_PYTHON)
  find_path(PYBIND_INCLUDE_DIR pybind11/pybind11.h PATHS ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/pybind11/include NO_DEFAULT_PATH)
    set(NG_INSTALL_PYBIND ON)
    if( NOT PYBIND_INCLUDE_DIR )
      # if the pybind submodule is missing, try to initialize and update all submodules
      execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      find_path(PYBIND_INCLUDE_DIR pybind11/pybind11.h PATHS ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/pybind11/include NO_DEFAULT_PATH)
    endif( NOT PYBIND_INCLUDE_DIR )
    if( PYBIND_INCLUDE_DIR )
        message("-- Found Pybind11: ${PYBIND_INCLUDE_DIR}")
    else( PYBIND_INCLUDE_DIR )
        message(FATAL_ERROR "Could NOT find pybind11!")
    endif( PYBIND_INCLUDE_DIR )
    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)

    set_vars(NETGEN_CMAKE_ARGS
      PYTHON_INCLUDE_DIRS
      PYTHON_LIBRARIES
      PYTHON_EXECUTABLE
      PYTHON_VERSION
      PYBIND_INCLUDE_DIR
      NG_INSTALL_PYBIND
      )
endif (USE_PYTHON)

#######################################################################

if(USE_GUI)
  include(cmake/external_projects/tcltk.cmake)
endif(USE_GUI)

if(USE_CGNS)
  include(cmake/external_projects/cgns.cmake)
endif(USE_CGNS)

#######################################################################
if(USE_MPI)
  if(UNIX)
    if (METIS_DIR)
      message(STATUS "Using external METIS at: ${METIS_DIR}")
    else (METIS_DIR)
      message(STATUS "Looking for system METIS")
      find_package(METIS QUIET)
      if(NOT METIS_FOUND)
	message(WARNING "Could not find METIS, it will be built from source (this might conflict with NGSolve MUMPS)!")
	include(cmake/external_projects/metis.cmake)
      endif(NOT METIS_FOUND)
    endif(METIS_DIR)
  else(UNIX)
    find_package(METIS REQUIRED)
  endif(UNIX)
endif(USE_MPI)


#######################################################################
# propagate cmake variables to Netgen subproject
set_vars( NETGEN_CMAKE_ARGS
  CMAKE_SHARED_LINKER_FLAGS
  CMAKE_SHARED_LINKER_FLAGS_RELEASE
  CMAKE_CXX_FLAGS
  CMAKE_CXX_FLAGS_RELEASE

  USE_GUI
  USE_PYTHON
  USE_MPI
  USE_VT
  USE_VTUNE
  USE_NUMA
  USE_CCACHE
  USE_NATIVE_ARCH
  USE_OCC
  USE_MPEG
  USE_JPEG
  USE_CGNS
  USE_INTERNAL_TCL
  INSTALL_PROFILES
  INTEL_MIC
  CMAKE_INSTALL_PREFIX
  ENABLE_UNIT_TESTS
  ENABLE_CPP_CORE_GUIDELINES_CHECK
  USE_SPDLOG
  DEBUG_LOG
  CHECK_RANGE
  TRACE_MEMORY
  BUILD_STUB_FILES
  BUILD_FOR_CONDA
  NG_COMPILE_FLAGS
  OpenCascade_ROOT
  ZLIB_INCLUDE_DIRS
  ZLIB_LIBRARIES

  NGLIB_LIBRARY_TYPE
  NGCORE_LIBRARY_TYPE
  NGGUI_LIBRARY_TYPE
  )

# propagate all variables set on the command line using cmake -DFOO=BAR
# to Netgen subproject
get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(CACHE_VAR ${CACHE_VARS})
  get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
  if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line." AND NOT CACHE_VAR STREQUAL "CMAKE_OSX_ARCHITECTURES")
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    string(REPLACE ";" "|" varvalue "${${CACHE_VAR}}" )
    set(NETGEN_CMAKE_ARGS ${NETGEN_CMAKE_ARGS};-D${CACHE_VAR}:${CACHE_VAR_TYPE}=${varvalue} CACHE INTERNAL "")
  endif()
endforeach()

if(${CMAKE_GENERATOR} STREQUAL "Unix Makefiles")
  set(NETGEN_BUILD_COMMAND $(MAKE) --silent )
else()
  set(NETGEN_BUILD_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/netgen --config ${CMAKE_BUILD_TYPE})
endif()


string(REPLACE ";" "|" NETGEN_CMAKE_PREFIX_PATH_ALT_SEP "${NETGEN_CMAKE_PREFIX_PATH}")
ExternalProject_Add (netgen
  ${SUBPROJECT_ARGS}
  DEPENDS ${NETGEN_DEPENDENCIES}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  CMAKE_ARGS
      -DUSE_SUPERBUILD=OFF
      ${NETGEN_CMAKE_ARGS}
      ${SUBPROJECT_CMAKE_ARGS}
      -DCMAKE_PREFIX_PATH=${NETGEN_CMAKE_PREFIX_PATH_ALT_SEP}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/netgen
  BUILD_COMMAND ${NETGEN_BUILD_COMMAND}
  STEP_TARGETS build
)

# Check if the git submodules (i.e. pybind11) are up to date
# in case, something is wrong, emit a warning but continue
 ExternalProject_Add_Step(netgen check_submodules
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_submodules.cmake
   DEPENDERS install # Steps on which this step depends
   )

# Due to 'ALWAYS 1', this step is always run which also forces a build of
# the Netgen subproject
 ExternalProject_Add_Step(netgen check_submodules1
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_submodules.cmake
   DEPENDEES configure # Steps on which this step depends
   DEPENDERS build     # Steps that depend on this step
   ALWAYS 1            # No stamp file, step always runs
   )


install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build . --target install --config ${CMAKE_BUILD_TYPE} WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}/netgen\")")

add_custom_target(test_netgen
  ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/netgen
                   --target test
                   --config ${CMAKE_BUILD_TYPE}
                   )

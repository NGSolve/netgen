include (ExternalProject)

set_property (DIRECTORY PROPERTY EP_PREFIX dependencies)

set (NETGEN_DEPENDENCIES)
set (LAPACK_DEPENDENCIES)
set (NETGEN_CMAKE_ARGS "" CACHE INTERNAL "")

macro(set_vars VAR_OUT)
  foreach(varname ${ARGN})
    if(NOT "${${varname}}" STREQUAL "")
      string(REPLACE ";" "$<SEMICOLON>" varvalue "${${varname}}" )
      set(${VAR_OUT} ${${VAR_OUT}};-D${varname}=${varvalue} CACHE INTERNAL "")
    endif()
  endforeach()
endmacro()
#######################################################################
if(WIN32)
  set (DEPS_DOWNLOAD_URL "https://github.com/NGSolve/ngsolve_dependencies/releases/download/v1.0.0" CACHE STRING INTERNAL)
  set (OCC_DOWNLOAD_URL_WIN "${DEPS_DOWNLOAD_URL}/occ_win64.zip" CACHE STRING INTERNAL)
  set (TCLTK_DOWNLOAD_URL_WIN "${DEPS_DOWNLOAD_URL}/tcltk_win64.zip" CACHE STRING INTERNAL)
  set (ZLIB_DOWNLOAD_URL_WIN "${DEPS_DOWNLOAD_URL}/zlib_win64.zip" CACHE STRING INTERNAL)
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_CXX_FLAGS_NEW ${CMAKE_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_NEW} CACHE STRING "compile flags" FORCE)
    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_CXX_FLAGS_NEW ${CMAKE_CXX_FLAGS_RELEASE})
    set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_NEW} CACHE STRING "compile flags" FORCE)

    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_SHARED_LINKER_FLAGS_NEW ${CMAKE_SHARED_LINKER_FLAGS})
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_NEW} /IGNORE:4217,4049" CACHE STRING "compile flags" FORCE)
    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_EXE_LINKER_FLAGS_NEW ${CMAKE_EXE_LINKER_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS_NEW}/IGNORE:4217,4049" CACHE STRING "compile flags" FORCE)

  endif(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
endif(WIN32)

if(UNIX)
  message("Checking for write permissions in install directory...")
  execute_process(COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX})
  execute_process(COMMAND test -w ${CMAKE_INSTALL_PREFIX} RESULT_VARIABLE res)
  if(res)
    message(WARNING "No write access at install directory, please set correct permissions")
  endif()
endif(UNIX)

if(NOT WIN32)
    find_package(ZLIB REQUIRED)
    set_vars(NETGEN_CMAKE_ARGS ZLIB_INCLUDE_DIRS ZLIB_LIBRARIES)
endif(NOT WIN32)

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

if(USE_OCC AND WIN32 AND NOT OCC_INCLUDE_DIR)
    ExternalProject_Add(win_download_occ
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tcl
      URL ${OCC_DOWNLOAD_URL_WIN}
      UPDATE_COMMAND "" # Disable update
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}
      LOG_DOWNLOAD 1
      )
    list(APPEND NETGEN_DEPENDENCIES win_download_occ)
endif(USE_OCC AND WIN32 AND NOT OCC_INCLUDE_DIR)

#######################################################################

include(cmake/external_projects/zlib.cmake)
if(USE_GUI)
  include(cmake/external_projects/tcltk.cmake)
endif(USE_GUI)

#######################################################################
if(USE_MPI)
  if(UNIX)
    find_package(METIS QUIET)
    if(NOT METIS_FOUND)
      message(STATUS "Could not find METIS, it will be built from source")
      include(cmake/external_projects/metis.cmake)
    endif()
  else(UNIX)
    find_package(METIS REQUIRED)
  endif(UNIX)
endif(USE_MPI)


#######################################################################
# propagate cmake variables to Netgen subproject
set_vars( NETGEN_CMAKE_ARGS
  CMAKE_CXX_COMPILER
  CMAKE_BUILD_TYPE
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
  USE_INTERNAL_TCL
  INSTALL_PROFILES
  INTEL_MIC
  CMAKE_PREFIX_PATH
  CMAKE_INSTALL_PREFIX
  )

# propagate all variables set on the command line using cmake -DFOO=BAR
# to Netgen subproject
get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(CACHE_VAR ${CACHE_VARS})
  get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
  if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    set(NETGEN_CMAKE_ARGS ${NETGEN_CMAKE_ARGS};-D${CACHE_VAR}:${CACHE_VAR_TYPE}=${${CACHE_VAR}} CACHE INTERNAL "")
  endif()
endforeach()

if(${CMAKE_GENERATOR} STREQUAL "Unix Makefiles")
  set(NETGEN_BUILD_COMMAND $(MAKE) --silent )
else()
  set(NETGEN_BUILD_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/netgen --config ${CMAKE_BUILD_TYPE})
endif()

ExternalProject_Add (netgen
  DEPENDS ${NETGEN_DEPENDENCIES}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  CMAKE_ARGS -DUSE_SUPERBUILD=OFF ${NETGEN_CMAKE_ARGS}
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


install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build . --target install --config ${CMAKE_BUILD_TYPE} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netgen)")

add_custom_target(test_netgen
  ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/netgen
                   --target test
                   --config ${CMAKE_BUILD_TYPE}
                   )

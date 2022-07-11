if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING INTERNAL)
endif(NOT CMAKE_BUILD_TYPE)

cmake_minimum_required(VERSION 3.13)
cmake_policy(VERSION 3.13)

option( USE_NATIVE_ARCH  "build for native cpu architecture" ON)

option( USE_GUI     "build with GUI" ON )
option( USE_PYTHON  "build with python interface" ON )
option( USE_MPI     "enable mpi parallelization" OFF )
option( USE_MPI4PY  "enable mpi4py interface" ON )
option( USE_OCC     "build with OpenCascade geometry kernel interface" OFF)
option( USE_STLGEOM     "build with STL geometry support" ON)
option( USE_CSG     "build with CSG kernel" ON)
option( USE_INTERFACE     "build nginterface" ON)
option( USE_GEOM2D     "build 2d geometry kernels" ON)
option( USE_JPEG    "enable snapshots using library libjpeg" OFF )
option( USE_MPEG    "enable video recording with FFmpeg, uses libavcodec" OFF )
option( USE_CGNS    "enable CGNS file read/write support" OFF )
option( USE_NUMA         "compile with NUMA-aware code")
option( INTEL_MIC        "cross compile for intel xeon phi")
option( INSTALL_PROFILES "install environment variable settings to /etc/profile.d" OFF )
option( USE_CCACHE       "use ccache")
option( USE_INTERNAL_TCL "Compile tcl files into the code and don't install them" ON)
option( ENABLE_UNIT_TESTS "Enable Catch unit tests")
option( ENABLE_CPP_CORE_GUIDELINES_CHECK "Enable cpp core guideline checks on ngcore" OFF)
option( USE_SPDLOG "Enable spd log logging" OFF)
option( DEBUG_LOG "Enable more debug output (may increase computation time) - only works with USE_SPDLOG=ON" OFF)
option( CHECK_RANGE "Check array range access, automatically enabled if built in debug mode" OFF)
option( BUILD_STUB_FILES "Build stub files for better autocompletion" ON)
option( BUILD_FOR_CONDA "Link python libraries only to executables" OFF)

option( USE_SUPERBUILD       "build dependencies automatically" ON)
option( TRACE_MEMORY     "Enable memory tracing" OFF)

set(NG_COMPILE_FLAGS "" CACHE STRING "Additional compile flags")

set(NGLIB_LIBRARY_TYPE SHARED CACHE STRING "nglib library type")
set(NGCORE_LIBRARY_TYPE SHARED CACHE STRING "ngcore library type")
set(NGGUI_LIBRARY_TYPE SHARED CACHE STRING "nggui library type")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_modules")

if(APPLE)
  set(INSTALL_DIR_DEFAULT /Applications/Netgen.app)
else(APPLE)
  if(WIN32)
    set(INSTALL_DIR_DEFAULT "C:/netgen")
  else(WIN32)
    set(INSTALL_DIR_DEFAULT /opt/netgen)
  endif(WIN32)
endif(APPLE)

if(INSTALL_DIR)
  message(WARNING "INSTALL_DIR is deprecated, use CMAKE_INSTALL_PREFIX instead")
  set(INSTALL_DIR_DEFAULT ${INSTALL_DIR})
endif(INSTALL_DIR)

if(UNIX AND USE_SUPERBUILD)
  message("Checking for write permissions in install directory...")
  execute_process(COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX})
  execute_process(COMMAND test -w ${CMAKE_INSTALL_PREFIX} RESULT_VARIABLE res)
  if(res)
    message(WARNING "No write access at install directory, please set correct permissions")
  endif()
endif(UNIX AND USE_SUPERBUILD)

if (USE_SUPERBUILD)
  project (SUPERBUILD)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR_DEFAULT}" CACHE PATH "Install directory" FORCE)
  endif()

  # execute the superbuild (this script will be invoked again without the
  # USE_SUPERBUILD option this time)
  include (cmake/SuperBuild.cmake)
  return()                      # stop processing this file further
else()
  project(Netgen)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR_DEFAULT}" CACHE PATH "Install directory" FORCE)
  endif()
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include (${CMAKE_CURRENT_LIST_DIR}/cmake/generate_version_file.cmake)
set(CPACK_PACKAGE_VERSION "${NETGEN_VERSION}")


#######################################################################
if(USE_CCACHE)
  find_program(CCACHE_FOUND NAMES ccache ccache.bat)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_FOUND})
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

#######################################################################
if(INTEL_MIC)
    set(MKL_ARCH "mic")
    include(cmake/mic.cmake)
else(INTEL_MIC)
    set(MKL_ARCH "intel64")
endif(INTEL_MIC)

#######################################################################
# Append install paths of software in non-standard paths (e.g. openmpi, metis, intel mkl, ...)
# cmake -DUSE_MPI=ON -DCMAKE_PREFIX_PATH="/opt/openmpi165;/opt/metis51" ../
set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, separated by ';'")
if (ADDITIONAL_PATHS)
  set(CMAKE_PREFIX_PATH ${ADDITIONAL_PATHS})
endif (ADDITIONAL_PATHS)

#######################################################################
# build options
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_SOURCE_DIR}/libsrc")
include_directories ("${PROJECT_SOURCE_DIR}/libsrc/include")
include_directories ("${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(USE_PYTHON)
    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)

    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE PYTHON_PACKAGES_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    file(TO_CMAKE_PATH ${PYTHON_PACKAGES_INSTALL_DIR} PYTHON_PACKAGES_INSTALL_DIR)
endif(USE_PYTHON)

set(NG_INSTALL_SUFFIX netgen CACHE STRING "Suffix appended to install directories (project name)")

if(APPLE)
  set(NG_INSTALL_DIR_BIN_DEFAULT Contents/MacOS)
  set(NG_INSTALL_DIR_LIB_DEFAULT Contents/MacOS)
  set(NG_INSTALL_DIR_CMAKE_DEFAULT Contents/Resources/CMake)
  set(NG_INSTALL_DIR_PYTHON_DEFAULT Contents/Resources/${PYTHON_PACKAGES_INSTALL_DIR})
  set(NG_INSTALL_DIR_RES_DEFAULT Contents/Resources/share)
  set(NG_INSTALL_DIR_INCLUDE_DEFAULT Contents/Resources/include)

  set(NG_RPATH_TOKEN "@loader_path")
else(APPLE)
  set(NG_INSTALL_DIR_BIN_DEFAULT bin)
  set(NG_INSTALL_DIR_LIB_DEFAULT lib)
  if(WIN32)
    set(NG_INSTALL_DIR_CMAKE_DEFAULT cmake)
  else(WIN32)
    set(NG_INSTALL_DIR_CMAKE_DEFAULT lib/cmake/${NG_INSTALL_SUFFIX})
  endif(WIN32)
  set(NG_INSTALL_DIR_PYTHON_DEFAULT ${PYTHON_PACKAGES_INSTALL_DIR})
  set(NG_INSTALL_DIR_RES_DEFAULT share)
  set(NG_INSTALL_DIR_INCLUDE_DEFAULT include)

  set(NG_RPATH_TOKEN "\$ORIGIN")
endif(APPLE)

set(NG_INSTALL_DIR_PYTHON ${NG_INSTALL_DIR_PYTHON_DEFAULT} CACHE STRING "Install directory for Python files")
set(NG_INSTALL_DIR_BIN ${NG_INSTALL_DIR_BIN_DEFAULT} CACHE STRING "Install directory for executables")
set(NG_INSTALL_DIR_LIB ${NG_INSTALL_DIR_LIB_DEFAULT} CACHE STRING "Install directory for libraries")
set(NG_INSTALL_DIR_INCLUDE ${NG_INSTALL_DIR_INCLUDE_DEFAULT} CACHE STRING "Install directory for header files")
set(NG_INSTALL_DIR_CMAKE ${NG_INSTALL_DIR_CMAKE_DEFAULT} CACHE STRING "Install directory for CMake files")
set(NG_INSTALL_DIR_RES ${NG_INSTALL_DIR_RES_DEFAULT} CACHE STRING "Install directory for resources")

get_filename_component(NETGEN_CMAKE_DIR_ABSOLUTE ${NG_INSTALL_DIR_CMAKE} ABSOLUTE BASE_DIR ${CMAKE_INSTALL_PREFIX})
get_filename_component(NETGEN_BINARY_DIR_ABSOLUTE ${NG_INSTALL_DIR_BIN} ABSOLUTE BASE_DIR ${CMAKE_INSTALL_PREFIX})
get_filename_component(NETGEN_LIBRARY_DIR_ABSOLUTE ${NG_INSTALL_DIR_LIB} ABSOLUTE BASE_DIR ${CMAKE_INSTALL_PREFIX})
get_filename_component(NETGEN_INCLUDE_DIR_ABSOLUTE ${NG_INSTALL_DIR_INCLUDE} ABSOLUTE BASE_DIR ${CMAKE_INSTALL_PREFIX})
get_filename_component(NETGEN_RESOURCE_DIR_ABSOLUTE ${NG_INSTALL_DIR_RES} ABSOLUTE BASE_DIR ${CMAKE_INSTALL_PREFIX})

file(RELATIVE_PATH NETGEN_DIR ${NETGEN_CMAKE_DIR_ABSOLUTE} ${CMAKE_INSTALL_PREFIX})
file(RELATIVE_PATH NETGEN_BINARY_DIR ${NETGEN_CMAKE_DIR_ABSOLUTE} ${NETGEN_BINARY_DIR_ABSOLUTE})
file(RELATIVE_PATH NETGEN_LIBRARY_DIR ${NETGEN_CMAKE_DIR_ABSOLUTE} ${NETGEN_LIBRARY_DIR_ABSOLUTE})
file(RELATIVE_PATH NETGEN_INCLUDE_DIR ${NETGEN_CMAKE_DIR_ABSOLUTE} ${NETGEN_INCLUDE_DIR_ABSOLUTE})
file(RELATIVE_PATH NETGEN_RESOURCE_DIR ${NETGEN_CMAKE_DIR_ABSOLUTE} ${NETGEN_RESOURCE_DIR_ABSOLUTE})

file(RELATIVE_PATH NETGEN_RPATH ${NETGEN_BINARY_DIR_ABSOLUTE} ${NETGEN_LIBRARY_DIR_ABSOLUTE})

if(USE_PYTHON)
    get_filename_component(NETGEN_PYTHON_DIR_ABSOLUTE ${NG_INSTALL_DIR_PYTHON} ABSOLUTE BASE_DIR ${CMAKE_INSTALL_PREFIX})
    file(RELATIVE_PATH NETGEN_PYTHON_DIR ${NETGEN_CMAKE_DIR_ABSOLUTE} ${NETGEN_PYTHON_DIR_ABSOLUTE})
    file(RELATIVE_PATH NETGEN_PYTHON_RPATH_BIN ${NETGEN_PYTHON_DIR_ABSOLUTE} ${NETGEN_BINARY_DIR_ABSOLUTE})
    file(RELATIVE_PATH NETGEN_PYTHON_RPATH ${NETGEN_PYTHON_DIR_ABSOLUTE} ${NETGEN_LIBRARY_DIR_ABSOLUTE})
    if(WIN32)
        set(NETGEN_PYTHON_RPATH ${NETGEN_PYTHON_RPATH_BIN})
    endif(WIN32)
endif(USE_PYTHON)

set(NG_INSTALL_DIR EXPORT netgen-targets RUNTIME DESTINATION ${NG_INSTALL_DIR_BIN} COMPONENT netgen LIBRARY DESTINATION ${NG_INSTALL_DIR_LIB} COMPONENT netgen_devel ARCHIVE DESTINATION ${NG_INSTALL_DIR_LIB} COMPONENT netgen_devel)
install(EXPORT netgen-targets DESTINATION ${NG_INSTALL_DIR_CMAKE} COMPONENT netgen_devel)

set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${NG_RPATH_TOKEN};${NG_RPATH_TOKEN}/${NETGEN_RPATH}")
if(BUILD_FOR_CONDA)
    file(RELATIVE_PATH py_rpath "/bin" "/${NG_INSTALL_DIR_LIB}")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${py_rpath}")
endif(BUILD_FOR_CONDA)

include (CheckIncludeFiles)
check_include_files (dlfcn.h HAVE_DLFCN_H)
if(HAVE_DLFCN_H)
  add_definitions(-DHAVE_DLFCN_H)
endif()

include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

#######################################################################
# platform specific definitions and flags
macro(get_WIN32_WINNT version)
    if (WIN32 AND CMAKE_SYSTEM_VERSION)
        if("${CMAKE_SYSTEM_VERSION}" MATCHES "^([0-9]+)\\.([0-9]+)")
            math(EXPR ver "${CMAKE_MATCH_1}*100 + ${CMAKE_MATCH_2}")
        endif()
        set(${version} "0x${ver}")
    endif()
endmacro()

macro(get_dll_from_lib dll_path lib_path)
        get_filename_component(parent_lib_path ${lib} DIRECTORY)
        get_filename_component(lib_name ${lib} name)
endmacro()

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_STANDARD 17)
if(WIN32)
    set(CMAKE_MFC_FLAG 0)
endif(WIN32)
if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif(APPLE)

#######################################################################
add_library(nglib ${NGLIB_LIBRARY_TYPE})
if(USE_GUI)
    add_library(nggui ${NGGUI_LIBRARY_TYPE})
    if(WIN32)
        set_target_properties( nggui PROPERTIES OUTPUT_NAME "libnggui")
    endif(WIN32)
endif(USE_GUI)

#######################################################################
if(NOT ZLIB_INCLUDE_DIRS)
  find_package(ZLIB REQUIRED)
endif(NOT ZLIB_INCLUDE_DIRS)
target_include_directories(nglib PRIVATE ${ZLIB_INCLUDE_DIRS})
if(USE_GUI)
    target_include_directories(nggui PRIVATE ${ZLIB_INCLUDE_DIRS})
endif(USE_GUI)
target_link_libraries(nglib PRIVATE ${ZLIB_LIBRARIES})

#######################################################################
if(WIN32)
    add_library(netgen_gui INTERFACE IMPORTED)
else()
    add_library(netgen_gui INTERFACE)
endif()

if (USE_GUI)
    find_package(TCL 8.5 REQUIRED)
    find_package(TclStub 8.5 REQUIRED)
    find_package(Threads REQUIRED)
    if(APPLE)
          set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework AppKit")
          target_link_libraries(netgen_gui INTERFACE "-framework AppKit")
    else(APPLE)
          find_package(X11 REQUIRED) 
          target_link_libraries( netgen_gui INTERFACE ${X11_Xmu_LIB} ${X11_X11_LIB})
    endif(APPLE)
    find_package(OpenGL REQUIRED)

    target_compile_definitions(netgen_gui INTERFACE -DTCL -DOPENGL -DUSE_TOGL_2 -DUSE_TCL_STUBS -DUSE_TK_STUBS)
    target_include_directories(netgen_gui INTERFACE ${TCL_INCLUDE_PATH} ${TK_INCLUDE_PATH})
    target_link_libraries(netgen_gui INTERFACE ${TCL_STUB_LIBRARY} ${TK_STUB_LIBRARY})
    if(NOT EXISTS ${TK_INCLUDE_PATH}/tkWin.h AND EXISTS ${TK_INCLUDE_PATH}/../win/tkWin.h)
       target_include_directories(netgen_gui INTERFACE ${TK_INCLUDE_PATH}/../win)
    endif()
    if(NOT EXISTS ${TK_INCLUDE_PATH}/x11/Xlib.h AND EXISTS ${TK_INCLUDE_PATH}/../xlib/X11/Xlib.h)
       target_include_directories(netgen_gui INTERFACE ${TK_INCLUDE_PATH}/../xlib)
    endif()

    target_link_libraries(nggui PUBLIC nglib togl PRIVATE "$<BUILD_INTERFACE:netgen_python>" )

    if(WIN32)
      target_compile_definitions(netgen_gui INTERFACE -DTOGL_WGL)
    endif()
    if(APPLE)
      target_compile_definitions(netgen_gui INTERFACE -DTOGL_NSOPENGL)
    endif()
    if(UNIX AND NOT APPLE)
      target_compile_definitions(netgen_gui INTERFACE -DTOGL_X11)
    endif()

endif (USE_GUI)

#######################################################################
if(WIN32)
    add_library(netgen_python INTERFACE IMPORTED)
else()
    add_library(netgen_python INTERFACE)
endif()

if (USE_PYTHON)
    add_subdirectory(external_dependencies/pybind11)
    find_path(PYBIND_INCLUDE_DIR pybind11/pybind11.h HINTS ${PYTHON_INCLUDE_DIR})
    if( PYBIND_INCLUDE_DIR )
        message(STATUS "Found Pybind11: ${PYBIND_INCLUDE_DIR}")
    else( PYBIND_INCLUDE_DIR )
        message(FATAL_ERROR "Could NOT find pybind11!")
    endif( PYBIND_INCLUDE_DIR )

    target_include_directories(netgen_python INTERFACE ${PYBIND_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
    target_include_directories(nglib PRIVATE ${PYBIND_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
    if(NOT ${BUILD_FOR_CONDA} OR WIN32)
        # Don't link python libraries in conda environments
        target_link_libraries(netgen_python INTERFACE ${PYTHON_LIBRARIES})
    endif()

    if(NG_INSTALL_PYBIND)
      install(DIRECTORY ${PYBIND_INCLUDE_DIR}/pybind11 DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)
      install(FILES ${PYBIND_INCLUDE_DIR}/../LICENSE DESTINATION ${NG_INSTALL_DIR_INCLUDE}/pybind11 COMPONENT netgen_devel)
    endif(NG_INSTALL_PYBIND)
endif (USE_PYTHON)

#######################################################################
add_library(netgen_mpi INTERFACE)
add_library(netgen_metis INTERFACE)
if (USE_MPI)
    find_package(MPI REQUIRED)
    target_include_directories(netgen_mpi INTERFACE ${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(netgen_mpi INTERFACE ${MPI_mpi_LIBRARY} ${MPI_CXX_LIBRARIES} )
    target_compile_definitions(netgen_mpi INTERFACE PARALLEL )

    find_package(METIS REQUIRED)
    target_include_directories(netgen_metis INTERFACE ${METIS_INCLUDE_DIR})
    target_link_libraries(netgen_metis INTERFACE ${METIS_LIBRARY} )
    target_compile_definitions(netgen_metis INTERFACE METIS )

    if(USE_MPI4PY AND USE_PYTHON)
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import mpi4py;print(mpi4py.get_include())" OUTPUT_VARIABLE mpi4py_path OUTPUT_STRIP_TRAILING_WHITESPACE)
      find_path(MPI4PY_INCLUDE_DIR mpi4py.h HINTS ${mpi4py_path}/mpi4py NO_DEFAULT_PATH REQUIRED)
      target_include_directories(netgen_metis INTERFACE ${MPI4PY_INCLUDE_DIR})
      target_compile_definitions(netgen_metis INTERFACE NG_MPI4PY )
      message(STATUS "Found mpi4py: ${MPI4PY_INCLUDE_DIR}")
    endif(USE_MPI4PY AND USE_PYTHON)
endif (USE_MPI)
install(TARGETS netgen_mpi netgen_metis ${NG_INSTALL_DIR})

#######################################################################
if (USE_OCC)
    find_package(OpenCascade NAMES OpenCASCADE opencascade REQUIRED CMAKE_FIND_ROOT_PATH_BOTH)
    add_definitions(-DOCCGEOMETRY)
    set(OCC_LIBRARIES
      TKBO
      TKBRep
      TKBool
      TKCAF
      TKCDF
      TKFillet
      TKG2d
      TKG3d
      TKGeomAlgo
      TKGeomBase
      TKHLR
      TKIGES
      TKLCAF
      TKMath
      TKMesh
      TKOffset
      TKPrim
      TKSTEP
      TKSTEP209
      TKSTEPAttr
      TKSTEPBase
      TKSTL
      TKService
      TKShHealing
      TKTopAlgo
      TKV3d
      TKVCAF
      TKXCAF
      TKXDEIGES
      TKXDESTEP
      TKXSBase
      TKernel
    )
    include_directories(${OpenCASCADE_INCLUDE_DIR})
    if(NOT OpenCASCADE_BUILD_SHARED_LIBS)
      if(OpenCASCADE_WITH_FREETYPE)
        find_library( FREETYPE NAMES freetype HINTS ${OpenCASCADE_INSTALL_PREFIX}/lib)
        list(APPEND OCC_LIBRARIES ${FREETYPE})
        if(UNIX AND NOT APPLE)
          find_package(Fontconfig REQUIRED)
          list(APPEND OCC_LIBRARIES ${Fontconfig_LIBRARIES})
        endif()
      endif(OpenCASCADE_WITH_FREETYPE)
      if(UNIX AND NOT APPLE)
        set(THREADS_PREFER_PTHREAD_FLAG ON)
        find_package(Threads REQUIRED)
        list(APPEND OCC_LIBRARIES Threads::Threads)
        list(PREPEND OCC_LIBRARIES -Wl,--start-group)
        list(APPEND OCC_LIBRARIES -Wl,--end-group)
      endif()
      if(WIN32)
          list(APPEND OCC_LIBRARIES Ws2_32.lib)
      endif()
    endif()
    message(STATUS "OCC DIRS ${OpenCASCADE_INCLUDE_DIR}")
    if(WIN32 AND USE_GUI)
      target_link_libraries(nggui PRIVATE ${OCC_LIBRARIES})
    endif(WIN32 AND USE_GUI)
endif (USE_OCC)

#######################################################################
if (USE_JPEG)
    find_package(JPEG REQUIRED)
    add_definitions(-DJPEGLIB)
    include_directories(${JPEG_INCLUDE_DIR})
endif (USE_JPEG)

#######################################################################
if (USE_MPEG)
    find_package(FFMPEG REQUIRED)
    add_definitions(-DFFMPEG -D__STDC_CONSTANT_MACROS)
    include_directories(${FFMPEG_INCLUDE_DIR})
endif (USE_MPEG)

#######################################################################
add_custom_target(ng_generate_version_file
  ${CMAKE_COMMAND}
    -DBDIR=${CMAKE_CURRENT_BINARY_DIR}
    -P ${CMAKE_CURRENT_LIST_DIR}/cmake/generate_version_file.cmake
    )
#######################################################################
if(INSTALL_PROFILES)
    file(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "#!/bin/sh\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export PATH=${CMAKE_INSTALL_PREFIX}/bin:$PATH\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export NETGENDIR=${CMAKE_INSTALL_PREFIX}/bin\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${PYTHON_PACKAGES_INSTALL_DIR}:.\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:.\n")
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh DESTINATION /etc/profile.d COMPONENT netgen)

    string(ASCII 27 Esc)
    file(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/postinst "#!/bin/sh\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/postinst "echo \"${Esc}[0;31mPlease log out and in again or do 'source /etc/profile.d/netgen.sh' to load the correct environment variables!${Esc}[m\"")
    set( CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_BINARY_DIR}/postinst")
endif(INSTALL_PROFILES)

#######################################################################
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/fixup.cmake "\
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/netgen_fixup.cmake)
set(APP ${CMAKE_INSTALL_PREFIX})
message(\${APP})
set(BU_CHMOD_BUNDLE_ITEMS ON)

file(GLOB libs ${CMAKE_INSTALL_PREFIX}/${NG_INSTALL_DIR_LIB}/*.dylib ${CMAKE_INSTALL_PREFIX}/${NG_INSTALL_DIR_LIB}/*.so)
message(\"\${libs}\")
netgen_fixup_bundle( \${APP}/Contents/MacOS/netgen \"\${libs}\" ${CMAKE_INSTALL_PREFIX}/${NG_INSTALL_DIR_LIB} )
execute_process(COMMAND ln -s /Applications ${CMAKE_INSTALL_PREFIX}/../Applications)
set (bundle_filename \$ENV{NETGEN_BUNDLE_NAME})
if(NOT bundle_filename)
  set(bundle_filename netgen)
endif(NOT bundle_filename)
execute_process(COMMAND hdiutil create -volname Netgen -srcfolder ${CMAKE_INSTALL_PREFIX} -ov -format UDZO \${bundle_filename}-${PACKAGE_VERSION}.dmg)
")
  add_custom_target(bundle COMMAND ${CMAKE_COMMAND} "-P" "${CMAKE_CURRENT_BINARY_DIR}/fixup.cmake")
#######################################################################
# CTest
enable_testing()
include(CTest)

if(ENABLE_UNIT_TESTS)
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/external_projects/catch.cmake)
endif(ENABLE_UNIT_TESTS)


#######################################################################

if(USE_SPDLOG)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/external_projects/spdlog.cmake)
  include_directories(${SPDLOG_INCLUDE_DIR})
endif(USE_SPDLOG)

if(ENABLE_CPP_CORE_GUIDELINES_CHECK)
  find_program(
    CLANG_TIDY_EXE
    NAMES "clang-tidy"
    DOC "Path to clang-tidy executable"
    )
  if(NOT CLANG_TIDY_EXE)
    message(WARNING "clang-tidy not found.")
  else()
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-header-filter=libsrc/core/")
  endif()
endif(ENABLE_CPP_CORE_GUIDELINES_CHECK)

add_library(netgen_cgns INTERFACE)
if(USE_CGNS)
    find_library( CGNS_LIBRARY NAMES cgns cgnsdll )
    find_path( CGNS_INCLUDE_DIR cgnslib.h )
    target_compile_definitions(netgen_cgns INTERFACE NG_CGNS)
    target_include_directories(netgen_cgns INTERFACE ${CGNS_INCLUDE_DIR})
    target_link_libraries(netgen_cgns INTERFACE ${CGNS_LIBRARY})
    if(NOT WIN32 AND NOT APPLE) # hdf5 is statically linked into cgns in Windows amd MacOS binaries
        find_library(HDF5_LIBRARY NAMES hdf5 hdf5_serial)
        target_link_libraries(netgen_cgns INTERFACE ${HDF5_LIBRARY})
    endif(NOT WIN32 AND NOT APPLE)
endif(USE_CGNS)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/netgen_version.hpp ${CMAKE_CURRENT_BINARY_DIR}/netgen_config.hpp DESTINATION ${NG_INSTALL_DIR_INCLUDE}/include COMPONENT netgen_devel)

# include instead of add_subdirectory to recognize the generated source files properly
include(rules/CMakeLists.txt)

add_subdirectory(windows)
add_subdirectory(libsrc)
add_subdirectory(ng)
add_subdirectory(tutorials)
add_subdirectory(py_tutorials)
add_subdirectory(doc)
add_subdirectory(nglib)
add_subdirectory(tests)

#######################################################################
if(USE_NATIVE_ARCH)
    if(WIN32)
        include(CheckCXXSourceRuns)
        check_cxx_source_runs("
                #include <immintrin.h>
                int main()
                {
                    __m256d a{1.,2.,3.,4.};
                    __m256d b{2.,0.,3.,5.};
                    __m256d c = _mm256_mul_pd(a,b);
                    return 0;
                } " NG_HAVE_AVX)
        check_cxx_source_runs("
                #include <immintrin.h>
                int main()
                {
                    __m256i a{1,2,3,4};
                    __m256i b{2,0,3,5};
                    __m256i c = _mm256_cmpgt_epi64 (a,b);
                    return 0;
                } " NG_HAVE_AVX2)
        check_cxx_source_runs("
                #include <immintrin.h>
                int main()
                {
                    __m512d a{1.,2.,3.,4.};
                    __m512d b{5.,6.,7.,8.};
                    __m512d c = _mm512_mul_pd(a,b);
                    return 0;
                } " NG_HAVE_AVX512)

        if(NG_HAVE_AVX512)
            target_compile_options(ngcore PUBLIC "/arch:AVX512")
            message(STATUS "Build for AVX512 CPU")
        elseif(NG_HAVE_AVX2)
            target_compile_options(ngcore PUBLIC "/arch:AVX2")
            message(STATUS "Build for AVX2 CPU")
        elseif(NG_HAVE_AVX)
            target_compile_options(ngcore PUBLIC "/arch:AVX")
            message(STATUS "Build for AVX CPU")
        else()
            message(STATUS "Build for generic CPU")
        endif()
    elseif(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
	# no flag necessary/available on Apple M1
    else()
        target_compile_options(ngcore PUBLIC "-march=native")
    endif(WIN32)
    if(APPLE)
        # work-around for bug in Xcode 11.3: https://forums.developer.apple.com/thread/121887
        target_compile_options(ngcore PUBLIC "-fno-stack-check")
    endif(APPLE)
endif(USE_NATIVE_ARCH)

if (USE_PYTHON)
  add_subdirectory(python)
endif (USE_PYTHON)

#######################################################################
# Debian packager

if(UNIX)
    set(CPACK_SOURCE_GENERATOR "TGZ")
    set(CPACK_SOURCE_IGNORE_FILES "/cmake/;/build/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")
    set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION} )

    set(CPACK_PACKAGE_NAME netgen)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "automatic 3d tetrahedral mesh generator")
    set(CPACK_PACKAGE_DESCRIPTION "NETGEN is an automatic 3d tetrahedral mesh generator. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STL file format. The connection to a geometry kernel allows the handling of IGES and STEP files. NETGEN contains modules for mesh optimization and hierarchical mesh refinement. Netgen is open source based on the LGPL license. It is available for Unix/Linux and Windows.")

    set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/netgen")
    execute_process(COMMAND grep CODENAME /etc/lsb-release OUTPUT_VARIABLE temp OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(temp)
        set(CPACK_GENERATOR "DEB")
        string(SUBSTRING ${temp} 17 -1 UBUNTU_VERSION)
        message("ubuntu version: ${UBUNTU_VERSION}")

        set(CPACK_DEBIAN_PACKAGE_DEPENDS "python3, libtk8.5, libtcl8.5, tix, libxmu6")
        execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Matthias Hochsteger <matthias.hochsteger@tuwien.ac.at>")
        if(USE_MPI)
            set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libmetis5, openmpi-bin")
            set(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}_mpi")
        endif(USE_MPI)
        if(USE_OCC)
            set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, liboce-ocaf-dev")
        endif(USE_OCC)
        set(CPACK_DEBIAN_PACKAGE_SECTION Science)
        set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
        set(CPACK_PACKAGE_FILE_NAME "netgen-${PACKAGE_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")   
    endif(temp)
endif(UNIX)

if(APPLE AND NOT SKBUILD)
    # create some auxiliary files
    set(mac_startup ${CMAKE_CURRENT_BINARY_DIR}/startup.sh)
    file(WRITE ${mac_startup} "\
#!/bin/sh
Netgen_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/MacOS\\/startup.sh//'`\"
Netgen_MACOS=\"$Netgen_BUNDLE/Contents/MacOS\"
export NETGENDIR=$Netgen_MACOS  
export DYLD_LIBRARY_PATH=$Netgen_MACOS:$DYLD_LIBRARY_PATH
# export TIX_LIBRARY=$Netgen_MACOS/library
# export TCLLIBPATH=$Netgen_MACOS:$TCLLIBPATH
export PYTHONPATH=$Netgen_BUNDLE/Contents/Resources/${PYTHON_PACKAGES_INSTALL_DIR}:$PYTHONPATH
cd $Netgen_MACOS
$Netgen_MACOS/netgen
")
    install(PROGRAMS ${mac_startup} DESTINATION ${NG_INSTALL_DIR_BIN})

    set(mac_ngsuite ${CMAKE_CURRENT_BINARY_DIR}/ngsuite.sh)
    file(WRITE ${mac_ngsuite} "\
#!/bin/sh
Netgen_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/MacOS\\/Netgen1//'`\"
Netgen_MACOS=\"$Netgen_BUNDLE/Contents/MacOS\"
open -a /Applications/Utilities/Terminal.app $Netgen_MACOS/startup.sh || open -a /System/Applications/Utilities/Terminal.app $Netgen_MACOS/startup.sh
")
    install(PROGRAMS ${mac_ngsuite} DESTINATION ${NG_INSTALL_DIR_BIN} RENAME Netgen1)

    set(mac_plist ${CMAKE_CURRENT_BINARY_DIR}/Info.plist)
    file(WRITE ${mac_plist} "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<!DOCTYPE plist PUBLIC \"-//Apple Computer//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">
<plist version=\"1.0\">
<dict>
	<key>CFBundleDevelopmentRegion</key>
	<string>English</string>
	<key>CFBundleExecutable</key>
	<string>Netgen1</string>
	<key>CFBundleIconFile</key>
	<string>Netgen.icns</string>
        <key>NSHighResolutionCapable</key>
        <string>True</string>
</dict>
</plist>
")
    install(FILES ${mac_plist} DESTINATION ${NG_INSTALL_DIR_BIN}/../)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/netgen.icns DESTINATION ${NG_INSTALL_DIR_RES}/../ RENAME Netgen.icns)

endif(APPLE AND NOT SKBUILD)

if(NOT APPLE)
  include(CPack) 
endif()

#######################################################################
# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

#######################################################################
# Generate package config file
get_directory_property(NETGEN_COMPILE_DEFINITIONS COMPILE_DEFINITIONS)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NetgenConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/NetgenConfig.cmake @ONLY ESCAPE_QUOTES)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NetgenConfig.cmake DESTINATION ${NG_INSTALL_DIR_CMAKE} COMPONENT netgen_devel)

#######################################################################
# Configure message
# TODO: other message in case of failure
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)
set(flags "${${name}} ${CMAKE_CXX_FLAGS}")
message("
------------------------------------------------------------------------
  ${PROJECT_NAME} ${PACKAGE_VERSION}:  Automatic configuration OK.

  Build type: ${CMAKE_BUILD_TYPE}
  Flags: ${flags}

  Enabled functionality:

    OCC: ............... ${USE_OCC}
    JPEGlib: ........... ${USE_JPEG}
    FFMPEG: ............ ${USE_MPEG}
    GUI: ............... ${USE_GUI}
    MPI: ............... ${USE_MPI}
    PYTHON: ............ ${USE_PYTHON}


  Building: ")
if(WIN32)
message("
    Open ${CMAKE_BINARY_DIR}/Netgen.sln and build solution to compile ${PROJECT_NAME}.

    Build \"INSTALL\" to install ${PROJECT_NAME}.
")
else(WIN32)
message("
    Type 'make' to compile ${PROJECT_NAME}.

    Type 'make install' to install ${PROJECT_NAME}.
")
endif(WIN32)
message("
  Install directory:

    ${CMAKE_INSTALL_PREFIX}

  Please set the following environment variables:

    NETGENDIR=${CMAKE_INSTALL_PREFIX}/bin")
if(USE_PYTHON)
message("    PYTHONPATH=.:${CMAKE_INSTALL_PREFIX}/${PYTHON_PACKAGES_INSTALL_DIR}")
endif(USE_PYTHON)
message("
------------------------------------------------------------------------
")
if (ADDITIONAL_PATHS)
  message(WARNING "The use of ADDITIONAL_PATHS is deprecated, use CMAKE_PREFIX_PATH instead")
endif (ADDITIONAL_PATHS)


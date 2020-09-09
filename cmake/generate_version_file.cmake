if(NOT BDIR)
  set(BDIR ${CMAKE_CURRENT_BINARY_DIR})
endif()

find_package(Git REQUIRED)
execute_process(COMMAND git describe --tags --match "v[0-9]*" --long --dirty WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR} OUTPUT_VARIABLE git_version_string RESULT_VARIABLE status ERROR_QUIET)

if(status AND NOT status EQUAL 0)
  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/../version.txt)
    # for source package files (generated for ubuntu builds on launchpad) read the version from version.txt
    if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/../version.txt)
      file(READ ${CMAKE_CURRENT_LIST_DIR}/../version.txt git_version_string )
    else()
      get_filename_component(git_version_string ${CMAKE_CURRENT_LIST_DIR}/.. NAME)
      string(REGEX REPLACE "^netgen(.*)" "\\1" git_version_string "${git_version_string}")
    endif()
  else()
    MESSAGE(WARNING "Could not determine git-version from source code - assuming 6.2.0.0")
    set(git_version_string "v6.2.0.0")
  endif()
endif()

string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" NETGEN_VERSION_MAJOR "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" NETGEN_VERSION_MINOR "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" NETGEN_VERSION_PATCH "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+\\-([0-9]+).*" "\\1" NETGEN_VERSION_TWEAK "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+\\-[0-9]+\\-([0-9a-z]+).*" "\\1" NETGEN_VERSION_HASH "${git_version_string}")

set(NETGEN_VERSION_SHORT ${NETGEN_VERSION_MAJOR}.${NETGEN_VERSION_MINOR}.${NETGEN_VERSION_PATCH})
set(NETGEN_VERSION_LONG ${NETGEN_VERSION_SHORT}-${NETGEN_VERSION_TWEAK}-${NETGEN_VERSION_HASH})

if(NETGEN_VERSION_TWEAK)
  # no release version - nightly build
  set(NETGEN_VERSION ${NETGEN_VERSION_LONG})
else()
  # TWEAK is 0 -> current version has a tag assigned
  set(NETGEN_VERSION ${NETGEN_VERSION_SHORT})
endif()

set(NETGEN_VERSION_LONG ${NETGEN_VERSION_SHORT}-${NETGEN_VERSION_TWEAK}-${NETGEN_VERSION_HASH})

set(version_file ${BDIR}/netgen_version.hpp)
set(new_version_file_string "\
#ifndef NETGEN_VERSION_HPP_INCLUDED
#define NETGEN_VERSION_HPP_INCLUDED
#define NETGEN_VERSION \"${NETGEN_VERSION}\"
#define NETGEN_VERSION_MAJOR ${NETGEN_VERSION_MAJOR}
#define NETGEN_VERSION_MINOR ${NETGEN_VERSION_MINOR}
#define NETGEN_VERSION_PATCH ${NETGEN_VERSION_PATCH}
#define NETGEN_VERSION_TWEAK ${NETGEN_VERSION_TWEAK}
#define NETGEN_VERSION_HASH \"${NETGEN_VERSION_HASH}\"
#endif // NETGEN_VERSION_HPP_INCLUDED
")
if(EXISTS ${version_file})
  file(READ ${version_file} old_version_file_string )
  if(${old_version_file_string} STREQUAL ${new_version_file_string})
  else()
    file(WRITE ${BDIR}/netgen_version.hpp ${new_version_file_string})
  endif()
else()
    file(WRITE ${BDIR}/netgen_version.hpp ${new_version_file_string})
endif()

file(GENERATE OUTPUT netgen_config.hpp CONTENT
"\
#ifndef NETGEN_CONFIG_HPP_INCLUDED___
#define NETGEN_CONFIG_HPP_INCLUDED___

#define NETGEN_USE_NATIVE_ARCH          $<BOOL:${USE_NATIVE_ARCH}>
#define NETGEN_USE_GUI                  $<BOOL:${USE_GUI}>
#define NETGEN_USE_PYTHON               $<BOOL:${USE_PYTHON}>
#define NETGEN_USE_MPI                  $<BOOL:${USE_MPI}}>
#define NETGEN_USE_MPI4PY               $<BOOL:${USE_MPI4PY}>
#define NETGEN_USE_OCC                  $<BOOL:${USE_OCC}}>
#define NETGEN_USE_JPEG                 $<BOOL:${USE_JPEG}}>
#define NETGEN_USE_MPEG                 $<BOOL:${USE_MPEG}}>
#define NETGEN_USE_CGNS                 $<BOOL:${USE_CGNS}}>
#define NETGEN_USE_NUMA                 $<BOOL:${USE_NUMA}}>
#define NETGEN_INTEL_MIC                $<BOOL:${USE_INTEL_MIC}}>
#define NETGEN_INSTALL_PROFILES         $<BOOL:${INSTALL_PROFILES}>
#define NETGEN_USE_CCACHE               $<BOOL:${USE_CCACHE}}>
#define NETGEN_USE_INTERNAL_TCL         $<BOOL:${USE_INTERNAL_TCL}>
#define NETGEN_ENABLE_UNIT_TESTS        $<BOOL:${ENABLE_UNIT_TESTS}>
#define NETGEN_ENABLE_CPP_CORE_GUIDELINES_CHECK $<BOOL:${ENABLE_CPP_CORE_GUIDELINES_CHECK}>
#define NETGEN_USE_SPDLOG               $<BOOL:${USE_SPDLOG}>
#define NETGEN_DEBUG_LOG                $<BOOL:${DEBUG_LOG}>
#define NETGEN_USE_CHECK_RANGE          $<BOOL:${CHECK_RANGE}>
#define NETGEN_BUILD_STUB_FILES         $<BOOL:${BUILD_STUB_FILES}>
#define NETGEN_BUILD_FOR_CONDA          $<BOOL:${BUILD_FOR_CONDA}>

#endif // NETGEN_CONFIG_HPP_INCLUDED___
")

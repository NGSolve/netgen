if(NOT BDIR)
  set(BDIR ${CMAKE_CURRENT_BINARY_DIR})
endif()

find_package(Git REQUIRED)

if(GIT_FOUND AND IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../.git)
  execute_process(COMMAND git describe --tags --match "v[0-9]*" --long --dirty WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR} OUTPUT_VARIABLE git_version_string)
else()
  # for source package files (generated for ubuntu builds on launchpad) read the version from version.txt
  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/../version.txt)
    file(READ ${CMAKE_CURRENT_LIST_DIR}/../version.txt git_version_string )
  else()
    get_filename_component(git_version_string ${CMAKE_CURRENT_LIST_DIR}/.. NAME)
    string(REGEX REPLACE "^netgen(.*)" "\\1" git_version_string "${git_version_string}")
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


if(UNIX AND NOT APPLE)
  set (LINUX TRUE)
endif()
if(LINUX)
    find_package(TclStub 8.5 REQUIRED)
else(LINUX)
if(SKBUILD)
# we are building a pip package - download the tcl/tk sources matching the tkinter version (for private headers not shipped with python)

execute_process(COMMAND ${Python3_EXECUTABLE} -c
"import tkinter;print(tkinter.Tcl().eval('info patchlevel').replace('.','-'))"
OUTPUT_VARIABLE PYTHON_TCL_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

set(TCL_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_tcl)
set(TK_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_tk)

ExternalProject_Add(project_tcl
  URL "https://github.com/tcltk/tcl/archive/refs/tags/core-${PYTHON_TCL_VERSION}.zip"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  ${SUBPROJECT_ARGS}
  DOWNLOAD_DIR download_tcl
)
ExternalProject_Add(project_tk
  URL "https://github.com/tcltk/tk/archive/refs/tags/core-${PYTHON_TCL_VERSION}.zip"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory macosx generic
  ${SUBPROJECT_ARGS}
  DOWNLOAD_DIR download_tk
  BUILD_IN_SOURCE 1
)

set(TCL_INCLUDE_PATH ${TCL_DIR}/generic)
set(TK_INCLUDE_PATH ${TK_DIR}/generic)
list(APPEND NETGEN_DEPENDENCIES project_tcl project_tk)

if(APPLE OR WIN32)
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.prefix)" OUTPUT_VARIABLE PYTHON_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    file(TO_CMAKE_PATH ${PYTHON_PREFIX} PYTHON_PREFIX)

    set(tcl_find_args
        REQUIRED
        NO_DEFAULT_PATH
        NO_PACKAGE_ROOT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH
        HINTS
        ${PYTHON_PREFIX}/lib
        ${PYTHON_PREFIX}/tcl
        ${PYTHON_PREFIX}/Frameworks
        ${PYTHON_PREFIX}/Frameworks/Tcl.framework
        ${PYTHON_PREFIX}/Frameworks/Tk.framework
        )
    find_library(TCL_STUB_LIBRARY NAMES tclstub85 tclstub8.5 tclstub86 tclstub8.6 ${tcl_find_args})
    find_library(TK_STUB_LIBRARY NAMES tkstub85 tkstub8.5 tkstub86 tkstub8.6 ${tcl_find_args})
    find_library(TCL_LIBRARY NAMES tcl85 tcl8.5 tcl86 tcl8.6 tcl86t Tcl ${tcl_find_args})
    find_library(TK_LIBRARY NAMES tk85 tk8.5 tk86 tk8.6 tk86t Tk ${tcl_find_args})
else()
    # use system tcl/tk on linux
    find_package(TclStub REQUIRED)
endif()

else(SKBUILD)
if(APPLE)
  set(tcl_prefix ${CMAKE_INSTALL_PREFIX})
  # URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.9/tcl8.6.9-src.tar.gz"
  # URL_MD5 aa0a121d95a0e7b73a036f26028538d4
  ExternalProject_Add(project_tcl
    URL "https://github.com/NGSolve/tcl/archive/7769161.zip"
    URL_MD5 1131f188dd26944df557913c475d43b4
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ../project_tcl/macosx/configure --enable-threads --enable-framework --prefix=${tcl_prefix} --libdir=${tcl_prefix}/Contents/Frameworks --bindir=${tcl_prefix}/Contents/Frameworks/Tcl.framework/bin
    BUILD_COMMAND make -j4 binaries libraries
    INSTALL_COMMAND make install-binaries install-headers install-libraries install-private-headers
    ${SUBPROJECT_ARGS}
    )

  # URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.9/tk8.6.9.1-src.tar.gz"
  # URL_MD5 9efe3976468352dc894dae0c4e785a8e
  ExternalProject_Add(project_tk
    DEPENDS project_tcl
    URL "https://github.com/NGSolve/tk/archive/e7c2bc7.zip"
    URL_MD5 94044140d4826069c22f1c60cedb6e59
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ../project_tk/macosx/configure --enable-aqua=yes --enable-threads --enable-framework --prefix=${tcl_prefix} --libdir=${tcl_prefix}/Contents/Frameworks --bindir=${tcl_prefix}/Contents/Frameworks/Tcl.framework/bin --with-tcl=${tcl_prefix}/Contents/Frameworks/Tcl.framework
    BUILD_COMMAND make -j4 binaries libraries
    INSTALL_COMMAND make install-binaries install-headers install-libraries install-private-headers
    ${SUBPROJECT_ARGS}
    )

  ExternalProject_Add(project_tkdnd
    URL "https://src.fedoraproject.org/repo/pkgs/tkdnd/tkdnd2.8-src.tar.gz/a6d47a996ea957416469b12965d4db91/tkdnd2.8-src.tar.gz"
    URL_MD5 a6d47a996ea957416469b12965d4db91
    DEPENDS project_tcl project_tk
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    PATCH_COMMAND  patch < ${CMAKE_CURRENT_LIST_DIR}/tkdnd_macosx.patch
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CMAKE_ARGS
           -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/Contents/MacOS
	   -DTCL_INCLUDE_PATH=${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tcl.framework/Headers
	   -DTK_INCLUDE_PATH=${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tk.framework/Headers
    ${SUBPROJECT_ARGS}
  )

  list(APPEND NETGEN_DEPENDENCIES project_tcl project_tk project_tkdnd)
  list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX}/Contents/Frameworks)
  set(TCL_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tcl.framework/Headers)
  set(TCL_LIBRARY ${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tcl.framework)
  set(TK_LIBRARY ${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tk.framework)
  set(TK_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tk.framework/Headers)

  set(TCL_STUB_LIBRARY ${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tcl.framework/libtclstub8.6.a)
  set(TK_STUB_LIBRARY ${CMAKE_INSTALL_PREFIX}/Contents/Frameworks/Tk.framework/libtkstub8.6.a)

#   # use system tcl/tk
#   if((${PYTHON_VERSION_STRING} VERSION_EQUAL "3.7") OR (${PYTHON_VERSION_STRING} VERSION_GREATER "3.7"))
#     # fetch tcl/tk sources to match the one used in Python 3.7
#     ExternalProject_Add(project_tcl
#       URL "https://prdownloads.sourceforge.net/tcl/tcl8.6.8-src.tar.gz"
#       URL_MD5 81656d3367af032e0ae6157eff134f89
#       DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
#       UPDATE_COMMAND "" # Disable update
#       CONFIGURE_COMMAND ""
#       BUILD_COMMAND ""
#       INSTALL_COMMAND ""
#       )
#     ExternalProject_Add(project_tk
#       URL "https://prdownloads.sourceforge.net/tcl/tk8.6.8-src.tar.gz"
#       URL_MD5 5e0faecba458ee1386078fb228d008ba
#       DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
#       UPDATE_COMMAND "" # Disable update
#       CONFIGURE_COMMAND ""
#       BUILD_COMMAND ""
#       INSTALL_COMMAND ""
#       )
# 
#     get_filename_component(PYTHON_LIB_DIR ${PYTHON_LIBRARY} DIRECTORY)
#     find_library(TCL_LIBRARY libtcl8.6.dylib PATHS ${PYTHON_LIB_DIR} NO_DEFAULT_PATH)
#     find_library(TK_LIBRARY libtk8.6.dylib PATHS ${PYTHON_LIB_DIR} NO_DEFAULT_PATH)
# 
#     set(TCL_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_tcl)
#     set(TK_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_tk)
#     set(TCL_INCLUDE_PATH "${TCL_DIR}/generic;${TCL_DIR}/macosx")
#     set(TK_INCLUDE_PATH "${TK_DIR}/generic;${TK_DIR}/macosx;${TK_DIR}/xlib")
#     string(REPLACE ";" "$<SEMICOLON>" TCL_INC "${TCL_INCLUDE_PATH}")
#     string(REPLACE ";" "$<SEMICOLON>" TK_INC "${TK_INCLUDE_PATH}")
# 
#     ExternalProject_Add(project_tkdnd
#       URL "http://sourceforge.net/projects/tkdnd/files/TkDND/TkDND%202.8/tkdnd2.8-src.tar.gz"
#       URL_MD5 a6d47a996ea957416469b12965d4db91
#       DEPENDS project_tcl project_tk
#       DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
#       PATCH_COMMAND  patch < ${CMAKE_CURRENT_LIST_DIR}/tkdnd_macosx.patch
#       UPDATE_COMMAND "" # Disable update
#       BUILD_IN_SOURCE 1
#       CMAKE_ARGS
# 	      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/Contents/MacOS
# 	      -DTCL_INCLUDE_PATH=${TCL_INC}
# 	      -DTK_INCLUDE_PATH=${TK_INC}
# 	      -DTK_LIBRARY=${TK_LIBRARY}
# 	      -DTCL_LIBRARY=${TCL_LIBRARY}
#       LOG_DOWNLOAD 1
#       LOG_CONFIGURE 1
#       LOG_BUILD 1
#       LOG_INSTALL 1
#     )
#       
# list(APPEND NETGEN_DEPENDENCIES project_tkdnd)
#   else()
#     find_package(TCL 8.5 REQUIRED)
#   endif()

elseif(WIN32)

  ExternalProject_Add(project_win_tcltk
    URL ${TCLTK_DOWNLOAD_URL_WIN}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory lib ${CMAKE_INSTALL_PREFIX}/${NG_INSTALL_DIR_LIB}
	    COMMAND ${CMAKE_COMMAND} -E copy_directory bin ${CMAKE_INSTALL_PREFIX}/${NG_INSTALL_DIR_BIN}
	    COMMAND ${CMAKE_COMMAND} -E copy_directory include ${CMAKE_INSTALL_PREFIX}/${NG_INSTALL_DIR_INCLUDE}

    ${SUBPROJECT_ARGS}
    )

  set (TK_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include)
  set (TCL_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include)
  set (TCL_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tcl86t.lib)
  set (TK_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tk86t.lib)
  set (TCL_STUB_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tclstub86.lib)
  set (TK_STUB_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tkstub86.lib)

  list(APPEND NETGEN_DEPENDENCIES project_win_tcltk)
else(WIN32)
    find_package(TCL 8.5 REQUIRED)
#     ExternalProject_Add(project_tkdnd
#       GIT_REPOSITORY https://github.com/petasis/tkdnd.git
#       GIT_TAG d7cfd96087b248255da5349086ef70cc4bbfb619
#       PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tkdnd
#       CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/lib
#       UPDATE_COMMAND ""
#       LOG_DOWNLOAD 1
#       LOG_BUILD 1
#       LOG_INSTALL 1
# )
# list(APPEND NETGEN_DEPENDENCIES project_tkdnd)
endif(APPLE)
endif(SKBUILD)
endif(LINUX)

# Propagate settings to Netgen subproject
set_vars(NETGEN_CMAKE_ARGS TCL_INCLUDE_PATH TCL_STUB_LIBRARY TCL_LIBRARY TK_STUB_LIBRARY TK_LIBRARY TK_INCLUDE_PATH TCL_TCLSH TK_WISH)

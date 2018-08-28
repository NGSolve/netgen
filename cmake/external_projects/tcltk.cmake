if(APPLE)
  # use system tcl/tk
  if(${PYTHON_VERSION_STRING} STREQUAL "3.7")
    # fetch tcl/tk sources to match the one used in Python 3.7
    ExternalProject_Add(project_tcl
      URL "https://prdownloads.sourceforge.net/tcl/tcl8.6.8-src.tar.gz"
      URL_MD5 81656d3367af032e0ae6157eff134f89
      DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
      UPDATE_COMMAND "" # Disable update
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      )
    ExternalProject_Add(project_tk
      URL "https://prdownloads.sourceforge.net/tcl/tk8.6.8-src.tar.gz"
      URL_MD5 5e0faecba458ee1386078fb228d008ba
      DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
      UPDATE_COMMAND "" # Disable update
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      )

    get_filename_component(PYTHON_LIB_DIR ${PYTHON_LIBRARY} DIRECTORY)
    find_library(TCL_LIBRARY libtcl8.6.dylib PATHS ${PYTHON_LIB_DIR} NO_DEFAULT_PATH)
    find_library(TK_LIBRARY libtk8.6.dylib PATHS ${PYTHON_LIB_DIR} NO_DEFAULT_PATH)

    set(TCL_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_tcl)
    set(TK_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_tk)
    set(TCL_INCLUDE_PATH ${TCL_DIR}/generic;${TCL_DIR}/macosx)
    set(TK_INCLUDE_PATH ${TK_DIR}/generic;${TK_DIR}/macosx;${TK_DIR}/xlib)
    find_package(TCL 8.6 REQUIRED)
    list(APPEND NETGEN_DEPENDENCIES project_tcl project_tk)

  else()
    find_package(TCL 8.5 REQUIRED)
  endif()
#   set(HOME $ENV{HOME})
#   set(tcl_prefix ${CMAKE_INSTALL_PREFIX}/../../)
#   ExternalProject_Add(project_tcl
#     URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tcl8.6.4-src.tar.gz"
#     URL_MD5 d7cbb91f1ded1919370a30edd1534304
#     DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
#     UPDATE_COMMAND "" # Disable update
#     CONFIGURE_COMMAND ../project_tcl/macosx/configure --enable-threads --enable-framework --prefix=${tcl_prefix} --libdir=${tcl_prefix}/Contents/Frameworks --bindir=${tcl_prefix}/Contents/Frameworks/Tcl.framework/bin
#     BUILD_COMMAND make -j4 binaries libraries
#     INSTALL_COMMAND make install-binaries install-headers install-libraries install-private-headers
#     LOG_DOWNLOAD 1
#     LOG_BUILD 1
#     LOG_CONFIGURE 1
#     LOG_INSTALL 1
#     )
# 
#   ExternalProject_Add(project_tk
#     DEPENDS project_tcl
#     URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tk8.6.4-src.tar.gz"
#     URL_MD5 261754d7dc2a582f00e35547777e1fea
#     DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
#     UPDATE_COMMAND "" # Disable update
#     CONFIGURE_COMMAND ../project_tk/macosx/configure --enable-aqua=yes --enable-threads --enable-framework --prefix=${tcl_prefix} --libdir=${tcl_prefix}/Contents/Frameworks --bindir=${tcl_prefix}/Contents/Frameworks/Tcl.framework/bin --with-tcl=${tcl_prefix}/Contents/Frameworks/Tcl.framework
#     BUILD_COMMAND make -j4 binaries libraries
#     INSTALL_COMMAND make install-binaries install-headers install-libraries install-private-headers
#     LOG_DOWNLOAD 1
#     LOG_BUILD 1
#     LOG_CONFIGURE 1
#     LOG_INSTALL 1
#     )
# 
#   ExternalProject_Add(project_tkdnd
#     URL "https://sourceforge.net/projects/tkdnd/files/OS%20X%20Binaries/TkDND%202.8/tkdnd2.8-OSX-MountainLion.tar.gz"
#     URL_MD5 2dbb471b1d66c5f391f3c3c5b71548fb
#     DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
#     BUILD_IN_SOURCE 1
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND ""
#     INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}/../MacOS
#     LOG_DOWNLOAD 1
#     LOG_CONFIGURE 1
#     LOG_BUILD 1
#     LOG_INSTALL 1
#     )
#  
#   list(APPEND NETGEN_DEPENDENCIES project_tcl project_tk project_tkdnd)
#   list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX}../Frameworks)
#   set(TCL_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/../Frameworks/Tcl.framework/Headers)
#   set(TCL_LIBRARY ${CMAKE_INSTALL_PREFIX}/../Frameworks/Tcl.framework)
#   set(TK_LIBRARY ${CMAKE_INSTALL_PREFIX}/../Frameworks/Tk.framework)
#   set(TK_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/../Frameworks/Tk.framework/Headers)
# 



ExternalProject_Add(project_tkdnd
  URL "http://sourceforge.net/projects/tkdnd/files/TkDND/TkDND%202.8/tkdnd2.8-src.tar.gz"
  URL_MD5 a6d47a996ea957416469b12965d4db91
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  PATCH_COMMAND  patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/tkdnd_macosx.patch
  UPDATE_COMMAND "" # Disable update
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX}/Contents/MacOS --libdir=${CMAKE_INSTALL_PREFIX}/Contents/MacOS
  BUILD_COMMAND make
  INSTALL_COMMAND make install
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)
list(APPEND NETGEN_DEPENDENCIES project_tkdnd)

elseif(WIN32)

  ExternalProject_Add(project_win_tcltk
    URL ${TCLTK_DOWNLOAD_URL_WIN}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}
    LOG_DOWNLOAD 1
    )

  set (TK_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include)
  set (TCL_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include)
  set (TCL_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tcl86t.lib)
  set (TK_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tk86t.lib)

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

# Propagate settings to Netgen subproject
set_vars(NETGEN_CMAKE_ARGS TCL_INCLUDE_PATH TCL_LIBRARY TK_LIBRARY TK_INCLUDE_PATH TCL_TCLSH TK_WISH)

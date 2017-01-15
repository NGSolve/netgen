if(APPLE)
  set(HOME $ENV{HOME})
  ExternalProject_Add(tcl
    URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tcl8.6.4-src.tar.gz"
    URL_MD5 d7cbb91f1ded1919370a30edd1534304
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C macosx install INSTALL_ROOT=${CMAKE_INSTALL_PREFIX}/../ INSTALL_PATH=/Frameworks #NATIVE_TCLSH=${HOME}/usr/local/bin/tclsh
    INSTALL_COMMAND ""
    LOG_DOWNLOAD 1
    LOG_BUILD 1
    LOG_INSTALL 1
    )

  ExternalProject_Add(tk
    DEPENDS tcl
    URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tk8.6.4-src.tar.gz"
    URL_MD5 261754d7dc2a582f00e35547777e1fea
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C macosx install INSTALL_ROOT=${CMAKE_INSTALL_PREFIX}/../ INSTALL_PATH=/Frameworks 
    INSTALL_COMMAND ""#make -C macosx install
    LOG_DOWNLOAD 1
    LOG_BUILD 1
    LOG_INSTALL 1
    )

  ExternalProject_Add(tkdnd
    DEPENDS tcl tk
    URL "http://sourceforge.net/projects/tkdnd/files/TkDND/TkDND%202.8/tkdnd2.8-src.tar.gz"
    URL_MD5 a6d47a996ea957416469b12965d4db91
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
    PATCH_COMMAND  patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/tkdnd_macosx.patch
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./configure --with-tcl=${CMAKE_INSTALL_PREFIX}/../Frameworks/Tcl.framework --with-tk=${CMAKE_INSTALL_PREFIX}/../Frameworks/Tk.framework --prefix=${CMAKE_INSTALL_PREFIX}/../MacOS --libdir=${CMAKE_INSTALL_PREFIX}/../MacOS
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
    )

  list(APPEND NETGEN_DEPENDENCIES tcl tk tkdnd)
  set(TCL_INCLUDE_PATH ${HOME}/Library/Frameworks/Tcl.framework/Headers)
  set(TCL_LIBRARY ${HOME}/Library/Frameworks/Tcl.framework)
  set(TK_LIBRARY ${HOME}/Library/Frameworks/Tk.framework)
  set(TK_INCLUDE_PATH ${HOME}/Library/Frameworks/Tk.framework/Headers)
  set(TCL_TCLSH ${HOME}/usr/local/bin/tclsh)
  set(TK_WISH ${HOME}/usr/local/bin/wish)

elseif(WIN32)

  ExternalProject_Add(win_extlibs
    URL ${EXT_LIBS_DOWNLOAD_URL_WIN}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${INSTALL_DIR}
    LOG_DOWNLOAD 1
    )

  list(APPEND NETGEN_DEPENDENCIES win_extlibs)
else(WIN32)
    find_package(TCL 8.5 REQUIRED)
endif(APPLE)

# Propagate settings to Netgen subproject
set_vars(NETGEN_CMAKE_ARGS TCL_INCLUDE_PATH TCL_LIBRARY TK_LIBRARY TK_INCLUDE_PATH TCL_TCLSH TK_WISH)

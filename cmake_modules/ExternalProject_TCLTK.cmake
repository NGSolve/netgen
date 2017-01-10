if(APPLE)
  set(HOME $ENV{HOME})
  ExternalProject_Add(tcl
    URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tcl8.6.4-src.tar.gz"
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C macosx install-embedded INSTALL_ROOT=${INSTALL_DIR}/../../ INSTALL_PATH=/Frameworks #NATIVE_TCLSH=${HOME}/usr/local/bin/tclsh
    INSTALL_COMMAND ""
    LOG_DOWNLOAD 1
    LOG_BUILD 1
    LOG_INSTALL 1
    )

  ExternalProject_Add(tk
    DEPENDS tcl
    URL "http://sourceforge.net/projects/tcl/files/Tcl/8.6.4/tk8.6.4-src.tar.gz"
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C macosx install-embedded INSTALL_ROOT=${INSTALL_DIR}/../../ INSTALL_PATH=/Frameworks 
    INSTALL_COMMAND ""#make -C macosx install
    LOG_DOWNLOAD 1
    LOG_BUILD 1
    LOG_INSTALL 1
    )

  ExternalProject_Add(tkdnd
    DEPENDS tcl tk
    URL "http://sourceforge.net/projects/tkdnd/files/TkDND/TkDND%202.8/tkdnd2.8-src.tar.gz"
    PATCH_COMMAND  patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/tkdnd_macosx.patch
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./configure --with-tcl=${INSTALL_DIR}/../../Frameworks/Tcl.framework --with-tk=${INSTALL_DIR}/../../Frameworks/Tk.framework --prefix=${INSTALL_DIR}
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

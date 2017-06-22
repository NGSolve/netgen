if(APPLE)
  # use system tcl/tk
  find_package(TCL 8.5 REQUIRED)
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
elseif(WIN32)

  ExternalProject_Add(project_win_extlibs
    URL ${EXT_LIBS_DOWNLOAD_URL_WIN}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}
    LOG_DOWNLOAD 1
    )

  set (TK_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include)
  set (TCL_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include)
  set (TCL_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tcl86.lib)
  set (TK_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/tk86.lib)

  list(APPEND NETGEN_DEPENDENCIES project_win_extlibs)
else(WIN32)
    find_package(TCL 8.5 REQUIRED)
endif(APPLE)

# Propagate settings to Netgen subproject
set_vars(NETGEN_CMAKE_ARGS TCL_INCLUDE_PATH TCL_LIBRARY TK_LIBRARY TK_INCLUDE_PATH TCL_TCLSH TK_WISH)

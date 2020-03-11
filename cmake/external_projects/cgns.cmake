if(WIN32)

  ExternalProject_Add(project_win_cgns
    URL ${CGNS_DOWNLOAD_URL_WIN}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}
    LOG_DOWNLOAD 1
    )

  list(APPEND NETGEN_DEPENDENCIES project_win_cgns)
endif(WIN32)

if(APPLE)
  ExternalProject_Add(project_mac_cgns
    URL ${CGNS_DOWNLOAD_URL_MAC}
    UPDATE_COMMAND "" # Disable update
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}
    LOG_DOWNLOAD 1
    )

  list(APPEND NETGEN_DEPENDENCIES project_mac_cgns)
  list(APPEND NETGEN_CMAKE_ARGS "-DCGNS_INCLUDE_DIR=${CMAKE_INSTALL_PREFIX}/Contents/Resources/include")
  list(APPEND NETGEN_CMAKE_ARGS "-DCGNS_LIBRARY=${CMAKE_INSTALL_PREFIX}/Contents/MacOS/libcgns.dylib")
endif(APPLE)

if(USE_INTERNAL_TCL)
  add_definitions(-DINTERNAL_TCL_DEFAULT=1)
else()
  add_definitions(-DINTERNAL_TCL_DEFAULT=0)
endif()

set(netgen_sources ngappinit.cpp onetcl.cpp)
if(WIN32)
    # add icon to netgen executable
    enable_language(RC)
    set(netgen_sources ${netgen_sources} ../windows/netgen.rc)
    # Don't use ccache here due to incompatibility with the resource compiler
    set_directory_properties(PROPERTIES RULE_LAUNCH_COMPILE "")
endif(WIN32)

if(USE_GUI)

    add_library(gui SHARED
        gui.cpp ngpkg.cpp demoview.cpp parallelfunc.cpp
        ../libsrc/stlgeom/stlpkg.cpp ../libsrc/visualization/visualpkg.cpp
        ../libsrc/csg/csgpkg.cpp ../libsrc/geom2d/geom2dpkg.cpp
        ../libsrc/occ/occpkg.cpp ../libsrc/occ/vsocc.cpp
    )

    add_executable(netgen ngappinit.cpp onetcl.cpp)

    target_link_libraries( gui PUBLIC nglib )
    target_link_libraries( gui PRIVATE ${LIBTOGL} ${ZLIB_LIBRARIES} ${JPEG_LIBRARIES} ${FFMPEG_LIBRARIES} ${X11_Xmu_LIB} ${X11_X11_LIB} ${OCC_LIBRARIES} )
    if(NOT APPLE)
      target_link_libraries( gui PRIVATE ${TCL_LIBRARY} ${TK_LIBRARY})
    endif(NOT APPLE)

    target_link_libraries( netgen nglib gui ${LIBTOGL} ${ZLIB_LIBRARIES} ${JPEG_LIBRARIES} ${FFMPEG_LIBRARIES} ${X11_Xmu_LIB} ${X11_X11_LIB} ${OCC_LIBRARIES} ${TK_LIBRARY} ${TCL_LIBRARY})

    if(NOT WIN32)
      target_link_libraries( netgen mesh stlvis stl geom2dvis interface geom2d csg stl visual csgvis )
      target_link_libraries( gui PUBLIC mesh stlvis stl geom2dvis interface geom2d csg stl visual csgvis )
    endif(NOT WIN32)

    install(TARGETS netgen ${NG_INSTALL_DIR})
    install(TARGETS gui ${NG_INSTALL_DIR})

    if(APPLE)
        set_target_properties(netgen PROPERTIES OUTPUT_NAME netgen)
    endif(APPLE)
    if(WIN32)
        set_target_properties( gui PROPERTIES OUTPUT_NAME libgui )
    endif(WIN32)

endif(USE_GUI)

if(USE_PYTHON)
    add_library(ngpy SHARED netgenpy.cpp)
    target_link_libraries( ngpy PUBLIC nglib PRIVATE ${PYTHON_LIBRARIES})
    if(APPLE)
        set_target_properties( ngpy PROPERTIES SUFFIX ".so")
    elseif(WIN32)
        set_target_properties( ngpy PROPERTIES SUFFIX ".pyd")
        set_target_properties( ngpy PROPERTIES OUTPUT_NAME "libngpy")
    endif()
    MESSAGE(STATUS "######################## CHECKING RPATH #######################################")
    MESSAGE(STATUS "${NG_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    set_target_properties(ngpy PROPERTIES INSTALL_RPATH "${NG_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    MESSAGE(STATUS "###############################################################################")
    install(TARGETS ngpy DESTINATION ${NG_INSTALL_DIR_PYTHON}/${NG_INSTALL_SUFFIX} COMPONENT netgen)
endif(USE_PYTHON)

if(USE_GUI)
  if(NOT USE_INTERNAL_TCL)
    install(FILES
        dialog.tcl menustat.tcl ngicon.tcl ng.tcl
        ngvisual.tcl sockets.tcl nghelp.tcl ngshell.tcl
        ngtesting.tcl parameters.tcl variables.tcl csgeom.tcl stlgeom.tcl
        occgeom.tcl acisgeom.tcl netgen.ocf drawing.tcl
        DESTINATION ${NG_INSTALL_DIR_BIN} COMPONENT netgen)
  endif()

  add_subdirectory(Togl2.1)
endif(USE_GUI)


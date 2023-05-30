if(USE_INTERNAL_TCL)
  add_definitions(-DINTERNAL_TCL_DEFAULT=1)
else()
  add_definitions(-DINTERNAL_TCL_DEFAULT=0)
endif()

if(WIN32)
  # add icon and version info to netgen executable
  enable_language(RC)
  # Don't use ccache here due to incompatibility with the resource compiler
  set_directory_properties(PROPERTIES RULE_LAUNCH_COMPILE "")
endif(WIN32)

target_sources(nglib PRIVATE onetcl.cpp)

if(USE_GUI)
    target_sources(nggui PRIVATE
        gui.cpp ngpkg.cpp demoview.cpp parallelfunc.cpp ngtcl.cpp
    )

    if(NOT BUILD_FOR_CONDA)
      add_executable(netgen ngappinit.cpp)
    if(WIN32)
      target_sources(netgen PRIVATE ../windows/netgen.rc)
    endif(WIN32)
    target_link_libraries( netgen nglib nggui netgen_python netgen_gui)
    install(TARGETS netgen ${NG_INSTALL_DIR})
    if(APPLE)
        set_target_properties(netgen PROPERTIES OUTPUT_NAME netgen)
    endif(APPLE)
    target_link_libraries( netgen ${PYTHON_LIBRARIES} ${TCL_LIBRARY} ${TK_LIBRARY})
    endif(NOT BUILD_FOR_CONDA)

    install(TARGETS nggui ${NG_INSTALL_DIR})
endif(USE_GUI)

if(USE_PYTHON)
    add_library(ngpy SHARED netgenpy.cpp)
    target_link_libraries( ngpy PUBLIC nglib PRIVATE "$<BUILD_INTERFACE:netgen_python>" )
    if(APPLE)
        set_target_properties( ngpy PROPERTIES SUFFIX ".so")
    elseif(WIN32)
        set_target_properties( ngpy PROPERTIES SUFFIX ".pyd")
        set_target_properties( ngpy PROPERTIES OUTPUT_NAME "libngpy")
    endif()
    set_target_properties(ngpy PROPERTIES INSTALL_RPATH "${NG_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    install(TARGETS ngpy DESTINATION ${NG_INSTALL_DIR_PYTHON}/${NG_INSTALL_SUFFIX} COMPONENT netgen)

    if(USE_GUI)
        add_library(ngguipy SHARED ngguipy.cpp)
        target_link_libraries( ngguipy PUBLIC nglib nggui PRIVATE "$<BUILD_INTERFACE:netgen_python>" $<BUILD_INTERFACE:netgen_gui>)
        if(APPLE)
            set_target_properties( ngguipy PROPERTIES SUFFIX ".so")
        elseif(WIN32)
            set_target_properties( ngguipy PROPERTIES SUFFIX ".pyd")
            set_target_properties( ngguipy PROPERTIES OUTPUT_NAME "libngguipy")
        endif()
        set_target_properties(ngguipy PROPERTIES INSTALL_RPATH "${NG_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
        install(TARGETS ngguipy DESTINATION ${NG_INSTALL_DIR_PYTHON}/${NG_INSTALL_SUFFIX} COMPONENT netgen)
    endif(USE_GUI)
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


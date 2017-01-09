add_library(stl ${NG_LIB_TYPE} 
        meshstlsurface.cpp stlgeom.cpp stlgeomchart.cpp 
        stlgeommesh.cpp stlline.cpp stltool.cpp stltopology.cpp python_stl.cpp
        )

if(NOT WIN32)
    target_link_libraries( stl mesh  ${PYTHON_LIBRARIES})
    target_link_libraries( stl ${PYTHON_LIBRARIES})
    install( TARGETS stl ${ng_install_dir})
endif(NOT WIN32)

if(USE_GUI)
    add_library(stlvis ${NG_LIB_TYPE}
            vsstl.cpp
            )
    if(NOT WIN32)
        target_link_libraries( stlvis stl )
        install( TARGETS stlvis ${ng_install_dir})
    endif(NOT WIN32)
endif(USE_GUI)


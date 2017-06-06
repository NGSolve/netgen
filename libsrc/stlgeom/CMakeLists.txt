add_library(stl ${NG_LIB_TYPE} 
        meshstlsurface.cpp stlgeom.cpp stlgeomchart.cpp 
        stlgeommesh.cpp stlline.cpp stltool.cpp stltopology.cpp python_stl.cpp
        )

if(NOT WIN32)
    target_link_libraries( stl mesh  ${PYTHON_LIBRARIES})
    target_link_libraries( stl ${PYTHON_LIBRARIES})
    install( TARGETS stl ${NG_INSTALL_DIR})
endif(NOT WIN32)

if(USE_GUI)
    add_library(stlvis ${NG_LIB_TYPE}
            vsstl.cpp
            )
    if(NOT WIN32)
        target_link_libraries( stlvis stl )
        install( TARGETS stlvis ${NG_INSTALL_DIR})
    endif(NOT WIN32)
endif(USE_GUI)

install(FILES
  meshstlsurface.hpp stlgeom.hpp stlline.hpp
  stltool.hpp stltopology.hpp vsstl.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/stlgeom COMPONENT netgen_devel
)

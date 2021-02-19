add_definitions(-DNGINTERFACE_EXPORTS)
add_library(stl ${NG_LIB_TYPE} 
        meshstlsurface.cpp stlgeom.cpp stlgeomchart.cpp 
        stlgeommesh.cpp stlline.cpp stltool.cpp stltopology.cpp python_stl.cpp
        )

if(NOT WIN32)
    target_link_libraries( stl PUBLIC mesh )
    install( TARGETS stl ${NG_INSTALL_DIR})
endif(NOT WIN32)

target_link_libraries( stl PUBLIC ngcore PRIVATE "$<BUILD_INTERFACE:netgen_python>" )

if(USE_GUI)
    add_library(stlvis ${NG_LIB_TYPE} vsstl.cpp)
    target_link_libraries(stlvis PRIVATE "$<BUILD_INTERFACE:netgen_python>" PUBLIC ngcore)
    if(NOT WIN32)
        target_link_libraries( stlvis PUBLIC stl )
        install( TARGETS stlvis ${NG_INSTALL_DIR})
    endif(NOT WIN32)
endif(USE_GUI)

install(FILES
  meshstlsurface.hpp stlgeom.hpp stlline.hpp
  stltool.hpp stltopology.hpp vsstl.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/stlgeom COMPONENT netgen_devel
)

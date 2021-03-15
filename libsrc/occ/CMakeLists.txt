add_definitions(-DNGINTERFACE_EXPORTS)
add_library(occ ${NG_LIB_TYPE}
        Partition_Inter2d.cxx Partition_Inter3d.cxx 
        Partition_Loop.cxx Partition_Loop2d.cxx Partition_Loop3d.cxx Partition_Spliter.cxx 
        occgenmesh.cpp occgeom.cpp occmeshsurf.cpp python_occ.cpp
        )
if(USE_GUI)
        add_library(occvis ${NG_LIB_TYPE} vsocc.cpp)
        target_link_libraries(occvis PUBLIC ngcore)
endif(USE_GUI)

target_link_libraries(occ PUBLIC ngcore PRIVATE "$<BUILD_INTERFACE:netgen_python>")

if(NOT WIN32)
    target_link_libraries( occ PRIVATE ${OCC_LIBRARIES} )
    if(USE_OCC AND APPLE)
        # Link AppKit in case OCE was built as static libraries
        find_library(AppKit AppKit)
        target_link_libraries( occ PRIVATE ${AppKit} )
    endif(USE_OCC AND APPLE)
    install( TARGETS occ ${NG_INSTALL_DIR})
    if (USE_GUI)
        target_link_libraries( occvis PUBLIC occ )
        install( TARGETS occvis ${NG_INSTALL_DIR})
    endif(USE_GUI)
endif(NOT WIN32)

install(FILES
  occgeom.hpp occmeshsurf.hpp vsocc.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/occ COMPONENT netgen_devel
)

add_definitions(-DNGINTERFACE_EXPORTS)
add_library(occ ${NG_LIB_TYPE}
        Partition_Inter2d.cxx Partition_Inter3d.cxx 
        Partition_Loop.cxx Partition_Loop2d.cxx Partition_Loop3d.cxx Partition_Spliter.cxx 
        occconstruction.cpp occgenmesh.cpp occgeom.cpp occmeshsurf.cpp python_occ.cpp
        )

add_library(occvis ${NG_LIB_TYPE} vsocc.cpp)

if(NOT WIN32)
    target_link_libraries( occ ${OCC_LIBRARIES} ${PYTHON_LIBRARIES})
    target_link_libraries( occvis occ )
    install( TARGETS occ occvis ${NG_INSTALL_DIR})
endif(NOT WIN32)

install(FILES
  occgeom.hpp occmeshsurf.hpp vsocc.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/occ COMPONENT netgen_devel
)

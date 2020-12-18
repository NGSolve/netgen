add_definitions(-DNGINTERFACE_EXPORTS)
install(FILES soldata.hpp DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel )

if(USE_GUI)
    set( LIB_VISUAL_SOURCES meshdoc.cpp mvdraw.cpp vsfieldlines.cpp vsmesh.cpp vssolution.cpp importsolution.cpp )
else(USE_GUI)
    set( LIB_VISUAL_SOURCES visual_dummy.cpp )
endif(USE_GUI)

add_library(visual ${NG_LIB_TYPE} ${LIB_VISUAL_SOURCES})

target_link_libraries( visual PUBLIC ngcore PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${MPI_CXX_LIBRARIES} ${OPENGL_LIBRARIES} )
install( TARGETS visual ${NG_INSTALL_DIR})

install(FILES
  meshdoc.hpp mvdraw.hpp
  vispar.hpp visual.hpp vssolution.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/visualization COMPONENT netgen_devel
)

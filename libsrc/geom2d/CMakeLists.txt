add_definitions(-DNGLIB_EXPORTS)
add_library(geom2d ${NG_LIB_TYPE} genmesh2d.cpp geom2dmesh.cpp geometry2d.cpp python_geom2d.cpp )
if(APPLE)
    set_target_properties( geom2d PROPERTIES SUFFIX ".so")
endif(APPLE)

if(NOT WIN32)
    target_link_libraries(geom2d mesh ${PYTHON_LIBRARIES})
    install( TARGETS geom2d ${NG_INSTALL_DIR})
endif(NOT WIN32)

if(USE_GUI)
    add_library(geom2dvis ${NG_LIB_TYPE} vsgeom2d.cpp)
    if(NOT WIN32)
        target_link_libraries(geom2dvis geom2d)
        install( TARGETS geom2dvis ${NG_INSTALL_DIR})
    endif(NOT WIN32)
endif(USE_GUI)

install(FILES
  geom2dmesh.hpp geometry2d.hpp spline2d.hpp
  splinegeometry2.hpp vsgeom2d.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/geom2d COMPONENT netgen_devel
)

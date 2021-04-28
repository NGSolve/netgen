add_definitions(-DNGINTERFACE_EXPORTS)
add_library(geom2d ${NG_LIB_TYPE} csg2d.cpp genmesh2d.cpp geometry2d.cpp python_geom2d.cpp )
if(APPLE)
    set_target_properties( geom2d PROPERTIES SUFFIX ".so")
endif(APPLE)

target_link_libraries(geom2d PUBLIC ngcore mesh PRIVATE "$<BUILD_INTERFACE:netgen_python>")
if(NOT WIN32)
  install( TARGETS geom2d ${NG_INSTALL_DIR})
endif(NOT WIN32)

if(USE_GUI)
    add_library(geom2dvis ${NG_LIB_TYPE} vsgeom2d.cpp)
    target_link_libraries(geom2dvis PUBLIC ngcore)
    if(NOT WIN32)
        target_link_libraries(geom2dvis PUBLIC geom2d)
        install( TARGETS geom2dvis ${NG_INSTALL_DIR})
    endif(NOT WIN32)
endif(USE_GUI)

install(FILES
  geometry2d.hpp spline2d.hpp
  vsgeom2d.hpp csg2d.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/geom2d COMPONENT netgen_devel
)

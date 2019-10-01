add_definitions(-DNGINTERFACE_EXPORTS)
add_library(gprim INTERFACE)
set(sdir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(gprim INTERFACE
        ${sdir}/adtree.cpp ${sdir}/geom2d.cpp ${sdir}/geom3d.cpp ${sdir}/geomfuncs.cpp
        ${sdir}/geomtest3d.cpp ${sdir}/transform3d.cpp ${sdir}/spline.cpp ${sdir}/splinegeometry.cpp
        )

install(FILES
  adtree.hpp geom2d.hpp geom3d.hpp geomfuncs.hpp
  geomobjects.hpp geomops2.hpp geomops.hpp geomtest3d.hpp gprim.hpp
  splinegeometry.hpp spline.hpp transform3d.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/gprim COMPONENT netgen_devel
)

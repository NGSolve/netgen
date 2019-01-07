add_definitions(-DNGINTERFACE_EXPORTS)
add_library(gprim OBJECT
        adtree.cpp geom2d.cpp geom3d.cpp geomfuncs.cpp 
        geomtest3d.cpp transform3d.cpp spline.cpp splinegeometry.cpp
        )

set_target_properties( gprim PROPERTIES POSITION_INDEPENDENT_CODE ON )

install(FILES
  adtree.hpp geom2d.hpp geom3d.hpp geomfuncs.hpp geomobjects2.hpp
  geomobjects.hpp geomops2.hpp geomops.hpp geomtest3d.hpp gprim.hpp
  splinegeometry.hpp spline.hpp transform3d.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/gprim COMPONENT netgen_devel
)

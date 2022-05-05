target_sources(nglib PRIVATE
    adtree.cpp
    geom2d.cpp
    geom3d.cpp
    geomfuncs.cpp
    geomtest3d.cpp
    spline.cpp
    splinegeometry.cpp
    transform3d.cpp
)

install(FILES
  adtree.hpp geom2d.hpp geom3d.hpp geomfuncs.hpp
  geomobjects.hpp geomops2.hpp geomops.hpp geomtest3d.hpp gprim.hpp
  splinegeometry.hpp spline.hpp transform3d.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/gprim COMPONENT netgen_devel
)

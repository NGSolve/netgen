target_sources(nglib PRIVATE
        algprim.cpp brick.cpp   
        bspline2d.cpp csgeom.cpp csgparser.cpp curve2d.cpp edgeflw.cpp	       
        explicitcurve2d.cpp extrusion.cpp gencyl.cpp genmesh.cpp identify.cpp  
        manifold.cpp meshsurf.cpp polyhedra.cpp revolution.cpp singularref.cpp 
        solid.cpp specpoin.cpp spline3d.cpp surface.cpp triapprox.cpp zrefine.cpp 
        python_csg.cpp splinesurface.cpp
)

if(USE_GUI)
    target_sources(nggui PRIVATE vscsg.cpp csgpkg.cpp)
endif(USE_GUI)

install(FILES
  algprim.hpp brick.hpp csgeom.hpp csg.hpp csgparser.hpp
  curve2d.hpp edgeflw.hpp explicitcurve2d.hpp extrusion.hpp
  gencyl.hpp geoml.hpp identify.hpp manifold.hpp meshsurf.hpp
  polyhedra.hpp revolution.hpp singularref.hpp solid.hpp
  specpoin.hpp spline3d.hpp splinesurface.hpp surface.hpp
  triapprox.hpp vscsg.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/csg COMPONENT netgen_devel
)

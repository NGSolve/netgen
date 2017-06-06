add_definitions(-DNGINTERFACE_EXPORTS)
add_library(csg ${NG_LIB_TYPE}
        algprim.cpp brick.cpp   
        bspline2d.cpp csgeom.cpp csgparser.cpp curve2d.cpp edgeflw.cpp	       
        explicitcurve2d.cpp extrusion.cpp gencyl.cpp genmesh.cpp identify.cpp  
        manifold.cpp meshsurf.cpp polyhedra.cpp revolution.cpp singularref.cpp 
        solid.cpp specpoin.cpp spline3d.cpp surface.cpp triapprox.cpp zrefine.cpp 
        python_csg.cpp splinesurface.cpp
        )
if(APPLE)
    set_target_properties( csg PROPERTIES SUFFIX ".so")
endif(APPLE)

if(NOT WIN32)
    target_link_libraries(csg mesh ${PYTHON_LIBRARIES})
    target_link_libraries(csg ${PYTHON_LIBRARIES})
    install( TARGETS csg ${NG_INSTALL_DIR})
endif(NOT WIN32)


if(USE_GUI)
    add_library(csgvis ${NG_LIB_TYPE} vscsg.cpp )
    if(NOT WIN32)
        target_link_libraries(csgvis csg visual)
        if(APPLE)
            set_target_properties( csgvis PROPERTIES SUFFIX ".so")
        endif(APPLE)
        install( TARGETS csgvis ${NG_INSTALL_DIR})
    endif(NOT WIN32)
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

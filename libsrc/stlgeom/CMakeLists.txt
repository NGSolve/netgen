target_sources(nglib PRIVATE
        meshstlsurface.cpp stlgeom.cpp stlgeomchart.cpp 
        stlgeommesh.cpp stlline.cpp stltool.cpp stltopology.cpp python_stl.cpp
)

if(USE_GUI)
    target_sources(nggui PRIVATE vsstl.cpp stlpkg.cpp)
endif(USE_GUI)

install(FILES
  meshstlsurface.hpp stlgeom.hpp stlline.hpp
  stltool.hpp stltopology.hpp vsstl.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/stlgeom COMPONENT netgen_devel
)

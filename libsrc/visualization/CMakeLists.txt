target_sources(nggui PRIVATE
    importsolution.cpp
    meshdoc.cpp
    mvdraw.cpp
    vsfieldlines.cpp
    vsmesh.cpp
    vssolution.cpp
    visualpkg.cpp
)
target_link_libraries( nggui PUBLIC "$<BUILD_INTERFACE:netgen_python>" ${OPENGL_LIBRARIES} nglib)

install(FILES
  meshdoc.hpp mvdraw.hpp visual_api.hpp
  vispar.hpp visual.hpp vssolution.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/visualization COMPONENT netgen_devel
)

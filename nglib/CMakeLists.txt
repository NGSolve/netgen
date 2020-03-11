add_definitions(-DNGLIB_EXPORTS)

if(WIN32)
    set(nglib_objects 
          $<TARGET_OBJECTS:mesh>
          $<TARGET_OBJECTS:stl>
          $<TARGET_OBJECTS:interface>
          $<TARGET_OBJECTS:geom2d>
          $<TARGET_OBJECTS:csg>

          $<TARGET_OBJECTS:visual>
          $<TARGET_OBJECTS:occ>
       )
       if(USE_GUI)
           set(nglib_objects ${nglib_objects}
	       $<TARGET_OBJECTS:stlvis>
	       $<TARGET_OBJECTS:geom2dvis>
               $<TARGET_OBJECTS:csgvis>
           )
       endif(USE_GUI)
endif(WIN32)

add_library(nglib SHARED nglib.cpp ${nglib_objects})
if(NOT WIN32)
      target_link_libraries( nglib PUBLIC mesh interface geom2d csg stl visual)
      if(USE_GUI)
          target_link_libraries( nglib PUBLIC stlvis geom2dvis csgvis )
      endif(USE_GUI)
endif(NOT WIN32)

# target_link_libraries(nglib PRIVATE gen la gprim PUBLIC ngcore)
target_link_libraries(nglib PUBLIC ngcore)

target_link_libraries( nglib PRIVATE ${OCC_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OPENGL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${X11_Xmu_LIB} ${JPEG_LIBRARIES} ${MKL_LIBRARIES} ${ZLIB_LIBRARIES} ${OCC_LIBRARIES} netgen_cgns )

if(USE_OCC AND NOT WIN32)
    target_link_libraries(nglib PUBLIC occ)
endif(USE_OCC AND NOT WIN32)

if(USE_PYTHON)
    target_link_libraries(nglib PRIVATE ${PYTHON_LIBRARIES})
endif(USE_PYTHON)

install(TARGETS nglib ${NG_INSTALL_DIR})
install(FILES nglib.h DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)

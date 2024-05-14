target_sources(nglib PRIVATE nglib.cpp)

if(USE_OCC)
    target_sources(nglib PRIVATE nglib_occ.cpp)
    install(FILES nglib_occ.h DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)
endif(USE_OCC)

if(EMSCRIPTEN)
  target_compile_options(nglib PUBLIC $<BUILD_INTERFACE:$<TARGET_PROPERTY:ngcore,INTERFACE_COMPILE_OPTIONS>>)
  target_compile_definitions(nglib PUBLIC $<TARGET_PROPERTY:ngcore,INTERFACE_COMPILE_DEFINITIONS>)
  target_include_directories(nglib PUBLIC $<TARGET_PROPERTY:ngcore,INTERFACE_INCLUDE_DIRECTORIES>)
else(EMSCRIPTEN)
  target_link_libraries(nglib PUBLIC ngcore)
  target_link_libraries( nglib PRIVATE ${CMAKE_THREAD_LIBS_INIT} ${JPEG_LIBRARIES} ${MKL_LIBRARIES} ${ZLIB_LIBRARIES} occ_libs netgen_cgns )
endif(EMSCRIPTEN)

install(TARGETS nglib netgen_cgns ${NG_INSTALL_DIR})
install(FILES nglib.h DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)

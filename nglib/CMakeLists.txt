target_sources(nglib PRIVATE nglib.cpp)

if(USE_OCC)
    target_sources(nglib PRIVATE nglib_occ.cpp)
    install(FILES nglib_occ.h DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)
endif(USE_OCC)

target_link_libraries(nglib PUBLIC ngcore PRIVATE ${ZLIB_LIBRARIES})

if(EMSCRIPTEN)
  target_include_directories(nglib PRIVATE $<TARGET_PROPERTY:occ_libs,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(nglib PRIVATE $<BUILD_INTERFACE:$<TARGET_PROPERTY:occ_libs,INTERFACE_LINK_LIBRARIES>>)
  target_link_options(nglib PUBLIC -sALLOW_MEMORY_GROWTH -sENVIRONMENT=web -sWASM_BIGINT -fwasm-exceptions)
else(EMSCRIPTEN)
  target_link_libraries( nglib PRIVATE ${CMAKE_THREAD_LIBS_INIT} ${JPEG_LIBRARIES} ${MKL_LIBRARIES} occ_libs netgen_cgns )
endif(EMSCRIPTEN)

install(TARGETS nglib netgen_cgns ${NG_INSTALL_DIR})
install(FILES nglib.h DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)

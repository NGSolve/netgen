target_sources(nglib PRIVATE nglib.cpp)

if(USE_OCC)
    target_sources(nglib PRIVATE nglib_occ.cpp)
    install(FILES nglib_occ.h DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)
endif(USE_OCC)

target_link_libraries(nglib PRIVATE ${ZLIB_LIBRARIES})
target_link_libraries( nglib PUBLIC ngcore PRIVATE ${CMAKE_THREAD_LIBS_INIT} ${JPEG_LIBRARIES} ${MKL_LIBRARIES} occ_libs netgen_cgns )

install(TARGETS nglib netgen_cgns ${NG_INSTALL_DIR})
install(FILES nglib.h DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)

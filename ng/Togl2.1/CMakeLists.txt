
if(APPLE)
    set(CMAKE_C_COMPILER "/usr/bin/gcc")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x objective-c")
endif(APPLE)

if(WIN32)
    set(TOGL_LIBRARY_TYPE SHARED)
else(WIN32)
    set(TOGL_LIBRARY_TYPE STATIC)
endif(WIN32)

add_library(togl ${TOGL_LIBRARY_TYPE} togl.c toglProcAddr.c toglStubInit.c)
target_link_libraries( togl PUBLIC $<BUILD_INTERFACE:netgen_gui>)

target_compile_definitions(togl PRIVATE -DBUILD_togl=1 -DSTDC_HEADERS=1)

if(WIN32)
  target_compile_definitions(togl PRIVATE -DUNICODE -D_UNICODE -DTOGL_USE_FONTS=0 -DSTDC_HEADER)
else(WIN32)
  target_compile_options(togl PRIVATE -fomit-frame-pointer -Wno-implicit-int)
  target_compile_definitions(togl PRIVATE -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=0 -DHAVE_LIMITS_H=1 -DHAVE_SYS_PARAM_H=1 -DUSE_THREAD_ALLOC=1 -D_REENTRANT=1 -D_THREAD_SAFE=1 -DTCL_THREADS=1 -D_LARGEFILE64_SOURCE=1 -DTCL_WIDE_INT_IS_LONG=1) 

  include_directories(BEFORE "${TCL_INCLUDE_PATH}/tcl-private/generic" "${TCL_INCLUDE_PATH}/tcl-private/unix")
  include_directories(BEFORE "${TK_INCLUDE_PATH}/tk-private/generic" "${TK_INCLUDE_PATH}/tk-private/unix" "${TK_INCLUDE_PATH}/tk-private")
  include_directories(BEFORE "${TCL_INCLUDE_PATH}/tk-private/generic/ttk") 
  include_directories(BEFORE "${TK_INCLUDE_PATH}/../PrivateHeaders") 
  include_directories(BEFORE "${TCL_INCLUDE_PATH}") 
  include_directories(BEFORE "${TK_INCLUDE_PATH}") 

endif(WIN32)

target_include_directories(togl PUBLIC ${OPENGL_INCLUDE_DIR})
target_link_libraries(togl PUBLIC ${OPENGL_LIBRARY})
set_target_properties(togl PROPERTIES POSITION_INDEPENDENT_CODE ON )
install(TARGETS togl DESTINATION ${NG_INSTALL_DIR} COMPONENT netgen)

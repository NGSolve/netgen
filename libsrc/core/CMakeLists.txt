
add_library(ngcore ${NGCORE_LIBRARY_TYPE}
  archive.cpp
  bitarray.cpp
  exception.cpp
  localheap.cpp
  logging.cpp
  flags.cpp
  paje_trace.cpp
  profiler.cpp
  table.cpp
  taskmanager.cpp
  utils.cpp
  version.cpp
  )

target_compile_options(ngcore PUBLIC "${NG_COMPILE_FLAGS}")

# Pybind11 2.3 Issue https://github.com/pybind/pybind11/issues/1604
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  target_compile_options(ngcore PUBLIC -fsized-deallocation -faligned-allocation)
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(ngcore PUBLIC stdc++fs)
endif()

if(USE_PYTHON)
  target_sources(ngcore PRIVATE python_ngcore.cpp)
  target_compile_definitions(ngcore PUBLIC NETGEN_PYTHON NG_PYTHON)
endif(USE_PYTHON)

if(WIN32)
  target_compile_options(ngcore PUBLIC /bigobj /MP /W1 /wd4068)
  get_WIN32_WINNT(ver)
  target_compile_definitions(ngcore PUBLIC _WIN32_WINNT=${ver} WNT WNT_WINDOW NOMINMAX MSVC_EXPRESS _CRT_SECURE_NO_WARNINGS HAVE_STRUCT_TIMESPEC WIN32)
  target_link_options(ngcore PUBLIC /ignore:4273 /ignore:4217 /ignore:4049)
endif(WIN32)

target_compile_definitions(ngcore PRIVATE NGCORE_EXPORTS)
target_include_directories(ngcore INTERFACE $<INSTALL_INTERFACE:${NG_INSTALL_DIR_INCLUDE}> $<INSTALL_INTERFACE:${NG_INSTALL_DIR_INCLUDE}/include>)

if(CHECK_RANGE OR CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  target_compile_definitions(ngcore PUBLIC NETGEN_ENABLE_CHECK_RANGE)
endif(CHECK_RANGE OR CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")

if(TRACE_MEMORY)
  target_compile_definitions(ngcore PUBLIC NETGEN_TRACE_MEMORY)
endif(TRACE_MEMORY)


if(USE_SPDLOG)
  include_directories(${SPDLOG_INCLUDE_DIR})
  install(DIRECTORY ${SPDLOG_INCLUDE_DIR}
    DESTINATION ${NG_INSTALL_DIR_INCLUDE}
    )
  add_dependencies(ngcore project_spdlog)
  target_compile_definitions(ngcore PUBLIC NETGEN_USE_SPDLOG)
  if(DEBUG_LOG)
    target_compile_definitions(ngcore PUBLIC NETGEN_LOG_DEBUG)
  endif(DEBUG_LOG)
endif(USE_SPDLOG)

if(USE_NUMA)
    find_library(NUMA_LIBRARY libnuma.so)
    target_compile_definitions(ngcore PUBLIC USE_NUMA)
    target_link_libraries(ngcore PRIVATE ${NUMA_LIBRARY})
endif(USE_NUMA)

install(TARGETS ngcore DESTINATION ${NG_INSTALL_DIR} COMPONENT netgen)

target_link_libraries(ngcore PUBLIC netgen_mpi PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${CMAKE_THREAD_LIBS_INIT})

install(FILES ngcore.hpp archive.hpp type_traits.hpp version.hpp ngcore_api.hpp logging.hpp memtracer.hpp
  exception.hpp symboltable.hpp paje_trace.hpp utils.hpp profiler.hpp mpi_wrapper.hpp
  array.hpp taskmanager.hpp concurrentqueue.h localheap.hpp python_ngcore.hpp flags.hpp
  xbool.hpp signal.hpp bitarray.hpp table.hpp hashtable.hpp ranges.hpp ngstream.hpp
  simd.hpp simd_avx.hpp simd_avx512.hpp simd_generic.hpp simd_sse.hpp simd_arm64.hpp
  register_archive.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/core COMPONENT netgen_devel)

if(ENABLE_CPP_CORE_GUIDELINES_CHECK)
  set_target_properties(ngcore PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif(ENABLE_CPP_CORE_GUIDELINES_CHECK)

add_dependencies(ngcore ng_generate_version_file)

if(USE_PYTHON)
  pybind11_add_module(pyngcore SHARED python_ngcore_export.cpp)
  target_link_libraries(pyngcore PUBLIC ngcore netgen_python)
  set_target_properties(pyngcore PROPERTIES INSTALL_RPATH "${NG_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
  install(TARGETS pyngcore DESTINATION ${NG_INSTALL_DIR_PYTHON}/pyngcore COMPONENT netgen)
endif(USE_PYTHON)


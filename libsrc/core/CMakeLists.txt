
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
  ng_mpi_wrapper.cpp
  )

string(REPLACE "|" ";" ng_compile_flags_replace_sep "${NG_COMPILE_FLAGS}")
target_compile_options(ngcore PUBLIC ${ng_compile_flags_replace_sep})

if(EMSCRIPTEN)
  set(PYTHON_MODULE_EXTENSION ".so")
  target_link_options(ngcore PUBLIC -sALLOW_MEMORY_GROWTH -sENVIRONMENT=web)
  target_compile_options(ngcore PUBLIC -sNO_DISABLE_EXCEPTION_CATCHING)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Python packages on Linux are compiled with the old ABI,
  # make sure that the same ABI is used in plugins aswell
  try_run(
    ret_val can_compile
    ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/_get_glibcxx_use_cxx11_abi.cpp
    RUN_OUTPUT_VARIABLE use_glibcxx_cxx11_abi
    )
  target_compile_definitions(ngcore PUBLIC -D_GLIBCXX_USE_CXX11_ABI=${use_glibcxx_cxx11_abi})
  if(USE_PYTHON)
    try_run(
      ret_val can_compile
      ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/_get_gxx_abi.cpp
      RUN_OUTPUT_VARIABLE cxx_abi_version
      )
    if(${can_compile} AND (${ret_val} EQUAL 0))
      # Different python modules using pybind11 need to use the same C++ ABI version
      # for compatibility
      message(STATUS "GNU C++ ABI version: ${cxx_abi_version}")
      target_compile_options(ngcore PUBLIC "-fabi-version=${cxx_abi_version}")
    endif()
  endif(USE_PYTHON)
endif()

if(USE_PYTHON)
  target_sources(ngcore PRIVATE python_ngcore.cpp)
  target_compile_definitions(ngcore PUBLIC NETGEN_PYTHON NG_PYTHON PYBIND11_SIMPLE_GIL_MANAGEMENT)
endif(USE_PYTHON)

if(WIN32)
  target_compile_options(ngcore PUBLIC /bigobj $<BUILD_INTERFACE:/MP;/W1;/wd4068>)
  get_WIN32_WINNT(ver)
  target_compile_definitions(ngcore PUBLIC _WIN32_WINNT=${ver} WNT WNT_WINDOW NOMINMAX MSVC_EXPRESS _CRT_SECURE_NO_WARNINGS HAVE_STRUCT_TIMESPEC WIN32)
  target_link_options(ngcore PUBLIC /ignore:4273 /ignore:4217 /ignore:4049)
else(WIN32)
  target_link_libraries(ngcore PUBLIC dl)
endif(WIN32)

target_compile_definitions(ngcore PRIVATE NGCORE_EXPORTS)
target_include_directories(ngcore INTERFACE $<INSTALL_INTERFACE:${NG_INSTALL_DIR_INCLUDE}> $<INSTALL_INTERFACE:${NG_INSTALL_DIR_INCLUDE}/include>)

if(CHECK_RANGE)
  target_compile_definitions(ngcore PUBLIC NETGEN_ENABLE_CHECK_RANGE)
endif(CHECK_RANGE)

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  target_compile_definitions(ngcore PUBLIC _DEBUG NETGEN_ENABLE_CHECK_RANGE)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")

if(TRACE_MEMORY)
  target_compile_definitions(ngcore PUBLIC NETGEN_TRACE_MEMORY)
endif(TRACE_MEMORY)

if(USE_NUMA)
    find_library(NUMA_LIBRARY libnuma.so)
    target_compile_definitions(ngcore PUBLIC USE_NUMA)
    target_link_libraries(ngcore PRIVATE ${NUMA_LIBRARY})
endif(USE_NUMA)

install(TARGETS ngcore DESTINATION ${NG_INSTALL_DIR} COMPONENT netgen)

target_link_libraries(ngcore PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${CMAKE_THREAD_LIBS_INIT})

install(FILES ngcore.hpp archive.hpp type_traits.hpp version.hpp ngcore_api.hpp logging.hpp memtracer.hpp
  exception.hpp symboltable.hpp paje_trace.hpp utils.hpp profiler.hpp mpi_wrapper.hpp
  array.hpp taskmanager.hpp concurrentqueue.h localheap.hpp python_ngcore.hpp flags.hpp
  xbool.hpp signal.hpp bitarray.hpp table.hpp hashtable.hpp ranges.hpp ngstream.hpp
  simd.hpp simd_avx.hpp simd_avx512.hpp simd_generic.hpp simd_sse.hpp simd_arm64.hpp
  register_archive.hpp autodiff.hpp autodiffdiff.hpp
  ng_mpi.hpp ng_mpi_generated_declarations.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/core COMPONENT netgen_devel)

if(ENABLE_CPP_CORE_GUIDELINES_CHECK)
  set_target_properties(ngcore PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif(ENABLE_CPP_CORE_GUIDELINES_CHECK)

add_dependencies(ngcore ng_generate_version_file)

if(USE_PYTHON)
  pybind11_add_module(pyngcore MODULE python_ngcore_export.cpp)
  target_link_libraries(pyngcore PUBLIC ngcore PRIVATE netgen_python)
  set_target_properties(pyngcore PROPERTIES INSTALL_RPATH "${NG_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
  if(EMSCRIPTEN)
    target_compile_definitions(pyngcore PRIVATE NGCORE_EXPORTS)
  endif(EMSCRIPTEN)
  install(TARGETS pyngcore DESTINATION ${NG_INSTALL_DIR_PYTHON}/pyngcore COMPONENT netgen)
endif(USE_PYTHON)

function (build_mpi_variant)
  set(target ng_${ARGV0})
  set(include_dir ${ARGV1})
  message("Building MPI variant: ${ARGV0} ${ARGV1}")
  add_library(${target} SHARED ng_mpi.cpp)
  target_link_libraries(${target} PUBLIC ngcore PRIVATE netgen_python)
  target_compile_definitions(${target} PUBLIC PARALLEL)
  target_include_directories(${target} PRIVATE ${include_dir})
  set_target_properties(${target} PROPERTIES PREFIX "")
  install(TARGETS ${target} RUNTIME DESTINATION ${NG_INSTALL_DIR_BIN} LIBRARY DESTINATION ${NG_INSTALL_DIR_LIB} COMPONENT netgen)
endfunction()

if(USE_MPI)
  target_compile_definitions(ngcore PUBLIC PARALLEL)

  message(STATUS "Found MPI version\n${MPI_C_LIBRARY_VERSION_STRING}")

  if(MPI_C_LIBRARY_VERSION_STRING MATCHES "Microsoft MPI.*")
    set(MICROSOFT_MPI_INCLUDE_DIR ${MPI_C_HEADER_DIR})
  endif()

  if(MPI_C_LIBRARY_VERSION_STRING MATCHES "Open MPI.*")
    set(OPENMPI_INCLUDE_DIR ${MPI_C_INCLUDE_PATH})
  endif()

  if(MPI_C_LIBRARY_VERSION_STRING MATCHES "MPICH.*")
    set(MPICH_INCLUDE_DIR ${MPI_C_INCLUDE_PATH})
  endif()

  if(MPI_C_LIBRARY_VERSION_STRING MATCHES "Intel.*")
    set(INTEL_MPI_INCLUDE_DIR ${MPI_C_INCLUDE_PATH})
  endif()

  if(OPENMPI_INCLUDE_DIR)
    build_mpi_variant(openmpi ${OPENMPI_INCLUDE_DIR})
  endif()
  if(MPICH_INCLUDE_DIR)
    build_mpi_variant(mpich ${MPICH_INCLUDE_DIR})
  endif()
  if(INTEL_MPI_INCLUDE_DIR)
    build_mpi_variant(intel_mpi ${INTEL_MPI_INCLUDE_DIR})
  endif()
  if(MICROSOFT_MPI_INCLUDE_DIR)
    build_mpi_variant(microsoft_mpi ${MICROSOFT_MPI_INCLUDE_DIR})
    target_link_libraries(ng_microsoft_mpi PUBLIC ${MPI_msmpi_LIBRARY})
  endif()

endif(USE_MPI)


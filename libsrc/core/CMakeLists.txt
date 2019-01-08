
add_library(ngcore SHARED archive.cpp logging.cpp paje_trace.cpp utils.cpp profiler.cpp)

target_compile_definitions(ngcore PRIVATE NGCORE_EXPORTS)
if(NOT WIN32)
  target_compile_options(ngcore PRIVATE -fvisibility=hidden)
endif(NOT WIN32)

target_compile_definitions(ngcore PUBLIC $<$<CONFIG:DEBUG>:NETGEN_ENABLE_CHECK_RANGE>)

if(CHECK_RANGE)
  target_compile_definitions(ngcore PUBLIC NETGEN_ENABLE_CHECK_RANGE)
endif(CHECK_RANGE)

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

install(TARGETS ngcore DESTINATION ${NG_INSTALL_DIR} COMPONENT netgen)

if(USE_PYTHON)
  target_compile_definitions(ngcore PUBLIC NETGEN_PYTHON)
  target_include_directories(ngcore PUBLIC ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(ngcore PUBLIC ${PYTHON_LIBRARIES})
endif(USE_PYTHON)

install(FILES ngcore.hpp archive.hpp type_traits.hpp version.hpp ngcore_api.hpp logging.hpp
  exception.hpp symboltable.hpp paje_trace.hpp utils.hpp profiler.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/core COMPONENT netgen_devel)

if(ENABLE_CPP_CORE_GUIDELINES_CHECK)
  set_target_properties(ngcore PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif(ENABLE_CPP_CORE_GUIDELINES_CHECK)

if(USE_PYTHON)
  pybind11_add_module(pyngcore SHARED python_ngcore.cpp)
  target_link_libraries(pyngcore PUBLIC ngcore ${PYTHON_LIBRARIES})
  set_target_properties(pyngcore PROPERTIES INSTALL_RPATH "${NG_RPATH_TOKEN}/${NETGEN_PYTHON_RPATH}")
  install(TARGETS pyngcore DESTINATION ${NG_INSTALL_DIR_PYTHON} COMPONENT netgen)
endif(USE_PYTHON)


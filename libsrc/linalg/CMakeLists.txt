add_library( la INTERFACE )
set(sdir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources( la INTERFACE
  ${sdir}/densemat.cpp ${sdir}/polynomial.cpp ${sdir}/bfgs.cpp ${sdir}/linopt.cpp ${sdir}/linsearch.cpp
)

install(FILES
  densemat.hpp linalg.hpp opti.hpp
  polynomial.hpp vector.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/linalg COMPONENT netgen_devel
)

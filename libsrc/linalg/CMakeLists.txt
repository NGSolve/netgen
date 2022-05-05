target_sources(nglib PRIVATE
  bfgs.cpp
  densemat.cpp
  linopt.cpp
  linsearch.cpp
  polynomial.cpp
)

install(FILES
  densemat.hpp linalg.hpp opti.hpp
  polynomial.hpp vector.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/linalg COMPONENT netgen_devel
)

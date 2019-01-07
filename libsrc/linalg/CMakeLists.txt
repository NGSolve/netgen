add_library( la OBJECT
        densemat.cpp polynomial.cpp bfgs.cpp linopt.cpp linsearch.cpp 
        )

set_target_properties(la PROPERTIES POSITION_INDEPENDENT_CODE ON )

install(FILES
  densemat.hpp linalg.hpp opti.hpp
  polynomial.hpp vector.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/linalg COMPONENT netgen_devel
)

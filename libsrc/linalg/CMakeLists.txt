add_library( la OBJECT
        densemat.cpp polynomial.cpp bfgs.cpp linopt.cpp linsearch.cpp 
        )

set_target_properties(la PROPERTIES POSITION_INDEPENDENT_CODE ON )

if(NOT WIN32)
  target_link_libraries(la ngcore)
endif(NOT WIN32)

install(FILES
  densemat.hpp linalg.hpp opti.hpp
  polynomial.hpp vector.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/linalg COMPONENT netgen_devel
)

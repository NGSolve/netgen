if(USE_PYTHON)
    add_test(NAME pytest COMMAND ${Python3_EXECUTABLE} -m pytest WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_custom_target(pytest ${Python3_EXECUTABLE} -m pytest WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties ( pytest PROPERTIES TIMEOUT 1800 )
    if(USE_MPI AND MPIEXEC_EXECUTABLE)
      # needs mpi4py in the test python, the tests skip themselves otherwise
      add_test(NAME pytest-mpi
               COMMAND ${MPIEXEC_EXECUTABLE} --allow-run-as-root --oversubscribe -np 3 ${Python3_EXECUTABLE} -m pytest
                       test_mpi_distribute.py test_mpi4py.py test_edgedescriptor_mpi.py
               WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      set_tests_properties ( pytest-mpi PROPERTIES TIMEOUT 600 )
    endif()
endif(USE_PYTHON)

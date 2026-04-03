add_subdirectory(catch)
add_subdirectory(pytest)

# this code goes here, because tests is the last add_subdirectory (otherwise it gets executed too early)
if(APPLE AND BUILD_FOR_CONDA)
    install(CODE "execute_process(COMMAND sh -c \"codesign --force -s - netgen/*.so netgen/*.dylib pyngcore/*.so\" WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX})")
endif(APPLE AND BUILD_FOR_CONDA)

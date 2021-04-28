configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
        meshing.py csg.py geom2d.py stl.py gui.py NgOCC.py occ.py read_gmsh.py
        DESTINATION ${NG_INSTALL_DIR_PYTHON}/${NG_INSTALL_SUFFIX}
        COMPONENT netgen
        )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import pybind11_stubgen; print(pybind11_stubgen.__file__)" OUTPUT_VARIABLE stubgen_path RESULT_VARIABLE pybind11_stubgen)
if(pybind11_stubgen AND NOT ${pybind11_stubgen} EQUAL 0)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
else()
  message("-- Found pybind11-stubgen: ${stubgen_path}")
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11_stubgen --no-setup-py netgen)")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/netgen-stubs/ DESTINATION ${NG_INSTALL_DIR_PYTHON}/netgen/ COMPONENT netgen)
endif()
endif(BUILD_STUB_FILES)

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
        meshing.py csg.py geom2d.py stl.py gui.py NgOCC.py occ.py read_gmsh.py
        DESTINATION ${NG_INSTALL_DIR_PYTHON}/${NG_INSTALL_SUFFIX}
        COMPONENT netgen
        )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
find_program(PYBIND11_STUBS NAMES pybind11-stubgen)
if(PYBIND11_STUBS)
  message("-- Found pybind11-stubgen: ${PYBIND11_STUBS}")
  install(CODE "execute_process(COMMAND ${PYBIND11_STUBS} --no-setup-py netgen)")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/netgen-stubs/ DESTINATION ${NG_INSTALL_DIR_PYTHON}/netgen/)
else(PYBIND11_STUBS)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
endif(PYBIND11_STUBS)
endif(BUILD_STUB_FILES)

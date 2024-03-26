get_target_property(ngcore_compile_definitions ngcore INTERFACE_COMPILE_DEFINITIONS)
get_property(have_options TARGET ngcore PROPERTY INTERFACE_COMPILE_OPTIONS SET)
if(have_options)
    get_target_property(ngcore_compile_options ngcore INTERFACE_COMPILE_OPTIONS)
endif(have_options)

configure_file(config/config_template.py ${CMAKE_CURRENT_BINARY_DIR}/config.py @ONLY)
configure_file(version_template.py ${CMAKE_CURRENT_BINARY_DIR}/version.py @ONLY)

install(FILES
        config/__init__.py
        config/__main__.py
        ${CMAKE_CURRENT_BINARY_DIR}/config.py
        DESTINATION ${NG_INSTALL_DIR_PYTHON}/${NG_INSTALL_SUFFIX}/config/
        COMPONENT netgen
        )

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/version.py
        __main__.py __init__.py
        meshing.py csg.py geom2d.py stl.py gui.py NgOCC.py occ.py
        read_gmsh.py read_meshio.py
        webgui.py
        DESTINATION ${NG_INSTALL_DIR_PYTHON}/${NG_INSTALL_SUFFIX}
        COMPONENT netgen
        )

install(FILES
        pyngcore/__init__.py
        DESTINATION ${NG_INSTALL_DIR_PYTHON}/pyngcore
        COMPONENT netgen
        )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
execute_process(COMMAND ${Python3_EXECUTABLE} -c "from importlib.metadata import version; print(version('pybind11-stubgen'))" OUTPUT_VARIABLE stubgen_version RESULT_VARIABLE pybind11_stubgen)
if(pybind11_stubgen AND NOT ${pybind11_stubgen} EQUAL 0)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
else()
  if(stubgen_version LESS "1.0")
    message(WARNING "pybind11-stubgen version is too old, if you want to create stub files for better autocompletion support upgrade it with pip.")
  else()
  message("-- Found pybind11-stubgen version: ${stubgen_version}")
  install(CODE "\
  set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/${NG_INSTALL_DIR_PYTHON})\n \
  execute_process(COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen --ignore-all-errors netgen)\n \
  execute_process(COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen --ignore-all-errors pyngcore)\n \
  ")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/ DESTINATION ${NG_INSTALL_DIR_PYTHON} COMPONENT netgen)
  endif()
endif()
endif(BUILD_STUB_FILES)

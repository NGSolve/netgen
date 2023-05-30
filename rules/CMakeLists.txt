# this file is included from the parent directory (otherwise generated source files are not recognized properly by cmake)

# generate .cpp files containing the string of the .rls meshing rule files
add_executable(makerls rules/makerlsfile.cpp)

set(rules
    hexrules
    prismrules2
    pyramidrules
    pyramidrules2
    quadrules
    tetrules
    triarules
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rules/)

foreach(rule ${rules})
    list(APPEND rules_sources ${CMAKE_CURRENT_BINARY_DIR}/rules/rule_${rule}.cpp)
    set(rule_file ${CMAKE_CURRENT_SOURCE_DIR}/rules/${rule}.rls)
    set(rule_cpp ${CMAKE_CURRENT_BINARY_DIR}/rules/rule_${rule}.cpp)

    add_custom_command(OUTPUT ${rule_cpp}
      COMMAND makerls ${rule_file} ${rule_cpp} ${rule}
      DEPENDS makerls ${rule_file}
    )
endforeach()

target_sources(nglib PRIVATE ${rules_sources})

set(UNIT_TESTS
    test_figure
    test_command_line_interface
    test_null_graphics
    test_default_graphics_regression
    test_default_graphics_robustness
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    if("${TEST}" STREQUAL "test_default_graphics_regression")
        target_compile_definitions(${TEST} PRIVATE
            ARIADNE_DRAWING_GOLDEN_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/golden\")
    endif()
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "drawing")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

# C++ graphics regression tests create and inspect their own PNGs.
# No Python interpreter, PNG comparison script or fixture ordering is needed.
set_tests_properties(test_default_graphics_regression test_default_graphics_robustness PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

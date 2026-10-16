set(UNIT_TESTS
    test_figure
    test_command_line_interface
    test_null_graphics
    test_backend_comparison
    test_default_graphics_robustness
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "drawing")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

# C++ graphics regression tests create and inspect their own PNGs.
# No Python interpreter, PNG comparison script or fixture ordering is needed.
set_tests_properties(test_backend_comparison test_default_graphics_robustness PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

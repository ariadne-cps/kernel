set(UNIT_TESTS
    test_figure
    test_command_line_interface
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "io")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

set(UNIT_TESTS
    test_expression
    test_expression_set
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "symbolic")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

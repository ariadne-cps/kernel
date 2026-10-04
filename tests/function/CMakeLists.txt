set(UNIT_TESTS
    test_function
    test_measurable_function
    test_multifunction
    test_formula
    test_procedure
    test_taylor_function
    test_taylor_model
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "function")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

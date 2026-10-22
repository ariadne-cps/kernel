set(UNIT_TESTS
    test_affine
    test_affine_model
    test_function_operations
    test_function_dispatch
    test_function_interface
    test_function_names
    test_function_patch_operations
    test_domain_scaling
    test_taylor_series
    test_function
    test_formula
    test_measurable_function
    test_multifunction
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

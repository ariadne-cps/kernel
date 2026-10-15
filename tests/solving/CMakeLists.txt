set(UNIT_TESTS
    test_constraint_solver
    test_linear_interior_point
    test_nonlinear_programming
    test_solvers
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "solving")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

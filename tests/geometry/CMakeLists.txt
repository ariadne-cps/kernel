set(UNIT_TESTS
    test_affine_sets
    test_binary_tree
    test_box
    test_point_curve
    test_constrained_image_set
    test_grid_paving
    test_paving
    test_measurable_set
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "geometry")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

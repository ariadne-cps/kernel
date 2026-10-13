set(UNIT_TESTS
    test_set
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME labelling_${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST labelling_${TEST} PROPERTY LABELS "labelling")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

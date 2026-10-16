set(UNIT_TESTS
    test_figure
    test_command_line_interface
    test_null_graphics
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "drawing")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

# Deterministic scene exported by the two implementations.
add_executable(test_backend_comparison test_backend_comparison.cpp)
target_link_libraries(test_backend_comparison ariadne-kernel kernel)
add_test(NAME drawing_backend_generate COMMAND $<TARGET_FILE:test_backend_comparison>)
set_tests_properties(drawing_backend_generate PROPERTIES LABELS "drawing"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
add_dependencies(tests-cpp test_backend_comparison)

if(ARIADNE_WITH_CAIRO)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME drawing_backend_compare
            COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/compare_backend_png.py"
                "comparison-cairo.png" "comparison-embedded.png")
        set_tests_properties(drawing_backend_compare PROPERTIES
            LABELS "drawing"
            DEPENDS drawing_backend_generate
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
endif()

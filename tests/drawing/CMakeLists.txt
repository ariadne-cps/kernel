set(UNIT_TESTS
    test_figure
    test_command_line_interface
    test_null_graphics
    test_backend_comparison
)
foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-kernel kernel)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "drawing")
endforeach()
add_dependencies(tests-cpp ${UNIT_TESTS})

# Compare matching raster fixtures rendered by both backends.
if(ARIADNE_WITH_CAIRO)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME compare_drawing_backends
            COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/compare_backend_png.py"
                "comparison-cairo.png" "comparison-embedded.png")
        set_tests_properties(compare_drawing_backends PROPERTIES
            DEPENDS test_backend_comparison
            LABELS "drawing")
    endif()
endif()

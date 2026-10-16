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

# The label legibility regression runs with either graphics configuration.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME test_embedded_labels
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_embedded_labels.py"
            "comparison-embedded.png")
    set_tests_properties(test_embedded_labels PROPERTIES
        DEPENDS test_backend_comparison LABELS "drawing")
endif()

# Compare matching raster fixtures rendered by both backends.
if(ARIADNE_WITH_CAIRO)
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

# Each scene is compared independently and emits a diagnostic image.
if(ARIADNE_WITH_CAIRO AND Python3_Interpreter_FOUND)
    foreach(SCENE RANGE 0 3)
        add_test(NAME compare_drawing_scene_${SCENE}
            COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/compare_backend_png.py"
                "comparison-scene${SCENE}-cairo.png"
                "comparison-scene${SCENE}-embedded.png"
                --diff "comparison-scene${SCENE}-diff.ppm")
        set_tests_properties(compare_drawing_scene_${SCENE} PROPERTIES
            DEPENDS test_backend_comparison LABELS "drawing")
    endforeach()
endif()

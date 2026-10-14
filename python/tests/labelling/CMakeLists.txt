add_custom_target(
    test_python_labelling
    COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/test_labelling.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_test(
    NAME test_python_labelling
    COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/test_labelling.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
set_property(TEST test_python_labelling PROPERTY
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:pyariadne>"
)
set_property(TEST test_python_labelling PROPERTY LABELS "python;labelling")
add_dependencies(test_python_labelling pyariadne)

execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pytest --version
    RESULT_VARIABLE PYTEST_result
    OUTPUT_QUIET
    ERROR_QUIET
)

if(PYTEST_result EQUAL 0)
    set(PYTHON_TESTS ON)
else()
    if(PYTHON_ENFORCE_TESTS)
        message(FATAL_ERROR "pytest is not available for ${Python3_EXECUTABLE}.")
    else()
        message(WARNING "pytest is not available; Python tests will not be performed.")
    endif()
endif()

if(PYTHON_TESTS)
    foreach(TESTNAME function calculus geometry solvers symbolic)
        add_custom_target(
            test_python_${TESTNAME}
            COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/test_${TESTNAME}.py"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        )
        add_test(
            NAME test_python_${TESTNAME}
            COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/test_${TESTNAME}.py"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        )
        set_property(TEST test_python_${TESTNAME} PROPERTY
            ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:pyariadne>"
        )
        set_property(TEST test_python_${TESTNAME} PROPERTY LABELS "python")
        add_dependencies(test_python_${TESTNAME} pyariadne)
    endforeach()
    add_dependencies(tests pyariadne)
endif()

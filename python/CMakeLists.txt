set(Python3_FIND_VIRTUALENV FIRST)
set(Python3_FIND_STRATEGY LOCATION)
set(PYBIND11_FINDPYTHON ON)

find_package(Python3 3.11 COMPONENTS Interpreter Development.Module)

option(PYTHON_ENFORCE_TESTS "Enforce Python tests" OFF)

if(Python3_FOUND)
    message(STATUS "Python executable: ${Python3_EXECUTABLE}")
    message(STATUS "Python version:    ${Python3_VERSION}")
    message(STATUS "Python includes:   ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python SOABI:      ${Python3_SOABI}")

    if(NOT TARGET pyariadne-algebra)
        message(FATAL_ERROR "The Algebra Python binding target is required to build pyariadne-kernel.")
    endif()

    if(NOT TARGET pybind11::module)
        get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)
        add_subdirectory(
            "${PYTHON_COMMON_SOURCE_DIR}/pybind11"
            "${CMAKE_CURRENT_BINARY_DIR}/pybind11"
        )
    endif()

    add_subdirectory(bindings)

    if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
        pybind11_add_module(pyariadne
            $<TARGET_OBJECTS:pyariadne-kernel-module-obj>
            $<TARGET_OBJECTS:KERNEL_FUNCTION_SRC>
            $<TARGET_OBJECTS:KERNEL_GEOMETRY_SRC>
            $<TARGET_OBJECTS:KERNEL_IO_SRC>
            $<TARGET_OBJECTS:KERNEL_SYMBOLIC_SRC>
            $<TARGET_OBJECTS:KERNEL_SOLVING_SRC>
            $<TARGET_OBJECTS:ALGEBRA_SRC>
            $<TARGET_OBJECTS:INTERVAL_SRC>
            $<TARGET_OBJECTS:NUMERIC_SRC>
            $<TARGET_OBJECTS:FOUNDATION_SRC>
            $<TARGET_OBJECTS:THREADING_SRC>
            $<TARGET_OBJECTS:LOGGING_SRC>
            $<TARGET_OBJECTS:UTILITY_SRC>
        )
        target_link_libraries(pyariadne PRIVATE
            pyariadne-kernel
            kernel
        )

        set(PYTHON_BINDINGS_INSTALL_DIR "${Python3_SITEARCH}"
            CACHE STRING "Install destination for Python targets")
        install(TARGETS pyariadne DESTINATION "${PYTHON_BINDINGS_INSTALL_DIR}")

        add_subdirectory(tests)
        add_dependencies(bindings pyariadne)
    endif()
else()
    if(PYTHON_ENFORCE_TESTS)
        message(FATAL_ERROR "Python >= 3.11 with Interpreter and Development.Module components is required when PYTHON_ENFORCE_TESTS=ON.")
    else()
        message(STATUS "Compatible Python development environment not found: Python bindings will not be built.")
    endif()
endif()

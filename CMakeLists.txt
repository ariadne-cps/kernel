cmake_minimum_required(VERSION 3.15)

project(Kernel VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/ss/configuration/cmake/ProjectOptions.cmake")
    message(FATAL_ERROR
        "Git submodules are not initialized. "
        "Run 'git submodule update --init --recursive' before configuring the project.")
endif()

include("${PROJECT_SOURCE_DIR}/ss/configuration/cmake/ProjectOptions.cmake")

require_initialized_submodules()
setup_project()


if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        REFERENCES
            "${PROJECT_SOURCE_DIR}|ss/configuration"
            "${PROJECT_SOURCE_DIR}/ss/algebra|ss/configuration"
            "${PROJECT_SOURCE_DIR}/ss/bdd|ss/configuration"
    )
    require_same_dependency_commit(Utility
        REFERENCES
            "${PROJECT_SOURCE_DIR}/ss/algebra/ss/interval/ss/numeric/ss/foundation|ss/utility"
            "${PROJECT_SOURCE_DIR}/ss/bdd/ss/threading/ss/logging|ss/utility"
    )
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(Cairo REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_subdirectory(ss)


add_library(kernel INTERFACE)
target_include_directories(kernel INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-kernel>
)
target_link_libraries(kernel INTERFACE algebra bdd Cairo::Cairo)


add_subdirectory(src)
foreach(target
    KERNEL_FUNCTION_SRC
    KERNEL_GEOMETRY_SRC
    KERNEL_IO_SRC
    KERNEL_LABELLED_SRC
    KERNEL_SOLVING_SRC
    KERNEL_SYMBOLIC_SRC
)
    target_link_libraries(${target} PRIVATE kernel)
endforeach()

ariadne_register_public_headers(
    TARGET kernel
    DIRECTORIES
        "${PROJECT_SOURCE_DIR}/include/function"
        "${PROJECT_SOURCE_DIR}/include/geometry"
        "${PROJECT_SOURCE_DIR}/include/io"
        "${PROJECT_SOURCE_DIR}/include/labelled"
        "${PROJECT_SOURCE_DIR}/include/solving"
        "${PROJECT_SOURCE_DIR}/include/symbolic"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-kernel.hpp"
)

setup_project_library(
    TARGET ariadne-kernel
    OBJECTS
        KERNEL_FUNCTION_SRC
        KERNEL_GEOMETRY_SRC
        KERNEL_IO_SRC
        KERNEL_LABELLED_SRC
        KERNEL_SYMBOLIC_SRC
        KERNEL_SOLVING_SRC
    LINK_LIBRARIES kernel
)

setup_standalone_project_tests(
    TARGET ariadne-kernel
    TEST_DIRECTORY tests
    TEST_TARGET tests
    MODULE_REPORTS
        function
        geometry
        io
        labelled
        solving
        symbolic
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET kernel
        DESTINATION include/ariadne-kernel
    )
endif()

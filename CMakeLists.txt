cmake_minimum_required(VERSION 3.15)

project(Kernel VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT TARGET ariadne-configuration)
    find_package(Git REQUIRED)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                    submodule update --init -- ss/configuration
                    RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize configuration submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/configuration/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        CHAINS
            "configuration"
            "algebra|configuration"
            "bdd|configuration"
    )
    require_same_dependency_commit(Utility
        CHAINS
            "algebra|interval|numeric|foundation|utility"
            "bdd|threading|logging|utility"
    )
endif()

setup_project()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
option(ARIADNE_WITH_CAIRO "Build optional Cairo reference backend" ON)
if(ARIADNE_WITH_CAIRO)
    find_package(Cairo QUIET)
    if(NOT Cairo_FOUND AND NOT CAIRO_FOUND)
        message(STATUS "Cairo unavailable: building Embedded backend only")
        set(ARIADNE_WITH_CAIRO OFF)
    endif()
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")

add_subdirectory(ss)

add_library(kernel INTERFACE)
target_include_directories(kernel INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-kernel>
)
target_link_libraries(kernel INTERFACE algebra bdd)
if(ARIADNE_WITH_CAIRO)
    target_link_libraries(kernel INTERFACE Cairo::Cairo)
    target_compile_definitions(kernel INTERFACE ARIADNE_HAS_CAIRO=1)
endif()

add_subdirectory(src)
foreach(target
    KERNEL_FUNCTION_SRC
    KERNEL_GEOMETRY_SRC
    KERNEL_DRAWING_SRC
    KERNEL_SOLVING_SRC
)
    target_link_libraries(${target} PRIVATE kernel)
endforeach()

ariadne_register_public_headers(
    TARGET kernel
    DIRECTORIES
        "${PROJECT_SOURCE_DIR}/include/function"
        "${PROJECT_SOURCE_DIR}/include/geometry"
        "${PROJECT_SOURCE_DIR}/include/drawing"
        "${PROJECT_SOURCE_DIR}/include/solving"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-kernel.hpp"
)

setup_project_library(
    TARGET ariadne-kernel
    OBJECTS
        KERNEL_FUNCTION_SRC
        KERNEL_GEOMETRY_SRC
        KERNEL_DRAWING_SRC
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
        drawing
        solving
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

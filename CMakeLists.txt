cmake_minimum_required(VERSION 3.15)

project(Kernel VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/submodules/configuration/cmake/ProjectOptions.cmake")
    message(FATAL_ERROR
        "Git submodules are not initialized. "
        "Run 'git submodule update --init --recursive' before configuring the project.")
endif()

include("${PROJECT_SOURCE_DIR}/submodules/configuration/cmake/ProjectOptions.cmake")

require_initialized_submodules()
setup_project()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        REFERENCES
            "${PROJECT_SOURCE_DIR}|submodules/configuration"
            "${PROJECT_SOURCE_DIR}/submodules/algebra|submodules/configuration"
            "${PROJECT_SOURCE_DIR}/submodules/threading|submodules/configuration"
    )
    require_same_dependency_commit(Utility
        REFERENCES
            "${PROJECT_SOURCE_DIR}/submodules/algebra/submodules/interval/submodules/numeric/submodules/foundation|submodules/utility"
            "${PROJECT_SOURCE_DIR}/submodules/threading/submodules/logging|submodules/utility"
    )
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(Threads REQUIRED)
find_package(Cairo REQUIRED)

include_directories(SYSTEM
    "${PROJECT_SOURCE_DIR}/include"
)

add_subdirectory(submodules)

if(MSVC)
    add_compile_options(/bigobj)
endif()

add_library(kernel INTERFACE)
target_include_directories(kernel INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-kernel>
)
target_link_libraries(kernel INTERFACE algebra threading Cairo::Cairo)

# Numeric is embedded into the monolithic kernel/Python binaries in this project.
# On Windows this disables dllimport/dllexport annotations inside those binaries.
target_compile_definitions(numeric INTERFACE ARIADNE_NUMERIC_EMBEDDED)


add_subdirectory(src)
foreach(target
    KERNEL_FUNCTION_SRC
    KERNEL_GEOMETRY_SRC
    KERNEL_IO_SRC
    KERNEL_SYMBOLIC_SRC
    KERNEL_SOLVING_SRC
)
    target_link_libraries(${target} PRIVATE kernel)
endforeach()

ariadne_register_public_headers(
    TARGET kernel
    DIRECTORIES
        "${PROJECT_SOURCE_DIR}/include/function"
        "${PROJECT_SOURCE_DIR}/include/geometry"
        "${PROJECT_SOURCE_DIR}/include/io"
        "${PROJECT_SOURCE_DIR}/include/symbolic"
        "${PROJECT_SOURCE_DIR}/include/solving"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-kernel.hpp"
)

if(WIN32)
    # CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS turns every global symbol in this
    # monolithic target into an import-library member. The aggregate exceeds
    # the 65535-member COFF library limit, so use a static archive on Windows.
    add_library(ariadne-kernel STATIC)
    foreach(_object_target
        KERNEL_FUNCTION_SRC
        KERNEL_GEOMETRY_SRC
        KERNEL_IO_SRC
        KERNEL_SYMBOLIC_SRC
        KERNEL_SOLVING_SRC
        ALGEBRA_SRC
        INTERVAL_SRC
        NUMERIC_SRC
        FOUNDATION_SRC
        THREADING_SRC
        LOGGING_SRC
        UTILITY_SRC
    )
        target_sources(ariadne-kernel PRIVATE $<TARGET_OBJECTS:${_object_target}>)
    endforeach()
    target_link_libraries(ariadne-kernel PUBLIC kernel)
    set_property(GLOBAL PROPERTY ARIADNE_PROJECT_LIBRARY_TARGET ariadne-kernel)
    install(TARGETS ariadne-kernel ARCHIVE DESTINATION lib)
else()
    setup_project_library(
        TARGET ariadne-kernel
        OBJECTS
            KERNEL_FUNCTION_SRC
            KERNEL_GEOMETRY_SRC
            KERNEL_IO_SRC
            KERNEL_SYMBOLIC_SRC
            KERNEL_SOLVING_SRC
            ALGEBRA_SRC
            INTERVAL_SRC
            NUMERIC_SRC
            FOUNDATION_SRC
            THREADING_SRC
            LOGGING_SRC
            UTILITY_SRC
        LINK_LIBRARIES kernel
    )
endif()

setup_standalone_project_tests(
    TARGET ariadne-kernel
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

add_subdirectory(tutorials)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET kernel
        DESTINATION include/ariadne-kernel
    )
endif()

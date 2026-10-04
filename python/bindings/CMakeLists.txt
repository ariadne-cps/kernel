include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-kernel-submodule-obj OBJECT
    utilities.cpp
    function_submodule.cpp
    function_bindings.cpp
    geometry_submodule.cpp
    geometry_bindings.cpp
    io_submodule.cpp
    io_bindings.cpp
    solving_submodule.cpp
    solving_bindings.cpp
    symbolic_submodule.cpp
)
target_include_directories(pyariadne-kernel-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-algebra,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-kernel-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    kernel
)
set_target_properties(pyariadne-kernel-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(pyariadne-kernel INTERFACE)
target_sources(pyariadne-kernel INTERFACE
    $<TARGET_OBJECTS:pyariadne-kernel-submodule-obj>
)
target_include_directories(pyariadne-kernel INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-algebra,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-kernel INTERFACE
    pyariadne-algebra
)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    add_library(pyariadne-kernel-module-obj OBJECT kernel_module.cpp)
    target_include_directories(pyariadne-kernel-module-obj PRIVATE
        "${PYTHON_COMMON_SOURCE_DIR}/include"
    )
    target_link_libraries(pyariadne-kernel-module-obj PRIVATE
        pyariadne-python-common
        pybind11::module
    )
    set_target_properties(pyariadne-kernel-module-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

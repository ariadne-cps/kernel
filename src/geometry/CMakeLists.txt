add_library(KERNEL_GEOMETRY_SRC OBJECT
    point.cpp
    box.cpp
    curve.cpp
    set.cpp
    function_set.cpp
    binary_tree.cpp
    grid.cpp
    grid_paving.cpp
    paver.cpp
    affine_set.cpp
    geometry.cpp
)

target_include_directories(KERNEL_GEOMETRY_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/geometry"
)

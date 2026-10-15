add_library(KERNEL_SOLVING_SRC OBJECT
    solver.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    interior_point_method.cpp
    nonlinear_programming.cpp
)

target_include_directories(KERNEL_SOLVING_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/solving"
)

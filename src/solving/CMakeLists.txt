add_library(KERNEL_SOLVING_SRC OBJECT
    equation_solver.cpp
    interval_newton_method.cpp
    krawczyk_method.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    interior_point_method.cpp
    nonlinear_programming.cpp
)

target_include_directories(KERNEL_SOLVING_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/solving"
)

add_library(KERNEL_SOLVING_SRC OBJECT
    solver.cpp
    smt_solver.cpp
    smt_boolean.cpp
    smt_theory.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    linear_programming.cpp
    nonlinear_programming.cpp
)

target_include_directories(KERNEL_SOLVING_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/solving"
)

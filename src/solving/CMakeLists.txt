add_library(KERNEL_SOLVING_SRC OBJECT
    solver.cpp
    smt_solver.cpp
    smt_boolean.cpp
    smt_theory.cpp
    integrator.cpp
    inclusion_integrator.cpp
    bounder.cpp
    constraint_solver.cpp
    simplex_algorithm.cpp
    linear_programming.cpp
    nonlinear_programming.cpp
    runge_kutta_integrator.cpp
)

target_include_directories(KERNEL_SOLVING_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/solving"
)

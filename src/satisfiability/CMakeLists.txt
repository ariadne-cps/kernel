add_library(KERNEL_SATISFIABILITY_SRC OBJECT
    smt_solver.cpp
    smt_boolean.cpp
    smt_theory.cpp
)

target_include_directories(KERNEL_SATISFIABILITY_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/satisfiability"
)

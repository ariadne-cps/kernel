add_library(KERNEL_SYMBOLIC_SRC OBJECT
    space.cpp
    expression.cpp
    expression_set.cpp
    formula.cpp
)

target_include_directories(KERNEL_SYMBOLIC_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/symbolic"
)

add_library(KERNEL_LABELLED_SRC OBJECT
    point.cpp
    curve.cpp
    figure.cpp
)

target_include_directories(KERNEL_LABELLED_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/labelled"
)

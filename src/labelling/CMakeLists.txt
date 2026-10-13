add_library(KERNEL_LABELLING_SRC OBJECT
    point.cpp
    curve.cpp
    figure.cpp
)

target_include_directories(KERNEL_LABELLING_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/labelled"
)

add_library(KERNEL_DRAWING_SRC OBJECT
    figure.cpp
    default_backend.cpp
    $<$<BOOL:${ARIADNE_WITH_CAIRO}>:cairo.cpp>
    graphics_manager.cpp
    drawer.cpp
    command_line_interface.cpp
)

target_include_directories(KERNEL_DRAWING_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/drawing"
)

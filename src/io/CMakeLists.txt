add_library(KERNEL_IO_SRC OBJECT
    figure.cpp
    cairo.cpp
    graphics_manager.cpp
    drawer.cpp
    command_line_interface.cpp
)

target_include_directories(KERNEL_IO_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/io"
)

cmake_minimum_required(VERSION 3.15)
project(AriadneKernelIOTutorial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(AriadneKernel REQUIRED)

add_executable(tutorial_io tutorial_io.cpp)
target_link_libraries(tutorial_io PRIVATE AriadneKernel::ariadne-kernel)

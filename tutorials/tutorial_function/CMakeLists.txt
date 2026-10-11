cmake_minimum_required(VERSION 3.15)
project(TutorialFunction VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(AriadneKernel REQUIRED)

add_executable(tutorial_function tutorial_function.cpp)
target_link_libraries(tutorial_function PRIVATE AriadneKernel::ariadne-kernel)

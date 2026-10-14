cmake_minimum_required(VERSION 3.15)
project(AriadneKernelSatisfiabilityTutorial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(AriadneKernel REQUIRED)

add_executable(tutorial_satisfiability tutorial_satisfiability.cpp)
target_link_libraries(tutorial_satisfiability PRIVATE AriadneKernel::ariadne-kernel)

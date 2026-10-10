include(CTest)

add_custom_target(tests)
add_custom_target(tests-cpp)
add_dependencies(tests tests-cpp)

add_subdirectory(function)
add_subdirectory(geometry)
add_subdirectory(io)
add_subdirectory(solving)
add_subdirectory(symbolic)

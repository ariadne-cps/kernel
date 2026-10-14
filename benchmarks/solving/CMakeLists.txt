set(SOLVING_BENCHMARKS
    benchmark_constraint_propagation
    benchmark_constraint_shaving
)
foreach(BENCHMARK ${SOLVING_BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne-kernel kernel)
endforeach()
add_dependencies(benchmarks ${SOLVING_BENCHMARKS})

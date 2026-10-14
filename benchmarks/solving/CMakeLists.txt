set(SOLVING_BENCHMARKS
    benchmark_constraint_propagation
    benchmark_constraint_shaving
    benchmark_smt_interval_newton
)
foreach(BENCHMARK ${SOLVING_BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne-kernel kernel)
endforeach()
add_dependencies(benchmarks ${SOLVING_BENCHMARKS})

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/data/smt_barr3_full64.bin")
    add_executable(benchmark_smt_barr3_verification benchmark_smt_barr3_verification.cpp)
    target_link_libraries(benchmark_smt_barr3_verification PRIVATE ariadne-kernel kernel)
    target_compile_definitions(benchmark_smt_barr3_verification PRIVATE
        ARIADNE_SMT_BARR3_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/smt_barr3_full64.bin")
    add_dependencies(benchmarks benchmark_smt_barr3_verification)
endif()

add_executable(fock_bench bench_library.cpp)
target_link_libraries(fock_bench fock benchmark::benchmark)

target_include_directories(fock_bench PRIVATE ${benchmark_SOURCE_DIR}/include)
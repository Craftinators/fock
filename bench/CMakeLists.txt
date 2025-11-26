cmake_minimum_required(VERSION 4.0)

include(FetchContent)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.9.4
)
FetchContent_MakeAvailable(googlebenchmark)

macro(add_fock_benchmark TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} PRIVATE fock benchmark::benchmark_main)
    target_compile_features(${TARGET_NAME} PRIVATE cxx_std_20)
endmacro()

add_fock_benchmark(benchmark_build_hamiltonian benchmark_build_hamiltonian.cpp)
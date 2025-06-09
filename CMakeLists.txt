cmake_minimum_required(VERSION 3.31)
project(fock LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create static library from source
add_library(fock STATIC src/library.cpp)

# Publicly expose the include directory
target_include_directories(fock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Use benchmarking
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)

include(FetchContent)

FetchContent_Declare(
        benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG        eddb0241389718a23a42db6af5f0164b6e0139af # Release v1.9.4
)

FetchContent_MakeAvailable(benchmark)

# Make benchmark headers available globally
target_include_directories(benchmark PUBLIC
        $<BUILD_INTERFACE:${benchmark_SOURCE_DIR}/include>
)

add_subdirectory(benchmarks)
add_subdirectory(examples)
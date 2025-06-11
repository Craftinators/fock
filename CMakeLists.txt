cmake_minimum_required(VERSION 3.31)
project(fock LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# mini-gmp target
add_library(mini-gmp STATIC external/mini-gmp/src/mini-gmp.c)

target_include_directories(mini-gmp PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/mini-gmp/include>
)

add_library(fock STATIC src/library.cpp)

target_include_directories(fock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

include(FetchContent)

# Disable google tests
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)

# Disable Eigen tests and examples to speed up build
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

# Eigen
FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0  # Latest stable version
)

# Google benchmark
FetchContent_Declare(
        benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG        eddb0241389718a23a42db6af5f0164b6e0139af # Release v1.9.4
)

FetchContent_MakeAvailable(Eigen3 benchmark)

target_link_libraries(fock
        PRIVATE mini-gmp
        PUBLIC Eigen3::Eigen
)

# Make benchmark headers available globally
target_include_directories(benchmark PUBLIC
        $<BUILD_INTERFACE:${benchmark_SOURCE_DIR}/include>
)

add_subdirectory(benchmarks)
add_subdirectory(examples)
cmake_minimum_required(VERSION 4.0)
project(fock VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(FOCK_VERSION ${CMAKE_PROJECT_VERSION})

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type." FORCE)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(FOCK_DEBUG_BUILD)
endif()

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    set(IS_ROOT_PROJECT ON)
else ()
    set(IS_ROOT_PROJECT OFF)
endif ()

option(FOCK_BUILD_EXAMPLES "Build examples" ${IS_ROOT_PROJECT})
option(FOCK_BUILD_BENCH "Build benchmarks" ${IS_ROOT_PROJECT})

include(FetchContent)

FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.0
)
FetchContent_MakeAvailable(Eigen3)

add_library(fock STATIC src/library.cpp)

target_include_directories(fock
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(fock PUBLIC cxx_std_20)
set_target_properties(fock PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME fock)

target_link_libraries(fock PUBLIC Eigen3::Eigen)

if (FOCK_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (FOCK_BUILD_BENCH)
    add_subdirectory(bench)
endif ()
cmake_minimum_required(VERSION 4.0)
<<<<<<< HEAD
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

find_package(Eigen3 REQUIRED)
find_package(Spectra REQUIRED)
find_package(Boost REQUIRED)

add_library(fock STATIC src/bitmask.cpp src/matrix.cpp)

target_include_directories(fock
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(fock PUBLIC cxx_std_20)
set_target_properties(fock PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME fock)

target_link_libraries(fock PUBLIC Eigen3::Eigen Boost::headers)

if (FOCK_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (FOCK_BUILD_BENCH)
    add_subdirectory(bench)
endif ()
=======
project(fock VERSION 1.0.0
    DESCRIPTION "A minimal quantum computational library for hard-core bosonic systems"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 5.0 REQUIRED)

add_library(${CMAKE_PROJECT_NAME} STATIC src/library.cpp)
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(fock PUBLIC Eigen3::Eigen)

add_subdirectory(examples)
>>>>>>> remapped

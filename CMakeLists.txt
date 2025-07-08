cmake_minimum_required(VERSION 3.31)
project(fock VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 20)
set(FOCK_VERSION ${CMAKE_PROJECT_VERSION})

include(FetchContent)

find_package(Eigen3 QUIET CONFIG)

if (NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found, fetching via FetchContent")
    FetchContent_Declare(Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG origin/master
    )
    FetchContent_MakeAvailable(Eigen3)
endif ()

add_library(fock STATIC src/single_particle.cpp)

target_include_directories(fock PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(fock PUBLIC Eigen3::Eigen)

add_subdirectory(examples)
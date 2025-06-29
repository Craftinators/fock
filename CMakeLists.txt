cmake_minimum_required(VERSION 3.31)
project(fock)

set(CMAKE_CXX_STANDARD 20)

add_library(fock STATIC src/library.cpp)

target_include_directories(fock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

include(FetchContent)

FetchContent_Declare(Eigen3
        GIT_REPOSITORY
        https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
)

FetchContent_Declare(
        Boost
        GIT_REPOSITORY https://github.com/boostorg/boost.git
        GIT_TAG        boost-1.88.0
)

FetchContent_MakeAvailable(Eigen3 Boost)

target_link_libraries(fock
        PUBLIC Eigen3::Eigen Boost::multiprecision
)

add_subdirectory(app)
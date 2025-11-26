cmake_minimum_required(VERSION 4.0)
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

cmake_minimum_required(VERSION 4.0)
project(fock VERSION 1.0.0
    DESCRIPTION "A minimal quantum computational library for hard-core bosonic systems"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EIGEN_REQUIRED_VERSION 5.0.0)
find_package(Eigen3 ${EIGEN_REQUIRED_VERSION} EXACT QUIET CONFIG)
if (NOT Eigen3_FOUND)
    message(STATUS "Eigen3 (version ${EIGEN_REQUIRED_VERSION}) not found in current environment. Fetching with FetchContent")
    include(FetchContent)
    FetchContent_Declare(eigen GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG ${EIGEN_REQUIRED_VERSION})
    FetchContent_MakeAvailable(eigen)
endif ()

add_library(${CMAKE_PROJECT_NAME} STATIC src/library.cpp)
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(fock PUBLIC Eigen3::Eigen)

add_subdirectory(examples)

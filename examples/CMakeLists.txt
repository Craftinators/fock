cmake_minimum_required(VERSION 4.0)

<<<<<<< HEAD
include(FetchContent)

FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG v1.2.2
)
FetchContent_MakeAvailable(matplotplusplus)

FetchContent_Declare(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog
        GIT_TAG v1.15.3
)
FetchContent_MakeAvailable(spdlog)

macro(add_fock_example TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} PRIVATE fock matplot spdlog::spdlog)
    target_compile_features(${TARGET_NAME} PRIVATE cxx_std_20)
endmacro()

add_subdirectory(minimal)
=======
macro(add_example TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} PRIVATE ${CMAKE_PROJECT_NAME})
endmacro()

add_example(bitmask_manipulation bitmask_manipulation.cpp)
>>>>>>> remapped

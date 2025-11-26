cmake_minimum_required(VERSION 4.0)

macro(add_example TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} PRIVATE ${CMAKE_PROJECT_NAME})
endmacro()

add_example(bitmask_manipulation bitmask_manipulation.cpp)

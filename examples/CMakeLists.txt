macro(add_fock_example TARGET_NAME)
    include(FetchContent)

    find_package(Matplot++ QUIET)
    if(NOT Matplot++_FOUND)
        message(STATUS "Matplot++ not found, fetching via FetchContent")
        FetchContent_Declare(matplotplusplus
                GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
                GIT_TAG origin/master
        )
        FetchContent_MakeAvailable(matplotplusplus)
    endif()

    find_package(spdlog QUIET)
    if(NOT spdlog_FOUND)
        message(STATUS "spdlog not found, fetching via FetchContent")
        include(FetchContent)
        FetchContent_Declare(spdlog
                GIT_REPOSITORY https://github.com/gabime/spdlog.git
                GIT_TAG origin/v1.x
        )
        FetchContent_MakeAvailable(spdlog)
    endif()

    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} PRIVATE fock matplot spdlog)
endmacro()

add_subdirectory(single_particle)
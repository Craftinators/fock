macro(add_fock_example TARGET_NAME)
    include(FetchContent)

    find_package(Matplot++ QUIET)
    if(NOT Matplot++_FOUND)
        message(STATUS "Matplot++ not found, fetching via FetchContent")
        FetchContent_Declare(matplotplusplus
                GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
                GIT_TAG fe308e13f6c3ac7c5bd2b2f4717961d1d56d07ea
        )
        FetchContent_MakeAvailable(matplotplusplus)
    endif()

    find_package(spdlog QUIET)
    if(NOT spdlog_FOUND)
        message(STATUS "spdlog not found, fetching via FetchContent")
        include(FetchContent)
        FetchContent_Declare(spdlog
                GIT_REPOSITORY https://github.com/gabime/spdlog.git
                GIT_TAG 6fa36017cfd5731d617e1a934f0e5ea9c4445b13
        )
        FetchContent_MakeAvailable(spdlog)
    endif()

    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} PRIVATE fock matplot spdlog)
endmacro()

add_subdirectory(minimal)
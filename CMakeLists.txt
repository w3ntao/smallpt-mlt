cmake_minimum_required(VERSION 3.10)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()

set(PROJ_NAME "smallpt-cpu")

project(${PROJ_NAME})

set(CMAKE_CXX_STANDARD 17)

find_package(Threads)

add_executable(${PROJ_NAME}
        src/main.cpp
        src/lodepng/lodepng.cpp
)

target_include_directories(
        ${PROJ_NAME} PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(${PROJ_NAME} Threads::Threads)

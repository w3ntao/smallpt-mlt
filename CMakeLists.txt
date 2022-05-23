cmake_minimum_required(VERSION 3.10)
project(smallpt)

set(CMAKE_CXX_STANDARD 17)

find_package (Threads)

add_executable(smallpt main.cpp)

target_link_libraries (smallpt Threads::Threads)

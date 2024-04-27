cmake_minimum_required(VERSION 3.27)
project(SiWiR2_EX1)

set(CMAKE_CXX_STANDARD 14)

add_compile_options(-O3 -Wall -Winline -Wshadow -std=c++17)

add_executable(SiWiR2_EX1
        mgsolve.cpp
        grid.cpp
        grid.hpp)

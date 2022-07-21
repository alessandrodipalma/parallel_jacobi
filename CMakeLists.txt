cmake_minimum_required(VERSION 3.22)
project(parallel_jacobi)

set(CMAKE_CXX_STANDARD 14)

add_executable(parallel_jacobi
        src/seq/sequential.cpp src/common/common.hpp src/common/tests.cpp)

cmake_minimum_required(VERSION 3.16)
project(parallel_jacobi)
#
#add_compile_definitions(PRINT_ITER)
#add_compile_definitions(CHECK_ONES)
add_compile_definitions(MEASURE_ITERATES)

include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE -O3)
#set(CMAKE_CXX_FLAGS_RELEASE -ftree-vectorize)
#set(CMAKE_CXX_FLAGS_RELEASE -fopt-info-vec-missed-all-note)
set(THREADS_PREFER_PTHREAD_FLAG ON)


find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# Fetch argparse from GitHub to parse arguments from command line
# in a modern C++ way.
message(STATUS "Fetching 'FastFlow' from GitHub...")
FetchContent_Declare(
        fastflow
        GIT_REPOSITORY https://github.com/fastflow/fastflow.git
)
FetchContent_MakeAvailable(fastflow)

# FastFlow is header-only
include_directories(${fastflow_SOURCE_DIR})

add_executable(${PROJECT_NAME}
        include/common.hpp
        include/solvers.h
        src/main.cpp
        include/timer.hpp
        src/jacobi_sequential.cpp
        src/jafarm_stdcpp.cpp
        src/jafarm_ff.cpp
        src/jafarm_omp.cpp
        src/jared_omp.cpp
        src/jacobi_double_reduce.cpp
        src/jared_ff.cpp src/jacobi_double_reduce.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC Threads::Threads  OpenMP::OpenMP_CXX)
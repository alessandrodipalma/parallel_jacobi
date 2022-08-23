cmake_minimum_required(VERSION 3.16)
project(parallel_jacobi)

include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_C_FLAGS -pthread)
set(CMAKE_CXX_FLAGS_RELEASE -O3)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

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
        src/main.cpp

        include/timer.hpp
        src/sequential.cpp
        src/parallel_jacobi.cpp
        )

target_link_libraries(${PROJECT_NAME} Threads::Threads)
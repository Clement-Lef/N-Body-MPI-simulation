cmake_minimum_required(VERSION 3.4)
project(NBodyMPI)



find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -std=c++11")

set(SOURCE_FILES main.cpp Body.cpp Body.h Node.cpp Node.h Time.cpp Time.h)
add_executable(NBodyMPI ${SOURCE_FILES})

target_link_libraries(NBodyMPI ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(NBodyMPI PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(NBodyMPI PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
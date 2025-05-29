cmake_minimum_required(VERSION 4.0.0)

project(kursovaya)

find_program(CMAKE_CXX_COMPILER NAMES g++ clang++) 
if(CMAKE_CXX_COMPILER)
   message(STATUS "Using compiler: ${CMAKE_CXX_COMPILER}")
else()
   message(WARNING "No C++ compiler found. Using default.")
endif()

include_directories(include)
add_executable(generate gridGeneration.cpp "common_includes.h")

add_executable(solve hexFEM.cpp input.cpp "common_includes.h" "input.h")
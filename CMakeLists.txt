cmake_minimum_required(VERSION 4.0.0)

project(kursovaya)

find_program(CMAKE_CXX_COMPILER NAMES g++ clang++) 
if(CMAKE_CXX_COMPILER)
   message(STATUS "Using compiler: ${CMAKE_CXX_COMPILER}")
else()
   message(WARNING "No C++ compiler found. Using default.")
endif()

include_directories(include)
add_executable(kursovaya input.cpp gridGeneration.cpp hexFEM.cpp "common_includes.h")
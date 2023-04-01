
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
project(testing)
add_subdirectory(pybind11)

pybind11_add_module(testing testing.cpp)

add_executable(testing_2 testing_2.cpp)

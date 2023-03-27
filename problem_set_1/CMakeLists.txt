
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
project(fdtd)
add_subdirectory(pybind11)

pybind11_add_module(fdtd fdtd.cpp)

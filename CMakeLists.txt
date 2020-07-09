cmake_minimum_required(VERSION 3.15)
project(Battery_Project_Threading)
set(CMAKE_CXX_STANDARD 14)
add_executable(Battery_Project_Threading main.cpp load_data.cpp load_data.h helpers.cpp helpers.h implicitRK4.cpp implicitRK4.h static_simulate.cpp static_simulate.h)
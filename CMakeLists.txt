cmake_minimum_required(VERSION 3.19)
project(MCMC)

set(CMAKE_CXX_STANDARD 14)

add_executable(MCMC main.cpp Solver.cpp Solver.h)
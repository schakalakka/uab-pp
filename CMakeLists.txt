cmake_minimum_required(VERSION 3.3)
project(pp)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c11 -g -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99 -lm -Ofast -fopenmp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")


add_executable(laplace-jacobi-iteration Lab-Ex-1/laplace-jacobi-iteration_3.c)
add_executable(laplace-jacobi-iteration-1 Lab-Ex-1/laplace-jacobi-iteration_1.c)
add_executable(laplace-jacobi-iteration-pointer Lab-Ex-1/laplace2d_5.c)

add_executable(mergesort Lab-Ex-1/mergesort.c)
add_executable(finite-differences Lab-Ex-1/finite-differences.c)
add_executable(finite-differences-3 Lab-Ex-1/finite-differences_3.c)

add_executable(Diff c_programming_performance_lab/diff.c)
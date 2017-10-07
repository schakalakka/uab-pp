cmake_minimum_required(VERSION 3.3)
project(pp)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c11 -g -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99 -lm -Ofast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")


add_executable(laplace-jacobi-iteration Lab-Ex-1/laplace-jacobi-iteration.c)

add_executable(L2Dg_withOfast Lab-Ex-1/ex1opt/laplace2d.c)
add_executable(L2Dg_1 Lab-Ex-1/ex1opt/laplace2d_1.c)
add_executable(L2Dg_2 Lab-Ex-1/ex1opt/laplace2d_2.c)
add_executable(L2Dg_3 Lab-Ex-1/ex1opt/laplace2d_3.c)
add_executable(L2Dg_4 Lab-Ex-1/ex1opt/laplace2d_4.c)
add_executable(L2Dg_5 Lab-Ex-1/ex1opt/laplace2d_5.c)


add_executable(mergesort Lab-Ex-1/mergesort.c)
add_executable(finite-differences Lab-Ex-1/finite-differences.c)
cmake_minimum_required(VERSION 3.16)

project(SimpleSparseMat)
set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

enable_testing()
add_test(SparseMatTest testMain)
add_executable(testMain test.cpp)

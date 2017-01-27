project(idrs_mtl_test)

cmake_minimum_required(VERSION 2.8)

find_package(MTL REQUIRED)
include_directories(${MTL_INCLUDE_DIRS})
add_definitions(${MTL_CXX_DEFINITIONS})

include_directories(../include)
add_executable(ex1 ex1.cpp)

cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0048 "NEW")

project(delaunay LANGUAGES CXX VERSION 1.0)

add_library(${PROJECT_NAME} STATIC ${PROJECT_SOURCE_DIR}/delabella.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

cmake_minimum_required(VERSION 3.16)

project(delabella LANGUAGES CXX VERSION 2.0)

add_library(${PROJECT_NAME} STATIC ${PROJECT_SOURCE_DIR}/delabella.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

add_library(delabella::delabella ALIAS delabella)

install(FILES ${PROJECT_SOURCE_DIR}/delabella.h DESTINATION include)
install(TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION lib)

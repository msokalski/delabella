cmake_minimum_required(VERSION 3.16)

project(delabella LANGUAGES CXX VERSION 2.0)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/delabella.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

add_library(delabella::delabella ALIAS delabella)

install(FILES ${PROJECT_SOURCE_DIR}/delabella.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

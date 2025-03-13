cmake_policy(SET CMP0048 NEW)
project(EDDI VERSION 1.0.0)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/io
    ${CMAKE_CURRENT_SOURCE_DIR}/include/density_field
    ${CMAKE_CURRENT_SOURCE_DIR}/include/molecule
)

option(BUILD_TESTS "Build the tests" ON)

add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()





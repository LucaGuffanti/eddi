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
    ${CMAKE_CURRENT_SOURCE_DIR}/include/util

)

option(BUILD_TESTS "Build the tests" ON)

add_subdirectory(src)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(EDDI PRIVATE -Wall -Wextra -DEDDI_DEBUG -g -O0)
    add_compile_definitions(EDDI_DEBUG)
    message(STATUS "Debug mode activated")
else(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(EDDI PRIVATE -O2 -DEDDI_DEBUG -Wall -Wextra -g)
    add_compile_definitions(EDDI_DEBUG)

    message(STATUS "Release mode activated")
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()





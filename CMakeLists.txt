cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0079 NEW)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/include/density
)

option(BUILD_TESTS "Build the tests" ON)

add_subdirectory(src)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug mode activated")

    add_compile_definitions(EDDI_DEBUG)

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -Wall -Wextra -DEDDI_DEBUG -g -O0 -static-libasan -fsanitize=address")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
    message(STATUS "Address sanitizer enabled")

    set(ENV{ASAN_OPTIONS} "detect_leaks=1")
    message(STATUS "Address sanitizer enabled to detect memory leaks")

else(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_definitions(EDDI_DEBUG)
    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O3")

    message(STATUS "Release mode activated")
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()





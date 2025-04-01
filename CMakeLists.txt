cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0079 NEW)

# Set project version
project(EDDI VERSION 1.0.0)

# Set minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set C standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Eddi Library (CPU-Based with OpenMP)
add_subdirectory(eddi)

# Eddi-g Library (GPU Accelerated with CUDA)
add_subdirectory(eddi-g)

# Add tests if BUILD_TESTS is ON
if(BUILD_TESTS)
    message(STATUS "=== Configuring Tests ===")
    enable_testing()
    add_subdirectory(test)
    message(STATUS "=== Tests configured ===")
endif()
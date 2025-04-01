cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0104 NEW)

# Set project version
project(EDDI VERSION 1.0.0 LANGUAGES C)

# Set minimum required version of CMake

# Set C standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find CUDA
find_package(CUDA QUIET)

# Eddi Library (CPU-Based with OpenMP)
add_subdirectory(eddi)

# Eddi-g Library (GPU Accelerated with CUDA)
if (NOT CUDA_FOUND)
    message(WARNING "CUDA not found. Please install CUDA to build EDDI-G.")
else()
    enable_language(CUDA)
    message(STATUS "CUDA found. Building EDDI-G.")
    add_subdirectory(eddi-g)
endif()

# Add tests if BUILD_TESTS is ON
if(BUILD_TESTS)
    message(STATUS "=== Configuring Tests ===")
    enable_testing()
    add_subdirectory(test)
    message(STATUS "=== Tests configured ===")
endif()
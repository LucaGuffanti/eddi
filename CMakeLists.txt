cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0079 NEW)

# Set project version
project(EDDI VERSION 1.0.0)

# Set minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set C standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Include directories
set(INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/io
    ${CMAKE_CURRENT_SOURCE_DIR}/include/density_field
    ${CMAKE_CURRENT_SOURCE_DIR}/include/molecule
    ${CMAKE_CURRENT_SOURCE_DIR}/include/util
    ${CMAKE_CURRENT_SOURCE_DIR}/include/density
)

# Option to build tests
option(BUILD_TESTS "Build the tests" ON)
option(BUILD_SHARED_LIB "Build EDDI as a shared library" ON)
option(ENABLE_OPENMP "Enable OpenMP support" ON)



# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Debug, Release, Profile)." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Profile")
endif()

# Debug mode configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug mode activated")

    # Enable debug definitions
    add_compile_definitions(EDDI_DEBUG)

    # AddressSanitizer settings
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer -O0 -g -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")

    # Set environment variable to catch memory leaks
    set(ENV{ASAN_OPTIONS} "detect_leaks=1")

    message(STATUS "AddressSanitizer enabled for memory leak detection")

elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
    message(STATUS "Profiling mode activated")

    # Profiling mode configuration
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3 -Wall -Wextra -g")

else()
    message(STATUS "Release mode activated")
    
    # Release mode configuration
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O3")
    
endif()

# Add the source directory
add_subdirectory(src)

# Add tests if BUILD_TESTS is ON
if(BUILD_TESTS)
    message(STATUS "=== Configuring Tests ===")
    enable_testing()
    add_subdirectory(test)
    message(STATUS "=== Tests configured ===")
    
endif()


# Print compilation and linking flags
message(STATUS "Compilation flags: ${CMAKE_C_FLAGS}")
message(STATUS "C++ Compilation flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Linker flags: ${CMAKE_LINKER_FLAGS}")
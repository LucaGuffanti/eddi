message(STATUS "=== Configuring EDDI ===")

# Collect all source files
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

# Create shared or static library
if (BUILD_SHARED_LIB)
    message(STATUS "Configuring EDDI as a shared library")
    add_library(EDDI SHARED ${SOURCES})
else()
    message(STATUS "Configuring EDDI as a static library")
    add_library(EDDI STATIC ${SOURCES})
endif()

# Enable OpenMP if requested
if (ENABLE_OPENMP)
    message(STATUS "Configuring EDDI with OpenMP")
    find_package(OpenMP REQUIRED)
    target_link_libraries(EDDI PUBLIC OpenMP::OpenMP_C)
endif()

# Debug mode settings
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Define debug macro
    message(STATUS "Configuring EDDI in DEBUG mode")
    target_compile_definitions(EDDI PRIVATE EDDI_DEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
    # Profile mode settings
    message(STATUS "Configuring EDDI in PROFILE mode")
else()
    # Release mode settings
    message(STATUS "Configuring EDDI in RELEASE mode")
endif()

# Include directories
target_include_directories(EDDI PUBLIC ${INCLUDE_DIRECTORIES})

# Link math library
target_link_libraries(EDDI PRIVATE m)

message(STATUS "=== EDDI configured ===")

# Include directories
set(EDDI_INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/io
    ${CMAKE_CURRENT_SOURCE_DIR}/include/density_field
    ${CMAKE_CURRENT_SOURCE_DIR}/include/molecule
    ${CMAKE_CURRENT_SOURCE_DIR}/include/util
    ${CMAKE_CURRENT_SOURCE_DIR}/include/density
)

message(STATUS "Include directories: ${EDDI_INCLUDE_DIRECTORIES}")


add_subdirectory(src)

# Print compilation and linking flags
message(STATUS "Compilation flags: ${CMAKE_C_FLAGS}")
message(STATUS "C++ Compilation flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Linker flags: ${CMAKE_LINKER_FLAGS}")
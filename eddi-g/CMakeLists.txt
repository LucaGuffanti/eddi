set(EDDI_G_INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/include    
)

message(STATUS "EDDI-G Include directories: ${EDDI_G_INCLUDE_DIRECTORIES}")

add_subdirectory(src)

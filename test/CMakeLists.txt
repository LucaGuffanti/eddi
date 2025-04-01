# Tests for the Eddi Library 
add_subdirectory(eddi)

# Tests for the Eddi-g Library
if (NOT CUDA_FOUND)
    message(STATUS "Skipping EDDI-G tests.")
else()
    message(STATUS "CUDA found. Building EDDI-G tests.")
    add_subdirectory(eddi-g)
endif()
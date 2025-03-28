#
# The test directory contains tests for the EDDI library.
#

# Collect all test source files
file(GLOB_RECURSE TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

# Macro to configure individual tests
macro(eddi_test_prepare TEST_SOURCE)
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    # Link the test with the EDDI library and math library
    target_link_libraries(${TEST_NAME} EDDI m)

    # Include necessary directories
    target_include_directories(${TEST_NAME} PUBLIC ${INCLUDE_DIRECTORIES})

    # Apply AddressSanitizer flags for Debug builds (only for GCC/Clang)
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        target_compile_definitions(${TEST_NAME} PRIVATE EDDI_DEBUG)
    endif()
endmacro()

# Macro to add a basic test
macro(add_eddi_test TEST_SOURCE)
    message(STATUS "Adding test ${TEST_SOURCE}")
    eddi_test_prepare(${TEST_SOURCE})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${ARGN})
endmacro()

# Macro to add a test with an input file
macro(add_eddi_test_with_input_file TEST_SOURCE INPUT_FILE)
    eddi_test_prepare(${TEST_SOURCE})

    message(STATUS "Adding test ${TEST_SOURCE} with input file ${INPUT_FILE}")

    if (EXISTS ${INPUT_FILE})
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${INPUT_FILE} ${ARGN})
    else()
        message(WARNING "Input file ${INPUT_FILE} does not exist. Skipping test ${TEST_NAME}.")
    endif()
endmacro()

message(STATUS "-> Tests in ${CMAKE_CURRENT_SOURCE_DIR}")

# Add tests
add_eddi_test(test_include.c)
add_eddi_test(test_h2.c h2.cube h2.bin)
add_eddi_test(test_gaussian_io.c h2.cube h22.cube)
add_eddi_test(test_binary_writer.c h2.cube)
add_eddi_test(test_binary_reader.c h2.cube h21.cube h22.bin h22.cube)
add_eddi_test(test_symbol_to_number.c)
add_eddi_test(test_he2.c he2.cube he2.bin)

# Add tests that require input files to be checked
add_eddi_test_with_input_file(test_h2o.c ${CMAKE_CURRENT_BINARY_DIR}/../../data/h2o.pdb h2o.cube h2o.bin)
add_eddi_test_with_input_file(test_pdb_reader.c ${CMAKE_CURRENT_BINARY_DIR}/../../data/8z9k.pdb)
add_eddi_test_with_input_file(test_mol.c ${CMAKE_CURRENT_BINARY_DIR}/../../data/ADP_ideal.sdf)
add_eddi_test_with_input_file(test_adenosine.c ${CMAKE_CURRENT_BINARY_DIR}/../../data/ADP_ideal.sdf adp.cube adp.bin)
# add_eddi_test_with_input_file(test_8z9k.c ${CMAKE_CURRENT_BINARY_DIR}/../../data/8z9k.pdb 8z9k.cube 8z9k.bin)
# add_eddi_test_with_input_file(test_9jt3.c ${CMAKE_CURRENT_BINARY_DIR}/../../data/9jt3.pdb 9jt3.cube 9jt3.bin)



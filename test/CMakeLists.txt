#
# The test directory contains tests for the EDDI library.
#

file(GLOB_RECURSE TEST_SOURCES "*.c")

macro(eddi_test_prepare TEST_SOURCE)
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} EDDI)
    target_include_directories(${TEST_NAME} PUBLIC ${INCLUDE_DIRECTORIES})
endmacro()

macro(add_eddi_test TEST_SOURCE)
    eddi_test_prepare(${TEST_SOURCE})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${ARGN})
endmacro()

macro(add_eddi_test_with_input_file TEST_SOURCE INPUT_FILE)
    eddi_test_prepare(${TEST_SOURCE})
    if (EXISTS ${INPUT_FILE})
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${INPUT_FILE} ${ARGN})
    else ()
        message(WARNING "Input file ${INPUT_FILE} does not exist. Skipping test ${TEST_NAME}.")
        return()
    endif()
endmacro()



message(STATUS "Configuring tests ${CMAKE_CURRENT_SOURCE_DIR}")

add_eddi_test(test_include.c)
add_eddi_test(test_h2.c h2.cube h2.bin)
add_eddi_test(test_gaussian_io.c h2.cube h22.cube)
add_eddi_test(test_binary_writer.c h2.cube)
add_eddi_test(test_binary_reader.c h2.cube h21.cube h22.bin h22.cube)
add_eddi_test_with_input_file(test_pdb_reader.c ${CMAKE_CURRENT_BINARY_DIR}/../../data/8z9k.pdb)
add_eddi_test(test_symbol_to_number.c)
add_eddi_test(test_he2.c he2.cube he2.bin)
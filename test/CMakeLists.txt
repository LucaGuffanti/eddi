#
# The test directory contains tests for the EDDI library.
#

file(GLOB_RECURSE TEST_SOURCES "*.c")
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} EDDI)
    target_include_directories(${TEST_NAME} PUBLIC ${INCLUDE_DIRECTORIES})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    
endforeach()
    

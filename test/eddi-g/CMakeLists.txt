

macro(eddi_g_test TEST_FILE)
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  add_executable(${TEST_NAME} ${TEST_FILE})
  
  target_link_libraries(${TEST_NAME} PRIVATE EDDI-G)
  
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
  
  target_include_directories(${TEST_NAME} PRIVATE ${EDDI_G_INCLUDE_DIRECTORIES})
  
  message(STATUS "Adding test ${TEST_NAME} from ${TEST_FILE}")
endmacro()

message(STATUS "-> Tests in ${CMAKE_CURRENT_SOURCE_DIR}")

# Add tests

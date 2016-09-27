include_directories (${TEST_SOURCE_DIR}/src)
add_executable (test_integration test_integration.c)
file(COPY table.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

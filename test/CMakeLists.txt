include_directories (${TEST_SOURCE_DIR}/src)

add_executable (test_linking test_linking.c)
target_link_libraries (test_linking LINK_PUBLIC liknorm m)

add_executable (test_integration test_integration.c)
target_link_libraries (test_integration LINK_PUBLIC liknorm m)
file(COPY table.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

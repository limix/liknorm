include_directories (${TEST_SOURCE_DIR}/src)

add_executable (test_linking test_linking.c)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_linking LINK_PUBLIC liknorm)
else()
  target_link_libraries (test_linking LINK_PUBLIC liknorm m)
endif()

add_executable (test_integration test_integration.c)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_integration LINK_PUBLIC liknorm)
else()
  target_link_libraries (test_integration LINK_PUBLIC liknorm m)
endif()

file(COPY table.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

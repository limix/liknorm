include_directories (${TEST_SOURCE_DIR}/src)

add_executable (test_linking test_linking.c)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_linking liknorm)
else()
  target_link_libraries (test_linking liknorm m)
endif()

add_executable (test_integration test_integration.c)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_integration liknorm)
else()
  target_link_libraries (test_integration liknorm m)
endif()

add_executable (test_specific test_specific.c)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_specific liknorm)
else()
  target_link_libraries (test_specific liknorm m)
endif()

configure_file(table.csv ${CMAKE_CURRENT_BINARY_DIR}/DESTINATION COPYONLY)

include_directories (${TEST_SOURCE_DIR}/src)


set(LLIBS liknorm)
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(LLIBS ${LLIBS} m)
endif()

if(REQUIRES_HCEPHES)
    include_directories(${HCEPHES_INCLUDE_DIRS})
    set(LLIBS ${LLIBS} ${HCEPHES_LIBRARIES})
endif()

add_executable (test_empty test_empty.c)
target_link_libraries (test_empty ${LLIBS})

add_executable (test_linking test_linking.c)
target_link_libraries (test_linking ${LLIBS})

add_executable (test_integration test_integration.c)
target_link_libraries (test_integration ${LLIBS})

add_executable (test_specific test_specific.c)
target_link_libraries (test_specific ${LLIBS})

configure_file(table.csv ${CMAKE_CURRENT_BINARY_DIR}/DESTINATION COPYONLY)

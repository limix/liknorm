include_directories (${TEST_SOURCE_DIR}/src)

set(LLIBS liknorm)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if(MSVC90)
        set(LLIBS ${LLIBS} hcephes)
    endif()
else()
    set(LLIBS ${LLIBS} m)
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

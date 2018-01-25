include_directories("${TEST_SOURCE_DIR}/src")

limix_add_test(test_empty ${PROJECT_NAME} "test_empty.c")
limix_add_test(test_linking ${PROJECT_NAME} "test_linking.c")
limix_add_test(test_integration ${PROJECT_NAME} "test_integration.c")
limix_add_test(test_specific ${PROJECT_NAME} "test_specific.c")

list(APPEND DATA_FILES data/table.csv)

foreach(DATA_FILE ${DATA_FILES})
    set(ORIGIN "data/${DATA_FILE}")
    file(COPY ${ORIGIN} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/data")
endforeach()

include_directories("${TEST_SOURCE_DIR}/src")

limix_add_test(empty ${PROJECT_NAME} "empty.c")
limix_add_test(linking ${PROJECT_NAME} "linking.c")
limix_add_test(integration ${PROJECT_NAME} "integration.c")
limix_add_test(specific ${PROJECT_NAME} "specific.c")
limix_add_test(logprod ${PROJECT_NAME} "logprod.c")

list(APPEND DATA_FILES data/table.csv)

foreach(DATA_FILE ${DATA_FILES})
    set(ORIGIN "${DATA_FILE}")
    file(COPY ${ORIGIN} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/data")
endforeach()

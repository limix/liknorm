include_directories (${TEST_SOURCE_DIR}/src)
add_executable (test_binomial test_binomial.c)
target_link_libraries (test_binomial liknorm)

include_directories (${TEST_SOURCE_DIR}/src)
add_executable (test_binomial_normal test_binomial_normal.c)
target_link_libraries (test_binomial_normal liknorm)

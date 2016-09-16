cmake_minimum_required (VERSION 3.4)
project (liknorm)
set(CMAKE_BUILD_TYPE Release)
set_property(TARGET liknorm PROPERTY C_STANDARD 11)

add_subdirectory (src)
add_subdirectory (test)

enable_testing()
add_test( NAME BinomialNormal COMMAND test_binomial_normal )

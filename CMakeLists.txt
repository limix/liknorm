cmake_minimum_required (VERSION 3.4)
project (liknorm)
set(CMAKE_BUILD_TYPE Release)

include_directories(include)
include_directories(lib)
add_subdirectory (src)
add_subdirectory (test)

enable_testing()
add_test( NAME BinomialNormal COMMAND test_binomial_normal )

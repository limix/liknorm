cmake_minimum_required (VERSION 3.4)
project (liknorm)

set(CMAKE_BUILD_TYPE Release)
SET(CMAKE_C_FLAGS "-Ofast")

include_directories(.)
add_subdirectory(test)

enable_testing()
add_test(NAME TestIntegration COMMAND test_integration)

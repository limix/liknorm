cmake_minimum_required (VERSION 3.2.2)

project(liknorm C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

include_directories(.)

add_subdirectory(liknorm)
add_subdirectory(test)

enable_testing()
add_test(NAME TestLinking COMMAND test_linking)
add_test(NAME TestIntegration COMMAND test_integration)

install (TARGETS liknorm
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib)

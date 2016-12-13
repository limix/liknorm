cmake_minimum_required (VERSION 3.2.2)

project(Liknorm C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

file(STRINGS "VERSION" LIKNORM_VERSION)

include_directories(.)

add_subdirectory(liknorm)

add_library(liknorm SHARED $<TARGET_OBJECTS:liknorm_objs>
            $<TARGET_OBJECTS:optimizer_objs> $<TARGET_OBJECTS:partition_objs>)

add_library(liknorm_static STATIC $<TARGET_OBJECTS:liknorm_objs>
            $<TARGET_OBJECTS:optimizer_objs> $<TARGET_OBJECTS:partition_objs>)

set_target_properties(liknorm
  PROPERTIES
    VERSION ${LIKNORM_VERSION}
    PUBLIC_HEADER "liknorm/liknorm.h")

install(TARGETS liknorm
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/liknorm)

install(TARGETS liknorm_static
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/liknorm)

add_subdirectory(test)

enable_testing()
add_test(NAME TestLinking COMMAND test_linking)
add_test(NAME TestIntegration COMMAND test_integration)

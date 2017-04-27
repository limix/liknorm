cmake_minimum_required (VERSION 2.6.4)

project(Liknorm C)

set(CMAKE_C_FLAGS "-std=c99")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

file(STRINGS "VERSION" LIKNORM_VERSION)

include_directories(.)

add_subdirectory(liknorm)

add_library(liknorm SHARED liknorm/normal.c liknorm/liknorm.c liknorm/interval.c liknorm/gfunc.c
            liknorm/integrate.c liknorm/optimizer/bracket.c liknorm/optimizer/brent.c liknorm/optimizer/zero.c
            liknorm/partition/binomial.c liknorm/partition/exponential.c liknorm/partition/geometric.c
            liknorm/partition/gamma.c liknorm/partition/poisson.c)

add_library(liknorm_static STATIC liknorm/normal.c liknorm/liknorm.c liknorm/interval.c liknorm/gfunc.c
            liknorm/integrate.c liknorm/optimizer/bracket.c liknorm/optimizer/brent.c liknorm/optimizer/zero.c
            liknorm/partition/binomial.c liknorm/partition/exponential.c liknorm/partition/geometric.c
            liknorm/partition/gamma.c liknorm/partition/poisson.c)

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
add_test(TestLinking test/test_linking)
add_test(TestIntegration test/test_integration)
add_test(TestSpecific test/test_specific)

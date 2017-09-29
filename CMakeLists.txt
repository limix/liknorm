cmake_minimum_required (VERSION 2.6.4)

project(Liknorm C)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)

  if (GCC_VERSION VERSION_GREATER 4.9)
    set (CMAKE_C_STANDARD 11)
  else ()
    set (CMAKE_C_FLAGS "--std=gnu99 ${CMAKE_C_FLAGS}")
  endif ()
endif ()

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-Dinline=__inline)
endif()

file(STRINGS "VERSION" LIKNORM_VERSION)

include_directories( include/ )
include_directories( liknorm/ )
include_directories(.)

add_subdirectory(liknorm)

file(GLOB LIKNORM_SRC
    "liknorm/*.c"
)

file(GLOB OPTIMIZER_SRC
    "liknorm/optimizer/*.c"
)

file(GLOB PARTITION_SRC
    "liknorm/partition/*.c"
)

add_library(liknorm
    SHARED ${LIKNORM_SRC} ${OPTIMIZER_SRC} ${PARTITION_SRC})

add_library(liknorm_static
    STATIC ${LIKNORM_SRC} ${OPTIMIZER_SRC} ${PARTITION_SRC})

set_target_properties(liknorm
  PROPERTIES
    VERSION ${LIKNORM_VERSION}
PUBLIC_HEADER include/liknorm/liknorm.h)

install(TARGETS liknorm
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION lib
        PUBLIC_HEADER DESTINATION include/liknorm)

install(TARGETS liknorm_static
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/liknorm)

add_subdirectory(test)

enable_testing()
add_test(test_linking test/test_linking)
add_test(test_integration test/test_integration)
add_test(test_specific test/test_specific)

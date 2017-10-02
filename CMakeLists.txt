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
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(REQUIRES_HCEPHES 0)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-Dinline=__inline)

  if(MSVC90)
      set(REQUIRES_HCEPHES 1)
  endif()
endif()

file(STRINGS "VERSION" LIKNORM_VERSION)

if(REQUIRES_HCEPHES)
    find_package( HCEPHES REQUIRED )
    if(NOT HCEPHES_FOUND)
        message(STATUS "Could not find hcephes library.
                        If you have it installed, please, define
                        HCEPHES_LIBRARYDIR and HCEPHES_INCLUDEDIR
                        appropriately to fix it.")
    endif()
    include_directories( ${HCEPHES_INCLUDE_DIRS} )
endif()

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
if(REQUIRES_HCEPHES)
    target_link_libraries( liknorm ${HCEPHES_LIBRARIES} )
endif()

add_library(liknorm_static
    STATIC ${LIKNORM_SRC} ${OPTIMIZER_SRC} ${PARTITION_SRC})
if(REQUIRES_HCEPHES)
    target_link_libraries( liknorm_static ${HCEPHES_LIBRARIES} )
endif()

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
add_test(test_empty test/test_empty)
add_test(test_linking test/test_linking)
add_test(test_integration test/test_integration)
add_test(test_specific test/test_specific)

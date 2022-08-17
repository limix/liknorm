cmake_minimum_required(VERSION 3.16.3 FATAL_ERROR)
project(
  liknorm
  VERSION 1.5.7
  LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(LIKNORM_MAIN_PROJECT OFF)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(LIKNORM_MAIN_PROJECT ON)
endif()

if(LIKNORM_MAIN_PROJECT)
  set(LIKNORM_BUILD_TESTS_DEFAULT ON)
else()
  set(LIKNORM_BUILD_TESTS_DEFAULT OFF)
endif()

option(LIKNORM_BUILD_TESTS "Build the unit tests"
       ${LIKNORM_BUILD_TESTS_DEFAULT})
message(STATUS "LIKNORM_MAIN_PROJECT: " ${LIKNORM_MAIN_PROJECT})
message(STATUS "LIKNORM_BUILD_TESTS: " ${LIKNORM_BUILD_TESTS})

include(cmake/compiler-options.cmake)
include(cmake/sanitizers.cmake)
include(cmake/CPM.cmake)

cpmaddpackage("gh:horta/logaddexp@2.1.9")

add_library(
  liknorm
  src/gfunc.c
  src/integrate.c
  src/interval.c
  src/liknorm.c
  src/normal.c
  src/optimizer/bracket.c
  src/optimizer/brent.c
  src/optimizer/zero.c
  src/partition/binomial.c
  src/partition/exponential.c
  src/partition/gamma.c
  src/partition/geometric.c
  src/partition/nbinomial.c
  src/partition/poisson.c
  src/report.c)
add_library(LIKNORM::liknorm ALIAS liknorm)

# Write out mytools_export.h to the current binary directory Must include
# ${CMAKE_CURRENT_BINARY_DIR} in the target's search path.
include(GenerateExportHeader)
generate_export_header(
  liknorm
  BASE_NAME
  LIKNORM
  INCLUDE_GUARD_NAME
  LIKNORM_EXPORT_H
  EXPORT_MACRO_NAME
  LIKNORM_API
  EXPORT_FILE_NAME
  ${CMAKE_CURRENT_BINARY_DIR}/liknorm/export.h)

target_include_directories(
  liknorm
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(liknorm PUBLIC LOGADDEXP::logaddexp)
set_property(TARGET liknorm PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(liknorm PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(liknorm PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(liknorm PROPERTIES DEFINE_SYMBOL LIKNORM_EXPORTS)
target_compile_features(liknorm PRIVATE c_std_11)
target_compile_definitions(liknorm
                           PRIVATE $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>)
target_compile_definitions(
  liknorm PRIVATE $<$<BOOL:${WIN32}>:_CRT_NONSTDC_NO_DEPRECATE>)
target_link_libraries(liknorm INTERFACE $<$<BOOL:${UNIX}>:m>)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  include(CheckIPOSupported)
  check_ipo_supported()
  set_target_properties(liknorm PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
endif()

include(CheckCSourceCompiles)
check_c_source_compiles(
  "
        static void f(char const* err, ...) __attribute__((format(printf, 1, 2)));
        int main(void) {return 0;}
    "
  HAVE_ATTR_FORMAT)
target_compile_definitions(
  liknorm PRIVATE $<$<BOOL:${HAVE_ATTR_FORMAT}>:HAVE_ATTR_FORMAT>)

install(
  TARGETS liknorm
  EXPORT liknorm-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include)
install(DIRECTORY include/ DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/liknorm/export.h
        DESTINATION include/liknorm)
install(
  EXPORT liknorm-targets
  FILE liknorm-targets.cmake
  NAMESPACE LIKNORM::
  DESTINATION lib/cmake/liknorm)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  liknorm-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/liknorm-config.cmake
  INSTALL_DESTINATION lib/cmake/liknorm
  NO_SET_AND_CHECK_MACRO)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/liknorm-config.cmake
        DESTINATION lib/cmake/liknorm)

if(LIKNORM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

set(CPACK_PACKAGE_NAME liknorm)
set(CPACK_PACKAGE_VENDOR "Danilo Horta")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Moments of the product of an exponential-family likelihood with a Normal distribution"
)
set(CPACK_PACKAGE_CONTACT "Danilo Horta")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md)
set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
set(CPACK_OUTPUT_FILE_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/package)
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_PACKAGE_RELOCATABLE YES)
set(CPACK_MONOLITHIC_INSTALL YES)
include(CPack)

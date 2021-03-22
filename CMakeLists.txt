cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(liknorm VERSION 1.5.5 LANGUAGES C)

# Generate compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Set RPATH only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# merge request fix: https://gitlab.kitware.com/cmake/cmake/-/merge_requests/1799
if(CMAKE_VERSION VERSION_GREATER 3.12.0)
    # Hide things by default
    set(CMAKE_C_VISIBILITY_PRESET hidden)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif()

add_library(liknorm
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
    src/report.c
)
add_library(LIKNORM::liknorm ALIAS liknorm)

# Write out mytools_export.h to the current binary directory
# Must include ${CMAKE_CURRENT_BINARY_DIR} in the target's
# search path.
include(GenerateExportHeader)
generate_export_header(liknorm BASE_NAME LIKNORM
    INCLUDE_GUARD_NAME LIKNORM_EXPORT_H
    EXPORT_MACRO_NAME LIKNORM_API
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/liknorm/export.h)

if(${CMAKE_C_COMPILER_ID} MATCHES "MSVC")
    set(WARNING_FLAGS /W3)
else()
    set(WARNING_FLAGS -Wall -Wextra -pedantic -Wstrict-prototypes -Wshadow
        -Wconversion -Wmissing-prototypes -Wno-unused-parameter -Wsign-conversion
        -Wno-unused-function # FIXME: added only to please clangd for now
        -Wno-empty-translation-unit # FIXME: added only to please clangd for now
    )
endif()

target_include_directories(liknorm
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
set_target_properties(liknorm PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(liknorm PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(liknorm PROPERTIES C_STANDARD 99)
set_target_properties(liknorm PROPERTIES DEFINE_SYMBOL LIKNORM_EXPORTS)
target_compile_definitions(liknorm PRIVATE $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>)
target_compile_definitions(liknorm PRIVATE $<$<BOOL:${WIN32}>:_CRT_NONSTDC_NO_DEPRECATE>)
target_link_libraries(liknorm INTERFACE $<$<BOOL:${UNIX}>:m>)

include(CheckCSourceCompiles)
check_c_source_compiles(
    "
        static void f(char const* err, ...) __attribute__((format(printf, 1, 2)));
        int main(void) {return 0;}
    "
    HAVE_ATTR_FORMAT
)
target_compile_definitions(liknorm PRIVATE $<$<BOOL:${HAVE_ATTR_FORMAT}>:HAVE_ATTR_FORMAT>)

install(TARGETS liknorm EXPORT liknorm-targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)
install(DIRECTORY include/ DESTINATION include)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/liknorm/export.h
    DESTINATION include/liknorm
)
install(EXPORT liknorm-targets
    FILE liknorm-targets.cmake
    NAMESPACE LIKNORM::
    DESTINATION lib/cmake/liknorm
)
include(CMakePackageConfigHelpers)
configure_package_config_file(liknorm-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/liknorm-config.cmake
    INSTALL_DESTINATION lib/cmake/liknorm
    NO_SET_AND_CHECK_MACRO
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/liknorm-config.cmake
    DESTINATION lib/cmake/liknorm
)

enable_testing()
add_subdirectory(test)

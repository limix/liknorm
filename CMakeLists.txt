cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
project(liknorm VERSION 1.5.2 LANGUAGES C)

# Generate compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Set RPATH only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
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

find_package(logaddexp REQUIRED)

target_include_directories(liknorm
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(liknorm PRIVATE LOGADDEXP::logaddexp)
set_target_properties(liknorm PROPERTIES PUBLIC_HEADER include/liknorm.h)
set_target_properties(liknorm PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(liknorm PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(liknorm PROPERTIES C_STANDARD 99)
set_target_properties(liknorm PROPERTIES DEFINE_SYMBOL "LIKNORM_EXPORTS;_CRT_SECURE_NO_WARNINGS")
target_compile_definitions(liknorm PRIVATE $<$<BOOL:${WIN32}>:_CRT_NONSTDC_NO_DEPRECATE>)
target_compile_definitions(liknorm PRIVATE $<$<BOOL:${WIN32}>:__extension__>)
target_link_libraries(liknorm INTERFACE $<$<BOOL:${UNIX}>:m>)
target_compile_options(liknorm PRIVATE
    $<$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>:/W3 /WX>
    $<$<NOT:$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>>:-Wall -Wextra -Wstrict-prototypes
        -Wshadow -Wconversion -Wmissing-prototypes -Wno-unused-parameter
        -Wno-empty-translation-unit>
)

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

cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
include(cmake/limix.cmake)
limix_initialise()
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION} LANGUAGES C)
limix_config()

set(CMAKE_C_STANDARD 11)

if (NOT c_restrict IN_LIST CMAKE_C_COMPILE_FEATURES)
    add_compile_definitions(restrict= )
endif()

if(NOT MSVC)
    add_compile_options(-Wall -Wextra -Wno-pedantic -Werror -Wno-missing-field-initializers -Wno-deprecated-declarations)
endif(NOT MSVC)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(LIBS "")
else()
    set(LIBS "m")
endif()

limix_process_default_dirs()

limix_add_library(${PROJECT_NAME} ${PROJECT_VERSION} "${SOURCES}" "${HEADERS}"
    "${PRIVATE_HEADERS}" "${LIBS}")

add_subdirectory(test)

cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
include(cmake/limix.cmake)
limix_initialise()
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION} LANGUAGES C)
limix_config()

set(CMAKE_C_STANDARD 11)

find_package(HCEPHES REQUIRED)
include_directories("${HCEPHES_INCLUDE_DIRS}")

limix_process_default_dirs()

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(LIBS "")
else()
    set(LIBS "m")
endif()

limix_add_library(${PROJECT_NAME} ${PROJECT_VERSION} "${SOURCES}" "${HEADERS}"
    "${HCEPHES_LIBRARIES};${LIBS}")

add_subdirectory(test)

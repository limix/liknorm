# FindHCEPHES
# ---------
#
# Find HCEPHES include dirs and libraries
#
# This module reads hints about search locations from variables::
#
#   HCEPHES_INCLUDEDIR      - Preferred include directory e.g. <prefix>/include
#   HCEPHES_LIBRARYDIR      - Preferred library directory e.g. <prefix>/lib

set(_hcephes_INCLUDE_SEARCH_DIRS "/usr/include" "/usr/local/include")
if(HCEPHES_INCLUDEDIR)
    list(APPEND _hcephes_INCLUDE_SEARCH_DIRS ${HCEPHES_INCLUDEDIR})
endif()

find_path(
    HCEPHES_INCLUDE_DIR
    NAMES hcephes/hcephes.h
    HINTS ${_hcephes_INCLUDE_SEARCH_DIRS}
)


set(_hcephes_LIBRARY_SEARCH_DIRS "/usr/lib" "/usr/local/lib")
if(HCEPHES_LIBRARYDIR)
    list(APPEND _hcephes_LIBRARY_SEARCH_DIRS ${HCEPHES_LIBRARYDIR})
endif()

find_library(
    HCEPHES_LIBRARY
    NAMES hcephes libhcephes
    HINTS ${_hcephes_LIBRARY_SEARCH_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HCEPHES DEFAULT_MSG
HCEPHES_LIBRARY HCEPHES_INCLUDE_DIR)

mark_as_advanced(HCEPHES_INCLUDE_DIR HCEPHES_LIBRARY)

set(HCEPHES_LIBRARIES ${HCEPHES_LIBRARY})
set(HCEPHES_INCLUDE_DIRS ${HCEPHES_INCLUDE_DIR})

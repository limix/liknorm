add_library (liknorm normal.c liknorm.c interval.c gfunc.c integrate.c)
add_subdirectory(optimizer)
add_subdirectory(partition)
target_link_libraries (liknorm LINK_PUBLIC optimizer partition)

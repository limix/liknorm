add_subdirectory(optimizer)
add_subdirectory(partition)

add_library(liknorm_objs OBJECT normal.c liknorm.c interval.c gfunc.c
            integrate.c $<TARGET_OBJECTS:optimizer_objs>
            $<TARGET_OBJECTS:partition_objs>)

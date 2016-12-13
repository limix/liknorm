add_subdirectory(optimizer)
add_subdirectory(partition)

add_library(liknorm_objs OBJECT normal.c liknorm.c interval.c gfunc.c
            integrate.c $<TARGET_OBJECTS:optimizer>
            $<TARGET_OBJECTS:partition>)

add_subdirectory(optimizer)
add_subdirectory(partition)

add_library (liknorm normal.c liknorm.c interval.c gfunc.c integrate.c
             optimizer/bracket.c optimizer/brent.c optimizer/zero.c
             partition/binomial.c partition/exponential.c partition/geometric.c
             partition/gamma.c partition/poisson.c)

set_target_properties(liknorm PROPERTIES PUBLIC_HEADER "liknorm.h")

install (TARGETS liknorm
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         PUBLIC_HEADER DESTINATION include/liknorm)

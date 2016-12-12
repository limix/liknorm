add_library (liknorm normal.c liknorm.c interval.c gfunc.c partition/binomial.c
             optimizer/bracket.c optimizer/brent.c)
add_subdirectory(partition)

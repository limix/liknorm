cmake_minimum_required (VERSION 3.4)
project (liknorm)

add_library (liknorm src/main.c)
set_property(TARGET tgt PROPERTY C_STANDARD 11)

enable_testing()
add_test( binomial exec_name )

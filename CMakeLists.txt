cmake_minimum_required (VERSION 3.4)
project (liknorm)

set(CMAKE_BUILD_TYPE Release)

file(GLOB SOURCES "src/*.cpp")

add_library (liknorm SHARED ${SOURCES})

set_property(TARGET liknorm PROPERTY C_STANDARD 11)

# enable_testing()
# add_test( binomial exec_name )

# install(TARGETS testStudent DESTINATION /usr/lib)

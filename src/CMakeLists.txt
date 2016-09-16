file(GLOB SOURCES "*.c")
add_library (liknorm SHARED ${SOURCES})
set_property(TARGET liknorm PROPERTY C_STANDARD 11)

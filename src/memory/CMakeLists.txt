FILE(GLOB sourcefiles "*.f90")
add_library(memory STATIC ${sourcefiles})
add_dependencies(memory base)

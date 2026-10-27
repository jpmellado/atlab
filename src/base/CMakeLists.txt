FILE(GLOB sourcefiles "*.f90")
add_library(base STATIC ${sourcefiles})

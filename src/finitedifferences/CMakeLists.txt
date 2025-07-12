FILE(GLOB sourcefiles "*.f90")
add_library(finitedifferences STATIC ${sourcefiles})
add_dependencies(finitedifferences utils base)
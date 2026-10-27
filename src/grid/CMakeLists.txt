FILE(GLOB sourcefiles "*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_MPI" IN_LIST defs )
else()
  FILE(GLOB to_remove "tlab_mpi*.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()
add_library(grid STATIC ${sourcefiles})
add_dependencies(grid memory base)

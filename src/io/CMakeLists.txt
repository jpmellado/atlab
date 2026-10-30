FILE(GLOB sourcefiles "*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_NETCDF" IN_LIST defs )
else()
  FILE(GLOB to_remove "*netcdf*.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()
add_library(io STATIC ${sourcefiles})
add_dependencies(io grid memory base)

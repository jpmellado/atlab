FILE(GLOB sourcefiles "*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_MPI" IN_LIST defs )
else()
  FILE(GLOB to_remove "*mpi*.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()
add_library(finitedifferences STATIC ${sourcefiles})
add_dependencies(finitedifferences utils base)

FILE(GLOB sourcefiles "fdm_derivative_X1.f90")
add_executable(vfdm1.x ${sourcefiles})
target_link_libraries(vfdm1.x finitedifferences utils base ${LIBS})

# FILE(GLOB sourcefiles "fdm_derivative_X2.f90")
# add_executable(vfdm2.x ${sourcefiles})
# target_link_libraries(vfdm2.x finitedifferences utils base ${LIBS})

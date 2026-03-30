if ( NOT BUILD_TYPE )
    set(BUILD_TYPE LITTLE)
endif()
message( STATUS "Build Type: " ${BUILD_TYPE} )

set(USER_Fortran_FLAGS "-fconvert=little-endian -cpp -ffree-form -ffree-line-length-none -fno-automatic -fallow-argument-mismatch")
set(USER_Fortran_FLAGS_RELEASE "-O3 -ffpe-summary=none -ffast-math -mtune=native -march=native")
set(USER_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow") # ,underflow,precision,denormal")

set(CMAKE_Fortran_COMPILER gfortran)

if (BUILD_TYPE STREQUAL "PARALLEL" )    # compiler for parallel build
    set(PARALLEL MPI_ONLY)
    add_definitions(-DUSE_MPI -DUSE_MPI_IO)
    # set(CMAKE_BUILD_TYPE DEBUG)
    
else()                                  # compiler for serial build
    if    (BUILD_TYPE  STREQUAL "BIG" )
        set(USER_Fortran_FLAGS_RELEASE "-fconvert=big-endian ${USER_Fortran_FLAGS_RELEASE}")
        
    elseif(BUILD_TYPE  STREQUAL "LITTLE" )
        
    elseif(BUILD_TYPE  STREQUAL "PROFILE" )
        set(USER_Fortran_FLAGS_DEBUG "-pg ${USER_Fortran_FLAGS_DEBUG}")
        set(CMAKE_BUILD_TYPE DEBUG)
        
    else()
        add_definitions(-D_DEBUG)
        set(CMAKE_BUILD_TYPE DEBUG)
        
    endif()
    
endif()

if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

add_definitions(-DUSE_FFTW)
#set(FFTW_INCLUDE_DIR   "/usr/local/include")
#set(FFTW_LIB           "/usr/local/lib/libfftw3.a")
set(FFTW_LIB           "-lfftw3")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

add_definitions(-DUSE_NETCDF)
set(NC_INCLUDE_DIR     "/usr/include")
set(NC_LIB             "-L/usr/lib -lnetcdff -lnetcdf")
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${NC_INCLUDE_DIR})
set(LIBS ${LIBS} ${NC_LIB})

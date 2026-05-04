if(NOT BUILD_TYPE)
    message(WARNING "Setting CMAKE_BUILD_TYPE to default value.")
    set(BUILD_TYPE LITTLE)
endif()

set(USER_Fortran_FLAGS "-cpp -ffree-form -ffree-line-length-none -fconvert=little-endian -fno-automatic -fallow-argument-mismatch")
set(USER_Fortran_FLAGS_RELEASE "-O3 -ffpe-summary=none -ffast-math -mtune=native -march=native")
set(USER_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow") # ,underflow,precision,denormal")

set(CMAKE_Fortran_COMPILER gfortran)

# set(USER_Fortran_FLAGS  " -fpp ${USER_profile_FLAGS} -nbs -save-temps -heap-arrays -simd -vec-threshold50 -unroll-aggressive ${USER_omp_FLAGS} " )
# set(USER_Fortran_FLAGS_RELEASE  " -march=core-avx2 -mtune=core-avx2 -qopt-prefetch -O3 -ipo" )
# set(USER_Fortran_FLAGS_DEBUG    " -g -traceback -debug all ")

# set(CMAKE_Fortran_COMPILER ifx)

if(${BUILD_TYPE} STREQUAL "PARALLEL")
    set(PARALLEL MPI_ONLY)
    add_definitions(-DUSE_MPI)
    # set(CMAKE_BUILD_TYPE DEBUG)

else()
    if(${BUILD_TYPE} STREQUAL "BIG")
        set(USER_Fortran_FLAGS_RELEASE "-fconvert=big-endian ${USER_Fortran_FLAGS_RELEASE}")

    elseif(${BUILD_TYPE} STREQUAL "LITTLE")

    elseif(${BUILD_TYPE} STREQUAL "PROFILE")
        set(USER_Fortran_FLAGS_DEBUG "-pg ${USER_Fortran_FLAGS_DEBUG}")
        set(CMAKE_BUILD_TYPE DEBUG)

    else(${BUILD_TYPE} STREQUAL "DEBUG")
        set(CMAKE_BUILD_TYPE DEBUG)
        add_definitions(-D_DEBUG)

    endif()
    
endif()

if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

add_definitions(-DUSE_FFTW)

# set(FFTW_INCLUDE_DIR   "/usr/local/include")
# set(FFTW_LIB           "/usr/local/lib/libfftw3.a")
set(FFTW_LIB "-lfftw3")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

add_definitions(-DUSE_NETCDF)
set(NC_INCLUDE_DIR "/sw/buster-x64/io/netcdf-c-4.7.4-fortran-4.5.2-cxx4-4.3.1-cxx-4.2-gccsys/include")
set(NC_LIB "-L/sw/buster-x64/io/netcdf-c-4.7.4-fortran-4.5.2-cxx4-4.3.1-cxx-4.2-gccsys/lib -Wl,-rpath -Wl,/sw/buster-x64/io/netcdf-c-4.7.4-fortran-4.5.2-cxx4-4.3.1-cxx-4.2-gccsys/lib -lnetcdff")
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${NC_INCLUDE_DIR})
set(LIBS ${LIBS} ${NC_LIB})

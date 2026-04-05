if ( NOT BUILD_TYPE )
  set(BUILD_TYPE "PARALLEL")
endif()
message( STATUS  "Build Type: " ${BUILD_TYPE} )

if ( NOT HYBRID )
  set(HYBRID FALSE)
else()
  message(WARNING "Compiling for hybrid openMP/MPI usage")
endif()

set(USER_Fortran_FLAGS          " -fpp -nbs -save-temps -heap-arrays -unroll -vec-threshold50 " )
set(USER_Fortran_FLAGS_RELEASE  " -Ofast -march=skylake-avx512 -axcommon-avx512,SSE4.2 -qopt-streaming-stores=always -qopt-zmm-usage=high -qopt-prefetch")
# -Ofast = -O3, -ipo, -no-prec-div, -static, -xHost; maybe check -axCORE-AVX2
set(USER_Fortran_FLAGS_DEBUG    " -g -traceback -debug all ")

set(CMAKE_Fortran_COMPILER ifx)

# compiler for parallel build and hybrid flags
if (BUILD_TYPE STREQUAL "PARALLEL" )
#   set(MPI_Fortran_COMPILER mpiifx)
  set(PARALLEL MPI_ONLY)
  add_definitions(-DUSE_MPI)

  if (HYBRID STREQUAL "TRUE" )
  endif()

endif()

if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE Release)
endif()

add_definitions(-DUSE_BLAS -DUSE_MKL)
set(BLAS_LIB 	 "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(LIBS         ${BLAS_LIB})

add_definitions(-DUSE_FFTW)
set(FFTW_INCLUDE_DIR "/dss/lrzsys/sys/spack/release/sles15.7/24.6.0/opt/skylake_avx512/fftw/3.3.10-gcc-yac4bby/include")
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${FFTW_INCLUDE_DIR})
set(LIBS         ${LIBS}         $ENV{FFTW_SHLIB})

add_definitions(-DUSE_NETCDF)
set(NC_INCLUDE_DIR     "/dss/lrzsys/sys/spack/release/sles15.7/24.6.0/opt/skylake_avx512/netcdf-fortran/4.6.1-oneapi-gp3ic7v/include")
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${NC_INCLUDE_DIR})
set(LIBS         ${LIBS}         $ENV{NETCDF_FORTRAN_SHLIB})

#message(STATUS "Build definitions : " ${})
set(GNU_SED "gsed")

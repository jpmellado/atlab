# ATLab

Tools to simulate and analyze 2D and 3D atmospheric turbulent flows. The numerical schemes are based on compact finite differences and Runge-Kutta methods. Meshes are structured and can be nonuniform. It assumes periodicity in the first coordinates. Derived from project TLab.

Some examples of applications can be found in this [website](https://jpmellado.github.io/gallery.html).

See [`doc/manual.pdf`](./doc/manual.pdf) for more information.

## Install

In order to compile the code, run the following commands:

```shell
cd ${PATH_TO_TLAB}
mkdir build
cd build
cmake ../src -DSYST={mpipc,juqueen,...} -DBUILD_TYPE={BIG,LITTLE,PARALLEL,...}
make
```
Instead of mpipc or juqueen, you have to use the corresponding file from the directory ${PATH_TO_ATLAB}/config

You can also run `./configure.sh`, which would create the different build_* directories automatically for your system if appropriately set up.

To clean the tree, simply delete the directories build*

**Prerequisites**
* cmake
* fortran compiler
* [FFTW](http://www.fftw.org/)
* Optionally, [NetCDF](https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html) for statistical data.

## Check

In order to check the code, run the following commands:

```shell
cd ${PATH_TO_DNS}
cd examples
make check BINS_PATH=${bins_LITTLE,bins_BIG,...} PRJS=${Case01,Case02,...}
make check-mpi BINS_PATH=${bins_PARALLEL,...} PRJS=${Case01,Case02,...}
```

Instead of bins_LITTLE or similar, use the corresponding directory whose executables you want to check.

If the variable PRJS is not passed, all projects in the Makefile will be checked.

To clean the examples tree, run `make clean`.

Use valgrind to check for memory leaks.

## Run

See directory [`examples`](./examples/README.md).
module TLabMPI_VARS
    use TLab_Constants, only: wp, wi
    use mpi_f08, only: MPI_Comm, MPI_Datatype
    implicit none

    type :: mpi_axis_dt
        type(MPI_Comm) :: comm
        integer :: num_processors                   ! we could name it size, in analogy to spatial grid, but this might be clearer
        integer :: rank
    end type

    type, extends(mpi_axis_dt) :: mpi_grid_dt
        type(mpi_axis_dt) axes(3)
    end type

    type(mpi_grid_dt), target :: mpiGrid
    type(mpi_axis_dt), pointer :: xMpi => null()    ! For readability in the code
    type(mpi_axis_dt), pointer :: yMpi => null()
    type(mpi_axis_dt), pointer :: zMpi => null()

    integer :: ims_err

    real(wp) :: ims_time_min, ims_time_max, ims_time_trans      ! Profiling

    type(MPI_Datatype) :: TLAB_MPI_REAL_TYPE                    ! MPI Type control

end module TLabMPI_VARS

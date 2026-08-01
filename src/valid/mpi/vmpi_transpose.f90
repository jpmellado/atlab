program vMpi_Transpose
    use TLab_Constants, only: wp, wi
    use mpi_f08
    use TLabMPI_VARS, only: mpiGrid, xMpi, yMpi, ims_err, ims_time_trans
    use TLabMPI_Transpose
    implicit none

    ! integer(wi), parameter :: nx = 8         ! full size of each linear system
    ! integer(wi), parameter :: nlines = 4    ! number of linear systems to solve, 32*2048
    integer(wi), parameter :: nx = 8192         ! full size of each linear system
    integer(wi), parameter :: nlines = 65536    ! number of linear systems to solve, 32*2048

    real(wp), allocatable :: u(:, :), f(:, :)
    real(wp), allocatable :: u_transposed(:, :)
    type(tmpi_transpose_x_dt) :: tmpi_trp

    integer nxLoc, nlinesLoc

    integer, parameter :: num_iterations = 10   ! Number of iterations to obtain a more representative time
    integer it
    real(wp) time_loc_1, time_loc_2

    ! -------------------------------------------------------------------
    call MPI_INIT(ims_err)

    mpiGrid%comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(mpiGrid%comm, mpiGrid%num_processors, ims_err)
    call MPI_COMM_RANK(mpiGrid%comm, mpiGrid%rank, ims_err)

    ! -------------------------------------------------------------------
    nxLoc = nx/mpiGrid%num_processors     ! task-local number of grid points along X
    allocate (u(nxLoc, nlines))
    allocate (f(nxLoc, nlines))

    ! -------------------------------------------------------------------
    ! Initialize
    ! call random_number(u)
    u(:, :) = mpiGrid%rank

    ! -------------------------------------------------------------------
    nlinesLoc = nlines/mpiGrid%num_processors
    allocate (u_transposed(nlinesLoc, nx))
    xMpi => mpiGrid%mpi_axis_dt
    yMpi => mpiGrid%mpi_axis_dt ! I need it for the initialization
    ! USE_MPI_DERIVED_TYPES = .false.
    call TLabMPI_Trp_Initialize()

    call tmpi_trp%initialize(nxLoc, nlines)

    ims_time_trans = 0.0_wp
    time_loc_1 = MPI_WTIME()

    if (mpiGrid%num_processors > 1) then
        do it = 1, num_iterations
            ! print *, "org", mpiGrid%rank, u
            call tmpi_trp%forward(u(:, 1), u_transposed(:, 1))
            ! print *, "trp", mpiGrid%rank, u_transposed

            call tmpi_trp%backward(u_transposed(:, 1), f(:, 1))
            ! print *, mpiGrid%rank, maxval(abs(f - u))
        end do
    end if
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm.'
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
        ! print *, 'Communication time in processor with rank 0 (seconds): ', ims_time_trans
        ! print *, 'Serial time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1 - ims_time_trans
    end if

    ! -------------------------------------------------------------------
    call MPI_FINALIZE(ims_err)

    stop

end program

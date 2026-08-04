program vMpi_Transpose
    use TLab_Constants, only: wp, wi
    use mpi_f08
    use TLabMPI_VARS, only: mpiGrid, xMpi, yMpi, ims_err, ims_time_trans
    use TLab_Memory, only: imax, jmax, kmax
    use TLabMPI_Transpose
    implicit none

    ! integer(wi), parameter :: nx = 8         ! full size of each linear system
    ! integer(wi), parameter :: nlines = 4    ! number of linear systems to solve, 32*2048
    integer(wi), parameter :: nx = 8192         ! full size of each linear system
    integer(wi), parameter :: nlines = 65536    ! number of linear systems to solve, 32*2048

    real(wp), allocatable :: u(:, :), f(:, :)
    real(wp), allocatable :: u_transposed(:, :)
    complex(wp), allocatable :: uc(:, :), fc(:, :)
    complex(wp), allocatable :: uc_transposed(:, :)
    type(tmpi_transpose_x_dt) :: tmpi_trp
    type(tmpi_transpose_block_dt) :: tmpi_trp_block

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
    !
    allocate (uc(nxLoc/2, nlines))          ! complex case
    allocate (fc(nxLoc/2, nlines))

    ! -------------------------------------------------------------------
    ! Initialize
    call random_number(u)
    !
    uc = cmplx(u(:nxloc/2, :), u(nxloc/2 + 1:, :))

    ! -------------------------------------------------------------------
    nlinesLoc = nlines/mpiGrid%num_processors
    allocate (u_transposed(nlinesLoc, nx))
    allocate (uc_transposed(nlinesLoc, nx/2))
    xMpi => mpiGrid%mpi_axis_dt
    yMpi => mpiGrid%mpi_axis_dt ! I need it for the initialization
    !
    imax = nx
    jmax = nlines
    kmax = 1
    call TLabMPI_Trp_Initialize()

    ! -------------------------------------------------------------------
    call tmpi_trp%initialize(nxLoc, nlines)

    time_loc_1 = MPI_WTIME()
    if (mpiGrid%num_processors > 1) then
        do it = 1, num_iterations
            call tmpi_trp%forward(u(:, 1), u_transposed(:, 1))
            call tmpi_trp%backward(u_transposed(:, 1), f(:, 1))
        end do
    end if
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm with MPI derived types'
        print *, 'Number of processors: ', mpiGrid%num_processors
        print *, 'Error in processors with rank 0: ', maxval(abs(f - u))
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
    end if

    ! -------------------------------------------------------------------
    call tmpi_trp_block%initialize(nxLoc, nlines, mpiGrid%mpi_axis_dt, locType=MPI_REAL8)
    !
    call tmpi_trp_block%initialize(nxLoc/2, nlines, mpiGrid%mpi_axis_dt, locType=MPI_DOUBLE_COMPLEX)

    time_loc_1 = MPI_WTIME()
    if (mpiGrid%num_processors > 1) then
        do it = 1, num_iterations
            call tmpi_trp_block%in_out(u(:, 1), u_transposed(:, 1), 'forward')
            call tmpi_trp_block%out_in(u_transposed(:, 1), f(:, 1), 'backward')
            ! call tmpi_trp_block%in_out(uc(:, 1), uc_transposed(:, 1), 'forward')
            ! call tmpi_trp_block%out_in(uc_transposed(:, 1), fc(:, 1), 'backward')
        end do
    end if
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm with raw types.'
        print *, 'Number of processors: ', mpiGrid%num_processors
        print *, 'Error in processors with rank 0: ', maxval(abs(f - u))
        ! print *, 'Error in processors with rank 0: ', maxval(abs(fc - uc))
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
    end if

    ! -------------------------------------------------------------------
    call MPI_FINALIZE(ims_err)

    stop

end program

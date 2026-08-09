program vMpi_Transpose
    use TLab_Constants, only: wp, wi
    use mpi_f08
    use TLabMPI_VARS, only: mpiGrid, xMpi, yMpi, ims_err !, ims_time_trans
    use TLabMPI_Transpose_DerivedTypes
    use TLabMPI_Transpose_X
    implicit none

    ! integer(wi), parameter :: nx = 8         ! full size of each linear system
    ! integer(wi), parameter :: nlines = 4    ! number of linear systems to solve, 32*2048
    integer(wi), parameter :: nx = 8192         ! full size of each linear system
    integer(wi), parameter :: nlines = 65536    ! number of linear systems to solve, 32*2048

    real(wp), allocatable :: u(:, :), f(:, :)
    real(wp), allocatable :: u_transposed(:, :)
    real(wp), allocatable :: wrk3d(:)
    complex(wp), allocatable :: uc(:, :), fc(:, :)
    complex(wp), allocatable :: uc_transposed(:, :)
    complex(wp), allocatable :: c_wrk3d(:)
    type(tmpi_transpose_x_dt) :: tmpi_trp
    type(tmpi_transposeX_inner_dt) :: tmpi_trpX_x
    type(tmpi_transposeX_outer_dt) :: tmpi_trpX_y

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
    !
    allocate (wrk3d(nxLoc*nlines))
    allocate (c_wrk3d(nxLoc/2*nlines))

    ! -------------------------------------------------------------------
    ! Initialize
    call random_number(u)
    !
    uc = cmplx(u(:nxloc/2, :), u(nxloc/2 + 1:, :), wp)

    ! -------------------------------------------------------------------
    nlinesLoc = nlines/mpiGrid%num_processors
    allocate (u_transposed(nlinesLoc, nx))
    allocate (uc_transposed(nlinesLoc, nx/2))
    xMpi => mpiGrid%mpi_axis_dt
    yMpi => mpiGrid%mpi_axis_dt ! I need it for the initialization
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
    call tmpi_trpX_x%initialize(nxLoc, nlines, mpiGrid%mpi_axis_dt, locType=MPI_REAL8)

    time_loc_1 = MPI_WTIME()
    if (mpiGrid%num_processors > 1) then
        do it = 1, num_iterations
            call tmpi_trpX_x%forward(u(:, 1), u_transposed(:, 1), wrk3d)
            call tmpi_trpX_x%backward(u_transposed(:, 1), f(:, 1), wrk3d)
        end do
    end if
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm with raw real types.'
        print *, 'Number of processors: ', mpiGrid%num_processors
        print *, 'Error in processors with rank 0: ', maxval(abs(f - u))
        ! print *, 'Error in processors with rank 0: ', maxval(abs(fc - uc))
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
    end if

    ! -------------------------------------------------------------------
    call tmpi_trpX_x%initialize(nxLoc/2, nlines, mpiGrid%mpi_axis_dt, locType=MPI_DOUBLE_COMPLEX)

    time_loc_1 = MPI_WTIME()
    if (mpiGrid%num_processors > 1) then
        do it = 1, num_iterations
            call tmpi_trpX_x%forward(uc(:, 1), uc_transposed(:, 1), c_wrk3d)
            call tmpi_trpX_x%backward(uc_transposed(:, 1), fc(:, 1), c_wrk3d)
        end do
    end if
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm with raw complex types.'
        print *, 'Number of processors: ', mpiGrid%num_processors
        print *, 'Error in processors with rank 0: ', maxval(abs(fc - uc))
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
    end if

    ! check inplace routines
    time_loc_1 = MPI_WTIME()
    if (mpiGrid%num_processors > 1) then
        do it = 1, num_iterations
            fc = uc
            call tmpi_trpX_x%forward(fc(:, 1), c_wrk3d)
            call tmpi_trpX_x%backward(fc(:, 1), c_wrk3d)
        end do
    end if
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm with raw complex types (in-place).'
        print *, 'Number of processors: ', mpiGrid%num_processors
        print *, 'Error in processors with rank 0: ', maxval(abs(fc - uc))
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
    end if

    ! -------------------------------------------------------------------
    call tmpi_trpX_y%initialize(nxLoc/2, nlines, mpiGrid%mpi_axis_dt, locType=MPI_DOUBLE_COMPLEX)

    time_loc_1 = MPI_WTIME()
    if (mpiGrid%num_processors > 1) then
        do it = 1, num_iterations
            call tmpi_trpX_y%forward(uc(:, 1), uc_transposed(:, 1), c_wrk3d)
            call tmpi_trpX_y%backward(uc_transposed(:, 1), fc(:, 1), c_wrk3d)
        end do
    end if
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm with raw y-complex types.'
        print *, 'Number of processors: ', mpiGrid%num_processors
        print *, 'Error in processors with rank 0: ', maxval(abs(fc - uc))
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
    end if

    ! -------------------------------------------------------------------
    call MPI_FINALIZE(ims_err)

    stop

end program

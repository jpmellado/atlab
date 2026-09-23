program vMpi_Thomas3_Scaling
    use TLab_Constants, only: wp, wi
    use mpi_f08
    use TLabMPI_VARS, only: mpiGrid, xMpi, ims_err, ims_time_trans
    use TLabMPI_PROCS, only: TLabMPI_Halos_X
    use TLab_Grid, only: axis_dt
    use FDM_Derivative_1order
    use FDM_Derivative_MPISplit
    use FDM_Base, only: FDM_COM6_DIRECT
    use FDM, only: FDM_CreatePlan_Der1
    use TLab_Arrays, only: wrk2d

    implicit none

    ! integer(wi), parameter :: nx = 4096         ! full size of each linear system
    ! integer(wi), parameter :: nlines = 32768    ! number of linear systems to solve, 32*1024

    integer(wi), parameter :: nx = 8192         ! full size of each linear system
    integer(wi), parameter :: nlines = 65536    ! number of linear systems to solve, 32*2048

    ! integer(wi), parameter :: nx = 16384        ! full size of each linear system
    ! integer(wi), parameter :: nlines = 131072   ! number of linear systems to solve, 32*4096

    ! integer(wi), parameter :: nx = 32768        ! full size of each linear system
    ! integer(wi), parameter :: nlines = 262144   ! number of linear systems to solve, 32*8192

    integer(wi), parameter :: batchsize = 64     ! # of nlines that are solved together, to test cache
    ! integer(wi), parameter :: batchsize = nlines    ! no locality

    type(axis_dt) :: x
    real(wp), allocatable :: u(:, :, :)         ! numerical solution of A u = f
    real(wp), allocatable :: f(:, :, :)         ! forcing
    real(wp) :: halos(batchsize, 2*3)

    class(der_dt), allocatable :: fdm_der1
    type(der_periodic_mpisplit) :: fdm_der1_split

    integer k, np, it, ib, nxLoc

    integer, parameter :: num_iterations = 10   ! Number of iterations to obtain a more representative time
    real(wp) time_loc_1, time_loc_2

    ! -------------------------------------------------------------------
    call MPI_INIT(ims_err)

    mpiGrid%comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(mpiGrid%comm, mpiGrid%num_processors, ims_err)
    call MPI_COMM_RANK(mpiGrid%comm, mpiGrid%rank, ims_err)

    ! -------------------------------------------------------------------
    ! Initialize
    x%size = nx
    x%scale = 1.0_wp
    ! x%periodic = .false.
    x%periodic = .true.
    x%uniform = .true.
    allocate (x%nodes(nx))
    x%nodes = [(real(k - 1, wp)/real(nx, wp)*x%scale, k=1, nx)]

    call FDM_CreatePlan_Der1(x, fdm_der1, FDM_COM6_DIRECT)

    select type (fdm_der1)
    type is (der1_periodic)
        call fdm_der1_split%initialize(fdm_der1, mpiGrid%mpi_axis_dt)
    end select

    xMpi => mpiGrid%mpi_axis_dt ! I need it in TLabMPI_Halos_X
    np = size(fdm_der1_split%rhs, 2)/2

    ! -------------------------------------------------------------------
    nxLoc = nx/mpiGrid%num_processors     ! task-local number of grid points along X
    allocate (u(batchsize, nxLoc, nlines/batchsize))
    allocate (f(batchsize, nxLoc, nlines/batchsize))

    allocate (wrk2d(batchsize, 2))

    f = 1.0_wp            ! forcing

    if (mpiGrid%rank == 0) then
    print *, new_line('a'), 'Solving ', nlines, ' systems of size ', nx, 'in batches of ', batchsize, ' over ', mpiGrid%num_processors, ' processors.'
    end if

    ! -------------------------------------------------------------------
    ! Solve using splitting algorithm

    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    ims_time_trans = 0.0_wp
    time_loc_1 = MPI_WTIME()
    do it = 1, num_iterations
        do ib = 1, nlines/batchsize
            call TLabMPI_Halos_X(f(1:batchsize*nxLoc, 1, ib), batchsize, np, halos(:, 1), halos(:, np + 1))
            call fdm_der1_split%compute(batchsize, f(:, :, ib), halos(:, 1:np), halos(:, np + 1:np + np), u(:, :, ib))
        end do
    end do
    time_loc_2 = MPI_WTIME()

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Splitting algorithm.'
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
        print *, 'Communication time in processor with rank 0 (seconds): ', ims_time_trans
        print *, 'Serial time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1 - ims_time_trans
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    ! -------------------------------------------------------------------
    call MPI_FINALIZE(ims_err)

    stop

end program

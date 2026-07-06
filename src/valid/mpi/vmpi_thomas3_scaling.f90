program vMpi_Thomas3_Scaling
    use TLab_Constants, only: wp, wi
    use mpi_f08
    use TLabMPI_VARS, only: mpiGrid, xMpi, yMpi, ims_err, ims_time_trans
    use Thomas_Parallel
    use Thomas_Circulant
    use TLabMPI_Transpose
    implicit none
    
    integer(wi), parameter :: nd = 3            ! number of diagonals
    ! integer(wi), parameter :: nx = 4096         ! full size of each linear system
    ! integer(wi), parameter :: nlines = 32768    ! number of linear systems to solve, 32*1024
    
    integer(wi), parameter :: nx = 8192         ! full size of each linear system
    integer(wi), parameter :: nlines = 65536    ! number of linear systems to solve, 32*2048
    
    ! integer(wi), parameter :: nx = 16384        ! full size of each linear system
    ! integer(wi), parameter :: nlines = 131072   ! number of linear systems to solve, 32*4096
    
    ! integer(wi), parameter :: nx = 32768        ! full size of each linear system
    ! integer(wi), parameter :: nlines = 262144   ! number of linear systems to solve, 32*8192

    real(wp) :: lhs(nx, nd)                     ! Diagonals of system matrix A
    real(wp), allocatable :: u(:, :)            ! numerical solution of A u = f
    real(wp), allocatable :: f(:, :)            ! forcing
    real(wp) :: wrk2d(nlines, 2)
    
    real(wp), allocatable :: u_transposed(:, :)
    type(tmpi_transpose_y_dt) :: tmpi_trp
    
    type(thomas_parallel_dt) split_mpi
    type(thomas_circulant_dt) :: thomas_circulant1
    
    integer k, np, it, nxLoc, nlinesLoc
    
    ! integer :: nseed
    ! integer, allocatable :: seed(:)
    
    integer, parameter :: num_iterations = 10   ! Number of iterations to obtain a more representative time
    real(wp) time_loc_1, time_loc_2
    
    ! -------------------------------------------------------------------
    call MPI_INIT(ims_err)
    
    mpiGrid%comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(mpiGrid%comm, mpiGrid%num_processors, ims_err)
    call MPI_COMM_RANK(mpiGrid%comm, mpiGrid%rank, ims_err)
    
    ! -------------------------------------------------------------------
    ! random number initialization for reproducibility
    ! from https://masuday.github.io/fortran_tutorial/random.html
    ! call random_seed(size=nseed)
    ! allocate (seed(nseed))
    ! ! call random_seed(get=seed)
    ! ! print *, seed
    ! seed = 123456789    ! putting arbitrary seed to all elements
    ! call random_seed(put=seed)
    ! ! call random_seed(get=seed)
    ! ! print *, seed
    ! deallocate (seed)
    
    ! -------------------------------------------------------------------
    ! Initialize
    ! call random_number(lhs)     ! diagonals in matrix A
    ! lhs(:,1) = 2.0_wp/11.0_wp   ! second order derivative
    ! lhs(:,2) = 1.0_wp
    ! lhs(:,3) = 2.0_wp/11.0_wp
    lhs(:, 1) = 1.0_wp/3.0_wp   ! first order derivative
    lhs(:, 2) = 1.0_wp
    lhs(:, 3) = 1.0_wp/3.0_wp
    
    ! -------------------------------------------------------------------
    nxLoc = nx/mpiGrid%num_processors     ! task-local number of grid points along X
    allocate (u(nlines, nxLoc))
    allocate (f(nlines, nxLoc))
    
    ! call random_number(f)       ! forcing
    f(:, :) = 1.0_wp            ! forcing
    
    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Solving ', nlines, ' systems of size ', nx, ' over ', mpiGrid%num_processors, ' processors.'
    end if
    
    ! -------------------------------------------------------------------
    ! Solve using splitting algorithm
    np = mpiGrid%num_processors     ! for clarity below
    call split_mpi%initialize(lhs, &
    [(k, k=nx/np, nx - nx/np, nx/np)], &
    block_id=mpiGrid%rank + 1, &
    circulant=.true.)
    split_mpi%mpi = mpiGrid%mpi_axis_dt
    
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    
    ims_time_trans = 0.0_wp
    time_loc_1 = MPI_WTIME()
    do it = 1, num_iterations
        u(:, :) = f(:, :)
        call split_mpi%SolveL(u)
        call split_mpi%SolveU(u)
        call split_mpi%reduce(u, wrk2d(:, 1), wrk2d(:, 2))
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
    ! Solve using transpose algorithm
    nlinesLoc = nlines/mpiGrid%num_processors
    
    call thomas_circulant1%initialize(lhs(:, 1:nd))
    
    allocate (u_transposed(nlinesLoc, nx))
    yMpi => mpiGrid%mpi_axis_dt
    xMpi => mpiGrid%mpi_axis_dt ! I need it for the initialization
    call TLabMPI_Trp_Initialize()
    
    call tmpi_trp%initialize(nxLoc, nlines)
    
    ims_time_trans = 0.0_wp
    time_loc_1 = MPI_WTIME()
    do it = 1, num_iterations
        u(:, :) = f(:, :)
        call tmpi_trp%forward(u(:, 1), u_transposed(:, 1))
        
        call thomas_circulant1%solveL(u_transposed)
        call thomas_circulant1%solveU(u_transposed)
        call thomas_circulant1%reduce(u_transposed, wrk2d(:, 1))
        
        call tmpi_trp%backward(u_transposed(:, 1), u(:, 1))
    end do
    time_loc_2 = MPI_WTIME()
    
    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Transpose algorithm.'
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
        print *, 'Communication time in processor with rank 0 (seconds): ', ims_time_trans
        print *, 'Serial time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1 - ims_time_trans
    end if
    
    ! -------------------------------------------------------------------
    call MPI_FINALIZE(ims_err)
    
    stop
    
end program

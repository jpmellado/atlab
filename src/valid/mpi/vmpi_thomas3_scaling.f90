program vMpi_Thomas3_Scaling
    use TLab_Constants, only: wp, wi
    use mpi_f08
    use TLabMPI_VARS, only: mpiGrid, ims_err, ims_time_trans
    use Thomas_Split
    implicit none

    integer(wi), parameter :: nd = 3            ! number of diagonals
    ! integer(wi), parameter :: nx = 4096         ! full size of each linear system
    ! integer(wi), parameter :: nlines = 32768    ! number of linear systems to solve, 32*1024

    integer(wi), parameter :: nx = 8192         ! full size of each linear system
    integer(wi), parameter :: nlines = 65536    ! number of linear systems to solve, 32*2048

    ! integer(wi), parameter :: nx = 16382        ! full size of each linear system
    ! integer(wi), parameter :: nlines = 131072   ! number of linear systems to solve, 32*4096

    real(wp) :: lhs(nx, nd)                     ! Diagonals of system matrix A
    real(wp), allocatable :: u(:, :)            ! numerical solution of A u = f
    real(wp), allocatable :: f(:, :)            ! forcing
    ! real(wp), allocatable :: u_a(:, :)          ! analytical solution
    real(wp) :: wrk2d(nlines, 2)

    type(thomas_split_dt) split_mpi

    integer k, np, it, nxLoc

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

    np = mpiGrid%num_processors     ! for clarity below
    call split_mpi%initialize(lhs, &
                              [(k, k=nx/np, nx, nx/np)], &
                              block_id=mpiGrid%rank + 1, &
                              circulant=.true.)
    split_mpi%mpi = mpiGrid%mpi_axis_dt

    ! -------------------------------------------------------------------
    nxLoc = nx/mpiGrid%num_processors     ! task-local number of grid points along X
    ! allocate(u_a(nlines, nxLoc))
    allocate (u(nlines, nxLoc))
    allocate (f(nlines, nxLoc))

    ! call random_number(f)       ! forcing
    f(:, :) = 1.0_wp            ! forcing
    ! -------------------------------------------------------------------
    ! Solve and reduce
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
        print *, 'Solving ', nlines, ' systems of size ', nx, ' over ', mpiGrid%num_processors, ' processors.'
        print *, 'Elapsed time in processor with rank 0 (seconds): ', time_loc_2 - time_loc_1
        ! print *, 'Communication time in processor with rank 0 (seconds): ', ims_time_trans
    end if
    ! call check(u_a, u)

    call MPI_FINALIZE(ims_err)

    stop

    ! ###################################################################
! contains
!     subroutine check(u, u_ref, name)
!         real(wp), intent(in) :: u(:, :), u_ref(:, :)
!         character(len=*), optional :: name

!         real(wp) dummy, error_l2, error_max
!         integer(wi) i, l

!         if (present(name)) then
!             open (20, file=name)
!         end if
!         error_l2 = 0.0_wp
!         error_max = 0.0_wp
!         dummy = 0.0_wp
!         do i = 1, size(u, 2)
!             do l = 1, size(u, 1)
!                 if (present(name)) then
!                     write (20, 1000) u(l, i), u_ref(l, i), u(l, i) - u_ref(l, i)
!                 end if
!                 dummy = dummy + u_ref(l, i)*u_ref(l, i)
!                 error_l2 = error_l2 + (u_ref(l, i) - u(l, i))**2
!                 error_max = max(error_max, abs(u_ref(l, i) - u(l, i)))
!             end do
!         end do
!         if (present(name)) then
!             close (20)
!         end if

!         write (*, *) 'Solution L2-norm ...........:', sqrt(dummy)/real(nlines, wp)
!         if (dummy == 0.0_wp) return
!         write (*, *) 'Relative Error L2-norm .....:', sqrt(error_l2)/sqrt(dummy)
!         write (*, *) 'Relative Error Linf-norm ...:', error_max/abs(maxval(u_ref))

!         return
! 1000    format(5(1x, e12.5))
!     end subroutine check

end program

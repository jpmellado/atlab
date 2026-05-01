program vThomas3_Parallel
    use TLab_Constants, only: wp, wi, BCS_NONE
    use Thomas
    use Thomas_Circulant
    use Thomas_Parallel
    use TLab_Arrays, only: wrk2d
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: mpiGrid, ims_err
#endif
    implicit none

    integer(wi), parameter :: nlines = 32
    integer(wi), parameter :: nsize = 1024
    integer(wi), parameter :: nd = 3
    integer(wi) n

    real(wp) :: rhs(nsize, nd), lhs(nsize, nd)
    real(wp) :: u(nlines, nsize), u_loc(nlines, nsize), f(nlines, nsize)

    integer :: nseed
    integer, allocatable :: seed(:)

    logical, parameter :: circulant = .true.

    type(thomas_dt) :: thomas1
    type(thomas_circulant_dt) :: thomas_circulant1
#ifdef USE_MPI
    type(thomas_parallel_dt) split_mpi
#endif
    type(thomas_parallel_dt), allocatable :: thomas_parallel1(:)  ! for testing in serial
    type(data_dt), allocatable :: data(:)
    target u_loc

    integer(wi) k, kk, np

    ! -------------------------------------------------------------------
#ifdef USE_MPI
    call MPI_INIT(ims_err)

    mpiGrid%comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(mpiGrid%comm, mpiGrid%num_processors, ims_err)
    call MPI_COMM_RANK(mpiGrid%comm, mpiGrid%rank, ims_err)
#endif

#ifdef USE_MPI
    np = mpiGrid%num_processors     ! for clarity below
#else
    np = 4
#endif

    ! -------------------------------------------------------------------
    ! random number initialization for reproducibility
    ! from https://masuday.github.io/fortran_tutorial/random.html
    call random_seed(size=nseed)
    allocate (seed(nseed))
    ! call random_seed(get=seed)
    ! print *, seed
    seed = 123456789    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    ! call random_seed(get=seed)
    ! print *, seed
    deallocate (seed)

    ! -------------------------------------------------------------------
    ! generate random data for f = Au
    call random_number(lhs)     ! diagonals in matrix A
    ! lhs(:,1) = 2.0_wp/11.0_wp   ! second order derivative
    ! lhs(:,2) = 1.0_wp
    ! lhs(:,3) = 2.0_wp/11.0_wp
    ! lhs(:,1) = 1.0_wp/3.0_wp   ! first order derivative
    ! lhs(:,2) = 1.0_wp
    ! lhs(:,3) = 1.0_wp/3.0_wp
    call random_number(u)       ! solution u

    n = 1
    f(:, n) = lhs(n, 2)*u(:, n) + lhs(n, 3)*u(:, n + 1)
    if (circulant) f(:, n) = f(:, n) + lhs(n, 1)*u(:, nsize)
    do n = 2, nsize - 1
        f(:, n) = lhs(n, 1)*u(:, n - 1) + lhs(n, 2)*u(:, n) + lhs(n, 3)*u(:, n + 1)
    end do
    f(:, n) = lhs(n, 1)*u(:, n - 1) + lhs(n, 2)*u(:, n)
    if (circulant) f(:, n) = f(:, n) + lhs(n, 3)*u(:, 1)

    ! -------------------------------------------------------------------
#ifdef USE_MPI
    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Running in parallel. Processor 0 doing the serial version.'

#endif

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Standard Thomas algorithm'

        u_loc(:, :) = f(:, :)   ! f = Au
        if (circulant) then
            call thomas_circulant1%initialize(lhs(:, 1:nd))

            allocate (wrk2d(nlines, 1))     ! needed in thomas_circulant
            call thomas_circulant1%solveL(u_loc)
            call thomas_circulant1%solveU(u_loc)
            call thomas_circulant1%reduce(u_loc, wrk2d(:, 1))

        else
            call thomas1%initialize(lhs(:, 1:nd))

            call thomas1%solveL(u_loc)
            call thomas1%solveU(u_loc)

        end if

        call check(u_loc, u, 'linear.dat')

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Splitting Thomas algorithm'

        allocate (thomas_parallel1(np))
        allocate (data(np))
        do k = 1, np
            call thomas_parallel1(k)%initialize(lhs(:, 1:3), &
                                                [(kk, kk=nsize/np, nsize, nsize/np)], &
                                                block_id=k, &
                                                circulant=circulant)

            data(k)%p => u_loc(1:nlines, thomas_parallel1(k)%nmin:thomas_parallel1(k)%nmax)
        end do

        u_loc(:, :) = f(:, :)   ! f = Au

        do k = 1, np
            call thomas_parallel1(k)%solveL(data(k)%p(:, :))
            call thomas_parallel1(k)%solveU(data(k)%p(:, :))
        end do
        call ThomasSplit_3_Reduce_Serial(thomas_parallel1, data)

        call check(u_loc, u, 'linear.dat')

#ifdef USE_MPI
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Parallel version.'
        print *, new_line('a'), 'Splitting Thomas algorithm'
    end if

    call split_mpi%initialize(lhs, &
                              [(kk, kk=nsize/np, nsize, nsize/np)], &
                              block_id=mpiGrid%rank + 1, &
                              circulant=circulant)
    split_mpi%mpi = mpiGrid%mpi_axis_dt

    u_loc(:, :) = f(:, :)   ! Each processor will only see its part of the array

    ! Solve and reduce
    if (allocated(wrk2d)) deallocate (wrk2d)
    allocate (wrk2d(nlines, 2))
    call split_mpi%SolveL(u_loc(1:nlines, split_mpi%nmin:split_mpi%nmax))
    call split_mpi%SolveU(u_loc(1:nlines, split_mpi%nmin:split_mpi%nmax))
    call split_mpi%reduce(u_loc(1:nlines, split_mpi%nmin:split_mpi%nmax), wrk2d(:, 1), wrk2d(:, 2))

    ! each processor checks its part
    call check(u_loc(1:nlines, split_mpi%nmin:split_mpi%nmax), &
               u(1:nlines, split_mpi%nmin:split_mpi%nmax))

    call MPI_FINALIZE(ims_err)

#endif

    stop

    ! ###################################################################
contains
    subroutine check(u, u_ref, name)
        real(wp), intent(in) :: u(:, :), u_ref(:, :)
        character(len=*), optional :: name

        real(wp) dummy, error_l2, error_max
        integer(wi) i, l

        if (present(name)) then
            open (20, file=name)
        end if
        error_l2 = 0.0_wp
        error_max = 0.0_wp
        dummy = 0.0_wp
        do i = 1, size(u, 2)
            do l = 1, size(u, 1)
                if (present(name)) then
                    write (20, 1000) u(l, i), u_ref(l, i), u(l, i) - u_ref(l, i)
                end if
                dummy = dummy + u_ref(l, i)*u_ref(l, i)
                error_l2 = error_l2 + (u_ref(l, i) - u(l, i))**2
                error_max = max(error_max, abs(u_ref(l, i) - u(l, i)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Solution L2-norm ...........:', sqrt(dummy)/real(nlines, wp)
        if (dummy == 0.0_wp) return
        write (*, *) 'Relative Error L2-norm .....:', sqrt(error_l2)/sqrt(dummy)
        write (*, *) 'Relative Error Linf-norm ...:', error_max/abs(maxval(u_ref))

        return
1000    format(5(1x, e12.5))
    end subroutine check

end program vThomas3_Parallel

program vThomas3_Split
    use TLab_Constants, only: wp, wi, BCS_NONE
    use Thomas
    use Thomas_Circulant
    ! use Thomas_Split, only: thomas3_split_dt, Thomas_Split_3_Initialize
    use Thomas_Split_X !, only: thomas_split_dt
    use TLab_Arrays, only: wrk2d
#ifdef USE_MPI
    use mpi_f08
    use Thomas_Split, only: thomas3_split_dt, ThomasSplit_3_Reduce_MPI
    use TLabMPI_VARS, only: mpiGrid
#endif
    implicit none

    integer(wi), parameter :: nlines = 32
    integer(wi), parameter :: nsize = 1024
    integer(wi), parameter :: nd = 3
    integer(wi) n

    real(wp) :: rhs(nsize, nd), lhs(nsize, nd)
    real(wp) :: u(nlines, nsize), u_loc(nlines, nsize), f(nlines, nsize)

    integer(wi) k
    integer, parameter :: nblocks = 4       ! number of blocks
    integer, parameter :: points(1:nblocks) = [(k, k=nsize/nblocks, nsize, nsize/nblocks)]

    ! type(thomas3_split_dt) split(nblocks)
    type(thomas_split_dt) split_X(nblocks)
    type(data_dt) data(nblocks)
    target u_loc

    integer :: nseed
    integer, allocatable :: seed(:)

    logical, parameter :: periodic = .true.

    type(thomas_dt) :: thomas
    type(thomas_circulant_dt) :: thomas_circulant

#ifdef USE_MPI
    integer ims_err
    ! type(thomas3_split_dt) split_mpi
    type(thomas_split_dt) split_mpi_X
#endif

    ! -------------------------------------------------------------------
#ifdef USE_MPI
    call MPI_INIT(ims_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiGrid%num_processors, ims_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpiGrid%rank, ims_err)
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
    call random_number(rhs)     ! diagonals in matrix A
    call random_number(u)       ! solution u

    n = 1
    f(:, n) = rhs(n, 2)*u(:, n) + rhs(n, 3)*u(:, n + 1)
    if (periodic) f(:, n) = f(:, n) + rhs(n, 1)*u(:, nsize)
    do n = 2, nsize - 1
        f(:, n) = rhs(n, 1)*u(:, n - 1) + rhs(n, 2)*u(:, n) + rhs(n, 3)*u(:, n + 1)
    end do
    f(:, n) = rhs(n, 1)*u(:, n - 1) + rhs(n, 2)*u(:, n)
    if (periodic) f(:, n) = f(:, n) + rhs(n, 3)*u(:, 1)

    ! -------------------------------------------------------------------
#ifdef USE_MPI
    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Running in parallel. Processor 0 doing the serial version.'

        if (nblocks /= mpiGrid%num_processors) then
            print *, 'Number of blocks must equal number of processors.'
            call MPI_FINALIZE(ims_err)
            stop
        end if
#endif

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Standard Thomas algorithm'

        lhs(:, :) = rhs(:, :)
        u_loc(:, :) = f(:, :)   ! f = Au
        if (periodic) then
            call thomas_circulant%initialize(lhs(:, 1:nd))

            allocate (wrk2d(nlines, 1))     ! needed in thomas_circulant
            call thomas_circulant%solveL(u_loc)
            call thomas_circulant%solveU(u_loc)
            call thomas_circulant%reduce(u_loc)

        else
            call thomas%initialize(lhs(:, 1:nd))

            call thomas%solveL(u_loc)
            call thomas%solveU(u_loc)

        end if

        call check(u_loc, u, 'linear.dat')

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Splitting Thomas algorithm'

        do k = 1, nblocks
            ! split(k)%circulant = periodic
            ! split(k)%block_id = k

            ! lhs(:, :) = rhs(:, :)
            ! call Thomas_Split_3_Initialize(lhs(:, 1:1), lhs(:, 2:3), &
            !                                points, split(k))
            ! data(k)%p => u_loc(1:nlines, split(k)%nmin:split(k)%nmax)

            !
            split_X(k)%circulant = periodic
            split_X(k)%block_id = k

            lhs(:, :) = rhs(:, :)
            call split_X(k)%initialize(lhs(:, 1:3), points)

            data(k)%p => u_loc(1:nlines, split_X(k)%nmin:split_X(k)%nmax)
        end do

        u_loc(:, :) = f(:, :)   ! f = Au

        do k = 1, nblocks
            ! call Thomas3_SolveL(split(k)%lhs(:, 1:1), data(k)%p(:, :))
            ! call Thomas3_SolveU(split(k)%lhs(:, 2:3), data(k)%p(:, :))
            call split_X(k)%solveL(data(k)%p(:, :))
            call split_X(k)%solveU(data(k)%p(:, :))
        end do
        ! call ThomasSplit_3_Reduce_Serial(split, data)
        call ThomasSplit_3_Reduce_Serial(split_X, data)

        call check(u_loc, u, 'linear.dat')

#ifdef USE_MPI
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    if (mpiGrid%rank == 0) then
        print *, new_line('a'), 'Parallel version.'
        print *, new_line('a'), 'Splitting Thomas algorithm'
    end if

    split_mpi_X%circulant = periodic
    split_mpi_X%block_id = mpiGrid%rank + 1
    split_mpi_X%mpi%comm = mpiGrid%comm
    split_mpi_X%mpi%rank = mpiGrid%rank
    split_mpi_X%mpi%num_processors = mpiGrid%num_processors

    lhs(:, :) = rhs(:, :)
    call split_mpi_X%initialize(lhs, points)

    u_loc(:, :) = f(:, :)   ! Each processor will only see its part of the array

    ! Solve and reduce
    if (allocated(wrk2d)) deallocate (wrk2d)
    allocate (wrk2d(nlines, 2))
    call split_mpi_X%SolveL(u_loc(1:nlines, split_mpi_X%nmin:split_mpi_X%nmax))
    call split_mpi_X%SolveU(u_loc(1:nlines, split_mpi_X%nmin:split_mpi_X%nmax))
    call split_mpi_X%reduce(u_loc(1:nlines, split_mpi_X%nmin:split_mpi_X%nmax), wrk2d(:, 1), wrk2d(:, 2))

    ! each processor checks its part
    call check(u_loc(1:nlines, split_mpi_X%nmin:split_mpi_X%nmax), &
               u(1:nlines, split_mpi_X%nmin:split_mpi_X%nmax))

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

end program vThomas3_Split

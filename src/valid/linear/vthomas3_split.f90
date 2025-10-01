program vThomas3_Split
    use TLab_Constants, only: wp, wi, BCS_NONE
    use Thomas3
    use Thomas3_Split
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_pro, ims_npro
#endif
    implicit none

    integer(wi), parameter :: nlines = 32
    integer(wi), parameter :: nsize = 1024
    integer(wi), parameter :: nd = 3
    integer(wi) n

    real(wp) :: lhs(nsize, nd), lhs_loc(nsize, nd)
    real(wp) :: u(nlines, nsize), u_loc(nlines, nsize), f(nlines, nsize)
    real(wp) :: z(nsize), wrk(nlines)       ! for circulant case

    integer(wi) k
    integer, parameter :: nblocks = 4       ! number of blocks
    integer, parameter :: points(1:nblocks) = [(k, k=nsize/nblocks, nsize, nsize/nblocks)]

    type(thomas3_split_dt) split(nblocks)
    type(data_dt) data(nblocks)
    target u_loc

    integer :: nseed
    integer, allocatable :: seed(:)

    logical, parameter :: periodic = .true.

#ifdef USE_MPI
    integer ims_err
    type(thomas3_split_dt) split_mpi
    real(wp), allocatable :: wrk2d(:, :)
#endif

    ! -------------------------------------------------------------------
#ifdef USE_MPI
    call MPI_INIT(ims_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ims_npro, ims_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)
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
    ! generate system
    call random_number(lhs)

    call random_number(u)

    n = 1
    f(:, n) = lhs(n, 2)*u(:, n) + lhs(n, 3)*u(:, n + 1)
    if (periodic) f(:, n) = f(:, n) + lhs(n, 1)*u(:, nsize)
    do n = 2, nsize - 1
        f(:, n) = lhs(n, 1)*u(:, n - 1) + lhs(n, 2)*u(:, n) + lhs(n, 3)*u(:, n + 1)
    end do
    f(:, n) = lhs(n, 1)*u(:, n - 1) + lhs(n, 2)*u(:, n)
    if (periodic) f(:, n) = f(:, n) + lhs(n, 3)*u(:, 1)

    ! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_pro == 0) then
        print *, new_line('a'), 'Running in parallel. Processor 0 doing the serial version.'

        if (nblocks /= ims_npro) then
            print *, 'Number of blocks must equal number of processors.'
            call MPI_FINALIZE(ims_err)
            stop
        end if
#endif

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Standard Thomas algorithm'

        u_loc(:, :) = f(:, :)

        lhs_loc = lhs
        if (periodic) then
            call Thomas3C_SMW_LU(lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3), z)
            call Thomas3C_SMW_Solve(lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3), z, u_loc, wrk)
        else
            call Thomas3_FactorLU(nsize, lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3))
            call Thomas3_SolveLU(nsize, nlines, lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3), u_loc)
        end if

        call check(u_loc, u, 'linear.dat')

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Splitting Thomas algorithm'

        u_loc(:, :) = f(:, :)

        do k = 1, nblocks
            split(k)%circulant = periodic
            split(k)%block_id = k

            lhs_loc = lhs
            call Thomas3_Split_Initialize(lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3), &
                                          points, split(k))
            data(k)%p => u_loc(1:nlines, split(k)%nmin:split(k)%nmax)

        end do

        ! call Thomas3_Split_Solve_Serial_Old(split, data)
        call Thomas3_Split_Solve_Serial(split, data)

        call check(u_loc, u, 'linear.dat')

#ifdef USE_MPI
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    if (ims_pro == 0) then
        print *, new_line('a'), 'Parallel version.'
        print *, new_line('a'), 'Splitting Thomas algorithm'
    end if

    split_mpi%circulant = periodic
    split_mpi%block_id = ims_pro + 1
    split_mpi%communicator = MPI_COMM_WORLD
    split_mpi%rank = ims_pro
    split_mpi%n_ranks = ims_npro

    lhs_loc = lhs
    call Thomas3_Split_Initialize(lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3), &
                                  points, split_mpi)

    u_loc(:, :) = f(:, :)   ! Each processor will only see its part of the array

    ! Solve and reduce
    allocate (wrk2d(nlines, 2))
    call Thomas3_Split_Solve_MPI(split_mpi, u_loc(1:nlines, split_mpi%nmin:split_mpi%nmax), wrk2d(:, 1), wrk2d(:, 2))

    ! each processor checks its part
    call check(u_loc(1:nlines, split_mpi%nmin:split_mpi%nmax), u(1:nlines, split_mpi%nmin:split_mpi%nmax))

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

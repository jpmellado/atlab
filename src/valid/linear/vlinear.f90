program vLinear
    use TLab_Constants, only: wp, wi, BCS_NONE
    use Thomas3
    ! use Thomas5
    use FDM_MatMul

    implicit none

    integer(wi), parameter :: len = 1
    integer(wi), parameter :: nx = 128
    integer(wi), parameter :: nd = 3

    real(wp) :: lhs(nx, nd)
    real(wp) :: u(len, nx), f(len, nx)

    integer(wi) id !, n

    integer :: nseed
    integer, allocatable :: seed(:)

    ! -------------------------------------------------------------------
    ! random number initialization for reproducibility
    ! from https://masuday.github.io/fortran_tutorial/random.html
    call random_seed(size=nseed)
    allocate (seed(nseed))
    call random_seed(get=seed)
    print *, seed
    seed = 123456789    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    call random_seed(get=seed)
    print *, seed
    deallocate (seed)

    ! -------------------------------------------------------------------
    ! generate system
    call random_number(lhs)
    do id = 1, nd               ! 1. upper diagonal is 1; see matmul_3d
        lhs(:, id) = lhs(:, id)/lhs(:, 3)
    end do
    lhs(1, 1) = 0.0_wp
    lhs(nx, 3) = 0.0_wp

    call random_number(u)
    call MatMul_3d(lhs, u, f, ibc=BCS_NONE)
    ! n = 1; f(:, n) = lhs(n, 2)*u(:, n) + lhs(n, 3)*u(:, n + 1)
    ! do n = 2, nx - 1
    !     f(:, n) = lhs(n, 1)*u(:, n - 1) + lhs(n, 2)*u(:, n) + lhs(n, 3)*u(:, n + 1)
    ! end do
    ! f(:, n) = lhs(n, 1)*u(:, n - 1) + lhs(n, 2)*u(:, n)

    ! -------------------------------------------------------------------
    ! solve using standard thomas algorithm
    call Thomas3_LU(nx, lhs(:, 1), lhs(:, 2), lhs(:, 3))
    call Thomas3_Solve(nx, len, lhs(:, 1), lhs(:, 2), lhs(:, 3), f)

    print *, 'Standard Thomas algorithm'
    call check(u, f, 'linear.dat')

    ! -------------------------------------------------------------------
    ! solve using standard splitting algorithm
    ! call SplittingThomas3_Initialize()
    ! call SplittingThomas3_Solve()

    print *, 'Splitting Thomas algorithm'
    ! call check(u, f, 'linear.dat')

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
                error_l2 = error_l2 + (u_ref(l, i) - u(l, i))**2.0_wp
                error_max = max(error_max, abs(u_ref(l, i) - u(l, i)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Solution L2-norm ...........:', sqrt(dummy)/real(len, wp)
        if (dummy == 0.0_wp) return
        write (*, *) 'Relative Error L2-norm .....:', sqrt(error_l2)/sqrt(dummy)
        write (*, *) 'Relative Error Linf-norm ...:', error_max/abs(maxval(u_ref))

        return
1000    format(5(1x, e12.5))
    end subroutine check

end program vLinear

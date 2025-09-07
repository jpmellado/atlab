program vLinSys
    use TLab_Constants, only: wp, wi, BCS_NONE
    use Thomas3
    use Thomas5
    implicit none

    integer(wi), parameter :: nlines = 1
    integer(wi), parameter :: nsize = 16
    integer, parameter :: ndmax = 7
    integer nd

    real(wp) :: lhs(nsize, ndmax + 2), lhs_loc(nsize, ndmax + 2)
    real(wp) :: u(nlines, nsize), f(nlines, nsize)
    real(wp) :: wrk(nlines)       ! for circulant case

    integer :: nseed
    integer, allocatable :: seed(:)

    character(len=32) str
    logical, parameter :: periodic = .false.

    ! ###################################################################
    ! random number initialization for reproducibility
    ! from https://masuday.github.io/fortran_tutorial/random.html
    call random_seed(size=nseed)
    allocate (seed(nseed))
    ! call random_seed(get=seed)
    ! print *, seed
    seed = 123456789    ! putting arbitrary seed to all elements
    seed = 132456789    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    ! call random_seed(get=seed)
    ! print *, seed
    deallocate (seed)

    ! ###################################################################
    ! generate system
    call random_number(lhs)

    call random_number(u)

    ! ###################################################################
    nd = 3

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve biased triadiagonal system'

    call matmul(lhs(:, 1:nd), u, f, periodic=.false.)
    lhs_loc(:, 1:nd) = lhs(:, 1:nd)
    call Thomas3_FactorLU(nsize, lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3))
    call Thomas3_SolveLU(nsize, nlines, lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3), f)
    ! call Thomas3_FactorUL(nsize, lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3))
    ! call Thomas3_SolveUL(nsize, nlines, lhs_loc(:, 1), lhs_loc(:, 2), lhs_loc(:, 3), f)

    write (str, *) nd
    call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve periodic triadiagonal system'

    call matmul(lhs(:, 1:nd), u, f, periodic=.true.)
    lhs_loc(:, 1:nd) = lhs(:, 1:nd)
    call Thomas3C_SMW_LU(lhs_loc(:, 1), &
                         lhs_loc(:, 2), &
                         lhs_loc(:, 3), &
                         lhs_loc(:, 4))
    call Thomas3C_SMW_Solve(lhs_loc(:, 1), &
                            lhs_loc(:, 2), &
                            lhs_loc(:, 3), &
                            lhs_loc(:, 4), f, wrk)

    write (str, *) nd
    call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    ! ###################################################################
    nd = 5

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve biased pentadiagonal system'

    call matmul(lhs(:, 1:nd), u, f, periodic=.false.)
    lhs_loc(:, 1:nd) = lhs(:, 1:nd)
    call Thomas5_FactorLU(nsize, &
                          lhs_loc(:, 1), &
                          lhs_loc(:, 2), &
                          lhs_loc(:, 3), &
                          lhs_loc(:, 4), &
                          lhs_loc(:, 5))
    call Thomas5_SolveLU(nsize, nlines, &
                         lhs_loc(:, 1), &
                         lhs_loc(:, 2), &
                         lhs_loc(:, 3), &
                         lhs_loc(:, 4), &
                         lhs_loc(:, 5), f)

    write (str, *) nd
    call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve periodic pentadiagonal system'

    call matmul(lhs(:, 1:nd), u, f, periodic=.true.)
    lhs_loc(:, 1:nd) = lhs(:, 1:nd)
    call Thomas5C_SMW_LU(nsize, &
                         lhs_loc(:, 1), &
                         lhs_loc(:, 2), &
                         lhs_loc(:, 3), &
                         lhs_loc(:, 4), &
                         lhs_loc(:, 5), &
                         lhs_loc(:, 6), &
                         lhs_loc(:, 7))
    call Thomas5C_SMW_Solve(nsize, nlines, &
                            lhs_loc(:, 1), &
                            lhs_loc(:, 2), &
                            lhs_loc(:, 3), &
                            lhs_loc(:, 4), &
                            lhs_loc(:, 5), &
                            lhs_loc(:, 6), &
                            lhs_loc(:, 7), f)

    write (str, *) nd
    call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    stop

    ! ###################################################################
contains
    subroutine matmul(lhs, u, f, periodic)
        real(wp), intent(in) :: lhs(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        logical, intent(in) :: periodic

        integer(wi) nsize, n
        integer(wi) ndl, idl, ic

        nsize = size(lhs, 1)    ! size of the system
        ndl = size(lhs, 2)      ! # of diagonals
        idl = ndl/2 + 1         ! index of centerline diagonal

        ! lower boundary
        do n = 1, idl - 1
            f(:, n) = lhs(n, idl)*u(:, n)
            do ic = 1, idl - 1
                f(:, n) = f(:, n) + &
                          lhs(n, idl + ic)*u(:, n + ic)
                if (n - ic >= 1) then
                    f(:, n) = f(:, n) + &
                              lhs(n, idl - ic)*u(:, n - ic)
                else
                    if (periodic) then
                        print *, n, ic, idl - ic, mod(n - ic + nsize - 1, nsize) + 1
                        f(:, n) = f(:, n) + &
                                  lhs(n, idl - ic)*u(:, mod(n - ic + nsize - 1, nsize) + 1)
                    end if
                end if
            end do
        end do

        ! interior points
        do n = idl, nsize - idl + 1
            f(:, n) = lhs(n, idl)*u(:, n)
            do ic = 1, idl - 1
                f(:, n) = f(:, n) + &
                          lhs(n, idl - ic)*u(:, n - ic) + &
                          lhs(n, idl + ic)*u(:, n + ic)
            end do
        end do

        ! upper boundary
        do n = nsize - idl + 2, nsize
            f(:, n) = lhs(n, idl)*u(:, n)
            do ic = 1, idl - 1
                f(:, n) = f(:, n) + &
                          lhs(n, idl - ic)*u(:, n - ic)
                if (n + ic <= nsize) then
                    f(:, n) = f(:, n) + &
                              lhs(n, idl + ic)*u(:, n + ic)
                else
                    if (periodic) then
                        print *, n, ic, idl + ic, mod(n + ic, nsize)
                        f(:, n) = f(:, n) + &
                                  lhs(n, idl + ic)*u(:, mod(n + ic, nsize))
                    end if
                end if
            end do
        end do

        return
    end subroutine matmul

    ! ###################################################################
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

end program vLinSys

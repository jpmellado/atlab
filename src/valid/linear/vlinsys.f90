program vLinSys
    use TLab_Constants, only: wp, wi, BCS_NONE
    use Thomas
    use Thomas3
    use Thomas5
    use Thomas7
    implicit none

    interface matmul
        procedure matmul_biased, matmul_circulant
    end interface matmul

    integer(wi), parameter :: nlines = 1
    integer(wi), parameter :: nsize = 1024
    integer, parameter :: ndmax = 7
    integer nd

    real(wp) :: lhs(nsize, ndmax + 2), lhs_loc(nsize, ndmax + 2)
    real(wp) :: u(nlines, nsize), f(nlines, nsize)
    real(wp) :: wrk(nlines)       ! for circulant case

    integer :: nseed
    integer, allocatable :: seed(:)

    character(len=32) str
    logical, parameter :: circulant = .false.

    ! ###################################################################
    ! random number initialization for reproducibility
    ! from https://masuday.github.io/fortran_tutorial/random.html
    call random_seed(size=nseed)
    allocate (seed(nseed))
    ! call random_seed(get=seed)
    ! print *, seed
    seed = 123456789    ! putting arbitrary seed to all elements
    ! seed = 132456789    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    ! call random_seed(get=seed)
    ! print *, seed
    deallocate (seed)

    ! ###################################################################
    ! generate random data for f = Au
    call random_number(lhs)     ! diagonals in matrix A
    call random_number(u)       ! solution u

    ! ###################################################################
    nd = 3

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve biased triadiagonal system'

    ! compute forcing
    call matmul(lhs(:, 1:nd), u, f, &
                rhs_b=lhs(1:nd/2, 1:nd), &
                rhs_t=lhs(nsize - nd/2 + 1:nsize, 1:nd))

    lhs_loc(:, 1:nd) = lhs(:, 1:nd)
    ! call Thomas3_FactorLU(nsize, &
    !                       lhs_loc(:, 1), &
    !                       lhs_loc(:, 2), &
    !                       lhs_loc(:, 3))
    ! call Thomas3_FactorLU_InPlace(lhs_loc(:, 1:nd/2), &
    !                               lhs_loc(:, nd/2 + 1:nd))
    call Thomas_FactorLU_InPlace(lhs_loc(:, 1:nd/2), &
                                  lhs_loc(:, nd/2 + 1:nd))

    ! call Thomas3_SolveLU(nsize, nlines, &
    !                      lhs_loc(:, 1), &
    !                      lhs_loc(:, 2), &
    !                      lhs_loc(:, 3), f)
    ! call Thomas3_SolveLU_L(lhs_loc(:, 1:nd/2), f)
    call Thomas_SolveLU_L(lhs_loc(:, 1:nd/2), f)
    ! call Thomas3_SolveLU_U(lhs_loc(:, nd/2 + 1:nd), f)
    call Thomas_SolveLU_U(lhs_loc(:, nd/2 + 1:nd), f)

    ! call Thomas3_FactorUL(nsize, &
    !                       lhs_loc(:, 1), &
    !                       lhs_loc(:, 2), &
    !                       lhs_loc(:, 3))
    ! call Thomas3_SolveUL(nsize, nlines, &
    !                      lhs_loc(:, 1), &
    !                      lhs_loc(:, 2), &
    !                      lhs_loc(:, 3), f)

    write (str, *) nd
    call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve circulant triadiagonal system'

    ! compute forcing
    call matmul(lhs(:, 1:nd), u, f)

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

    ! compute forcing
    call matmul(lhs(:, 1:nd), u, f, &
                rhs_b=lhs(1:nd/2, 1:nd), &
                rhs_t=lhs(nsize - nd/2 + 1:nsize, 1:nd))

    lhs_loc(:, 1:nd) = lhs(:, 1:nd)
    ! call Thomas5_FactorLU(nsize, &
    !                       lhs_loc(:, 1), &
    !                       lhs_loc(:, 2), &
    !                       lhs_loc(:, 3), &
    !                       lhs_loc(:, 4), &
    !                       lhs_loc(:, 5))
    ! call Thomas5_FactorLU_InPlace(lhs_loc(:, 1:nd/2), &
    !                               lhs_loc(:, nd/2 + 1:nd))
    call Thomas_FactorLU_InPlace(lhs_loc(:, 1:nd/2), &
                                 lhs_loc(:, nd/2 + 1:nd))
    ! call Thomas5_SolveLU(nsize, nlines, &
    !                      lhs_loc(:, 1), &
    !                      lhs_loc(:, 2), &
    !                      lhs_loc(:, 3), &
    !                      lhs_loc(:, 4), &
    !                      lhs_loc(:, 5), f)
    ! call Thomas5_SolveLU_L(lhs_loc(:, 1:nd/2), f)
    call Thomas_SolveLU_L(lhs_loc(:, 1:nd/2), f)
    ! call Thomas5_SolveLU_U(lhs_loc(:, nd/2 + 1:nd), f)
    call Thomas_SolveLU_U(lhs_loc(:, nd/2 + 1:nd), f)

    ! call Thomas5_FactorUL(nsize, &
    !                       lhs_loc(:, 1), &
    !                       lhs_loc(:, 2), &
    !                       lhs_loc(:, 3), &
    !                       lhs_loc(:, 4), &
    !                       lhs_loc(:, 5))
    ! call Thomas5_SolveUL(nsize, nlines, &
    !                      lhs_loc(:, 1), &
    !                      lhs_loc(:, 2), &
    !                      lhs_loc(:, 3), &
    !                      lhs_loc(:, 4), &
    !                      lhs_loc(:, 5), f)

    write (str, *) nd
    call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve circulant pentadiagonal system'

    ! compute forcing
    ! lhs(:, 1) = 1.4629948364887945e-003 ! this comes from fdm1-penta and works. Why?
    ! lhs(:, 2) = 9.0361445783137314e-003
    ! lhs(:, 3) = 1.6135972461273688e-002
    ! lhs(:, 4) = 9.0361445783132474e-003
    ! lhs(:, 5) = 1.4629948364888149e-003

    call matmul(lhs(:, 1:nd), u, f)

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

    ! ###################################################################
    nd = 7

    ! -------------------------------------------------------------------
    print *, new_line('a'), 'Solve biased heptadiagonal system'

    ! compute forcing
    call matmul(lhs(:, 1:nd), u, f, &
                rhs_b=lhs(1:nd/2, 1:nd), &
                rhs_t=lhs(nsize - nd/2 + 1:nsize, 1:nd))

    lhs_loc(:, 1:nd) = lhs(:, 1:nd)
    ! call Thomas7_FactorLU(nsize, &
    !                       lhs_loc(:, 1), &
    !                       lhs_loc(:, 2), &
    !                       lhs_loc(:, 3), &
    !                       lhs_loc(:, 4), &
    !                       lhs_loc(:, 5), &
    !                       lhs_loc(:, 6), &
    !                       lhs_loc(:, 7))
    ! call Thomas7_FactorLU_InPlace(lhs_loc(:, 1:nd/2), &
    !                               lhs_loc(:, nd/2 + 1:nd))
    call Thomas_FactorLU_InPlace(lhs_loc(:, 1:nd/2), &
                                 lhs_loc(:, nd/2 + 1:nd))
    ! call Thomas7_SolveLU(nsize, nlines, &
    !                      lhs_loc(:, 1), &
    !                      lhs_loc(:, 2), &
    !                      lhs_loc(:, 3), &
    !                      lhs_loc(:, 4), &
    !                      lhs_loc(:, 5), &
    !                      lhs_loc(:, 6), &
    !                      lhs_loc(:, 7), f)
    !   call Thomas7_SolveLU_L(lhs_loc(:, 1:nd/2), f)
    call Thomas_SolveLU_L(lhs_loc(:, 1:nd/2), f)
    ! call Thomas7_SolveLU_U(lhs_loc(:, nd/2 + 1:nd), f)
    call Thomas_SolveLU_U(lhs_loc(:, nd/2 + 1:nd), f)

    write (str, *) nd
    call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    stop

contains
    ! ###################################################################
    ! ###################################################################
    ! calculate f = A u, where A is narrow banded with diagonals given by rhs
    ! Allowing for different band size at the bottom and at the top
    subroutine matmul_biased(rhs, u, f, rhs_b, rhs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)

        integer(wi) nx, ir, nmin, nmax
        integer(wi) ndr, idr, ic
        integer nx_b, ndr_b, idr_b, nx_t, ndr_t, idr_t

        ! ###################################################################
        nx = size(rhs, 1)       ! size of the system
        ndr = size(rhs, 2)      ! # of diagonals
        idr = ndr/2 + 1         ! index of centerline diagonal

        ! -------------------------------------------------------------------
        ! lower boundary
        nx_b = size(rhs_b, 1)
        ndr_b = size(rhs_b, 2)
        idr_b = ndr_b/2 + 1

        do ir = 1, ndr_b/2
            f(:, ir) = rhs_b(ir, idr_b - ir + 1)*u(:, 1)
            do ic = 2, ndr_b/2 + ir
                f(:, ir) = f(:, ir) + &
                           rhs_b(ir, idr_b - ir + ic)*u(:, ic)
            end do
        end do

        do ir = idr_b, nx_b
            f(:, ir) = rhs_b(ir, idr_b)*u(:, ir)
            do ic = 1, ndr_b/2
                f(:, ir) = f(:, ir) + &
                           rhs_b(ir, idr_b + ic)*u(:, ir + ic) + &
                           rhs_b(ir, idr_b - ic)*u(:, ir - ic)
            end do
        end do

        ! do ir = 1, nx_b
        !     f(:, ir) = rhs_b(ir, idr_b)*u(:, ir)
        !     do ic = 1, idr_b - 1
        !         f(:, ir) = f(:, ir) + &
        !                    rhs_b(ir, idr_b + ic)*u(:, ir + ic)
        !         if (ir - ic >= 1) then
        !             ! print *, ir, ic, idr - ic, ir - ic
        !             f(:, ir) = f(:, ir) + &
        !                        rhs_b(ir, idr_b - ic)*u(:, ir - ic)
        !         end if
        !     end do
        ! end do

        nmin = nx_b + 1

        ! -------------------------------------------------------------------
        ! upper boundary
        nx_t = size(rhs_t, 1)
        ndr_t = size(rhs_t, 2)
        idr_t = ndr_t/2 + 1

        do ir = 0, ndr_t/2 - 1
            f(:, nx - ir) = rhs_t(nx_t - ir, idr_t + ir)*u(:, nx)
            do ic = 1, ndr_t/2 + ir
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs_t(nx_t - ir, idr_t + ir - ic)*u(:, nx - ic)
            end do
        end do

        do ir = ndr_t/2, nx_t - 1
            f(:, nx - ir) = rhs_t(nx_t - ir, idr_t)*u(:, nx - ir)
            do ic = 1, ndr_t/2
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs_t(nx_t - ir, idr_t - ic)*u(:, nx - ir - ic) + &
                                rhs_t(nx_t - ir, idr_t + ic)*u(:, nx - ir + ic)
            end do
        end do

        ! do ir = 0, nx_t - 1
        !     f(:, nx - ir) = rhs_t(nx_t - ir, idr_t)*u(:, nx - ir)
        !     do ic = 1, idr_t - 1
        !         f(:, nx - ir) = f(:, nx - ir) + &
        !                         rhs_t(nx_t - ir, idr_t - ic)*u(:, nx - ir - ic)
        !         if (ic - ir <= 0) then
        !             ! print *, n, ic, idr + ic, n + ic
        !             f(:, nx - ir) = f(:, nx - ir) + &
        !                             rhs_t(nx_t - ir, idr_t + ic)*u(:, nx - ir + ic)
        !         end if
        !     end do
        ! end do

        nmax = nx - nx_t

        ! -------------------------------------------------------------------
        ! interior points
        do ir = nmin, nmax
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ic)*u(:, ir - ic) + &
                           rhs(ir, idr + ic)*u(:, ir + ic)
            end do
        end do

        return
    end subroutine matmul_biased

    ! ###################################################################
    ! ###################################################################
    subroutine matmul_circulant(rhs, u, f)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)

        integer(wi) nx, ir
        integer(wi) ndr, idr, ic

        ! ###################################################################
        ndr = size(rhs, 2)      ! # of diagonals
        idr = ndr/2 + 1         ! index of centerline diagonal
        nx = size(rhs, 1)       ! size of the system

        ! -------------------------------------------------------------------
        ! lower boundary
        do ir = 1, idr - 1
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr + ic)*u(:, ir + ic) + &
                           rhs(ir, idr - ic)*u(:, mod(ir - ic + nx - 1, nx) + 1)
            end do
        end do

        ! -------------------------------------------------------------------
        ! upper boundary
        do ir = 0, idr - 2
            f(:, nx - ir) = rhs(nx - ir, idr)*u(:, nx - ir)
            do ic = 1, idr - 1
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs(nx - ir, idr - ic)*u(:, nx - ir - ic) + &
                                rhs(nx - ir, idr + ic)*u(:, mod(nx - ir + ic, nx))
            end do
        end do

        ! -------------------------------------------------------------------
        ! interior points
        do ir = idr, nx - idr + 1
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ic)*u(:, ir - ic) + &
                           rhs(ir, idr + ic)*u(:, ir + ic)
            end do
        end do

        return
    end subroutine matmul_circulant

    ! ###################################################################
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

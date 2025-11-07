program vLinSys
    use TLab_Constants, only: wp, wi, BCS_NONE
    use Matmul_Halo
    use Matmul_Halo_Thomas
    use Thomas
    use Thomas_Circulant
    implicit none

    interface matmul
        procedure matmul_biased, matmul_biased_solveL
    end interface matmul

    integer(wi), parameter :: nlines = 1
    integer(wi), parameter :: nsize = 1024
    integer, parameter :: ndmax = 7
    integer nd

    real(wp) :: rhs(nsize, ndmax), lhs(nsize, ndmax + 2)
    real(wp) :: u(nlines, nsize), f(nlines, nsize)
    real(wp) :: wrk(nlines)       ! for circulant case

    integer :: nseed
    integer, allocatable :: seed(:)

    integer ic
    integer, parameter :: list_of_cases(1:3) = [3, 5, 7]

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
    call random_number(rhs)     ! diagonals in matrix A
    call random_number(u)       ! solution u

    ! ###################################################################
    do ic = 1, size(list_of_cases)
        nd = list_of_cases(ic)

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Solve biased system, bands ', nd

        lhs(:, 1:nd) = rhs(:, 1:nd)
        call Thomas_FactorLU_InPlace(lhs(:, 1:nd/2), &
                                     lhs(:, nd/2 + 1:nd))
        ! select case (nd)
        ! case (3)
        !     call Thomas3_FactorLU_InPlace(lhs(:, 1:nd/2), &
        !                                   lhs(:, nd/2 + 1:nd))
        ! case (5)
        !     call Thomas5_FactorLU_InPlace(lhs(:, 1:nd/2), &
        !                                   lhs(:, nd/2 + 1:nd))
        ! case (7)
        !     call Thomas7_FactorLU_InPlace(lhs(:, 1:nd/2), &
        !                                   lhs(:, nd/2 + 1:nd))
        ! end select

        ! compute forcing
        call matmul(rhs=rhs(:, 1:nd), &
                    rhs_b=rhs(1:nd/2, 1:nd), &
                    rhs_t=rhs(nsize - nd/2 + 1:nsize, 1:nd), &
                    u=u, &
                    f=f, &
                    L=lhs(:, 1:nd/2))

        ! call Thomas_SolveL(lhs(:, 1:nd/2), f)
        call Thomas_SolveU(lhs(:, nd/2 + 1:nd), f)

        ! select case (nd)
        ! case (3)
        !     call Thomas3_SolveL(lhs(:, 1:nd/2), f)
        !     call Thomas3_SolveU(lhs(:, nd/2 + 1:nd), f)
        ! case (5)
        !     call Thomas5_SolveL(lhs(:, 1:nd/2), f)
        !     call Thomas5_SolveU(lhs(:, nd/2 + 1:nd), f)
        ! case (7)
        !     call Thomas7_SolveL(lhs(:, 1:nd/2), f)
        !     call Thomas7_SolveU(lhs(:, nd/2 + 1:nd), f)
        ! end select

        write (str, *) nd
        call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

        ! -------------------------------------------------------------------
        print *, new_line('a'), 'Solve circulant system, bands', nd

        lhs(:, 1:nd) = rhs(:, 1:nd)
        select case (nd)
        case (3)
            call ThomasCirculant_3_Initialize(lhs(:, 1:nd/2), &
                                              lhs(:, nd/2 + 1:nd), &
                                              lhs(1, nd + 1))
        case (5)
            call ThomasCirculant_5_Initialize(lhs(:, 1:nd/2), &
                                              lhs(:, nd/2 + 1:nd), &
                                              lhs(1, nd + 1))
        case (7)
            cycle
        end select

        ! compute forcing
        ! call MatMul_Halo_X(rhs(:, 1:nd), u, u(:, nsize - nd/2 + 1:nsize), u(:, 1:nd/2), f)
        call MatMul_Halo_X_ThomasL_Y(rhs=rhs(:, 1:nd), &
                                     u=u, &
                                     u_halo_m=u(:, nsize - nd/2 + 1:nsize), &
                                     u_halo_p=u(:, 1:nd/2), &
                                     f=f, &
                                     L=lhs(:, 1:nd/2))

        select case (nd)
        case (3)
            ! call Thomas3_SolveL(lhs(:, 1:nd/2), f)
            call Thomas3_SolveU(lhs(:, nd/2 + 1:nd), f)
            call ThomasCirculant_3_Reduce(lhs(:, 1:nd/2), &
                                          lhs(:, nd/2 + 1:nd), &
                                          lhs(:, nd + 1), &
                                          f, wrk)
        case (5)
            ! call Thomas5_SolveL(lhs(:, 1:nd/2), f)
            call Thomas5_SolveU(lhs(:, nd/2 + 1:nd), f)
            call ThomasCirculant_5_Reduce(lhs(:, 1:nd/2), &
                                          lhs(:, nd/2 + 1:nd), &
                                          lhs(:, nd + 1), &
                                          f)
        case (7)
        end select

        write (str, *) nd
        call check(f, u, 'linsys-'//trim(adjustl(str))//'.dat')

    end do

    stop

contains
    ! ###################################################################
    ! ###################################################################
    ! calculate f = A u, where A is narrow banded with diagonals given by rhs
    ! Allowing for different band size at the bottom and at the top
    subroutine matmul_biased(rhs, rhs_b, rhs_t, u, f)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)

        integer(wi) nx, ir
        integer(wi) ndr, idr, ic
        integer nx_b, ndr_b, idr_b, nx_t, ndr_t, idr_t

        ! ###################################################################
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

        ! -------------------------------------------------------------------
        ! interior points
        nx = size(rhs, 1)       ! size of the system
        ndr = size(rhs, 2)      ! # of diagonals
        idr = ndr/2 + 1         ! index of centerline diagonal

        ! upper boundary
        nx_t = size(rhs_t, 1)
        ndr_t = size(rhs_t, 2)
        idr_t = ndr_t/2 + 1

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ic)*u(:, ir - ic) + &
                           rhs(ir, idr + ic)*u(:, ir + ic)
            end do
        end do

        ! -------------------------------------------------------------------
        ! upper boundary
        do ir = nx_t - 1, ndr_t/2, -1
            f(:, nx - ir) = rhs_t(nx_t - ir, idr_t)*u(:, nx - ir)
            do ic = 1, ndr_t/2
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs_t(nx_t - ir, idr_t - ic)*u(:, nx - ir - ic) + &
                                rhs_t(nx_t - ir, idr_t + ic)*u(:, nx - ir + ic)
            end do
        end do

        do ir = ndr_t/2 - 1, 0, -1
            f(:, nx - ir) = rhs_t(nx_t - ir, idr_t + ir)*u(:, nx)
            do ic = 1, ndr_t/2 + ir
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs_t(nx_t - ir, idr_t + ir - ic)*u(:, nx - ic)
            end do
        end do

        return
    end subroutine matmul_biased

    ! ###################################################################
    ! ###################################################################
    ! Assumes that ndl is less or equal than nx_b, nx_t
    subroutine matmul_biased_solveL(rhs, rhs_b, rhs_t, u, f, L)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, ir
        integer(wi) ndr, idr, ic
        integer nx_b, ndr_b, idr_b, nx_t, ndr_t, idr_t
        integer(wi) ndl

        ! ###################################################################
        ! interior points
        nx = size(rhs, 1)       ! size of the system
        ndr = size(rhs, 2)      ! # of diagonals
        idr = ndr/2 + 1         ! index of centerline diagonal

        ! lower boundary
        nx_b = size(rhs_b, 1)
        ndr_b = size(rhs_b, 2)
        idr_b = ndr_b/2 + 1

        ! upper boundary
        nx_t = size(rhs_t, 1)
        ndr_t = size(rhs_t, 2)
        idr_t = ndr_t/2 + 1

        ! array L
        ndl = size(L, 2)
        if (any([nx_b, nx_t] < ndl)) then
            print *, __FILE__//'Error'
        end if

        ! -------------------------------------------------------------------
        ! lower boundary
        do ir = 1, ndr_b/2
            f(:, ir) = rhs_b(ir, idr_b - ir + 1)*u(:, 1)
            do ic = 2, ndr_b/2 + ir
                f(:, ir) = f(:, ir) + &
                           rhs_b(ir, idr_b - ir + ic)*u(:, ic)
            end do
            do ic = 1, min(ir - 1, ndl)  ! solve L
                f(:, ir) = f(:, ir) + f(:, ir - ic)*L(ir, ndl - ic + 1)
            end do
        end do

        do ir = idr_b, nx_b
            f(:, ir) = rhs_b(ir, idr_b)*u(:, ir)
            do ic = 1, ndr_b/2
                f(:, ir) = f(:, ir) + &
                           rhs_b(ir, idr_b + ic)*u(:, ir + ic) + &
                           rhs_b(ir, idr_b - ic)*u(:, ir - ic)
            end do
            do ic = 1, min(ir - 1, ndl)  ! solve L
                f(:, ir) = f(:, ir) + f(:, ir - ic)*L(ir, ndl - ic + 1)
            end do
        end do

        ! -------------------------------------------------------------------
        ! interior points
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ic)*u(:, ir - ic) + &
                           rhs(ir, idr + ic)*u(:, ir + ic)
            end do
            do ic = 1, ndl      ! solve L
                f(:, ir) = f(:, ir) + f(:, ir - ic)*L(ir, ndl - ic + 1)
            end do
        end do

        ! -------------------------------------------------------------------
        ! upper boundary
        do ir = nx_t - 1, ndr_t/2, -1
            f(:, nx - ir) = rhs_t(nx_t - ir, idr_t)*u(:, nx - ir)
            do ic = 1, ndr_t/2
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs_t(nx_t - ir, idr_t - ic)*u(:, nx - ir - ic) + &
                                rhs_t(nx_t - ir, idr_t + ic)*u(:, nx - ir + ic)
            end do
            do ic = 1, ndl      ! solve L
                f(:, nx - ir) = f(:, nx - ir) + f(:, nx - ir - ic)*L(nx - ir, ndl - ic + 1)
            end do
        end do

        do ir = ndr_t/2 - 1, 0, -1
            f(:, nx - ir) = rhs_t(nx_t - ir, idr_t + ir)*u(:, nx)
            do ic = 1, ndr_t/2 + ir
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs_t(nx_t - ir, idr_t + ir - ic)*u(:, nx - ic)
            end do
            do ic = 1, ndl      ! solve L
                f(:, nx - ir) = f(:, nx - ir) + f(:, nx - ir - ic)*L(nx - ir, ndl - ic + 1)
            end do
        end do

        return
    end subroutine matmul_biased_solveL

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

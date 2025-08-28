program VPARTIAL
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Memory, only: imax, jmax, kmax, isize_field, isize_wrk1d, isize_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
    use TLab_Arrays, only: wrk2d, txc
    use TLab_Grid, only: grid_dt
    use Thomas3
    use Thomas5
    use FDM, only: fdm_dt, FDM_CreatePlan
    use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
    use FDM_Derivative, only: FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA, FDM_COM6_JACOBIAN_HYPER, FDM_COM4_DIRECT, FDM_COM6_DIRECT
    use FDM_ComX_Direct
    use FDM_Base
    use FDM_MatMul
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian

    implicit none

    integer(wi) :: i, l, nlines

    real(wp), dimension(:, :), pointer :: u
    real(wp), dimension(:, :), pointer :: du1_a, du1_b, du1_c, du1_n
    real(wp), dimension(:, :), pointer :: du2_a, du2_n1, du2_n2, du2_n3
    real(wp) :: wk, x_0, coef(5)!, dummy
    integer(wi) :: test_type, ibc, ip, ic, ndr, idr, ndl, idl, im, ib, ip2
    integer(wi) :: nmin, nmax, nsize
    real(wp) rhsr_b(5, 0:7), rhsr_t(0:4, 8)

    real(wp), allocatable :: w(:, :), z(:, :)         ! for case 5

    integer :: bcs_cases(4)
    integer :: fdm_cases(5) = [FDM_COM4_JACOBIAN, &
                               FDM_COM6_JACOBIAN, &
                               FDM_COM6_JACOBIAN_PENTA, &
                               FDM_COM4_DIRECT, &
                               FDM_COM6_DIRECT]
    character(len=32) :: fdm_names(5) = [character(len=32) :: &
                                         'Jacobian 4', &
                                         'Jacobian 6', &
                                         'Jacobian 6 penta-diagonal', &
                                         'Direct 4', &
                                         'Direct 6']

    type(grid_dt) :: x
    type(fdm_dt) g

    ! ###################################################################
    ! Initialize
    imax = 2
    jmax = 3
    kmax = 256
    nlines = imax*jmax

    x%size = kmax
    x%scale = 1.0_wp
    ! x%periodic = .true.
    allocate (x%nodes(kmax))

    isize_field = imax*jmax*kmax
    isize_txc_field = isize_field
    isize_wrk3d = isize_txc_field
    isize_wrk1d = kmax
    isize_wrk2d = max(imax*jmax, max(imax*kmax, jmax*kmax))

    inb_txc = 9

    call TLab_Initialize_Memory(__FILE__)

    u(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 1)

    du1_a(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 2)
    du1_b(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    du1_c(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    du1_n(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    du2_a(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    du2_n1(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 7)
    du2_n2(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 8)
    du2_n3(1:nlines, 1:kmax) => txc(1:imax*jmax*kmax, 9)

    allocate (w(1, kmax), z(1, kmax))

    print *, '1. First-order derivative.'
    print *, '2. Second-order derivative.'
    print *, '3. Reduction routines.'
    print *, '4. Boundary conditions.'
    print *, '5. Boundary conditions decomposition.'
    read (*, *) test_type

    !  ###################################################################
    if (x%periodic) then
        do i = 1, kmax
            x%nodes(i) = real(i - 1, wp)/real(kmax, wp)*x%scale
        end do

    else
        do i = 1, kmax
            x%nodes(i) = real(i - 1, wp)/real(kmax - 1, wp)*x%scale
        end do
        ! open (21, file='y.dat')
        ! do i = 1, kmax
        !     read (21, *) x%nodes(i)
        ! end do
        ! close (21)
    end if

    g%periodic = x%periodic
    g%der1%mode_fdm = FDM_COM6_JACOBIAN     ! default
    g%der2%mode_fdm = g%der1%mode_fdm
    call FDM_CreatePlan(x, g)
    ndr = g%der1%nb_diag(2)
    ndl = g%der1%nb_diag(1)

    ! ###################################################################
    ! Define the function and analytic derivatives
    x_0 = 0.75_wp
    wk = 1.0_wp

    do i = 1, kmax
        ! single-mode
        u(:, i) = 1.0_wp + sin(2.0_wp*pi_wp/g%scale*wk*(g%nodes(i) - x_0*x%scale)) ! + pi_wp/4.0_wp)
        du1_a(:, i) = (2.0_wp*pi_wp/g%scale*wk) &
                      *cos(2.0_wp*pi_wp/g%scale*wk*(g%nodes(i) - x_0*x%scale))! + pi_wp/4.0_wp)
        du2_a(:, i) = -(2.0_wp*pi_wp/g%scale*wk)**2.0_wp &
                      *sin(2.0_wp*pi_wp/g%scale*wk*(g%nodes(i) - x_0*x%scale))! + pi_wp/4.0_wp)
        ! ! Gaussian
        ! u(:, i) = exp(-(g%nodes(i) - x_0*g%scale)**2/(2.0_wp*(g%scale/wk)**2))
        ! du1_a(:, i) = -(g%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*u(:, i)
        ! du2_a(:, i) = -(g%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*du1_a(:, i) &
        !               - 1.0_wp/(g%scale/wk)**2*u(:, i)
        ! ! exponential
        ! u(:, i) = exp(-g%nodes(i)*wk)
        ! du1_a(:, i) = -wk*u(:, i)
        ! step
        ! u(:, i) = max(0.0_wp, (g%nodes(i) - g%nodes(kmax/2))*x_0)
        ! du1_a(:, i) = (1.0_wp + sign(1.0_wp, g%nodes(i) - g%nodes(kmax/2)))*0.5_wp*x_0
        ! ! tanh
        ! u(:, i) = x_0*log(1.0_wp + exp((g%nodes(i) - g%nodes(kmax/2))/x_0))
        ! du1_a(:, i) = 0.5_wp*(1.0_wp + tanh(0.5_wp*(g%nodes(i) - g%nodes(kmax/2))/x_0))
        ! ! Polynomial
        ! dummy = 4.0_wp
        ! u(:, i) = ((g%scale - g%nodes(i))/wk)**dummy
        ! du1_a(:, i) = -dummy*((g%scale - g%nodes(i))/wk)**(dummy - 1.0_wp)
        ! ! zero
        ! u(:, i) = 0.0_wp
        ! du1_a(:, i) = 0.0_wp
        ! ! delta-function
        ! u(:, i) = max(0.0_wp, 2.0_wp - real(i, wp))
        ! du1_a(:, i) = 0.0_wp
        ! du2_a(:, i) = 0.0_wp
    end do

    ! ###################################################################
    ! First-order derivative
    ! ###################################################################
    select case (test_type)
    case (1)
        do im = 1, size(fdm_cases)
            print *, new_line('a'), fdm_names(im)

            g%der1%mode_fdm = fdm_cases(im)
            call FDM_CreatePlan(x, g)

            call FDM_Der1_Solve(nlines, BCS_NONE, g%der1, g%der1%lu, u, du1_n, wrk2d)

            call check(u, du1_a, du1_n, 'partial.dat')

        end do

        ! ###################################################################
        ! Second-order derivative
        ! ###################################################################
    case (2)
        fdm_cases(3) = FDM_COM6_JACOBIAN_HYPER
        fdm_names(3) = 'Jacobian 6, hyper-diffusive'

        do im = 1, size(fdm_cases)
            print *, new_line('a'), fdm_names(im)

            g%der2%mode_fdm = fdm_cases(im)
            call FDM_CreatePlan(x, g)

            call FDM_Der1_Solve(nlines, BCS_NONE, g%der1, g%der1%lu, u, du1_n, wrk2d)  ! I need du1_n in Jacobian formulation
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, u, du2_n1, du1_n, wrk2d)

            call check(u, du2_a, du2_n1, 'partial.dat')

        end do

        ! ###################################################################
        !   Testing the reduction routines
        ! ###################################################################
    case (3)
        bcs_cases(1:3) = [BCS_MIN, BCS_MAX, BCS_BOTH]

        do ip = 1, 3
            ibc = bcs_cases(ip)
            print *, new_line('a'), 'Bcs case ', ibc

            nmin = 1
            nmax = g%size
            if (any([BCS_MIN, BCS_BOTH] == ibc)) then
                nmin = nmin + 1
            end if
            if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                nmax = nmax - 1
            end if
            nsize = nmax - nmin + 1

            do im = 1, 3
                g%der1%mode_fdm = fdm_cases(im)
                print *, fdm_names(im)

                select case (g%der1%mode_fdm)
                case (FDM_COM4_JACOBIAN)
                    call FDM_C1N4_Jacobian(kmax, g%jac, g%der1%lu, g%der1%rhs, g%der1%nb_diag, coef, g%periodic)

                case (FDM_COM6_JACOBIAN)
                    call FDM_C1N6_Jacobian(kmax, g%jac, g%der1%lu, g%der1%rhs, g%der1%nb_diag, coef, g%periodic)

                case (FDM_COM6_JACOBIAN_PENTA)
                    call FDM_C1N6_Jacobian_Penta(kmax, g%jac, g%der1%lu, g%der1%rhs, g%der1%nb_diag, coef, g%periodic)

                end select
                ndl = g%der1%nb_diag(1)
                idl = g%der1%nb_diag(1)/2 + 1
                ndr = g%der1%nb_diag(2)
                idr = g%der1%nb_diag(2)/2 + 1

                ! g%rhsr_b = 0.0_wp
                ! g%rhsr_t = 0.0_wp
                call FDM_Bcs_Reduce(ibc, g%der1%lu(:, 1:ndl), g%der1%rhs(:, 1:ndr), rhsr_b, rhsr_t)

                select case (g%der1%nb_diag(1))
                case (3)
                    call Thomas3_LU(nsize, g%der1%lu(nmin:nmax, 1), g%der1%lu(nmin:nmax, 2), g%der1%lu(nmin:nmax, 3))
                case (5)
   call Thomas5_LU(nsize, g%der1%lu(nmin:nmax, 1), g%der1%lu(nmin:nmax, 2), g%der1%lu(nmin:nmax, 3), g%der1%lu(nmin:nmax, 4), g%der1%lu(nmin:nmax, 5))
                end select

                du1_n(:, 1) = u(:, 1)           ! boundary condition
                du1_n(:, kmax) = u(:, kmax)
                select case (g%der1%nb_diag(2))
                case (3)
                    call MatMul_3d_antisym(g%der1%rhs(:, 1:3), u, du1_n, &
                                           ibc, rhs_b=rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                case (5)
                    call MatMul_5d_antisym(g%der1%rhs(:, 1:5), u, du1_n, &
                                           ibc, rhs_b=rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                case (7)
                    call MatMul_7d_antisym(g%der1%rhs(:, 1:7), u, du1_n, &
                                           ibc, rhs_b=rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                end select

                select case (g%der1%nb_diag(1))
                case (3)
                    call Thomas3_Solve(nsize, nlines, g%der1%lu(nmin:nmax, 1), g%der1%lu(nmin:nmax, 2), g%der1%lu(nmin:nmax, 3), du1_n(:, nmin:nmax))
                case (5)
                    call Thomas5_Solve(nsize, nlines, g%der1%lu(nmin:nmax, 1), g%der1%lu(nmin:nmax, 2), g%der1%lu(nmin:nmax, 3), g%der1%lu(nmin:nmax, 4), g%der1%lu(nmin:nmax, 5), du1_n(:, nmin:nmax))
                end select

                if (any([BCS_MIN, BCS_BOTH] == ibc)) then
                    du1_n(:, 1) = wrk2d(:, 1)
                    do ic = 1, idl - 1
                        du1_n(:, 1) = du1_n(:, 1) + g%der1%lu(1, idl + ic)*du1_n(:, 1 + ic)
                    end do
                    du1_n(:, 1) = du1_n(:, 1) + g%der1%lu(1, 1)*du1_n(:, 1 + ic)
                end if
                if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                    du1_n(:, kmax) = wrk2d(:, 2)
                    do ic = 1, idl - 1
                        du1_n(:, kmax) = du1_n(:, kmax) + g%der1%lu(kmax, idl - ic)*du1_n(:, kmax - ic)
                    end do
                    du1_n(:, kmax) = du1_n(:, kmax) + g%der1%lu(kmax, ndl)*du1_n(:, kmax - ic)
                end if
                call check(u, du1_a, du1_n, 'partial.dat')

            end do

        end do

        ! ###################################################################
        ! Boundary conditions
        ! ###################################################################
    case (4)
        bcs_cases(1:4) = [BCS_DD, BCS_ND, BCS_DN, BCS_NN]

#define bcs_hb(i) wrk2d(i,1)
#define bcs_ht(i) wrk2d(i,2)

        do im = 1, size(fdm_cases)
            g%der1%mode_fdm = fdm_cases(im)
            print *, new_line('a'), fdm_names(im)

            g%der1%mode_fdm = fdm_cases(im)
            call FDM_CreatePlan(x, g)

            do ib = 1, 4
                ibc = bcs_cases(ib)
                print *, new_line('a'), 'Bcs case ', ibc

                nmin = 1; nmax = g%size
                if (any([BCS_ND, BCS_NN] == ibc)) then
                    du1_n(:, 1) = du1_a(:, 1)
                    nmin = nmin + 1
                end if
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    du1_n(:, kmax) = du1_a(:, kmax)
                    nmax = nmax - 1
                end if
                nsize = nmax - nmin + 1

                ! -------------------------------------------------------------------
                ! Calculate RHS in system of equations A u' = B u
                call g%der1%matmul(g%der1%rhs, u, du1_n, ibc, g%der1%rhs_b, g%der1%rhs_t, bcs_hb(:), bcs_ht(:))

                ! -------------------------------------------------------------------
                ! Solve for u' in system of equations A u' = B u
                ip = ibc*5

                select case (g%der1%nb_diag(1))
                case (3)
                    call Thomas3_Solve(nsize, nlines, &
                                       g%der1%lu(nmin:nmax, ip + 1), &
                                       g%der1%lu(nmin:nmax, ip + 2), &
                                       g%der1%lu(nmin:nmax, ip + 3), &
                                       du1_n(:, nmin:nmax))
                case (5)
                    call Thomas5_Solve(nsize, nlines, &
                                       g%der1%lu(nmin:nmax, ip + 1), &
                                       g%der1%lu(nmin:nmax, ip + 2), &
                                       g%der1%lu(nmin:nmax, ip + 3), &
                                       g%der1%lu(nmin:nmax, ip + 4), &
                                       g%der1%lu(nmin:nmax, ip + 5), &
                                       du1_n(:, nmin:nmax))
                end select

                call check(u(:, nmin:nmax), du1_a(:, nmin:nmax), du1_n(:, nmin:nmax), 'partial.dat')

                idl = g%der1%nb_diag(1)/2 + 1
                if (any([BCS_ND, BCS_NN] == ibc)) then
                    do ic = 1, idl - 1
                        bcs_hb(1:nlines) = bcs_hb(1:nlines) + g%der1%lu(1, ip + idl + ic)*du1_n(:, 1 + ic)
                    end do
                    print *, u(:, 1)
                    print *, bcs_hb(1:nlines)
                end if
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    do ic = 1, idl - 1
                        bcs_ht(1:nlines) = bcs_ht(1:nlines) + g%der1%lu(kmax, ip + idl - ic)*du1_n(:, kmax - ic)
                    end do
                    print *, u(:, kmax)
                    print *, bcs_ht(1:nlines)
                end if

            end do

        end do

#undef bcs_hb
#undef bcs_ht

        ! ###################################################################
        ! Decomposition of Neumann boundary condition
        ! ###################################################################
    case (5)
        bcs_cases(1:3) = [BCS_MIN, BCS_MAX, BCS_BOTH]

        do im = 2, 2!size(fdm_cases)
            g%der1%mode_fdm = fdm_cases(im)
            print *, fdm_names(im)

            g%der1%mode_fdm = fdm_cases(im)
            call FDM_CreatePlan(x, g)

            do ip = 1, 3
                ibc = bcs_cases(ip)
                print *, new_line('a'), 'Bcs case ', ibc

                i = 1
                w = 0.0_wp; w(1, i) = 1.0_wp
                z = 0.0_wp; z(1, i) = 1.0_wp

                nmin = 1
                nmax = g%size
                if (any([BCS_ND, BCS_NN] == ibc)) then
                    z(1, 1) = w(1, 1)
                    nmin = nmin + 1
                end if
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    z(1, kmax) = w(1, kmax)
                    nmax = nmax - 1
                end if
                nsize = nmax - nmin + 1

                ! call FDM_Der1_Solve(1, ibc, g%der1, g%der1%lu, w, z, wrk2d)
                call g%der1%matmul(g%der1%rhs, w, z, ibc, g%der1%rhs_b, g%der1%rhs_t)

                ip2 = ibc*5
                select case (g%der1%nb_diag(1))
                case (3)
                    call Thomas3_Solve(nsize, 1, g%der1%lu(nmin:, ip2 + 1), g%der1%lu(nmin:, ip2 + 2), g%der1%lu(nmin:, ip2 + 3), &
                                       z(:, nmin:))
                case (5)
                    print *, 'undeveloped'
                    ! call Thomas5_Solve(nsize, nlines, lu1(nmin:, ip + 1), lu1(nmin:, ip + 2), lu1(nmin:, ip + 3), lu1(nmin:, ip + 4), lu1(nmin:, ip + 5), &
                    !                    result(:, nmin:))
                end select

                do i = 1, kmax !nmin, nmax
                    print *, i, z(1, i)
                end do

            end do
        end do

    end select

    stop

    ! ###################################################################
contains
    subroutine check(u, du_a, du_n, name)
        real(wp), intent(in) :: u(:, :), du_a(:, :), du_n(:, :)
        character(len=*), optional :: name

        real(wp) dummy, error_l2, error_max

        if (present(name)) then
            open (20, file=name)
        end if
        error_l2 = 0.0_wp
        error_max = 0.0_wp
        dummy = 0.0_wp
        do i = 1, size(u, 2)
            do l = 1, size(u, 1)
                if (present(name)) then
                    write (20, 1000) g%nodes(i), u(l, i), du_a(l, i), du_n(l, i), du_a(l, i) - du_n(l, i)
                end if
                dummy = dummy + du_a(l, i)*du_a(l, i)
                error_l2 = error_l2 + (du_a(l, i) - du_n(l, i))**2.0_wp
                error_max = max(error_max, abs(du_a(l, i) - du_n(l, i)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Solution L2-norm ...........:', sqrt(g%jac(1, 1)*dummy)/real(nlines, wp)
        if (dummy == 0.0_wp) return
        write (*, *) 'Relative Error L2-norm .....:', sqrt(g%jac(1, 1)*error_l2)/maxval(abs(du1_a))
        write (*, *) 'Relative Error Linf-norm ...:', error_max/maxval(abs(du1_a))

        return
1000    format(5(1x, e12.5))
    end subroutine check

end program VPARTIAL

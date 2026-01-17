program VPARTIAL
    use TLab_Constants, only: wp, wi, pi_wp, roundoff_wp
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Memory, only: imax, jmax, kmax, isize_field, isize_wrk1d, isize_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
    use TLab_Arrays, only: wrk1d, wrk2d, txc
    use TLab_Grid, only: grid_dt
    use Thomas
    use FDM, only: fdm_dt, FDM_CreatePlan, FDM_CreatePlan_Der1
    use FDM_Derivative
    use FDM_derivative_Neumann
    use FDM_ComX_Direct
    use FDM_Base
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use FDM_Derivative_1order_X

    implicit none

    integer(wi) :: i, l, nlines

    real(wp), dimension(:, :), pointer :: u
    real(wp), dimension(:, :), pointer :: du1_a, du1_b, du1_c, du1_n
    real(wp), dimension(:, :), pointer :: du2_a, du2_n1, du2_n2, du2_n3
    real(wp) :: wk, x_0, dummy
    integer(wi) :: test_type, ibc, ip, ic, ndr, idr, ndl, idl, im, ib
    integer(wi) :: nmin, nmax, nsize

    real(wp), allocatable :: c(:)         ! for case 5

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
    character(len=32) str
    type(grid_dt) :: x
    type(fdm_dt) g
    class(der_dt), allocatable :: fdm_der1

    ! ###################################################################
    ! Initialize
    imax = 2
    jmax = 3
    kmax = 768
    nlines = imax*jmax

    x%size = kmax
    x%scale = 1.0_wp
    x%periodic = .false.
    x%periodic = .true.
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

    allocate (c(kmax))

    print *, '1. First-order derivative.'
    print *, '2. Second-order derivative.'
    print *, '3. Boundary conditions.'
    print *, '4. Neumann decomposition.'
    print *, '5. Reduction routines.'
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
        open (21, file='z.dat')
        do i = 1, kmax
            read (21, *) x%nodes(i)
        end do
        ! wrk1d(1:kmax, 1) = x%nodes(1:kmax)  ! reverse
        ! do i = 1, kmax
        !     x%nodes(i) = x%nodes(kmax) - wrk1d(kmax - i + 1, 1)
        ! end do
        close (21)
    end if

    g%periodic = x%periodic
    g%der1%mode_fdm = FDM_COM6_JACOBIAN     ! default
    g%der2%mode_fdm = g%der1%mode_fdm
    call FDM_CreatePlan(x, g)
    ndr = g%der1%nb_diag(2)
    ndl = g%der1%nb_diag(1)

    call FDM_CreatePlan_Der1(x, fdm_der1, FDM_COM6_JACOBIAN)

    ! ###################################################################
    ! Define the function and analytic derivatives
    x_0 = 0.1_wp
    wk = 6.0_wp

    do i = 1, kmax
        ! single-mode
        u(:, i) = 1.0_wp + sin(2.0_wp*pi_wp/g%scale*wk*(x%nodes(i) - x_0*x%scale)) ! + pi_wp/4.0_wp)
        du1_a(:, i) = (2.0_wp*pi_wp/g%scale*wk) &
                      *cos(2.0_wp*pi_wp/g%scale*wk*(x%nodes(i) - x_0*x%scale))! + pi_wp/4.0_wp)
        du2_a(:, i) = -(2.0_wp*pi_wp/g%scale*wk)**2 &
                      *sin(2.0_wp*pi_wp/g%scale*wk*(x%nodes(i) - x_0*x%scale))! + pi_wp/4.0_wp)
        ! ! Gaussian
        ! u(:, i) = exp(-(x%nodes(i) - x_0*g%scale)**2/(2.0_wp*(g%scale/wk)**2))
        ! du1_a(:, i) = -(x%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*u(:, i)
        ! du2_a(:, i) = -(x%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*du1_a(:, i) &
        !               - 1.0_wp/(g%scale/wk)**2*u(:, i)
        ! ! exponential
        ! ! u(:, i) = exp(-x%nodes(i)*wk)
        ! ! du1_a(:, i) = -wk*u(:, i)
        ! ! du2_a(:, i) = wk**2*u(:, i)
        ! u(:, i) = exp(x%nodes(i)*wk)
        ! du1_a(:, i) = wk*u(:, i)
        ! du2_a(:, i) = wk**2*u(:, i)
        ! step
        ! u(:, i) = max(0.0_wp, (x%nodes(i) - x%nodes(kmax/2))*x_0)
        ! du1_a(:, i) = (1.0_wp + sign(1.0_wp, x%nodes(i) - x%nodes(kmax/2)))*0.5_wp*x_0
        ! ! tanh
        ! u(:, i) = x_0*log(1.0_wp + exp((x%nodes(i) - x%nodes(kmax/2))/x_0))
        ! du1_a(:, i) = 0.5_wp*(1.0_wp + tanh(0.5_wp*(x%nodes(i) - x%nodes(kmax/2))/x_0))
        ! ! Polynomial
        ! ! dummy = 5.0_wp
        ! ! u(:, i) = ((g%scale - x%nodes(i))*wk)**dummy
        ! ! du1_a(:, i) = -dummy*wk*((g%scale - x%nodes(i))*wk)**(dummy - 1.0_wp)
        ! dummy = 4.0_wp
        ! u(:, i) = 1.0_wp + (x%nodes(i)*wk)**dummy
        ! du1_a(:, i) = dummy*wk*(x%nodes(i)*wk)**(dummy - 1.0_wp)
        ! du2_a(:, i) = dummy*(dummy - 1.0_wp)*wk**2*(x%nodes(i)*wk)**(dummy - 2.0_wp)
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

            ! old formulation
            g%der1%mode_fdm = fdm_cases(im)
            call FDM_CreatePlan(x, g)

            call FDM_Der1_Solve(nlines, g%der1, g%der1%lu, u, du1_n, wrk2d)

            write (str, *) im
            call check(u, du1_a, du1_n, 'partial-old-'//trim(adjustl(str))//'.dat')
            call write_scheme(g%der1%lhs(:, 1:g%der1%nb_diag(1)), &
                              g%der1%rhs(:, 1:g%der1%nb_diag(2)), 'fdm1-old-'//trim(adjustl(str)))

            ! new formulation
            call FDM_CreatePlan_Der1(x, fdm_der1, fdm_cases(im))

            call fdm_der1%compute(nlines, u, du1_n)

            call check(u, du1_a, du1_n, 'partial-'//trim(adjustl(str))//'.dat')
            call write_scheme(fdm_der1%lhs, &
                              fdm_der1%rhs, 'fdm1-'//trim(adjustl(str)))

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

            call FDM_Der1_Solve(nlines, g%der1, g%der1%lu, u, du1_n, wrk2d)  ! I need du1_n in Jacobian formulation
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, u, du2_n1, du1_n, wrk2d)

            write (str, *) im
            call check(u, du2_a, du2_n1, 'partial-'//trim(adjustl(str))//'.dat')
            call write_scheme(g%der2%lhs(:, 1:g%der2%nb_diag(1)), &
                              g%der2%rhs(:, 1:g%der2%nb_diag(2)), 'fdm2-'//trim(adjustl(str)))

        end do

        ! ###################################################################
        ! Boundary conditions
        ! ###################################################################
    case (3)
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
                    bcs_hb(1:nlines) = du1_a(1:nlines, 1)
                    nmin = nmin + 1
                end if
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    du1_n(:, kmax) = du1_a(:, kmax)
                    bcs_ht(1:nlines) = du1_a(1:nlines, kmax)
                    nmax = nmax - 1
                end if
                nsize = nmax - nmin + 1

                ndl = g%der1%nb_diag(1)
                idl = ndl/2 + 1
                ndr = g%der1%nb_diag(2)
                idr = ndr/2 + 1

                ip = ibc*5

                ! -------------------------------------------------------------------
                ! Calculate RHS in system of equations A u' = B u
                select case (ibc)
                case (BCS_DD)
                    ! call g%der1%matmul(rhs=g%der1%rhs(:, 1:ndr), &
                    !                         rhs_b=g%der1%rhs(1:ndr/2, 1:ndr), &
                    !                         rhs_t=g%der1%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                    !                         u=u, &
                    !                         f=du1_n)
                    call g%der1%matmul_thomas(rhs=g%der1%rhs(:, 1:ndr), &
                                              rhs_b=g%der1%rhs(1:ndr/2, 1:ndr), &
                                              rhs_t=g%der1%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                              u=u, &
                                              L=g%der1%lu(:, ip + 1:ip + ndl/2), &
                                              f=du1_n)
                case (BCS_ND)
                    ! call g%der1%matmul(rhs=g%der1%rhs(:, 1:ndr), &
                    !                         rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         rhs_t=g%der1%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                    !                         u=u, &
                    !                         f=du1_n, bcs_b=bcs_hb(1:nlines))
                    call g%der1%matmul_thomas(rhs=g%der1%rhs(:, 1:ndr), &
                                              rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              rhs_t=g%der1%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                              u=u, &
                                              L=g%der1%lu(:, ip + 1:ip + ndl/2), &
                                              f=du1_n, bcs_b=bcs_hb(1:nlines))
                case (BCS_DN)
                    ! call g%der1%matmul(rhs=g%der1%rhs(:, 1:ndr), &
                    !                         rhs_b=g%der1%rhs(1:ndr/2, 1:ndr), &
                    !                         rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         u=u, &
                    !                         f=du1_n, bcs_t=bcs_ht(1:nlines))
                    call g%der1%matmul_thomas(rhs=g%der1%rhs(:, 1:ndr), &
                                              rhs_b=g%der1%rhs(1:ndr/2, 1:ndr), &
                                              rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              u=u, &
                                              L=g%der1%lu(:, ip + 1:ip + ndl/2), &
                                              f=du1_n, bcs_t=bcs_ht(1:nlines))
                case (BCS_NN)
                    ! call g%der1%matmul(rhs=g%der1%rhs(:, 1:ndr), &
                    !                         rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         u=u, &
                    !                         f=du1_n, bcs_b=bcs_hb(1:nlines), bcs_t=bcs_ht(1:nlines))
                    call g%der1%matmul_thomas(rhs=g%der1%rhs(:, 1:ndr), &
                                              rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              u=u, &
                                              L=g%der1%lu(:, ip + 1:ip + ndl/2), &
                                              f=du1_n, bcs_b=bcs_hb(1:nlines), bcs_t=bcs_ht(1:nlines))

                end select

                ! -------------------------------------------------------------------
                ! Solve for u' in system of equations A u' = B u
                ! call g%der1%thomasL(g%der1%lu(nmin:nmax, ip + 1:ip + ndl/2), du1_n(:, nmin:nmax))
                call g%der1%thomasU(g%der1%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), du1_n(:, nmin:nmax))

                write (str, *) im
                call check(u, du1_a, du1_n, 'partial-'//trim(adjustl(str))//'.dat')

                if (any([BCS_ND, BCS_NN] == ibc)) then
                    do ic = 1, idl - 1
                        bcs_hb(1:nlines) = bcs_hb(1:nlines) + g%der1%lu(1, ip + idl + ic)*du1_n(:, 1 + ic)
                    end do
                    bcs_hb(1:nlines) = bcs_hb(1:nlines)/g%der1%rhs(1, idr)
                    print *, u(:, 1)
                    print *, bcs_hb(1:nlines)
                    write (*, *) 'Boundary Relative Error Linf-norm ...:', &
                        maxval(abs(bcs_hb(1:nlines) - u(1:nlines, 1)))/maxval(abs(u(1:nlines, 1)))
                end if
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    do ic = 1, idl - 1
                        bcs_ht(1:nlines) = bcs_ht(1:nlines) + g%der1%lu(kmax, ip + idl - ic)*du1_n(:, kmax - ic)
                    end do
                    bcs_ht(1:nlines) = bcs_ht(1:nlines)/g%der1%rhs(kmax, idr)
                    print *, u(:, kmax)
                    print *, bcs_ht(1:nlines)
                    write (*, *) 'Boundary Relative Error Linf-norm ...:', &
                        maxval(abs(bcs_ht(1:nlines) - u(1:nlines, kmax)))/maxval(abs(u(1:nlines, kmax)))
                end if

            end do

        end do

#undef bcs_hb
#undef bcs_ht

        ! ###################################################################
        ! Decomposition of Neumann boundary condition
        ! ###################################################################
    case (4)
        bcs_cases(1:3) = [BCS_MIN, BCS_MAX, BCS_BOTH]

#define bcs_hb(i) wrk2d(i,1)
#define bcs_ht(i) wrk2d(i,2)

#define c_b(i) wrk1d(i,1)
#define c_t(i) wrk1d(i,2)

        do im = 1, size(fdm_cases)
            g%der1%mode_fdm = fdm_cases(im)
            print *, new_line('a'), fdm_names(im)

            g%der1%mode_fdm = fdm_cases(im)
            call FDM_CreatePlan(x, g)

            do ib = 1, 3
                ibc = bcs_cases(ib)
                print *, new_line('a'), 'Bcs case ', ibc

                ! truncated version
                call FDM_Der1_NeumannMin_Initialize(g%der1, c_b(:), wrk1d(1, 3), wrk1d(1, 4), nmax)
                ! print *, nmax
                call FDM_Der1_NeumannMax_Initialize(g%der1, c_t(:), wrk1d(1, 3), wrk1d(1, 4), nmax)
                ! print *, nmax

                if (any([BCS_ND, BCS_NN] == ibc)) then
                    bcs_hb(1:nlines) = du1_a(:, 1)*c_b(1)
                    do i = 2, nmax
                        bcs_hb(1:nlines) = bcs_hb(1:nlines) + c_b(i)*u(:, i)
                    end do
                    write (*, *) 'Boundary Relative Error Linf-norm ...:', &
                        maxval(abs(bcs_hb(1:nlines) - u(1:nlines, 1)))/maxval(abs(u(1:nlines, 1)))
                end if
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    bcs_ht(1:nlines) = du1_a(:, kmax)*c_t(kmax)
                    do i = g%size - 1, g%size - nmax + 1, -1
                        bcs_ht(1:nlines) = bcs_ht(1:nlines) + c_t(i)*u(:, i)
                    end do
                    write (*, *) 'Boundary Relative Error Linf-norm ...:', &
                        maxval(abs(bcs_ht(1:nlines) - u(1:nlines, kmax)))/maxval(abs(u(1:nlines, kmax)))
                end if

                ! ! full version
                ! call FDM_Der1_Neumann_Initialize(ibc, g%der1, c_b(:), c_t(:), wrk1d(1, 3), wrk1d(1, 4))

                ! nmin = 1; nmax = g%size
                ! if (any([BCS_ND, BCS_NN] == ibc)) then
                !     bcs_hb(1:nlines) = du1_a(:, 1)*c_b(1)
                !     bcs_ht(1:nlines) = du1_a(:, 1)*c_t(1)           ! only used in BCS_NN
                !     nmin = nmin + 1
                ! end if
                ! if (any([BCS_DN, BCS_NN] == ibc)) then
                !     bcs_hb(1:nlines) = du1_a(:, kmax)*c_b(kmax)     ! only used in BCS_NN
                !     bcs_ht(1:nlines) = du1_a(:, kmax)*c_t(kmax)
                !     nmax = nmax - 1
                ! end if

                ! if (any([BCS_ND, BCS_NN] == ibc)) then
                !     do i = nmin, nmax
                !         bcs_hb(1:nlines) = bcs_hb(1:nlines) + c_b(i)*u(:, i)
                !     end do
                !     write (*, *) 'Relative Error Linf-norm ...:', &
                !         maxval(abs(bcs_hb(1:nlines) - u(1:nlines, 1)))/maxval(abs(u(1:nlines, 1)))
                !     ! print *, u(:, 1)
                !     ! print *, bcs_hb(1:nlines)
                ! end if
                ! if (any([BCS_DN, BCS_NN] == ibc)) then
                !     do i = nmax, nmin, -1
                !         bcs_ht(1:nlines) = bcs_ht(1:nlines) + c_t(i)*u(:, i)
                !     end do
                !     write (*, *) 'Relative Error Linf-norm ...:', &
                !         maxval(abs(bcs_ht(1:nlines) - u(1:nlines, kmax)))/maxval(abs(u(1:nlines, kmax)))
                !     ! print *, u(:, kmax)
                !     ! print *, bcs_ht(1:nlines)
                ! end if

            end do
        end do

#undef c_b
#undef c_t

#undef bcs_hb
#undef bcs_ht

        ! ###################################################################
        !   Testing the reduction routines
        ! ###################################################################
    case (5)
        bcs_cases(1:3) = [BCS_MIN, BCS_MAX, BCS_BOTH]

#define bcs_hb(i) wrk2d(i,1)
#define bcs_ht(i) wrk2d(i,2)

        do im = 1, size(fdm_cases)
            g%der1%mode_fdm = fdm_cases(im)
            print *, new_line('a'), fdm_names(im)

            g%der1%mode_fdm = fdm_cases(im)
            call FDM_CreatePlan(x, g)

            do ib = 1, 3
                ibc = bcs_cases(ib)
                print *, new_line('a'), 'Bcs case ', ibc

                nmin = 1; nmax = g%size
                if (any([BCS_MIN, BCS_BOTH] == ibc)) then
                    ! du1_n(:, 1) = u(:, 1)
                    bcs_hb(1:nlines) = u(1:nlines, 1)           ! boundary condition
                    nmin = nmin + 1
                end if
                if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                    ! du1_n(:, kmax) = u(:, kmax)
                    bcs_ht(1:nlines) = u(1:nlines, kmax)
                    nmax = nmax - 1
                end if
                nsize = nmax - nmin + 1

                ndl = g%der1%nb_diag(1)
                idl = g%der1%nb_diag(1)/2 + 1
                ndr = g%der1%nb_diag(2)
                idr = g%der1%nb_diag(2)/2 + 1

                g%der1%rhs_b1 = 0.0_wp
                g%der1%rhs_t1 = 0.0_wp
                g%der1%lu(:, 1:ndl) = g%der1%lhs(:, 1:ndl)
                ! call FDM_Bcs_Reduce(ibc, g%der1%lu(:, 1:ndl), g%der1%rhs(:, 1:ndr), g%der1%rhs_b1, g%der1%rhs_t1(:,2:))
                call FDM_Bcs_Reduce(ibc, g%der1%lu(:, 1:ndl), g%der1%rhs(:, 1:ndr), g%der1%rhs_b1, g%der1%rhs_t1)

                ! new format; extending to ndr+2 diagonals
                ! longer stencil
                i = 1
                g%der1%rhs_b1(i, ndr + 2) = g%der1%rhs_b1(i, 2); g%der1%rhs_b1(i, 2) = 0.0_wp
                i = max(idl, idr + 1)
                g%der1%rhs_t1(i, 1) = g%der1%rhs_t1(i, ndr + 1); g%der1%rhs_t1(i, ndr + 1) = 0.0_wp

                call Thomas_FactorLU_InPlace(g%der1%lu(nmin:nmax, 1:ndl/2), &
                                             g%der1%lu(nmin:nmax, ndl/2 + 1:ndl))

                ! -------------------------------------------------------------------
                ! Calculate RHS in system of equations A u' = B u
                select case (ibc)
                case (BCS_MIN)
                    ! call g%der1%matmul(rhs=g%der1%rhs(:, 1:ndr), &
                    !                         rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         rhs_t=g%der1%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                    !                         u=u, &
                    !                         f=du1_n)!, bcs_b=bcs_hb(1:nlines))
                    call g%der1%matmul_thomas(rhs=g%der1%rhs(:, 1:ndr), &
                                              rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              rhs_t=g%der1%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                              u=u, &
                                              L=g%der1%lu(:, 1:ndl/2), &
                                              f=du1_n, bcs_b=bcs_hb(1:nlines))
                case (BCS_MAX)
                    ! call g%der1%matmul(rhs=g%der1%rhs(:, 1:ndr), &
                    !                         rhs_b=g%der1%rhs(1:ndr/2, 1:ndr), &
                    !                         rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         u=u, &
                    !                         f=du1_n)!, bcs_t=bcs_ht(:))
                    call g%der1%matmul_thomas(rhs=g%der1%rhs(:, 1:ndr), &
                                              rhs_b=g%der1%rhs(1:ndr/2, 1:ndr), &
                                              rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              u=u, &
                                              L=g%der1%lu(:, 1:ndl/2), &
                                              f=du1_n, bcs_t=bcs_ht(1:nlines))
                case (BCS_BOTH)
                    ! call g%der1%matmul(rhs=g%der1%rhs(:, 1:ndr), &
                    !                         rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                    !                         u=u, &
                    !                         f=du1_n)!, bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))
                    call g%der1%matmul_thomas(rhs=g%der1%rhs(:, 1:ndr), &
                                              rhs_b=g%der1%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              rhs_t=g%der1%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                              u=u, &
                                              L=g%der1%lu(:, 1:ndl/2), &
                                              f=du1_n, bcs_b=bcs_hb(1:nlines), bcs_t=bcs_ht(1:nlines))

                end select

                ! -------------------------------------------------------------------
                ! Solve for u' in system of equations A u' = B u
                ! call g%der1%thomasL(g%der1%lu(nmin:nmax, 1:ndl/2), du1_n(:, nmin:nmax))
                call g%der1%thomasU(g%der1%lu(nmin:nmax, ndl/2 + 1:ndl), du1_n(:, nmin:nmax))

                if (any([BCS_MIN, BCS_BOTH] == ibc)) then
                    du1_n(:, 1) = bcs_hb(1:nlines)
                    do ic = 1, idl - 1
                        du1_n(:, 1) = du1_n(:, 1) + g%der1%lu(1, idl + ic)*du1_n(:, 1 + ic)
                    end do
                    du1_n(:, 1) = du1_n(:, 1) + g%der1%lu(1, 1)*du1_n(:, 1 + ic)
                    du1_n(:, 1) = du1_n(:, 1)/g%der1%lu(1, idl)
                end if
                if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                    du1_n(:, kmax) = bcs_ht(1:nlines)
                    do ic = 1, idl - 1
                        du1_n(:, kmax) = du1_n(:, kmax) + g%der1%lu(kmax, idl - ic)*du1_n(:, kmax - ic)
                    end do
                    du1_n(:, kmax) = du1_n(:, kmax) + g%der1%lu(kmax, ndl)*du1_n(:, kmax - ic)
                    du1_n(:, kmax) = du1_n(:, kmax)/g%der1%lu(kmax, idl)
                end if

                write (str, *) im
                call check(u, du1_a, du1_n, 'partial-'//trim(adjustl(str))//'.dat')
                ! call write_scheme(g%der1%lhs(:, 1:g%der1%nb_diag(1)), &
                !                   g%der1%rhs(:, 1:g%der1%nb_diag(2)), 'fdm1-'//trim(adjustl(str)))

            end do

        end do

#undef bcs_hb
#undef bcs_ht

    end select

    stop

contains
    ! ###################################################################
    ! ###################################################################
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
                    write (20, 1000) x%nodes(i), u(l, i), du_a(l, i), du_n(l, i), du_a(l, i) - du_n(l, i)
                end if
                dummy = dummy + du_a(l, i)*du_a(l, i)
                error_l2 = error_l2 + (du_a(l, i) - du_n(l, i))**2
                error_max = max(error_max, abs(du_a(l, i) - du_n(l, i)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Solution L2-norm ...........:', sqrt(dummy)/real(size(u), wp)
        if (dummy == 0.0_wp) return
        write (*, *) 'Relative Error L2-norm .....:', sqrt(error_l2)/sqrt(dummy)
        write (*, *) 'Relative Error Linf-norm ...:', error_max/maxval(abs(du1_a))

        return
1000    format(5(1x, e12.5))
    end subroutine check

    subroutine write_scheme(lhs, rhs, name)
        real(wp), intent(in) :: lhs(:, :)
        real(wp), intent(in) :: rhs(:, :)
        character(len=*), intent(in) :: name

        integer nx, n

        nx = size(lhs, 1)

        open (21, file=trim(adjustl(name))//'-lhs.dat')
        do n = 1, nx
            write (21, *) lhs(n, :)
        end do
        close (21)

        open (22, file=trim(adjustl(name))//'-rhs.dat')
        do n = 1, nx
            write (22, *) rhs(n, :)
        end do
        close (23)

        return
    end subroutine write_scheme

end program VPARTIAL

program VELLIPTIC
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: ifile, gfile
    use TLab_Time, only: itime
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: imax, jmax, kmax, inb_txc
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_err, ims_pro, ims_offset_i, ims_offset_j
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use FDM, only: FDM_Initialize
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Grid
    use IO_Fields
    use OPR_Fourier
    use OPR_Partial
    use OPR_Elliptic
    use FI_VECTORCALCULUS
    use Averages
    use Integration
    implicit none

    real(wp), dimension(:, :), allocatable :: bcs_hb, bcs_ht
    real(wp), dimension(:, :, :), pointer :: a, b, c, d, e, f

    integer(wi) type_of_operator, type_of_problem
    integer(wi) i, idsp, jdsp
    integer ibc, ib

    integer :: bcs_cases(1:4) = [BCS_DD, BCS_DN, BCS_ND, BCS_NN]
    character(len=32) :: bcs_names(4) = [character(len=32) :: &
                                         'Dirichlet/Dirichlet', &
                                         'Dirichlet/Neumann', &
                                         'Neumann/Dirichlet', &
                                         'Neumann/Neumann']

    character(len=32) str
    real(wp) mean, lambda
    real(wp) wk, x_0
    real(wp) delta

! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)

    inb_txc = 8

    call TLab_Initialize_Memory(__FILE__)

    allocate (bcs_ht(imax, jmax), bcs_hb(imax, jmax))

    a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    d(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    e(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 7)
    f(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 8)

    call OPR_Partial_Initialize(ifile)
    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

    type_of_operator = 1    ! default
#ifndef USE_MPI
    print *, '1. Poisson routines'
    print *, '2. Helmholtz routines'
    read (*, *) type_of_operator

    if (type_of_operator == 2) then
        write (*, *) 'Eigenvalue (negative)?'
        read (*, *) lambda
        if (lambda > 0.0_wp) then
            print *, 'Eigenvalue for Helmholtz operator needs to be negative'
            call TLab_Stop(0)
        end if
    end if
#endif

    ! type_of_problem = 1     ! the forcing in the rhs is given
    type_of_problem = 2     ! the field in the lhs is given

    ! ###################################################################
    ! Define the function and analytic derivatives
    x_0 = 0.75_wp
    wk = 1.0_wp

#ifdef USE_MPI
    idsp = ims_offset_i; jdsp = ims_offset_j
#else
    idsp = 0; jdsp = 0
#endif

    do i = 1, imax
        ! single-mode
        a(i, :, :) = 1.0_wp !+ sin(2.0_wp*pi_wp/x%scale*wk*x%nodes(i + idsp)) ! + pi_wp/4.0_wp)
    end do

    do i = 1, jmax
        ! single-mode
        a(:, i, :) = a(:, i, :)*cos(2.0_wp*pi_wp/y%scale*wk*y%nodes(i + jdsp)) ! + pi_wp/4.0_wp)
    end do

    do i = 1, kmax
        ! ! exponential
        ! u(:, i) = exp(-z%nodes(i)*wk)
        ! du1_a(:, i) = -wk*u(:, i)
        ! step
        ! u(:, i) = max(0.0_wp, (z%nodes(i) - z%nodes(kmax/2))*x_0)
        ! du1_a(:, i) = (1.0_wp + sign(1.0_wp, z%nodes(i) - z%nodes(kmax/2)))*0.5_wp*x_0
        ! ! tanh
        ! u(:, i) = x_0*log(1.0_wp + exp((z%nodes(i) - z%nodes(kmax/2))/x_0))
        ! du1_a(:, i) = 0.5_wp*(1.0_wp + tanh(0.5_wp*(z%nodes(i) - z%nodes(kmax/2))/x_0))
        ! Polynomial
        lambda = 4.0_wp
        a(:, :, i) = a(:, :, i)*(z%nodes(i)*wk)**lambda
        ! ! zero
        ! u(:, i) = 0.0_wp
        ! du1_a(:, i) = 0.0_wp
        ! ! delta-function
        ! u(:, i) = max(0.0_wp, 2.0_wp - real(i, wp))
        ! du1_a(:, i) = 0.0_wp
        ! du2_a(:, i) = 0.0_wp
    end do

    ! call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, a, params)

    ! call random_seed()
    ! call random_number(a)

    ! -------------------------------------------------------------------
    ! DC level at lower/upper boundary set to zero
    ! mean = AVG_IK(imax, jmax, kmax, 1, a)
    mean = AVG_IK(imax, jmax, kmax, kmax, a)
    a = a - mean

    ! Just single precision, no header to visualize easily
    io_datatype = IO_TYPE_SINGLE

    ! ###################################################################
    select case (type_of_problem)
    case (1) ! The input field a is used as rhs in lap b = a and we solve for b
        f = a

        ! ! remove 2\Delta x wave
        ! call OPR_FILTER(imax, jmax, kmax, Dealiasing, a, txc)

        ibc = BCS_DD
        select case (ibc)
        case (BCS_DD)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
        case (BCS_DN)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
        case (BCS_ND)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
        case (BCS_NN)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
            ! For Neumann conditions, we need to satisfy the compatibility constraint dpdy_top-dpdy_bottom=int f
            mean = AVG_IK(imax, jmax, 1, 1, bcs_hb)
            call AVG_IK_V(imax, jmax, kmax, a, wrk1d(:, 1), wrk1d(:, 2))
            delta = mean + Int_Simpson(wrk1d(1:kmax, 1), z%nodes(1:kmax))
            mean = AVG_IK(imax, jmax, 1, 1, bcs_ht)
            bcs_ht = bcs_ht - mean + delta
        end select

        if (type_of_operator == 1) then
            call OPR_Poisson(imax, jmax, kmax, ibc, a, txc(:, 1), txc(:, 2), bcs_hb, bcs_ht)

        else if (type_of_operator == 2) then
            call OPR_Helmholtz(imax, jmax, kmax, ibc, lambda, a, txc(:, 1), txc(:, 2), bcs_hb, bcs_ht)

        end if

        ! -------------------------------------------------------------------
        ! With the calculated a, we calculate the b = lap a
        ! call OPR_Partial_X(OPR_P1, imax, jmax, kmax, a, c)
        ! call OPR_Partial_X(OPR_P1, imax, jmax, kmax, c, b)
        call OPR_Partial_X(OPR_P2_P1, imax, jmax, kmax, a, b, txc(:, 1))

        ! call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, a, c)
        ! call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, c, e)
        call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, a, e, txc(:, 1))
        b = b + e

        ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, a, c)
        ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, c, e)
        call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, a, e, txc(:, 1))
        b = b + e

        if (type_of_operator == 2) then
            b = b + lambda*a
        end if

        ! -------------------------------------------------------------------
        b(:, :, 1) = f(:, :, 1)         ! The boundary points do not satisfy the PDE
        b(:, :, kmax) = f(:, :, kmax)
        call IO_Write_Fields('field.out', imax, jmax, kmax, itime, 1, b)!, io_header_s(1:1))
        call check(f, b, txc(:, 1), 'field.dif')

! ###################################################################
    case (2) ! The input field a is used to construct the forcing term as lap a = f
        ! 1. option: calculate f = div grad a
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, a, b)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, a, c)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, a, d)
        call FI_INVARIANT_P(imax, jmax, kmax, b, c, d, f, txc(:, 1))
        f = -f

        ! ! 2. option: calculate f = lap a
        ! call OPR_Partial_X(OPR_P2_P1, imax, jmax, kmax, a, b, txc(:, 1))
        ! call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, a, c, txc(:, 1))
        ! call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, a, f, d)
        ! f = f + b + c

        if (type_of_operator == 2) then
            f = f + lambda*a
        end if

        do ib = 1, 4, 3 ! size(bcs_names)
            ibc = bcs_cases(ib)
            print *, new_line('a'), bcs_names(ib)

            select case (ibc)
            case (BCS_DD)
                bcs_hb(:, :) = a(:, :, 1); bcs_ht(:, :) = a(:, :, kmax)
                ! case (BCS_DN)         ! undeveloped
                !     bcs_hb(:, :) = a(:, :, 1); bcs_ht(:, :) = d(:, :, kmax)
                ! case (BCS_ND)
                !     bcs_hb(:, :) = d(:, :, 1); bcs_ht(:, :) = a(:, :, kmax)
            case (BCS_NN)
                bcs_hb(:, :) = d(:, :, 1); bcs_ht(:, :) = d(:, :, kmax)
            end select

            e = f       ! to save f for other cases in the loop
            if (type_of_operator == 1) then
                call OPR_Poisson(imax, jmax, kmax, ibc, e, txc(:, 1), txc(:, 2), bcs_hb, bcs_ht)

            else if (type_of_operator == 2) then
                call OPR_Helmholtz(imax, jmax, kmax, ibc, lambda, e, txc(:, 1), txc(:, 2), bcs_hb, bcs_ht)

            end if

            ! -------------------------------------------------------------------
            write (str, *) ib
            call check(a, e, txc(:, 1), 'elliptic-'//trim(adjustl(str)))

            call IO_Write_Fields('field.out', imax, jmax, kmax, itime, 1, e)!, io_header_s(1:1))

        end do

    end select

    call TLab_Stop(0)

! ###################################################################
contains
    subroutine check(a1, a2, dif, name)
        real(wp), intent(in) :: a1(*), a2(*)
        real(wp), intent(inout) :: dif(*)
        character(len=*), optional :: name

        real(wp) dummy, error
#ifdef USE_MPI
        real(wp) sum_mpi
#endif
        error = 0.0_wp
        dummy = 0.0_wp
        do i = 1, imax*jmax*kmax
            dif(i) = a2(i) - a1(i)
            error = error + dif(i)*dif(i)
            dummy = dummy + a1(i)*a1(i)
        end do
#ifdef USE_MPI
        sum_mpi = error
        call MPI_ALLREDUCE(sum_mpi, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        sum_mpi = dummy
        call MPI_ALLREDUCE(sum_mpi, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)

        if (ims_pro == 0) then
#endif
            write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
        end if
#endif

        if (present(name)) then
            call IO_Write_Fields(name, imax, jmax, kmax, itime, 1, dif)
        end if

        return
    end subroutine check

end program VELLIPTIC

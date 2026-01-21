#include "tlab_error.h"

! Split the routines into the ones that are initialized and the ones that not?
! If not initialized, you can enter with any kmax, but the periodic directions need to be the global ones because of OPR_Fourier.
module OPR_Elliptic
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: efile
    use TLab_Memory, only: TLab_Allocate_Real
    use TLab_Memory, only: imax, jmax, kmax, isize_txc_field
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk1d, wrk2d, wrk3d
    use TLab_Grid, only: x, y, z
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_j, ims_pro_i
#endif
    use FDM, only: fdm_der1_Z, FDM_CreatePlan_Der2
    use FDM_Derivative_2order_X
    use FDM_Integral
    use OPR_Fourier, only: OPR_Fourier_XY_Backward, OPR_Fourier_XY_Forward
    use OPR_ODES
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    public :: OPR_Elliptic_Initialize
    public :: OPR_Poisson
    public :: OPR_Helmholtz

    ! -----------------------------------------------------------------------
    procedure(OPR_Poisson_interface) :: OPR_Poisson_dt      ! Implicit pointer (Procedure type)
    abstract interface
        subroutine OPR_Poisson_interface(nx, ny, nz, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht)
            use TLab_Constants, only: wi, wp
            integer(wi), intent(in) :: nx, ny, nz
            integer, intent(in) :: ibc                                      ! Dirichlet/Neumman BCs at kmin/kmax: BCS_DD, BCS_ND, BCS_DN, BCS_NN
            real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
            real(wp), intent(inout), target :: tmp1(2*nz, nx/2 + 1, ny)
            real(wp), intent(inout), target :: tmp2(2*nz, nx/2 + 1, ny)
            real(wp), intent(in) :: bcs_hb(nx, ny), bcs_ht(nx, ny)          ! Boundary-condition fields
        end subroutine
    end interface
    procedure(OPR_Poisson_dt), pointer :: OPR_Poisson

    procedure(OPR_Helmholtz_interface) :: OPR_Helmholtz_dt  ! Implicit pointer (Procedure type)
    abstract interface
        subroutine OPR_Helmholtz_interface(nx, ny, nz, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
            use TLab_Constants, only: wi, wp
            integer(wi), intent(in) :: nx, ny, nz
            integer, intent(in) :: ibc                                      ! Dirichlet/Neumman BCs at kmin/kmax: BCS_DD, BCS_ND, BCS_DN, BCS_NN
            real(wp), intent(in) :: alpha
            real(wp), intent(inout) :: a(nx, ny, nz)                        ! Forcing term, and solution field a
            real(wp), intent(inout), target :: tmp1(2*nz, nx/2 + 1, ny)
            real(wp), intent(inout), target :: tmp2(2*nz, nx/2 + 1, ny)
            real(wp), intent(in) :: bcs_hb(nx, ny), bcs_ht(nx, ny)          ! Boundary-condition fields
        end subroutine
    end interface
    procedure(OPR_Helmholtz_dt), pointer :: OPR_Helmholtz

    real(wp) norm
    integer(wi) i_sing(2), j_sing(2)                                ! singular modes
    integer(wi) i, j, i_max, isize_line

    ! type(fdm_dt) fdm_loc                                            ! scheme used for the elliptic solvers
    class(der2_dt), allocatable :: fdm_der2

    type(fdm_integral_dt), allocatable :: fdm_int1(:, :, :)         ! factorized method
    real(wp), allocatable, target :: rhs_b(:, :), rhs_t(:, :)       ! rhs to free memory space
    type(fdm_integral_dt) :: fdm_int1_loc(2)

    type(fdm_integral_dt), allocatable :: fdm_int2(:, :)            ! direct method
    real(wp), allocatable, target :: rhs_d(:, :)                    ! rhs to free memory space
    type(fdm_integral_dt) :: fdm_int2_loc

    real(wp), allocatable :: lambda(:, :)

    complex(wp), pointer :: c_tmp1(:) => null(), c_tmp2(:) => null()
    real(wp), pointer :: p_wrk3d_loc(:, :, :) => null()

contains
    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Elliptic_Initialize(inifile)
        use FDM_Base, only: FDM_COM4_DIRECT, FDM_COM6_DIRECT
        use FDM, only: fdm_der1_X, fdm_der1_Y, fdm_der1_Z
        use FDM, only: fdm_der2_X, fdm_der2_Y, fdm_der2_Z
        use FDM_Derivative_1order_X, only: der1_periodic, der1_biased, FDM_Der1_ModifyWavenumbers
        use FDM_Derivative_2order_X, only: der2_periodic, der2_biased, FDM_Der2_ModifyWavenumbers

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer imode_elliptic
        integer, parameter :: TYPE_FACTORIZE = 1
        integer, parameter :: TYPE_DIRECT = 2

        integer(wi) :: ndl, ndr, nd
        character(len=32) bakfile, block
        character(len=512) sRes

        integer(wi) iglobal, jglobal
        integer(wi) fft_offset_i, fft_offset_j

        real(wp), allocatable :: mwn_x(:), mwn_y(:)

        ! ###################################################################
        ! Reading
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Space'

        imode_elliptic = TYPE_FACTORIZE     ! default is the finite-difference method used for the derivatives
        !                                   to impose zero divergence down to round-off error in the interior points
        !                                   If factorize type, it needs to be equal to the scheme used to calculate
        !                                   derivatives and thus Neumann boundary conditions,
        !                                   because of the compatibility constraint in pressure-Poisson equation

        call ScanFile_Char(bakfile, inifile, block, 'SchemeElliptic', 'void', sRes)
        ! call ScanFile_Char(bakfile, inifile, block, 'SchemeElliptic', 'compactdirect6', sRes)
        select case (trim(adjustl(sRes)))
        case ('void')
        case ('compactdirect4')
            imode_elliptic = TYPE_DIRECT
            call FDM_CreatePlan_Der2(z, fdm_der2, FDM_COM4_DIRECT, fdm_der1_Z)

        case ('compactdirect6')
            imode_elliptic = TYPE_DIRECT
            call FDM_CreatePlan_Der2(z, fdm_der2, FDM_COM6_DIRECT, fdm_der1_Z)

        case default
            call TLab_Write_ASCII(efile, __FILE__//'. Undeveloped SchemeElliptic.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end select

        ! ###################################################################
        ! Initializing
        isize_line = imax/2 + 1

        allocate (lambda(isize_line, jmax))
        norm = 1.0_wp/real(x%size*y%size, wp)

        select case (imode_elliptic)
        case (TYPE_FACTORIZE)
            OPR_Poisson => OPR_Poisson_FourierXZ_Factorize
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Factorize

            ndl = size(fdm_der1_Z%lhs, 2)
            ndr = size(fdm_der1_Z%rhs, 2)
            nd = ndl
            allocate (fdm_int1(2, isize_line, jmax))
            call TLab_Allocate_Real(__FILE__, rhs_b, [z%size, nd], 'rhs_b')
            call TLab_Allocate_Real(__FILE__, rhs_t, [z%size, nd], 'rhs_t')

            i_sing = [1, x%size/2 + 1]      ! global indexes, transformed below to task-local indexes.
            j_sing = [1, y%size/2 + 1]

            if (x%size > 1) &
                call FDM_Der1_ModifyWavenumbers(size(fdm_der1_X%lhs, 1), fdm_der1_X%lhs(1, :), fdm_der1_X%rhs(1, :), mwn_x)
            if (y%size > 1) &
                call FDM_Der1_ModifyWavenumbers(size(fdm_der1_Y%lhs, 1), fdm_der1_Y%lhs(1, :), fdm_der1_Y%rhs(1, :), mwn_y)

        case (TYPE_DIRECT)
            OPR_Poisson => OPR_Poisson_FourierXZ_Direct
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Direct

            ndl = size(fdm_der2%lhs, 2)
            ndr = size(fdm_der2%rhs, 2)
            nd = ndl
            allocate (fdm_int2(isize_line, jmax))
            call TLab_Allocate_Real(__FILE__, rhs_d, [z%size, nd], 'rhs_d')

            i_sing = [1, 1]                 ! 2nd order FDMs are non-zero at Nyquist
            j_sing = [1, 1]

            if (x%size > 1) &
                call FDM_Der2_ModifyWavenumbers(size(fdm_der2_X%lhs, 1), fdm_der2_X%lhs(1, :), fdm_der2_X%rhs(1, :), mwn_x)
            if (y%size > 1) &
                call FDM_Der2_ModifyWavenumbers(size(fdm_der2_Y%lhs, 1), fdm_der2_Y%lhs(1, :), fdm_der2_Y%rhs(1, :), mwn_y)

        end select

#ifdef USE_MPI
        fft_offset_i = ims_pro_i*isize_line
        fft_offset_j = ims_offset_j

#else
        fft_offset_i = 0
        fft_offset_j = 0
#endif

        i_sing = i_sing - [fft_offset_i, fft_offset_i]          ! Singular modes in task-local variables
        j_sing = j_sing - [fft_offset_j, fft_offset_j]
        i_max = min(x%size/2 + 1 - fft_offset_i, isize_line)    ! Maximum mode is x direction

        do i = 1, i_max
#ifdef USE_MPI

            iglobal = i + fft_offset_i
#else
            iglobal = i
#endif

            do j = 1, jmax
#ifdef USE_MPI
                jglobal = j + fft_offset_j
#else
                jglobal = j
#endif

                select case (imode_elliptic)
                case (TYPE_FACTORIZE)
                    ! Define \lambda based on modified wavenumbers (real)
                    if (y%size > 1) then
                        lambda(i, j) = mwn_x(iglobal)**2 + mwn_y(jglobal)**2
                    else
                        lambda(i, j) = mwn_x(iglobal)**2
                    end if

                    call FDM_Int1_Initialize(fdm_der1_Z, &
                                             sqrt(lambda(i, j)), BCS_MIN, fdm_int1(BCS_MIN, i, j))
                    ! ! multiply by the FFT normalization
                    ! idl = ndl/2 + 1
                    ! idr = ndr/2 + 1
                    ! fdm_int1(BCS_MIN, i, j)%rhs(:, :) = fdm_int1(BCS_MIN, i, j)%rhs(:, :)*norm
                    ! fdm_int1(BCS_MIN, i, j)%rhs(1, idl) = fdm_int1(BCS_MIN, i, j)%rhs(1, idl)/norm  ! not this one
                    ! fdm_int1(BCS_MIN, i, j)%rhs_b(:, :) = fdm_int1(BCS_MIN, i, j)%rhs_b(:, :)*norm
                    ! fdm_int1(BCS_MIN, i, j)%rhs_t(:, :) = fdm_int1(BCS_MIN, i, j)%rhs_t(:, :)*norm
                    ! fdm_int1(BCS_MIN, i, j)%lhs(1, idr) = fdm_int1(BCS_MIN, i, j)%lhs(1, idr)*norm

                    call FDM_Int1_Initialize(fdm_der1_Z, &
                                             -sqrt(lambda(i, j)), BCS_MAX, fdm_int1(BCS_MAX, i, j))
                    ! ! multiply by the FFT normalization
                    ! fdm_int1(BCS_MAX, i, j)%rhs(:, :) = fdm_int1(BCS_MAX, i, j)%rhs(:, :)*norm
                    ! fdm_int1(BCS_MAX, i, j)%rhs(z%size, idl) = fdm_int1(BCS_MAX, i, j)%rhs(z%size, idl)*norm ! not this one
                    ! fdm_int1(BCS_MAX, i, j)%rhs_b(:, :) = fdm_int1(BCS_MAX, i, j)%rhs_b(:, :)*norm
                    ! fdm_int1(BCS_MAX, i, j)%rhs_t(:, :) = fdm_int1(BCS_MAX, i, j)%rhs_t(:, :)*norm
                    ! fdm_int1(BCS_MAX, i, j)%lhs(z%size, idr) = fdm_int1(BCS_MAX, i, j)%lhs(z%size, idr)*norm

                    if (any(i_sing == i) .and. any(j_sing == j)) then
                    else                                        ! free memory that is independent of lambda
                        rhs_b(:, :) = fdm_int1(BCS_MIN, i, j)%rhs(:, :)
                        if (allocated(fdm_int1(BCS_MIN, i, j)%rhs)) deallocate (fdm_int1(BCS_MIN, i, j)%rhs)

                        rhs_t(:, :) = fdm_int1(BCS_MAX, i, j)%rhs(:, :)
                        if (allocated(fdm_int1(BCS_MAX, i, j)%rhs)) deallocate (fdm_int1(BCS_MAX, i, j)%rhs)

                    end if

                    ! idr = ndr/2 + 1
                    ! fdm_int1(BCS_MIN, i, j)%lhs(2:, idr) = fdm_int1(BCS_MIN, i, j)%lhs(2:, idr)*norm
                    ! fdm_int1(BCS_MAX, i, j)%lhs(:z%size - 1, idr) = fdm_int1(BCS_MAX, i, j)%lhs(:z%size - 1, idr)*norm

                case (TYPE_DIRECT)     ! only for case BCS_NN
                    ! Define \lambda based on modified wavenumbers (real)
                    if (y%size > 1) then
                        lambda(i, j) = mwn_x(iglobal) + mwn_y(jglobal)
                    else
                        lambda(i, j) = mwn_x(iglobal)
                    end if

                    ! Compatibility constraint. The reference value of p at the lower boundary will be set to zero
                    if (any(i_sing == i) .and. any(j_sing == j)) then
                        call FDM_Int2_Initialize(z%nodes(:), fdm_der2, lambda(i, j), BCS_DN, fdm_int2(i, j))
                    else
                        call FDM_Int2_Initialize(z%nodes(:), fdm_der2, lambda(i, j), BCS_NN, fdm_int2(i, j))
                    end if
                    ! multiply by the FFT normalization
                    fdm_int2(i, j)%rhs = fdm_int2(i, j)%rhs*norm
                    fdm_int2(i, j)%rhs_t1 = fdm_int2(i, j)%rhs_t1*norm
                    fdm_int2(i, j)%rhs_b1 = fdm_int2(i, j)%rhs_b1*norm

                    ! free memory that is independent of lambda
                    rhs_d(:, :) = fdm_int2(i, j)%rhs(:, :)
                    if (allocated(fdm_int2(i, j)%rhs)) deallocate (fdm_int2(i, j)%rhs)

                end select

            end do
        end do

        return
    end subroutine OPR_Elliptic_Initialize

    !########################################################################
    !#
    !# Solve Lap p = f using Fourier in xOy planes, to rewrite the problem as
    !#
    !#     \hat{p}''-\lambda \hat{p} = \hat{f}
    !#
    !# where \lambda = kx^2+ky^2
    !#
    !# The reference value of p at the lower boundary is set to zero
    !#
    !########################################################################
    subroutine OPR_Poisson_FourierXZ_Factorize(nx, ny, nz, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout), target :: tmp1(2*nz, nx/2 + 1, ny)
        real(wp), intent(inout), target :: tmp2(2*nz, nx/2 + 1, ny)
        real(wp), intent(in) :: bcs_hb(nx, ny), bcs_ht(nx, ny)          ! Boundary-condition fields

        ! -----------------------------------------------------------------------
        real(wp) bcs(2, 2)

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field/2])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field/2])
        p_wrk3d_loc(1:2*nz, 1:nx/2 + 1, 1:ny) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Construct forcing term in Fourier space, \hat{f}
        p(1:nx, 1:ny, 1) = bcs_hb(1:nx, 1:ny)               ! Add boundary conditions to forcing array
        p(1:nx, 1:ny, nz) = bcs_ht(1:nx, 1:ny)
        call OPR_Fourier_XY_Forward(p(:, 1, 1), c_tmp1, c_tmp2)

        ! ! multiply by the FFT normalization
        ! tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
#define f(k,i,j) tmp1(k,i,j)
#define u(k,i,j) tmp2(k,i,j)
#define v(k,i,j) p_wrk3d_loc(k,i,j)

        ! Solve for each (kx,ky) a system of 1 complex equation as 2 independent real equations
        do j = 1, ny
            do i = 1, i_max
                bcs(1:2, 1) = f(1:2, i, j)                  ! bottom boundary conditions
                bcs(1:2, 2) = f(2*nz - 1:2*nz, i, j)        ! top boundary conditions

                select case (ibc)
                case (BCS_NN)       ! Neumann & Neumann boundary conditions
                    if (any(i_sing == i) .and. any(j_sing == j)) then
                        call OPR_ODE2_Factorize_NN_Sing(2, fdm_int1(:, i, j), &
                                                        u(:, i, j), f(:, i, j), bcs, v(:, i, j), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_Factorize_NN(2, fdm_int1(:, i, j), rhs_b, rhs_t, &
                                                   u(:, i, j), f(:, i, j), bcs, v(:, i, j), wrk1d, wrk2d)
                    end if

                case (BCS_DD)       ! Dirichlet & Dirichlet boundary conditions
                    if (any(i_sing == i) .and. any(j_sing == j)) then
                        call OPR_ODE2_Factorize_DD_Sing(2, fdm_int1(:, i, j), &
                                                        u(:, i, j), f(:, i, j), bcs, v(:, i, j), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_Factorize_DD(2, fdm_int1(:, i, j), rhs_b, rhs_t, &
                                                   u(:, i, j), f(:, i, j), bcs, v(:, i, j), wrk1d, wrk2d)
                    end if

                end select

                ! multiply by the FFT normalization
                u(:, i, j) = u(:, i, j)*norm

            end do
        end do

        ! ###################################################################
        ! Transform solution to physical space
        call OPR_Fourier_XY_Backward(c_tmp2, p(:, 1, 1), c_tmp1)

        nullify (c_tmp1, c_tmp2, p_wrk3d_loc)
#undef f
#undef v
#undef u

        return
    end subroutine OPR_Poisson_FourierXZ_Factorize

    !########################################################################
    !########################################################################
    subroutine OPR_Poisson_FourierXZ_Direct(nx, ny, nz, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout), target :: tmp1(2*nz, nx/2 + 1, ny)
        real(wp), intent(inout), target :: tmp2(2*nz, nx/2 + 1, ny)
        real(wp), intent(in) :: bcs_hb(nx, ny), bcs_ht(nx, ny)          ! Boundary-condition fields

        ! -----------------------------------------------------------------------
        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field/2])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field/2])
        p_wrk3d_loc(1:2*nz, 1:nx/2 + 1, 1:ny) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Construct forcing term in Fourier space, \hat{f}
        p(1:nx, 1:ny, 1) = bcs_hb(1:nx, 1:ny)               ! Add boundary conditions to forcing array
        p(1:nx, 1:ny, nz) = bcs_ht(1:nx, 1:ny)
        call OPR_Fourier_XY_Forward(p(:, 1, 1), c_tmp1, c_tmp2)

        ! ###################################################################
        ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
#define f(k,i,j) tmp1(k,i,j)
#define u(k,i,j) tmp2(k,i,j)

        ! Solve for each (kx,ky) a system of 1 complex equation as 2 independent real equations
        do j = 1, ny
            do i = 1, i_max
                u(1:2, i, j) = f(1:2, i, j)                         ! bottom boundary conditions
                u(2*nz - 1:2*nz, i, j) = f(2*nz - 1:2*nz, i, j)     ! top boundary conditions

                select case (ibc)
                case (BCS_NN)           ! use precalculated LU factorization
                    ! Compatibility constraint for singular modes. The reference value of p at bottom is set to zero
                    if (any(i_sing == i) .and. any(j_sing == j)) u(1:2, i, j) = 0.0_wp

                    call FDM_Int2_Solve(2, fdm_int2(i, j), rhs_d, f(:, i, j), u(:, i, j), wrk2d)

                    ! multiply by the FFT normalization
                    ! already done at initialization

                case default            ! Need to calculate and factorize LHS
                    ! call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2, lambda(i, j), ibc, fdm_int2_loc)
                    call FDM_Int2_Initialize(z%nodes(:), fdm_der2, lambda(i, j), ibc, fdm_int2_loc)
                    call FDM_Int2_Solve(2, fdm_int2_loc, fdm_int2_loc%rhs, f(:, i, j), u(:, i, j), wrk2d)

                    ! multiply by the FFT normalization
                    u(:, i, j) = u(:, i, j)*norm

                end select

            end do
        end do

        ! ###################################################################
        ! Transform solution to physical space
        call OPR_Fourier_XY_Backward(c_tmp2, p(:, 1, 1), c_tmp1)

        nullify (c_tmp1, c_tmp2, p_wrk3d_loc)
#undef f
#undef u

        return
    end subroutine OPR_Poisson_FourierXZ_Direct

    !########################################################################
    !#
    !# Solve Lap a + alpha a = f using Fourier in xOy planes, to rewrite the problem as
    !#
    !#      \hat{a}''-(\lambda-alpha) \hat{a} = \hat{f}
    !#
    !# where \lambda = kx^2+ky^2
    !#
    !########################################################################
    subroutine OPR_Helmholtz_FourierXZ_Factorize(nx, ny, nz, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        use FDM, only: fdm_der1_Z
        use FDM_Derivative_1order_X, only: der1_biased
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                        ! Forcing term, and solution field
        real(wp), intent(inout), target :: tmp1(2*nz, nx/2 + 1, ny)
        real(wp), intent(inout), target :: tmp2(2*nz, nx/2 + 1, ny)
        real(wp), intent(in) :: bcs_hb(nx, ny), bcs_ht(nx, ny)          ! Boundary-condition fields

        ! -----------------------------------------------------------------------
        real(wp) bcs(2, 2)

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field/2])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field/2])
        p_wrk3d_loc(1:2*nz, 1:nx/2 + 1, 1:ny) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Construct forcing term in Fourier space, \hat{f}
        a(1:nx, 1:ny, 1) = bcs_hb(1:nx, 1:ny)               ! Add boundary conditions to forcing array
        a(1:nx, 1:ny, nz) = bcs_ht(1:nx, 1:ny)
        call OPR_Fourier_XY_Forward(a(:, 1, 1), c_tmp1, c_tmp2)

        ! ! multiply by the FFT normalization
        ! tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE (\hat{p}')'-(\lambda + alpha) \hat{p} = \hat{f}
#define f(k,i,j) tmp1(k,i,j)
#define u(k,i,j) tmp2(k,i,j)
#define v(k,i,j) p_wrk3d_loc(k,i,j)

        ! Solve for each (kx,ky) a system of 1 complex equation as 2 independent real equations
        do i = 1, i_max
            do j = 1, ny
                bcs(1:2, 1) = f(1:2, i, j)                  ! bottom boundary conditions
                bcs(1:2, 2) = f(2*nz - 1:2*nz, i, j)        ! top boundary conditions

                ! call FDM_Int1_Initialize(fdm_loc%der1, &
                call FDM_Int1_Initialize(fdm_der1_Z, &
                                         sqrt(lambda(i, j) - alpha), BCS_MIN, fdm_int1_loc(BCS_MIN))

                ! call FDM_Int1_Initialize(fdm_loc%der1, &
                call FDM_Int1_Initialize(fdm_der1_Z, &
                                         -sqrt(lambda(i, j) - alpha), BCS_MAX, fdm_int1_loc(BCS_MAX))

                select case (ibc)
                case (BCS_NN)
                    call OPR_ODE2_Factorize_NN(2, fdm_int1_loc, fdm_int1_loc(BCS_MIN)%rhs, fdm_int1_loc(BCS_MAX)%rhs, &
                                               u(:, i, j), f(:, i, j), bcs, v(:, i, j), wrk1d, wrk2d)

                case (BCS_DD)
                    call OPR_ODE2_Factorize_DD(2, fdm_int1_loc, fdm_int1_loc(BCS_MIN)%rhs, fdm_int1_loc(BCS_MAX)%rhs, &
                                               u(:, i, j), f(:, i, j), bcs, v(:, i, j), wrk1d, wrk2d)
                end select

                ! multiply by the FFT normalization
                u(:, i, j) = u(:, i, j)*norm

            end do
        end do

        ! ###################################################################
        ! Transform solution to physical space
        call OPR_Fourier_XY_Backward(c_tmp2, a(:, 1, 1), c_tmp1)

        nullify (c_tmp1, c_tmp2, p_wrk3d_loc)
#undef f
#undef v
#undef u

        return
    end subroutine OPR_Helmholtz_FourierXZ_Factorize

    !########################################################################
    !########################################################################
    subroutine OPR_Helmholtz_FourierXZ_Direct(nx, ny, nz, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                        ! Forcing term, and solution field
        real(wp), intent(inout), target :: tmp1(2*nz, nx/2 + 1, ny)
        real(wp), intent(inout), target :: tmp2(2*nz, nx/2 + 1, ny)
        real(wp), intent(in) :: bcs_hb(nx, ny), bcs_ht(nx, ny)          ! Boundary-condition fields

        ! -----------------------------------------------------------------------
        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field/2])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field/2])
        p_wrk3d_loc(1:2*nz, 1:nx/2 + 1, 1:ny) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Construct forcing term in Fourier space, \hat{f}
        a(1:nx, 1:ny, 1) = bcs_hb(1:nx, 1:ny)               ! Add boundary conditions to forcing array
        a(1:nx, 1:ny, nz) = bcs_ht(1:nx, 1:ny)
        call OPR_Fourier_XY_Forward(a(:, 1, 1), c_tmp1, c_tmp2)

        ! ! multiply by the FFT normalization
        ! tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE \hat{p}''-(\lambda + alpha) \hat{p} = \hat{f}
#define f(k,i,j) tmp1(k,i,j)
#define u(k,i,j) tmp2(k,i,j)

        ! Solve for each (kx,ky) a system of 1 complex equation as 2 independent real equations
        do i = 1, i_max
            do j = 1, ny
                u(1:2, i, j) = f(1:2, i, j)
                u(2*nz - 1:2*nz, i, j) = f(2*nz - 1:2*nz, i, j)

                ! call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2, lambda(i, j) - alpha, ibc, fdm_int2_loc)
                call FDM_Int2_Initialize(z%nodes(:), fdm_der2, lambda(i, j) - alpha, ibc, fdm_int2_loc)
                call FDM_Int2_Solve(2, fdm_int2_loc, fdm_int2_loc%rhs, f(:, i, j), u(:, i, j), wrk2d)

                ! multiply by the FFT normalization
                u(:, i, j) = u(:, i, j)*norm

            end do
        end do

        ! ###################################################################
        ! Transform solution to physical space
        call OPR_Fourier_XY_Backward(c_tmp2, a(:, 1, 1), c_tmp1)

        nullify (c_tmp1, c_tmp2, p_wrk3d_loc)
#undef f
#undef u

        return
    end subroutine OPR_Helmholtz_FourierXZ_Direct

end module OPR_Elliptic

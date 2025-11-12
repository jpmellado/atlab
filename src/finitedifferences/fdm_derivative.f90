#include "tlab_error.h"

module FDM_Derivative
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: lfile, efile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_NONE, BCS_PERIODIC, BCS_MIN, BCS_MAX
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas
    use Thomas_Circulant
    use MatmulDevel
    use Matmul_Thomas
    use Matmul_Halo
    use Matmul_Halo_Thomas
    use Preconditioning
    use FDM_Base
    use FDM_ComX_Direct
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use FDM_Com0_Jacobian
    implicit none
    private

    type, public :: fdm_derivative_dt
        sequence
        integer mode_fdm                            ! finite-difference method
        integer(wi) size                            ! # of grid points, for convenience in the code
        logical :: periodic = .false.
        logical :: need_1der = .false.              ! In nonuniform, Jacobian formulation, need 1. order derivative for the 2. order one
        integer nb_diag(2)                          ! # of left and right diagonals  (max 5/7)
        real(wp) :: rhs_b(4 + 1, 0:7), rhs_t(0:4, 7)    ! Neumann boundary conditions, max. # of diagonals is 7, # rows is 7/2+1
        real(wp), allocatable :: lhs(:, :)          ! memory space for LHS
        real(wp), allocatable :: rhs(:, :)          ! memory space for RHS
        real(wp), allocatable :: rhs_b1(:, :), rhs_t1(:, :)
        real(wp), allocatable :: mwn(:)             ! memory space for modified wavenumbers
        !
        real(wp), allocatable :: lu(:, :)           ! memory space for LU decomposition

        ! procedure(matmul_halo_ice), pointer, nopass :: matmul_halo => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul_halo_thomas => null()
        ! procedure(matmuldevel_ice), pointer, nopass :: matmuldevel => null()
        ! procedure(matmuldevel_add_ice), pointer, nopass :: matmuldevel_add => null()
        procedure(matmuldevel_thomas_ice), pointer, nopass :: matmuldevel_thomas => null()
        procedure(matmuldevel_add_thomas_ice), pointer, nopass :: matmuldevel_add_thomas => null()
        procedure(thomas_ice), pointer, nopass :: thomasL => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()

        ! type(fdm_bcs) :: bcs(0:3)           ! linear system for 4 different boundary conditions, 0 is the default
    end type fdm_derivative_dt

    ! type :: fdm_bcs
    !     integer bcs_type
    !     real(wp), allocatable :: lu(:, :)               ! memory space for LU decomposition
    !     real(wp) :: rhs_b(4 + 1, 0:7), rhs_t(0:4, 7)    ! Neumann boundary conditions, max. # of diagonals is 7, # rows is 7/2+1
    ! end type fdm_bcs

    public :: FDM_Der1_Initialize
    ! public :: FDM_Der1_CreateSystem
    public :: FDM_Der1_Solve

    public :: FDM_Der2_Initialize
    ! public :: FDM_Der2_CreateSystem
    public :: FDM_Der2_Solve

    integer, parameter, public :: FDM_NONE = 0

    integer, parameter, public :: FDM_COM4_JACOBIAN = 1
    integer, parameter, public :: FDM_COM6_JACOBIAN = 2
    integer, parameter, public :: FDM_COM6_JACOBIAN_HYPER = 3
    integer, parameter, public :: FDM_COM6_JACOBIAN_PENTA = 4

    integer, parameter, public :: FDM_COM4_DIRECT = 11
    integer, parameter, public :: FDM_COM6_DIRECT = 12
    integer, parameter, public :: FDM_COM6_DIRECT_HYPER = 13

    ! -----------------------------------------------------------------------
    ! abstract interface
    !     subroutine matmul_ice(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
    !         use TLab_Constants, only: wp
    !         real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
    !         real(wp), intent(in) :: u(:, :)                                 ! vector u
    !         real(wp), intent(out) :: f(:, :)                                ! vector f = B u
    !         integer, intent(in) :: ibc
    !         real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
    !         real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)
    !     end subroutine
    ! end interface

    abstract interface
        subroutine matmul_halo_ice(rhs, u, u_halo_m, u_halo_p, f)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:)              ! diagonals of B
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B u
        end subroutine
    end interface

    abstract interface
        subroutine thomas_ice(A, f)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: A(:, :)
            real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        end subroutine
    end interface

    abstract interface
        subroutine matmul_halo_thomas_ice(rhs, u, u_halo_m, u_halo_p, f, L)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:)              ! diagonals of B
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B u
            real(wp), intent(in) :: L(:, :)
        end subroutine
    end interface

    abstract interface
        subroutine matmuldevel_ice(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine
    end interface

    abstract interface
        subroutine matmuldevel_add_ice(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(in) :: rhs_add(:, :)
            real(wp), intent(in) :: u_add(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine
    end interface

    abstract interface
        subroutine matmuldevel_thomas_ice(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(in) :: L(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine
    end interface

    abstract interface
        subroutine matmuldevel_add_thomas_ice(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, L, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(in) :: rhs_add(:, :)
            real(wp), intent(in) :: u_add(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
            real(wp), intent(in) :: L(:, :)
        end subroutine
    end interface

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_Initialize(x, dx, g, periodic, bcs_cases)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(:)                   ! Jacobian
        type(fdm_derivative_dt), intent(inout) :: g
        logical, intent(in) :: periodic
        integer, intent(in) :: bcs_cases(:)

        ! -------------------------------------------------------------------
        integer ndl, ndr, idl, idr, i
        integer(wi) ib, ip
        integer(wi) nmin, nmax, nsize

        ! ###################################################################
        call FDM_Der1_CreateSystem(x, dx, g, periodic)

        ! -------------------------------------------------------------------
        ndl = g%nb_diag(1)                      ! number of diagonals in lhs
        idl = ndl/2 + 1
        ndr = g%nb_diag(2)                      ! number of diagonals in rhs
        idr = ndr/2 + 1

        ! Preconditioning
        call Precon_Rhs(g%lhs(:, 1:ndl), g%rhs(:, 1:ndr), periodic=periodic)

        ! LU decomposition
        if (allocated(g%lu)) deallocate (g%lu)
        if (g%periodic) then
            allocate (g%lu(g%size, ndl + 2))
        else
            allocate (g%lu(g%size, 5*size(bcs_cases)))
        end if
        g%lu(:, :) = 0.0_wp

        g%rhs_b(:, :) = 0.0_wp
        g%rhs_t(:, :) = 0.0_wp
        ! new format; extending to ndr+2 diagonals
        if (allocated(g%rhs_b1)) deallocate (g%rhs_b1)
        allocate (g%rhs_b1(max(idl, idr + 1), 1:ndr + 2))
        if (allocated(g%rhs_t1)) deallocate (g%rhs_t1)
        allocate (g%rhs_t1(max(idl, idr + 1), 1:ndr + 2))
        g%rhs_b1(:, :) = 0.0_wp
        g%rhs_t1(:, :) = 0.0_wp

        if (periodic) then
            g%lu(:, 1:ndl) = g%lhs(:, 1:ndl)

            select case (ndl)
            case (3)
                call ThomasCirculant_3_Initialize(g%lu(:, 1:ndl/2), &
                                                  g%lu(:, ndl/2 + 1:ndl), &
                                                  g%lu(1, ndl + 1))
            case (5)
                call ThomasCirculant_5_Initialize(g%lu(:, 1:ndl/2), &
                                                  g%lu(:, ndl/2 + 1:ndl), &
                                                  g%lu(1, ndl + 1))
            end select

        else                            ! biased,  different BCs
            do ib = 1, size(bcs_cases)
                ip = (ib - 1)*5

                call FDM_Der1_Neumann_Reduce(g%lhs(:, 1:ndl), g%rhs(:, 1:ndr), bcs_cases(ib), g%lu(:, ip + 1:ip + ndl), g%rhs_b, g%rhs_t)

                ! new format; extending to ndr+2 diagonals
                do i = 1, max(idl, idr + 1)
                    g%rhs_b1(i, 1:ndr + 1) = g%rhs_b(i, 0:ndr)
                    g%rhs_t1(i, 2:ndr + 2) = g%rhs_t(i - 1, 1:ndr + 1)
                end do
                ! longer stencil
                i = 1
                g%rhs_b1(i, ndr + 2) = g%rhs_b1(i, 2); g%rhs_b1(i, 2) = 0.0_wp
                i = max(idl, idr + 1)
                g%rhs_t1(i, 1) = g%rhs_t1(i, ndr + 1); g%rhs_t1(i, ndr + 1) = 0.0_wp

                nmin = 1; nmax = g%size
                if (any([BCS_ND, BCS_NN] == bcs_cases(ib))) nmin = nmin + 1
                if (any([BCS_DN, BCS_NN] == bcs_cases(ib))) nmax = nmax - 1
                nsize = nmax - nmin + 1

                call Thomas_FactorLU_InPlace(g%lu(nmin:nmax, ip + 1:ip + ndl/2), &
                                             g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl))
            end do

        end if

        ! -------------------------------------------------------------------
        ! Procedure pointers to linear solvers
        select case (ndl)
        case (3)
            g%thomasL => Thomas3_SolveL
            g%thomasU => Thomas3_SolveU
        case (5)
            g%thomasL => Thomas5_SolveL
            g%thomasU => Thomas5_SolveU
        case (7)
            g%thomasL => Thomas7_SolveL
            g%thomasU => Thomas7_SolveU
        end select

        ! -------------------------------------------------------------------
        ! Procedure pointers to matrix multiplication to calculate the right-hand side
        if (periodic) then
            select case (ndr)
            case (3)
                ! g%matmul_halo => MatMul_Halo_3d_antisym
                if (ndl == 3) g%matmul_halo_thomas => MatMul_Halo_3_antisym_ThomasL_3
            case (5)
                ! g%matmul_halo => MatMul_Halo_5d_antisym
                if (ndl == 3) g%matmul_halo_thomas => MatMul_Halo_5_antisym_ThomasL_3
                if (ndl == 5) g%matmul_halo_thomas => MatMul_Halo_5_antisym_ThomasL_5
            case (7)
                ! g%matmul_halo => MatMul_Halo_7d_antisym
                if (ndl == 3) g%matmul_halo_thomas => MatMul_Halo_7_antisym_ThomasL_3
                if (ndl == 5) g%matmul_halo_thomas => MatMul_Halo_7_antisym_ThomasL_5
            end select

        else
            if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm)) then
                select case (ndr)
                case (3)
                    ! g%matmuldevel => MatMul_3
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_3_ThomasL_3
                case (5)
                    ! g%matmuldevel => MatMul_5
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_5_ThomasL_3
                    if (ndl == 5) g%matmuldevel_thomas => MatMul_5_ThomasL_5
                end select
            else
                select case (ndr)
                case (3)
                    ! g%matmuldevel => MatMul_3_antisym
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_3_antisym_ThomasL_3
                case (5)
                    ! g%matmuldevel => MatMul_5_antisym
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_5_antisym_ThomasL_3
                    if (ndl == 5) g%matmuldevel_thomas => MatMul_5_antisym_ThomasL_5
                case (7)
                    ! g%matmuldevel => MatMul_7_antisym
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_7_antisym_ThomasL_3
                    if (ndl == 5) g%matmuldevel_thomas => MatMul_7_antisym_ThomasL_5
                end select
            end if

        end if

        return
    end subroutine FDM_Der1_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_CreateSystem(x, dx, g, periodic)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(:)                   ! Jacobian
        type(fdm_derivative_dt), intent(inout) :: g
        logical, intent(in) :: periodic

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx
        integer, parameter :: ndl_max = 5, ndr_max = 7

        ! ###################################################################
        g%size = size(x)                ! # grid points
        nx = g%size                     ! for code readability

        if (allocated(g%lhs)) deallocate (g%lhs)
        if (allocated(g%rhs)) deallocate (g%rhs)
        if (allocated(g%mwn)) deallocate (g%mwn)
        allocate (g%lhs(nx, ndl_max))
        allocate (g%rhs(nx, ndr_max))
        allocate (g%mwn(nx))
        g%lhs(:, :) = 0.0_wp
        g%rhs(:, :) = 0.0_wp

        g%periodic = periodic

        ! -------------------------------------------------------------------
        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C1N4_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)

        case (FDM_COM6_DIRECT)
            call FDM_C1N6_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (periodic) then
            nx = g%size

#define wn(i) g%mwn(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            g%mwn(:) = 2.0_wp*(coef(3)*sin(wn(:)) + coef(4)*sin(2.0_wp*wn(:)) + coef(5)*sin(3.0_wp*wn(:))) &
                       /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(wn(:)))

#undef wn

        end if

        return
    end subroutine FDM_Der1_CreateSystem

! #######################################################################
! #######################################################################
    subroutine FDM_Der1_Neumann_Reduce(lhs, rhs, ibc, r_lhs, r_rhs_b, r_rhs_t)
        real(wp), intent(in) :: lhs(:, :)
        real(wp), intent(in) :: rhs(:, :)
        integer, intent(in) :: ibc
        real(wp), intent(out) :: r_lhs(:, :)                        ! new, reduced lhs
        real(wp), intent(inout) :: r_rhs_b(:, 0:), r_rhs_t(0:, :)   ! new, reduced rhs

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, nx
        real(wp), allocatable :: aux(:, :)
        real(wp) locRhs_b(5, 0:7), locRhs_t(0:4, 7)

        ! ###################################################################
        ndl = size(lhs, 2)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = size(lhs, 1)           ! # grid points

        ! For A_22, we need idl >= idr -1
        if (idl < idr - 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. LHS array is too small.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        ! For b_21, we need idr >= idl
        if (idr < idl) then
            call TLab_Write_ASCII(efile, __FILE__//'. RHS array is too small.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        if (allocated(aux)) deallocate (aux)
        allocate (aux(1:nx, 1:ndr))

        ! -------------------------------------------------------------------
        r_lhs(:, 1:ndl) = lhs(:, 1:ndl)
        aux(1:nx, 1:ndr) = rhs(1:nx, 1:ndr)       ! array changed in FDM_Bcs_Reduce

        locRhs_b = 0.0_wp
        locRhs_t = 0.0_wp
        call FDM_Bcs_Reduce(ibc, aux, lhs(:, 1:ndl), locRhs_b, locRhs_t)

        ! reorganize data
        if (any([BCS_ND, BCS_NN] == ibc)) then
            r_rhs_b = 0.0_wp
            r_rhs_b(1:idr + 1, 1:ndr) = aux(1:idr + 1, 1:ndr)
            r_rhs_b(1, idr) = lhs(1, idl)       ! save a_11 for nonzero bc
            do ir = 1, idr - 1                  ! save -a^R_{21} for nonzero bc
                r_rhs_b(1 + ir, idr - ir) = -locRhs_b(1 + ir, idl - ir)
            end do

            r_lhs(2:idl + 1, 1:ndl) = locRhs_b(2:idl + 1, 1:ndl)
            r_lhs(1, idl) = rhs(1, idr)

        end if

        if (any([BCS_DN, BCS_NN] == ibc)) then
            r_rhs_t = 0.0_wp
            r_rhs_t(0:idr, 1:ndr) = aux(nx - idr:nx, 1:ndr)
            r_rhs_t(idr, idr) = lhs(nx, idl)
            do ir = 1, idr - 1              ! change sign in a^R_{21} for nonzero bc
                r_rhs_t(idr - ir, idr + ir) = -locRhs_t(idl - ir, idl + ir)
            end do

            r_lhs(nx - idl:nx - 1, 1:ndl) = locRhs_t(0:idl - 1, 1:ndl)
            r_lhs(nx, idl) = rhs(nx, idr)

        end if

        if (allocated(aux)) deallocate (aux)

        return
    end subroutine FDM_Der1_Neumann_Reduce

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_Solve(nlines, g, lu1, u, result, wrk2d, ibc)
        integer(wi), intent(in) :: nlines   ! # of lines to be solved
        type(fdm_derivative_dt), intent(in) :: g
        real(wp), intent(in) :: lu1(:, :)
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: wrk2d(nlines, 2)
        integer, intent(in), optional :: ibc          ! Boundary condition [BCS_DD=BCS_NONE, BCS_DN, BCS_ND, BCS_NN]

        ! -------------------------------------------------------------------
        integer(wi) nmin, nmax, nsize, ip
        integer ndl, ndr, idl, idr
        integer ibc_loc

        ! ###################################################################
        ndl = g%nb_diag(1)
        idl = ndl/2 + 1
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1

        if (g%periodic) then
            ! Calculate RHS in system of equations A u' = B u
            ! call g%matmul_halo(rhs=g%rhs(1, 1:ndr), &
            !                    u=u, &
            !                    u_halo_m=u(:, g%size - ndr/2 + 1:g%size), &
            !                    u_halo_p=u(:, 1:ndr/2), &
            !                    f=result)
            call g%matmul_halo_thomas(rhs=g%rhs(1, 1:ndr), &
                                      u=u, &
                                      u_halo_m=u(:, g%size - ndr/2 + 1:g%size), &
                                      u_halo_p=u(:, 1:ndr/2), &
                                      f=result, &
                                      L=lu1(:, 1:ndl/2))

            ! Solve for u' in system of equations A u' = B u
            ! call g%thomasL(lu1(:, 1:ndl/2), result)
            call g%thomasU(lu1(:, ndl/2 + 1:ndl), result)
            select case (g%nb_diag(1))
            case (3)
                call ThomasCirculant_3_Reduce(lu1(:, 1:ndl/2), &
                                              lu1(:, ndl/2 + 1:ndl), &
                                              lu1(:, ndl + 1), &
                                              result, wrk2d(:, 1))
            case (5)
                call ThomasCirculant_5_Reduce(lu1(:, 1:ndl/2), &
                                              lu1(:, ndl/2 + 1:ndl), &
                                              lu1(:, ndl + 1), &
                                              result)!, wrk2d)
            end select

        else    ! biased
            if (present(ibc)) then
                ibc_loc = ibc
            else
                ibc_loc = BCS_NONE
            end if

#define bcs_hb(i) wrk2d(i,1)
#define bcs_ht(i) wrk2d(i,2)

            nmin = 1; nmax = g%size
            if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                result(:, 1) = 0.0_wp                   ! homogeneous Neumann bcs
                bcs_hb(:) = 0.0_wp
                nmin = nmin + 1
            end if
            if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                result(:, g%size) = 0.0_wp              ! homogeneous Neumann bcs
                bcs_ht(:) = 0.0_wp
                nmax = nmax - 1
            end if
            nsize = nmax - nmin + 1

            ip = ibc_loc*5

            select case (ibc_loc)
            case (BCS_DD)
                ! call g%matmuldevel(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                    u=u, &
                !                    f=result)
                call g%matmuldevel_thomas(rhs=g%rhs(:, 1:ndr), &
                                          rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                                          rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                          u=u, &
                                          f=result, &
                                          L=lu1(:, ip + 1:ip + ndl/2))
            case (BCS_ND)
                ! call g%matmuldevel(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                    u=u, &
                !                    f=result, &
                !                    bcs_b=bcs_hb(:))
                call g%matmuldevel_thomas(rhs=g%rhs(:, 1:ndr), &
                                          rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                          rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                          u=u, &
                                          f=result, &
                                          L=lu1(:, ip + 1:ip + ndl/2), &
                                          bcs_b=bcs_hb(:))
            case (BCS_DN)
                ! call g%matmuldevel(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                    rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    u=u, &
                !                    f=result, &
                !                    bcs_t=bcs_ht(:))
                call g%matmuldevel_thomas(rhs=g%rhs(:, 1:ndr), &
                                          rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                                          rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                          u=u, &
                                          f=result, &
                                          L=lu1(:, ip + 1:ip + ndl/2), &
                                          bcs_t=bcs_ht(:))
            case (BCS_NN)
                ! call g%matmuldevel(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    u=u, &
                !                    f=result, &
                !                    bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))
                call g%matmuldevel_thomas(rhs=g%rhs(:, 1:ndr), &
                                          rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                          rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                          u=u, &
                                          f=result, &
                                          L=lu1(:, ip + 1:ip + ndl/2), &
                                          bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))

            end select

#undef bcs_hb
#undef bcs_ht

            ! Solve for u' in system of equations A u' = B u
            ! call g%thomasL(lu1(nmin:nmax, ip + 1:ip + ndl/2), result(:, nmin:nmax))
            call g%thomasU(lu1(nmin:nmax, ip + ndl/2 + 1:ip + ndl), result(:, nmin:nmax))

        end if

        return
    end subroutine FDM_Der1_Solve

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der2_Initialize(x, dx, g, periodic, uniform)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(:, :)                ! Jacobians
        type(fdm_derivative_dt), intent(inout) :: g     ! fdm plan for 2. order derivative
        logical, intent(in) :: periodic, uniform

        integer ndl, ndr

        ! ###################################################################
        call FDM_Der2_CreateSystem(x, dx, g, periodic, uniform)

        ! -------------------------------------------------------------------
        ! LU decomposition
        ndl = g%nb_diag(1)                              ! number of diagonals in lhs
        ndr = g%nb_diag(2)                              ! number of diagonals in rhs

        if (allocated(g%lu)) deallocate (g%lu)
        if (g%periodic) then
            allocate (g%lu(g%size, ndl + 2))
        else
            allocate (g%lu(g%size, ndl*1))              ! Only 1 bcs
        end if
        g%lu(:, :) = 0.0_wp

        g%lu(:, 1:ndl) = g%lhs(:, 1:ndl)
        if (g%periodic) then
            select case (ndl)
            case (3)
                call ThomasCirculant_3_Initialize(g%lu(:, 1:ndl/2), &
                                                  g%lu(:, ndl/2 + 1:ndl), &
                                                  g%lu(1, ndl + 1))
            end select

        else
            call Thomas_FactorLU_InPlace(g%lu(:, 1:ndl/2), &
                                         g%lu(:, ndl/2 + 1:ndl))

        end if

        ! -------------------------------------------------------------------
        ! Procedure pointers to linear solvers
        select case (ndl)
        case (3)
            g%thomasL => Thomas3_SolveL
            g%thomasU => Thomas3_SolveU
        case (5)
            g%thomasL => Thomas5_SolveL
            g%thomasU => Thomas5_SolveU
        case (7)
            g%thomasL => Thomas7_SolveL
            g%thomasU => Thomas7_SolveU
        end select

        ! -------------------------------------------------------------------
        ! Procedure pointers to matrix multiplication to calculate the right-hand side
        if (periodic) then
            select case (ndr)
            case (3)
                ! g%matmul_halo => MatMul_Halo_3d_sym
                if (ndl == 3) g%matmul_halo_thomas => MatMul_Halo_3_sym_ThomasL_3
            case (5)
                ! g%matmul_halo => MatMul_Halo_5d_sym
                if (ndl == 3) g%matmul_halo_thomas => MatMul_Halo_5_sym_ThomasL_3
                if (ndl == 5) g%matmul_halo_thomas => MatMul_Halo_5_sym_ThomasL_5
            case (7)
                ! g%matmul_halo => MatMul_Halo_7d_sym
                if (ndl == 3) g%matmul_halo_thomas => MatMul_Halo_7_sym_ThomasL_3
                if (ndl == 5) g%matmul_halo_thomas => MatMul_Halo_7_sym_ThomasL_5
            end select

        else
            if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT, FDM_COM6_DIRECT_HYPER] == g%mode_fdm)) then
                select case (ndr)
                case (5)
                    ! g%matmuldevel => MatMul_5
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_5_ThomasL_3
                    if (ndl == 5) g%matmuldevel_thomas => MatMul_5_ThomasL_5
                end select
            else
                select case (ndr)
                case (5)
                    ! g%matmuldevel => MatMul_5_sym
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_5_sym_ThomasL_3
                    if (ndl == 5) g%matmuldevel_thomas => MatMul_5_sym_ThomasL_5
                    ! g%matmuldevel_add => MatMul_5_sym_add_3
                    if (ndl == 3) g%matmuldevel_add_thomas => MatMul_5_sym_add_3_ThomasL_3
                    if (ndl == 5) g%matmuldevel_add_thomas => MatMul_5_sym_add_3_ThomasL_5
                case (7)
                    ! g%matmuldevel => MatMul_7_sym
                    if (ndl == 3) g%matmuldevel_thomas => MatMul_7_sym_ThomasL_3
                    if (ndl == 5) g%matmuldevel_thomas => MatMul_7_sym_ThomasL_5
                    ! g%matmuldevel_add => MatMul_7_sym_add_3
                    if (ndl == 3) g%matmuldevel_add_thomas => MatMul_7_sym_add_3_ThomasL_3
                    if (ndl == 5) g%matmuldevel_add_thomas => MatMul_7_sym_add_3_ThomasL_5
                end select
            end if

        end if

        return
    end subroutine FDM_Der2_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der2_CreateSystem(x, dx, g, periodic, uniform)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(:, :)                ! Jacobians
        type(fdm_derivative_dt), intent(inout) :: g     ! fdm plan for 2. order derivative
        logical, intent(in) :: periodic, uniform

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx
        integer, parameter :: ndl_max = 5, ndr_max = 7

        ! ###################################################################
        g%size = size(x)                ! # grid points
        nx = g%size                     ! for code readability

        if (allocated(g%lhs)) deallocate (g%lhs)
        if (allocated(g%rhs)) deallocate (g%rhs)
        if (allocated(g%mwn)) deallocate (g%mwn)
        allocate (g%lhs(nx, ndl_max))
        allocate (g%rhs(nx, ndr_max + ndl_max))     ! ndl_max is space for du correction in nonuniform case
        allocate (g%mwn(nx))
        g%lhs(:, :) = 0.0_wp
        g%rhs(:, :) = 0.0_wp

        g%periodic = periodic

        ! -------------------------------------------------------------------
        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)
            if (.not. uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN)
            call FDM_C2N6_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)
            if (.not. uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)
            if (.not. uniform) g%need_1der = .true.

        case (FDM_COM4_DIRECT)
            call FDM_C2N4_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)
            g%need_1der = .false.

        case (FDM_COM6_DIRECT)
            call FDM_C2N6_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)
            g%need_1der = .false.

        case (FDM_COM6_DIRECT_HYPER)
            call TLab_Write_ASCII(lfile, 'Direct, hyper-diffusive scheme undeveloped, use standard one.')
            call FDM_C2N6_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)
            g%need_1der = .false.

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (periodic) then

#define wn(i) g%mwn(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            g%mwn(:) = 2.0_wp*(coef(3)*(1.0_wp - cos(wn(:))) + coef(4)*(1.0_wp - cos(2.0_wp*wn(:))) + coef(5)*(1.0_wp - cos(3.0_wp*wn(:)))) &
                       /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(2.0_wp*wn(:)))

#undef wn

        end if

        return
    end subroutine FDM_Der2_CreateSystem

    ! ###################################################################################
    ! ###################################################################################
    subroutine FDM_Der2_Solve(nlines, g, lu, u, result, du, wrk2d)
        integer(wi), intent(in) :: nlines                   ! # of lines to be solved
        type(fdm_derivative_dt), intent(in) :: g            ! plan for 2. order derivative
        real(wp), intent(in) :: lu(:, :)
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(in) :: du(nlines, g%size)          ! 1. derivative for correction in case of Jacobian formulation
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(out) :: wrk2d(nlines)

        ! -------------------------------------------------------------------
        integer(wi) ip, ndl, ndr

        ! ###################################################################
        ndl = g%nb_diag(1)
        ndr = g%nb_diag(2)

        if (g%periodic) then
            ! Calculate RHS in system of equations A u' = B u
            ! call g%matmul_halo(rhs=g%rhs(1, 1:ndr), &
            !                    u=u, &
            !                    u_halo_m=u(:, g%size - ndr/2 + 1:g%size), &
            !                    u_halo_p=u(:, 1:ndr/2), &
            !                    f=result)
            call g%matmul_halo_thomas(rhs=g%rhs(1, 1:ndr), &
                                      u=u, &
                                      u_halo_m=u(:, g%size - ndr/2 + 1:g%size), &
                                      u_halo_p=u(:, 1:ndr/2), &
                                      f=result, &
                                      L=lu(:, 1:ndl/2))

            ! Solve for u' in system of equations A u' = B u
            ! call g%thomasL(lu(:, 1:ndl/2), result)
            call g%thomasU(lu(:, ndl/2 + 1:ndl), result)
            select case (g%nb_diag(1))
            case (3)
                call ThomasCirculant_3_Reduce(lu(:, 1:ndl/2), &
                                              lu(:, ndl/2 + 1:ndl), &
                                              lu(:, ndl + 1), &
                                              result, wrk2d)
            end select

        else    ! biased
            ! Calculate RHS in system of equations A u' = B u
            if (g%need_1der) then           ! add Jacobian correction A_2 dx2 du
                ip = g%nb_diag(2)           ! so far, only tridiagonal systems
                ! call g%matmuldevel_add(rhs=g%rhs(:, 1:ndr), &
                !                        rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                        rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                        u=u, &
                !                        rhs_add=g%rhs(:, ip + 1:ip + 3), &
                !                        u_add=du, &
                !                        f=result)
                call g%matmuldevel_add_thomas(rhs=g%rhs(:, 1:ndr), &
                                              rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                                              rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                              u=u, &
                                              rhs_add=g%rhs(:, ip + 1:ip + 3), &
                                              u_add=du, &
                                              f=result, &
                                              L=lu(:, 1:ndl/2))
            else
                ! call g%matmuldevel(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                    u=u, &
                !                    f=result)
                call g%matmuldevel_thomas(rhs=g%rhs(:, 1:ndr), &
                                          rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                                          rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                          u=u, &
                                          f=result, &
                                          L=lu(:, 1:ndl/2))
            end if

            ! Solve for u' in system of equations A u' = B u
            ! call g%thomasL(lu(:, 1:ndl/2), result)
            call g%thomasU(lu(:, ndl/2 + 1:ndl), result)

        end if

        return
    end subroutine FDM_Der2_Solve

end module FDM_Derivative

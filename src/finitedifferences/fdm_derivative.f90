#include "tlab_error.h"

module FDM_Derivative
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: lfile, efile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_PERIODIC, BCS_NONE
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas
    use Thomas_Circulant
    ! use MatMul
    use MatMul_Thomas
    ! use MatMul_Halo
    use MatMul_Halo_Thomas
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
        real(wp), allocatable :: lhs(:, :)                  ! memory space for LHS
        real(wp), allocatable :: rhs(:, :)                  ! memory space for RHS
        real(wp), allocatable :: mwn(:)                     ! memory space for modified wavenumbers
        !
        real(wp), allocatable :: lu(:, :)                   ! memory space for LU decomposition
        real(wp), allocatable :: rhs_b1(:, :), rhs_t1(:, :) ! memory space Neumann boundary conditions
        real(wp), allocatable :: rhs_d1(:, :)               ! memory space for 1. order derivative to 2. order one in Jacobian formulation

        ! procedure(matmul_halo_ice), pointer, nopass :: matmul_halo => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul_halo_thomas => null()
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        ! procedure(matmul_add_ice), pointer, nopass :: matmul_add => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul_thomas => null()
        procedure(matmul_add_thomas_ice), pointer, nopass :: matmul_add_thomas => null()
        ! procedure(thomas_ice), pointer, nopass :: thomasL => null()
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
    abstract interface
        subroutine thomas_ice(A, f)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: A(:, :)
            real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        end subroutine
    end interface

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

    ! abstract interface
    !     subroutine matmul_ice(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
    !         use TLab_Constants, only: wp
    !         real(wp), intent(in) :: rhs(:, :)
    !         real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
    !         real(wp), intent(in) :: u(:, :)
    !         real(wp), intent(out) :: f(:, :)
    !         real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
    !     end subroutine
    ! end interface

    ! abstract interface
    !     subroutine matmul_add_ice(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, bcs_b, bcs_t)
    !         use TLab_Constants, only: wp
    !         real(wp), intent(in) :: rhs(:, :)
    !         real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
    !         real(wp), intent(in) :: u(:, :)
    !         real(wp), intent(out) :: f(:, :)
    !         real(wp), intent(in) :: rhs_add(:, :)
    !         real(wp), intent(in) :: u_add(:, :)
    !         real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
    !     end subroutine
    ! end interface

    abstract interface
        subroutine matmul_thomas_ice(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
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
        subroutine matmul_add_thomas_ice(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, L, bcs_b, bcs_t)
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
        integer ndl, ndr, idl, idr
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
        call Precon_Rhs(g%lhs, g%rhs, periodic=periodic)

        ! LU decomposition
        if (allocated(g%lu)) deallocate (g%lu)
        if (g%periodic) then
            allocate (g%lu(g%size, ndl + ndl/2))
        else
            allocate (g%lu(g%size, 5*size(bcs_cases)))
        end if
        g%lu(:, :) = 0.0_wp

        ! extending rhs at the boundaries to ndr+2 diagonals
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

                call FDM_Der1_Neumann_Reduce(g%lhs(:, 1:ndl), &
                                             g%rhs(:, 1:ndr), &
                                             bcs_cases(ib), &
                                             g%lu(:, ip + 1:ip + ndl), &
                                             g%rhs_b1, &
                                             g%rhs_t1)

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
            ! g%thomasL => Thomas3_SolveL
            g%thomasU => Thomas3_SolveU
        case (5)
            ! g%thomasL => Thomas5_SolveL
            g%thomasU => Thomas5_SolveU
        case (7)
            ! g%thomasL => Thomas7_SolveL
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
                    ! g%matmul => MatMul_3
                    if (ndl == 3) g%matmul_thomas => MatMul_3_ThomasL_3
                case (5)
                    ! g%matmul => MatMul_5
                    if (ndl == 3) g%matmul_thomas => MatMul_5_ThomasL_3
                    if (ndl == 5) g%matmul_thomas => MatMul_5_ThomasL_5
                end select
            else
                select case (ndr)
                case (3)
                    ! g%matmul => MatMul_3_antisym
                    if (ndl == 3) g%matmul_thomas => MatMul_3_antisym_ThomasL_3
                case (5)
                    ! g%matmul => MatMul_5_antisym
                    if (ndl == 3) g%matmul_thomas => MatMul_5_antisym_ThomasL_3
                    if (ndl == 5) g%matmul_thomas => MatMul_5_antisym_ThomasL_5
                case (7)
                    ! g%matmul => MatMul_7_antisym
                    if (ndl == 3) g%matmul_thomas => MatMul_7_antisym_ThomasL_3
                    if (ndl == 5) g%matmul_thomas => MatMul_7_antisym_ThomasL_5
                end select
            end if

        end if

        return
    end subroutine FDM_Der1_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_CreateSystem(x, dx, g, periodic)
        use FDM_Base, only: MultiplyByDiagonal
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(:)                   ! Jacobian
        type(fdm_derivative_dt), intent(inout) :: g
        logical, intent(in) :: periodic

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx

        ! ###################################################################
        g%size = size(x)                ! # grid points
        nx = g%size                     ! for code readability
        g%periodic = periodic           ! flag for periodic direction

        ! -------------------------------------------------------------------
        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(g%size, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_Jacobian(g%size, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(g%size, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C1N4_Direct(x, g%lhs, g%rhs)

        case (FDM_COM6_DIRECT)
            call FDM_C1N6_Direct(x, g%lhs, g%rhs)

        end select

        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA)
            call MultiplyByDiagonal(g%lhs, dx)    ! multiply by the Jacobian
        end select

        ! For code readability later in the code
        g%nb_diag = [size(g%lhs, 2), size(g%rhs, 2)]

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (allocated(g%mwn)) deallocate (g%mwn)
        allocate (g%mwn(nx))
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
                       /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(2.0_wp*wn(:)))

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
        real(wp), intent(inout) :: r_rhs_b(:, :), r_rhs_t(:, :)     ! new, reduced rhs, extended diagonals

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, nx
        integer(wi) idr_t, ndr_t, idr_b, ndr_b, nx_t
        real(wp), allocatable :: aux(:, :)
        real(wp) locRhs_b(7, 8), locRhs_t(7, 8)

        ! ###################################################################
        nx = size(lhs, 1)           ! # grid points

        ndl = size(lhs, 2)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = ndr/2 + 1             ! center diagonal in rhs

        ndr_b = size(r_rhs_b, 2)    ! they can have a different number of diagonals than rhs
        idr_b = ndr_b/2 + 1
        ndr_t = size(r_rhs_t, 2)
        idr_t = ndr_t/2 + 1
        nx_t = max(idl, idr + 1)

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
        aux(1:nx, 1:ndr) = rhs(1:nx, 1:ndr)         ! array changed in FDM_Bcs_Reduce

        locRhs_b = 0.0_wp
        locRhs_t = 0.0_wp
        call FDM_Bcs_Reduce(ibc, aux, lhs, &
                            locRhs_b(1:max(idl, idr + 1), 1:ndr_b), &
                            locRhs_t(1:max(idl, idr + 1), 1:ndr_t))

        ! reorganize data
        if (any([BCS_ND, BCS_NN] == ibc)) then
            r_rhs_b(:, :) = 0.0_wp
            r_rhs_b(1:idr + 1, idr_b - ndr/2:idr_b + ndr/2) = aux(1:idr + 1, 1:ndr)
            r_rhs_b(1, idr_b) = lhs(1, idl)         ! save a_11 for nonzero bc
            do ir = 1, idr - 1                      ! save -a^R_{21} for nonzero bc
                r_rhs_b(1 + ir, idr_b - ir) = -locRhs_b(1 + ir, idr_b - ir)
            end do

            r_lhs(2:idl + 1, 1:ndl) = locRhs_b(2:idl + 1, idr_b - ndl/2:idr_b + ndl/2)
            r_lhs(1, idl) = rhs(1, idr)

            ! moving extended stencil in first element of old array to natural position
            r_rhs_b(1, idr_b + ndr/2 + 1) = r_rhs_b(1, idr_b - ndr/2)
            r_rhs_b(1, idr_b - ndr/2) = 0.0_wp

        end if

        if (any([BCS_DN, BCS_NN] == ibc)) then
            r_rhs_t(:, :) = 0.0_wp
            r_rhs_t(nx_t - idr:nx_t, idr_t - ndr/2:idr_t + ndr/2) = aux(nx - idr:nx, 1:ndr)
            r_rhs_t(nx_t, idr_t) = lhs(nx, idl)
            do ir = 1, idr - 1              ! change sign in a^R_{21} for nonzero bc
                r_rhs_t(nx_t - ir, idr_t + ir) = -locRhs_t(nx_t - ir, idr_t + ir)
            end do

            r_lhs(nx - idl:nx - 1, 1:ndl) = locRhs_t(nx_t - idl:nx_t - 1, idr_t - ndl/2:idr_t + ndl/2)
            r_lhs(nx, idl) = rhs(nx, idr)

            ! moving extended stencil in last element of old array to natural position
            r_rhs_t(nx_t, idr_t - ndr/2 - 1) = r_rhs_t(nx_t, idr_t + ndr/2)
            r_rhs_t(nx_t, idr_t + ndr/2) = 0.0_wp

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
                ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                    u=u, &
                !                    f=result)
                call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
                                     rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                                     rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                     u=u, &
                                     f=result, &
                                     L=lu1(:, ip + 1:ip + ndl/2))
            case (BCS_ND)
                ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                    u=u, &
                !                    f=result, &
                !                    bcs_b=bcs_hb(:))
                call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
                                     rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                     rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                     u=u, &
                                     f=result, &
                                     L=lu1(:, ip + 1:ip + ndl/2), &
                                     bcs_b=bcs_hb(:))
            case (BCS_DN)
                ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                    rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    u=u, &
                !                    f=result, &
                !                    bcs_t=bcs_ht(:))
                call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
                                     rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                                     rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                     u=u, &
                                     f=result, &
                                     L=lu1(:, ip + 1:ip + ndl/2), &
                                     bcs_t=bcs_ht(:))
            case (BCS_NN)
                ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                !                    u=u, &
                !                    f=result, &
                !                    bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))
                call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
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
        use FDM_Base, only: MultiplyByDiagonal
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

        ! Preconditioning
        if (g%need_1der) then
            call NormalizeByDiagonal(g%rhs, &
                                     1, &                           ! use 1. upper diagonal in rhs
                                     g%lhs, &
                                     g%rhs_d1, &
                                     switchAtBoundary=.not. (periodic))
        else
            call NormalizeByDiagonal(g%rhs, &
                                     1, &                           ! use 1. upper diagonal in rhs
                                     g%lhs, &
                                     switchAtBoundary=.not. (periodic))
        end if

        if (allocated(g%lu)) deallocate (g%lu)
        if (g%periodic) then
            allocate (g%lu(g%size, ndl + 1))
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
            ! g%thomasL => Thomas3_SolveL
            g%thomasU => Thomas3_SolveU
        case (5)
            ! g%thomasL => Thomas5_SolveL
            g%thomasU => Thomas5_SolveU
        case (7)
            ! g%thomasL => Thomas7_SolveL
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
                    ! g%matmul => MatMul_5
                    if (ndl == 3) g%matmul_thomas => MatMul_5_ThomasL_3
                    if (ndl == 5) g%matmul_thomas => MatMul_5_ThomasL_5
                end select
            else
                select case (ndr)
                case (5)
                    ! g%matmul => MatMul_5_sym
                    if (ndl == 3) g%matmul_thomas => MatMul_5_sym_ThomasL_3
                    if (ndl == 5) g%matmul_thomas => MatMul_5_sym_ThomasL_5
                    ! g%matmul_add => MatMul_5_sym_add_3
                    if (ndl == 3) g%matmul_add_thomas => MatMul_5_sym_add_3_ThomasL_3
                    if (ndl == 5) g%matmul_add_thomas => MatMul_5_sym_add_3_ThomasL_5
                case (7)
                    ! g%matmul => MatMul_7_sym
                    if (ndl == 3) g%matmul_thomas => MatMul_7_sym_ThomasL_3
                    if (ndl == 5) g%matmul_thomas => MatMul_7_sym_ThomasL_5
                    ! g%matmul_add => MatMul_7_sym_add_3
                    if (ndl == 3) g%matmul_add_thomas => MatMul_7_sym_add_3_ThomasL_3
                    if (ndl == 5) g%matmul_add_thomas => MatMul_7_sym_add_3_ThomasL_5
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

        ! ###################################################################
        g%size = size(x)                ! # grid points
        nx = g%size                     ! for code readability
        g%periodic = periodic           ! flag for periodic direction

        ! -------------------------------------------------------------------
        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(g%size, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN)
            call FDM_C2N6_Jacobian(g%size, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(g%size, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C2N4_Direct(x, g%lhs, g%rhs)

        case (FDM_COM6_DIRECT)
            call FDM_C2N6_Direct(x, g%lhs, g%rhs)

        end select

        g%need_1der = .false.       ! Just is case you overwrite type and it was .true. before
        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER)
            if (.not. uniform) then
                g%need_1der = .true.
                if (allocated(g%rhs_d1)) deallocate (g%rhs_d1)  ! Contribution from 1. order derivative in nonuniform grids
                allocate (g%rhs_d1, mold=g%lhs)
                g%rhs_d1 = -g%lhs
                call MultiplyByDiagonal(g%rhs_d1, dx(:, 2))
            end if
            call MultiplyByDiagonal(g%lhs, dx(:, 1))            ! multiply by the Jacobians
            call MultiplyByDiagonal(g%lhs, dx(:, 1))
        end select

        ! For code readability later in the code
        g%nb_diag = [size(g%lhs, 2), size(g%rhs, 2)]

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (allocated(g%mwn)) deallocate (g%mwn)
        allocate (g%mwn(nx))
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
        integer(wi) ndl, ndr

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
                ! call g%matmul_add(rhs=g%rhs(:, 1:ndr), &
                !                        rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                        rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                        u=u, &
                !                        rhs_add=g%rhs(:, ip + 1:ip + 3), &
                !                        u_add=du, &
                !                        f=result)
                call g%matmul_add_thomas(rhs=g%rhs(:, 1:ndr), &
                                         rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                                         rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                                         u=u, &
                                         rhs_add=g%rhs_d1, &
                                         u_add=du, &
                                         f=result, &
                                         L=lu(:, 1:ndl/2))
            else
                ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
                !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
                !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
                !                    u=u, &
                !                    f=result)
                call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
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

#include "tlab_error.h"

!########################################################################
! Boundary-value problems and integrals based on the compact schemes.
! Should we move this to OPR_ODE in operators?
!########################################################################

module FDM_Integral
    use TLab_Constants, only: wp, wi, efile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas
    use MatMul
    use MatMul_Thomas
    use Preconditioning
    use FDM_Derivative, only: fdm_derivative_dt
    use FDM_Base
    implicit none
    private

    type, public :: fdm_integral_dt
        sequence
        integer mode_fdm                                    ! original finite-difference method; only informative
        real(wp) :: lambda                                  ! constant of the equation
        integer :: bc                                       ! type of boundary condition, [ BCS_MIN, BCS_MAX ]
        real(wp), allocatable :: lhs(:, :)                  ! Often overwritten to LU decomposition.
        real(wp), allocatable :: rhs(:, :)
        real(wp), allocatable :: rhs_b1(:, :), rhs_t1(:, :) ! boundary conditions
        procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul_thomas => null()
        procedure(thomas_ice), pointer, nopass :: thomasL => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
    end type fdm_integral_dt
    ! This type is used in elliptic operators for different eigenvalues. This can lead to fragmented memory.
    ! One could use pointers instead of allocatable for lhs and rhs, and point the pointers to the
    ! corresponding memory space.

    public FDM_Int1_Initialize                      ! Prepare to solve u' +\lambda u = f
    public FDM_Int1_CreateSystem
    public FDM_Int1_Solve

    public FDM_Int2_Initialize                      ! Prepare to solve (u')' - \lamba^2 u = f
    ! public FDM_Int2_CreateSystem
    public FDM_Int2_Solve

    ! -----------------------------------------------------------------------
    abstract interface
        subroutine thomas_ice(A, f)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: A(:, :)
            real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        end subroutine
    end interface

    abstract interface
        subroutine matmul_ice(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine
    end interface

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

contains
    !########################################################################
    !#
    !#     u'_i + \lamba u_i = f_i  N-1 eqns
    !#     u_1 or u_N given         1   eqn
    !#     Au' = Bu                 N   eqns
    !#
    !# starting from generic diagonal matrices A (lhs) and B (rhs).
    !#
    !# The system of N-1 eqns:
    !#
    !#                    (B + \lambda A)u = Af
    !#
    !# is established in this routine (see notes).
    !#
    !# System normalized s.t. RHS diagonals are O(1), LHS diagonals O(h^2)
    !# System normalized s.t. 1. upper-diagonal in B is 1 (except at boundaries)
    !#
    !########################################################################
    subroutine FDM_Int1_Initialize(g, lambda, ibc, fdmi)
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, ndl, ndr

        !########################################################################
        call FDM_Int1_CreateSystem(g, lambda, ibc, fdmi)

        ! LU decomposition
        nx = size(fdmi%lhs, 1)              ! # of grid points
        ndl = size(fdmi%lhs, 2)             ! # of diagonals
        ndr = size(fdmi%rhs, 2)             ! # of diagonals

        call Thomas_FactorLU_InPlace(fdmi%lhs(2:nx - 1, 1:ndl/2), &
                                     fdmi%lhs(2:nx - 1, ndl/2 + 1:ndl))

        ! -------------------------------------------------------------------
        ! Procedure pointers to linear solvers
        select case (ndl)
        case (3)
            fdmi%thomasL => Thomas3_SolveL
            fdmi%thomasU => Thomas3_SolveU
        case (5)
            fdmi%thomasL => Thomas5_SolveL
            fdmi%thomasU => Thomas5_SolveU
        case (7)
            fdmi%thomasL => Thomas7_SolveL
            fdmi%thomasU => Thomas7_SolveU
        end select

        ! -------------------------------------------------------------------
        ! Procedure pointers to matrix multiplication to calculate the right-hand side
        select case (ndr)
        case (3)
            ! fdmi%matmul => MatMul_3
            if (ndl == 3) fdmi%matmul_thomas => MatMul_3_ThomasL_3
            if (ndl == 5) fdmi%matmul_thomas => MatMul_3_ThomasL_5
        case (5)
            ! fdmi%matmul => MatMul_5
            if (ndl == 3) fdmi%matmul_thomas => MatMul_5_ThomasL_3
            if (ndl == 5) fdmi%matmul_thomas => MatMul_5_ThomasL_5
            if (ndl == 7) fdmi%matmul_thomas => MatMul_5_ThomasL_7
        end select

        return
    end subroutine FDM_Int1_Initialize

    !########################################################################
    !########################################################################
    subroutine FDM_Int1_CreateSystem(g, lambda, ibc, fdmi)
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, ic, nx, nmin, nmax
        integer(wi) idr_t, ndr_t, idr_b, ndr_b, nx_t
        ! real(wp) locRhs_b(5, 0:7), locRhs_t(0:4, 8)
        real(wp) locRhs_b(7, 8), locRhs_t(7, 8)

        ! ###################################################################
        ndl = g%nb_diag(1)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = g%size                 ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%mode_fdm = g%mode_fdm
        fdmi%lambda = lambda
        fdmi%bc = ibc

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A), independent of lambda
        if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        allocate (fdmi%rhs(nx, ndl))
        fdmi%rhs(1:nx, 1:ndl) = g%lhs(1:nx, 1:ndl)

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C), dependent on lambda
        !                    | 0 a_12 |
        !   C = B + lamnda h | 0 A_22 |     for BCS_MIN
        !
        ! and correspondingly for BCS_MAX
        if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        allocate (fdmi%lhs(nx, ndr))
        fdmi%lhs(1:nx, 1:ndr) = g%rhs(1:nx, 1:ndr)

        select case (fdmi%bc)
        case (BCS_MIN)          ! first column of lambda h A is zero
            nmin = 2; nmax = nx
        case (BCS_MAX)          ! last column of lambda h A is zero
            nmin = 1; nmax = nx - 1
        end select

        fdmi%lhs(nmin:nmax, idr) = fdmi%lhs(nmin:nmax, idr) + lambda*g%lhs(nmin:nmax, idl)  ! center diagonal
        do ic = 1, idl - 1                                                                  ! off-center diagonals
            fdmi%lhs(nmin + ic:nx, idr - ic) = fdmi%lhs(nmin + ic:nx, idr - ic) + lambda*g%lhs(nmin + ic:nx, idl - ic)
            fdmi%lhs(1:nmax - ic, idr + ic) = fdmi%lhs(1:nmax - ic, idr + ic) + lambda*g%lhs(1:nmax - ic, idl + ic)
        end do

        ! -------------------------------------------------------------------
        ! Reduction; extending 2 diagonals the rhs at the boundaries
        if (allocated(fdmi%rhs_b1)) deallocate (fdmi%rhs_b1)
        if (allocated(fdmi%rhs_t1)) deallocate (fdmi%rhs_t1)
        ! allocate (fdmi%rhs_b1(max(idr, idl + 1), 1:max(ndl, ndr) + 2))  ! should be ndl+2
        ! allocate (fdmi%rhs_t1(max(idr, idl + 1), 1:max(ndl, ndr) + 2))
        ! allocate (fdmi%rhs_b1(max(idr, idl + 1), 1:ndl + 2))
        ! allocate (fdmi%rhs_t1(max(idr, idl + 1), 1:ndl + 2))
        allocate (fdmi%rhs_b1(max(idr, idl) + 1, 1:max(ndl, ndr) + 2))
        allocate (fdmi%rhs_t1(max(idr, idl) + 1, 1:max(ndl, ndr) + 2))
        fdmi%rhs_b1(:, :) = 0.0_wp
        fdmi%rhs_t1(:, :) = 0.0_wp

        ndr_b = size(fdmi%rhs_b1, 2)    ! they can have a different number of diagonals than rhs
        idr_b = ndr_b/2 + 1
        ndr_t = size(fdmi%rhs_t1, 2)
        idr_t = ndr_t/2 + 1

        nx_t = size(fdmi%rhs_t1, 1)

        locRhs_b = 0.0_wp
        locRhs_t = 0.0_wp
        ! call FDM_Bcs_Reduce_Old(fdmi%bc, fdmi%rhs, fdmi%lhs, locRhs_b, locRhs_t)
        call FDM_Bcs_Reduce(fdmi%bc, fdmi%rhs, fdmi%lhs, locRhs_b(1:nx_t, 1:ndr_b), locRhs_t(1:nx_t, 1:ndr_t))

        select case (fdmi%bc)
        case (BCS_MIN)
            ! fdmi%lhs(1:idr, 1:ndr) = locRhs_b(1:idr, 1:ndr)
            fdmi%lhs(1:nx_t, 1:ndr) = locRhs_b(1:nx_t, idr_b - ndr/2:idr_b + ndr/2)

            ! fdmi%rhs_b1(1:idl + 1, 2:ndl + 1) = fdmi%rhs(1:idl + 1, 1:ndl)
            fdmi%rhs_b1(1:nx_t, idr_b - ndl/2:idr_b + ndl/2) = fdmi%rhs(1:nx_t, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                ! fdmi%rhs_b1(1 + ir, 1 + idl - ir) = -locRhs_b(1 + ir, idr - ir)
                fdmi%rhs_b1(1 + ir, idr_b - ir) = -locRhs_b(1 + ir, idr_b - ir)
            end do

            ! reducing system in the opposite end to account for the case of extended stencils
            ! call FDM_Bcs_Reduce_Old(BCS_MAX, fdmi%lhs, fdmi%rhs, rhs_t=fdmi%rhs_t1(:, 2:))
            call FDM_Bcs_Reduce(BCS_MAX, fdmi%lhs, fdmi%rhs, rhs_t=fdmi%rhs_t1)

        case (BCS_MAX)
            ! fdmi%lhs(nx - idr + 1:nx, 1:ndr) = locRhs_t(1:idr, 1:ndr)
            fdmi%lhs(nx - nx_t + 1:nx, 1:ndr) = locRhs_t(1:nx_t, idr_t - ndr/2:idr_t + ndr/2)

            ! fdmi%rhs_t1(1:idl + 1, 2:ndl + 1) = fdmi%rhs(nx - idl:nx, 1:ndl)
            fdmi%rhs_t1(1:nx_t, idr_t - ndl/2:idr_t + ndl/2) = fdmi%rhs(nx - nx_t + 1:nx, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                ! fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir) = -locRhs_t(idr - ir, idr + ir)
                fdmi%rhs_t1(nx_t - ir, idr_t + ir) = -locRhs_t(nx_t - ir, idr_t + ir)
            end do

            ! reducing system in the opposite end to account for the case of extended stencils
            ! call FDM_Bcs_Reduce_Old(BCS_MIN, fdmi%lhs, fdmi%rhs, rhs_b=fdmi%rhs_b1)
            call FDM_Bcs_Reduce(BCS_MIN, fdmi%lhs, fdmi%rhs, rhs_b=fdmi%rhs_b1)

        end select

        ! moving extended stencil in last element of old array to natural position
        ir = 1
        ! fdmi%rhs_b1(ir, ndl + 2) = fdmi%rhs_b1(ir, 2); fdmi%rhs_b1(ir, 2) = 0.0_wp
        fdmi%rhs_b1(ir, idr_b + ndl/2 + 1) = fdmi%rhs_b1(ir, idr_b - ndl/2)
        fdmi%rhs_b1(ir, idr_b - ndl/2) = 0.0_wp
        ir = nx_t
        ! fdmi%rhs_t1(ir, 1) = fdmi%rhs_t1(ir, ndl + 1); fdmi%rhs_t1(ir, ndl + 1) = 0.0_wp
        fdmi%rhs_t1(ir, idr_t - ndl/2 - 1) = fdmi%rhs_t1(ir, idr_t + ndl/2)
        fdmi%rhs_t1(ir, idr_t + ndl/2) = 0.0_wp

        ! -------------------------------------------------------------------
        ! preconditioning
        !
        ! so far, based on the rhs
        ! use of lhs brings lambda to the rhs, in conflict with opr_elliptic.
        call Precon_Rhs(fdmi%lhs, fdmi%rhs, fdmi%rhs_b1, fdmi%rhs_t1)

        return
    end subroutine FDM_Int1_CreateSystem

    !########################################################################
    !########################################################################
    ! Allow to pass separate rhs because this part does not depend on lambda
    subroutine FDM_Int1_Solve(nlines, fdmi, rhsi, f, result, wrk2d, du_boundary)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi
        real(wp), intent(in) :: rhsi(:, :)
        real(wp), intent(in) :: f(nlines, size(fdmi%lhs, 1))
        real(wp), intent(inout) :: result(nlines, size(fdmi%lhs, 1))   ! contains bcs
        real(wp), intent(inout) :: wrk2d(nlines, 2)
        real(wp), intent(out), optional :: du_boundary(nlines)

        ! -------------------------------------------------------------------
        integer(wi) :: nx
        integer(wi) :: idl, ndl, idr, ndr, ic

        ! ###################################################################
        nx = size(fdmi%lhs, 1)

        ndl = size(fdmi%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(rhsi, 2)
        idr = ndr/2 + 1

#define bcs_hb(i) wrk2d(i,1)
#define bcs_ht(i) wrk2d(i,2)

        select case (fdmi%bc)
        case (BCS_MIN)
            bcs_hb(:) = result(:, 1)
            bcs_ht(:) = f(:, nx)
        case (BCS_MAX)
            bcs_hb(:) = f(:, 1)
            bcs_ht(:) = result(:, nx)
        end select

        ! call fdmi%matmul(rhs=rhsi(:, 1:ndr), &
        !                       rhs_b=fdmi%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
        !                       rhs_t=fdmi%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
        !                       u=f, &
        !                       f=result, &
        !                       bcs_b=bcs_hb(:), &
        !                       bcs_t=bcs_ht(:))
        call fdmi%matmul_thomas(rhs=rhsi(:, 1:ndr), &
                                ! rhs_b=fdmi%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                ! rhs_t=fdmi%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                rhs_b=fdmi%rhs_b1, &
                                rhs_t=fdmi%rhs_t1, &
                                u=f, &
                                f=result, &
                                L=fdmi%lhs(:, 1:ndl/2), &
                                bcs_b=bcs_hb(:), &
                                bcs_t=bcs_ht(:))

        ! call fdmi%thomasL(fdmi%lhs(2:nx - 1, 1:ndl/2), result(:, 2:nx - 1))
        call fdmi%thomasU(fdmi%lhs(2:nx - 1, ndl/2 + 1:ndl), result(:, 2:nx - 1))

        select case (fdmi%bc)
        case (BCS_MIN)
            result(:, nx) = bcs_ht(:)
            do ic = 1, idl - 1
                result(:, nx) = result(:, nx) + fdmi%lhs(nx, idl - ic)*result(:, nx - ic)
            end do
            result(:, nx) = result(:, nx) + fdmi%lhs(nx, ndl)*result(:, nx - ic)

            result(:, nx) = result(:, nx)/fdmi%lhs(nx, idl)

            if (present(du_boundary)) then      ! calculate u'1
                du_boundary(:) = result(:, 1)*fdmi%lhs(1, idl)
                do ic = 1, idl - 1
                    du_boundary(:) = du_boundary(:) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
                end do
                ic = idl                        ! longer stencil at the boundary
                du_boundary(:) = du_boundary(:) + fdmi%lhs(1, 1)*result(:, 1 + ic)

                do ic = 1, idr - 1
                    du_boundary(:) = du_boundary(:) + rhsi(1, idr + ic)*f(:, 1 + ic)
                end do

                du_boundary(:) = du_boundary(:)/rhsi(1, idr)

            end if

        case (BCS_MAX)
            result(:, 1) = bcs_hb(:)
            do ic = 1, idl - 1
                result(:, 1) = result(:, 1) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
            end do
            result(:, 1) = result(:, 1) + fdmi%lhs(1, 1)*result(:, 1 + ic)

            result(:, 1) = result(:, 1)/fdmi%lhs(1, idl)

            if (present(du_boundary)) then      ! calculate u'n
                du_boundary(:) = result(:, nx)*fdmi%lhs(nx, idl)
                do ic = 1, idl - 1
                    du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, idl - ic)*result(:, nx - ic)
                end do
                ic = idl                        ! longer stencil at the boundary
                du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, ndl)*result(:, nx - ic)

                do ic = 1, idr - 1
                    du_boundary(:) = du_boundary(:) + rhsi(nx, idr - ic)*f(:, nx - ic)
                end do

                du_boundary(:) = du_boundary(:)/rhsi(nx, idr)

            end if

        end select

#undef bcs_hb
#undef bcs_ht

        return
    end subroutine FDM_Int1_Solve

    !########################################################################
    !#
    !#     u''_i - \lamba^2 u_i = f_i  N-2 eqns
    !#     u_1 and u_N given           2   eqns
    !#     Au'' = Bu                   N   eqns
    !#
    !# starting from the matrices A (lhs, tridiagonal) and B (rhs, pentadiagonal, with unitary central diagonal).
    !#
    !# The system of N-2 eqns:
    !#
    !#                    (B - \lambda^2 A)u = Af = g
    !#
    !# is established in this routine, giving diagonals a-e and g (see notes).
    !#
    !# System normalized s.t. RHS diagonals are O(1), LHS diagonals O(h^2)
    !# System normalized s.t. 1. upper-diagonal in B is 1 (except at boundaries)
    !#
    !########################################################################
    subroutine FDM_Int2_Initialize(x, g, lambda2, ibc, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda2                 ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, ndl, ndr

        !########################################################################
        call FDM_Int2_CreateSystem(x, g, lambda2, ibc, fdmi)

        ! LU decomposition
        nx = size(fdmi%lhs, 1)              ! # of grid points
        ndl = size(fdmi%lhs, 2)             ! # of diagonals
        ndr = size(fdmi%rhs, 2)             ! # of diagonals

        call Thomas_FactorLU_InPlace(fdmi%lhs(2:nx - 1, 1:ndl/2), &
                                     fdmi%lhs(2:nx - 1, ndl/2 + 1:ndl))

        ! -------------------------------------------------------------------
        ! Procedure pointers to linear solvers
        select case (ndl)
        case (3)
            fdmi%thomasL => Thomas3_SolveL
            fdmi%thomasU => Thomas3_SolveU
        case (5)
            fdmi%thomasL => Thomas5_SolveL
            fdmi%thomasU => Thomas5_SolveU
        case (7)
            fdmi%thomasL => Thomas7_SolveL
            fdmi%thomasU => Thomas7_SolveU
        end select

        ! -------------------------------------------------------------------
        ! Procedure pointers to matrix multiplication to calculate the right-hand side
        select case (ndr)
        case (3)
            ! fdmi%matmul => MatMul_3
            if (ndl == 3) fdmi%matmul_thomas => MatMul_3_ThomasL_3
            if (ndl == 5) fdmi%matmul_thomas => MatMul_3_ThomasL_5
        case (5)
            ! fdmi%matmul => MatMul_5
            if (ndl == 3) fdmi%matmul_thomas => MatMul_5_ThomasL_3
            if (ndl == 5) fdmi%matmul_thomas => MatMul_5_ThomasL_5
            if (ndl == 7) fdmi%matmul_thomas => MatMul_5_ThomasL_7
        end select

        return
    end subroutine FDM_Int2_Initialize

    !########################################################################
    !########################################################################
    ! Follows FDM_Int1_CreateSystem very much (see notes)
    subroutine FDM_Int2_CreateSystem(x, g, lambda2, ibc, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda2                 ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, nx, ic
        real(wp) locRhs_b(5, 0:7), locRhs_t(0:4, 8)
        real(wp) coef(5)

        ! ###################################################################
        ndl = g%nb_diag(1)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = g%size                 ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%mode_fdm = g%mode_fdm
        fdmi%lambda = lambda2
        fdmi%bc = ibc

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A22R), independent of lambda
        if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        allocate (fdmi%rhs(nx, ndl))
        fdmi%rhs(1:nx, 1:ndl) = g%lhs(1:nx, 1:ndl)

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C22R); remember rhs center diagonal is not saved because it was 1
        if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        allocate (fdmi%lhs(nx, ndr))
        fdmi%lhs(1:nx, 1:ndr) = g%rhs(1:nx, 1:ndr)

        fdmi%lhs(:, idr) = fdmi%lhs(:, idr) - lambda2*g%lhs(:, idl)             ! center diagonal
        do ic = 1, idl - 1                                                      ! off-diagonals
            fdmi%lhs(1 + ic:nx, idr - ic) = fdmi%lhs(1 + ic:nx, idr - ic) - lambda2*g%lhs(1 + ic:nx, idl - ic)
            fdmi%lhs(1:nx - ic, idr + ic) = fdmi%lhs(1:nx - ic, idr + ic) - lambda2*g%lhs(1:nx - ic, idl + ic)
        end do

        ! -------------------------------------------------------------------
        ! Reduction; extending 2 diagonals the rhs at the boundaries
        if (allocated(fdmi%rhs_b1)) deallocate (fdmi%rhs_b1)
        if (allocated(fdmi%rhs_t1)) deallocate (fdmi%rhs_t1)
        allocate (fdmi%rhs_b1(max(idr, idl + 1), 1:max(ndl, ndr) + 2))
        allocate (fdmi%rhs_t1(max(idr, idl + 1), 1:max(ndl, ndr) + 2))
        fdmi%rhs_b1(:, :) = 0.0_wp
        fdmi%rhs_t1(:, :) = 0.0_wp

        locRhs_b = 0.0_wp
        locRhs_t = 0.0_wp
        call FDM_Bcs_Reduce_Old(BCS_BOTH, fdmi%rhs, g%rhs(:, 1:ndr), locRhs_b, locRhs_t)

        ! bcs min
        fdmi%rhs_b1(1:idl + 1, 2:ndl + 1) = fdmi%rhs(1:idl + 1, 1:ndl)
        do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
            fdmi%rhs_b1(1 + ir, 1 + idl - ir) = -locRhs_b(1 + ir, idr - ir)
        end do

        fdmi%lhs(2:idr, 1:ndr) = locRhs_b(2:idr, 1:ndr)
        do ir = 1, idr - 1
            fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) - lambda2*fdmi%rhs_b1(1 + ir, 2:ndl + 1)
        end do

        ! bcs min
        fdmi%rhs_t1(1:idl + 1, 2:ndl + 1) = fdmi%rhs(nx - idl:nx, 1:ndl)
        do ir = 1, idr - 1              ! change sign in b^R_{2n} for nonzero bc
            fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir) = -locRhs_t(idr - ir, idr + ir)
        end do

        fdmi%lhs(nx - idr + 1:nx - 1, 1:ndr) = locRhs_t(1:idr - 1, 1:ndr)
        do ir = 1, idr - 1
            fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) - lambda2*fdmi%rhs_t1(1 + idl - ir, 2:ndl + 1)
        end do

        ! -------------------------------------------------------------------
        ! Corrections to the BCS_DD to account for Neumann using third-order fdm for derivative at the boundary
        if (any([BCS_ND, BCS_NN] == fdmi%bc)) then
            ! Coefficients in FDM p'_1 = b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, 1)                  ! second-order
            ! coef(1:4) = coef_e1n3_biased(x, 1)                  ! third-order
            coef(1:5) = coef_c1n4_biased(x, 1)                  ! fourth-order

            ! Solve for p_1 (see notes)
            fdmi%lhs(1, :) = 0.0_wp
            fdmi%lhs(1, 1:3) = -coef(2:4)/coef(1)               ! vector d_2
            fdmi%rhs_b1(1, :) = 0.0_wp
            fdmi%rhs_b1(1, 1 + idl) = 1.0_wp/coef(1)                 ! coefficient d_1
            fdmi%rhs_b1(1, 1 + idl + 1) = -coef(5)/coef(1)           ! vector e_2, only 1 component

            ! Construct vector d + lambda^2h^2 e, e only 1 component
            fdmi%lhs(1, 1) = fdmi%lhs(1, 1) + lambda2*fdmi%rhs_b1(1, 1 + idl + 1)

            do ir = 1, idr - 1
                fdmi%lhs(1 + ir, idr - ir + 1:idr - ir + 1 + 2) = fdmi%lhs(1 + ir, idr - ir + 1:idr - ir + 1 + 2) &
                                                                  - fdmi%rhs_b1(1 + ir, 1 + idl - ir)*fdmi%lhs(1, 1:3)       ! in reduced C matrix

                fdmi%rhs_b1(1 + ir, 1 + idl - ir + 1) = fdmi%rhs_b1(1 + ir, 1 + idl - ir + 1) &
                                                        + fdmi%rhs_b1(1 + ir, 1 + idl - ir)*fdmi%rhs_b1(1, 1 + idl + 1)                ! in reduced A matrix

                fdmi%rhs_b1(1 + ir, 1 + idl - ir) = fdmi%rhs_b1(1 + ir, 1 + idl - ir)*fdmi%rhs_b1(1, 1 + idl)                          ! d_1 b^R_{21}
            end do

        end if

        if (any([BCS_DN, BCS_NN] == fdmi%bc)) then
            ! Coefficients in FDM p'_n = b_1 p_n + b_2 p_{n-1} + b_3 p_{n-2} +...
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, nx, backwards=.true.)
            ! coef(1:4) = coef_e1n3_biased(x, nx, backwards=.true.)
            coef(1:5) = coef_c1n4_biased(x, nx, backwards=.true.)

            ! Solve for p_n (see notes)
            fdmi%lhs(nx, :) = 0.0_wp
            fdmi%lhs(nx, ndr - 2:ndr) = -coef([4, 3, 2])/coef(1)  ! vector d_n-1
            fdmi%rhs_t1(1 + idl, :) = 0.0_wp
            fdmi%rhs_t1(1 + idl, 1 + idl) = 1.0_wp/coef(1)                 ! coefficient d_n
            fdmi%rhs_t1(1 + idl, 1 + idl - 1) = -coef(5)/coef(1)           ! vector e_n-1, only 1 component

            ! Construct vector d + lambda^2h^2 e, e only 1 component
            fdmi%lhs(nx, ndr) = fdmi%lhs(nx, ndr) + lambda2*fdmi%rhs_t1(1 + idl, 1 + idl - 1)

            do ir = 1, idr - 1
                fdmi%lhs(nx - ir, ir - 1 + 1:ir - 1 + 3) = fdmi%lhs(nx - ir, ir - 1 + 1:ir - 1 + 3) &
                                                           - fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir)*fdmi%lhs(nx, ndr - 2:ndr)              ! in reduced C matrix

                fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir - 1) = fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir - 1) &
                                                              + fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir)*fdmi%rhs_t1(1 + idl, 1 + idl - 1)      ! in reduced A matrix

                fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir) = fdmi%rhs_t1(1 + idl - ir, 1 + idl + ir)*fdmi%rhs_t1(1 + idl, 1 + idl)                ! d_n b^R_{2n}
            end do

        end if

        ! moving extended stencil in last element of old array to natural position
        ir = 1
        fdmi%rhs_b1(ir, ndl + 2) = fdmi%rhs_b1(ir, 2); fdmi%rhs_b1(ir, 2) = 0.0_wp
        ir = max(idr, idl + 1)
        fdmi%rhs_t1(ir, 1) = fdmi%rhs_t1(ir, ndl + 1); fdmi%rhs_t1(ir, ndl + 1) = 0.0_wp

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        ! First and last rows are not normalized
        call Precon_Rhs(fdmi%lhs(2:nx - 1, :), &
                        fdmi%rhs(2:nx - 1, :), &
                        fdmi%rhs_b1(2:max(idr, idl + 1), :), &
                        fdmi%rhs_t1(1:max(idr, idl + 1) - 1, :))

        return
    contains
        !########################################################################
        ! 1. derivatie of interpolation polynomial between equations (15) and (16)
        !    p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
        !
        ! Notation in Shukla and Zhong (2005), JCP, 204, 404–429 for the interpolation:
        !
        !       +                    I_n: set of points where the function and derivatives are given
        !   +---+---+---+---...
        !   +       +   +            I_m: set of points where only the function is given.
        !########################################################################
        function coef_c1n4_biased(x, i, backwards) result(coef)
            real(wp), intent(in) :: x(:)
            integer(wi), intent(in) :: i
            logical, intent(in), optional :: backwards
            real(wp) coef(5)

            real(wp) a2, b1, b2, b3, b4
            real(wp) dx1, dx3, dx4
            real(wp) D
            integer(wi) set_m(3), i1, i2, i3, i4

            i1 = i
            if (present(backwards)) then
                ! same as fowards, but changing the signs of the increments w.r.t. i
                ! To understand it, e.g., define a new variable k = -j, where k is the
                ! discrete variable moving around i
                i2 = i - 1
                i3 = i - 2
                i4 = i - 3
            else
                i2 = i + 1
                i3 = i + 2
                i4 = i + 3
            end if
            dx1 = x(i2) - x(i1)
            dx3 = x(i2) - x(i3)
            dx4 = x(i2) - x(i4)
            set_m = [i1, i3, i4]

            ! -------------------------------------------------------------------
            a2 = 0.5_wp*(Pi(x, i1, set_m) - dx1*Pi_p(x, i1, set_m))/Pi_p(x, i2, set_m)

            b2 = Pi_p(x, i1, set_m)*(2.0_wp*Pi_p(x, i2, set_m) + dx1*Pi_pp_3(x, i2, set_m)) &
                 - Pi(x, i1, set_m)*Pi_pp_3(x, i2, set_m)
            b2 = 0.5_wp*b2/Pi(x, i2, set_m)/Pi_p(x, i2, set_m)

            ! -------------------------------------------------------------------
            D = Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m)
            b1 = Lag(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + 2*dx1*Lag_p(x, i2, i1, set_m)) &
                 - dx1*Lag_p(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m))
            b1 = -b1/dx1/D

            D = Lag(x, i2, i3, set_m) + dx3*Lag_p(x, i2, i3, set_m)
            b3 = Lag(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + 2*dx1*Lag_p(x, i2, i3, set_m)) &
                 - dx1*Lag_p(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + dx1*Lag_p(x, i2, i3, set_m))
            b3 = -b3/dx3/D

            D = Lag(x, i2, i4, set_m) + dx4*Lag_p(x, i2, i4, set_m)
            b4 = Lag(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + 2*dx1*Lag_p(x, i2, i4, set_m)) &
                 - dx1*Lag_p(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + dx1*Lag_p(x, i2, i4, set_m))
            b4 = -b4/dx4/D

            coef = [b1, b2, b3, b4, a2]

            ! if uniform, we should have ( -29/6 54/6 -27/6 2/6 )/h and 3h
            ! print*, [b1, b2, b3, b4]*(x(2)-x(1))
            ! print*, a2/(x(2)-x(1))

            return
        end function

    end subroutine FDM_Int2_CreateSystem

    !########################################################################
    !########################################################################
    ! Allow to pass separate rhs because this part does not depend on lambda
    subroutine FDM_Int2_Solve(nlines, fdmi, rhsi, f, result, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi
        real(wp), intent(in) :: rhsi(:, :)
        real(wp), intent(in) :: f(nlines, size(fdmi%lhs, 1))
        real(wp), intent(inout) :: result(nlines, size(fdmi%lhs, 1))   ! contains bcs
        real(wp), intent(inout) :: wrk2d(nlines, 2)

        ! -------------------------------------------------------------------
        integer(wi) :: nx
        integer(wi) :: ndl, ndr, idl, idr, ic

        ! ###################################################################
        nx = size(fdmi%lhs, 1)

        ndl = size(fdmi%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(rhsi, 2)
        idr = ndr/2 + 1

#define bcs_hb(i) wrk2d(i,1)
#define bcs_ht(i) wrk2d(i,2)

        bcs_hb(:) = result(:, 1)
        bcs_ht(:) = result(:, nx)
        ! call fdmi%matmul(rhs=rhsi(:, 1:ndr), &
        !                       rhs_b=fdmi%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
        !                       rhs_t=fdmi%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
        !                       u=f, &
        !                       f=result, &
        !                       bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))
        call fdmi%matmul_thomas(rhs=rhsi(:, 1:ndr), &
                                rhs_b=fdmi%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
                                rhs_t=fdmi%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
                                u=f, &
                                f=result, &
                                L=fdmi%lhs(:, 1:ndl/2), &
                                bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))

        ! call fdmi%thomasL(fdmi%lhs(2:nx - 1, 1:ndl/2), result(:, 2:nx - 1))
        call fdmi%thomasU(fdmi%lhs(2:nx - 1, ndl/2 + 1:ndl), result(:, 2:nx - 1))

        !   Corrections to the BCS_DD to account for Neumann
        if (any([BCS_ND, BCS_NN] == fdmi%bc)) then
            result(:, 1) = bcs_hb(:)
            do ic = 1, ndl
                result(:, 1) = result(:, 1) + fdmi%lhs(1, ic)*result(:, 1 + ic)
            end do
        end if

        if (any([BCS_DN, BCS_NN] == fdmi%bc)) then
            result(:, nx) = bcs_ht(:)
            do ic = 1, ndl
                result(:, nx) = result(:, nx) + fdmi%lhs(nx, ndl - ic + 1)*result(:, nx - ic)
            end do
        end if

#undef bcs_hb
#undef bcs_ht

        return
    end subroutine FDM_Int2_Solve

end module FDM_Integral

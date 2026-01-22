#include "tlab_error.h"

!########################################################################
! Boundary-value problems and integrals based on the compact schemes.
! Should we move this to OPR_ODE in operators?
!########################################################################

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

module FDM_Integral_1
    use TLab_Constants, only: wp, wi, efile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas
    use MatMul
    use MatMul_Thomas
    use Preconditioning
    use FDM_Integral_Base
    use FDM_Base
    implicit none
    private

    public fdm_integral_dt                          ! so that I do not need to load FDM_Integral_Base
    public FDM_Int1_Initialize                      ! Prepare to solve u' +\lambda u = f
    public FDM_Int1_CreateSystem
    public FDM_Int1_Solve

contains
    !########################################################################
    !########################################################################
    subroutine FDM_Int1_Initialize(fdm_der, lambda, ibc, fdmi)
        use FDM_Derivative_1order, only: der_dt !, der1_biased
        class(der_dt), intent(in) :: fdm_der
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, ndl, ndr

        !########################################################################
        call FDM_Int1_CreateSystem(fdm_der, lambda, ibc, fdmi)

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
    subroutine FDM_Int1_CreateSystem(fdm_der, lambda, ibc, fdmi)
        use FDM_Derivative_1order, only: der_dt
        class(der_dt), intent(in) :: fdm_der            ! derivative plan to be inverted
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, ic, nx, nmin, nmax
        integer(wi) idr_t, ndr_t, idr_b, ndr_b, nx_t
        real(wp), allocatable :: locRhs_b(:, :), locRhs_t(:, :)

        ! ###################################################################
        ndl = size(fdm_der%lhs, 2)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = size(fdm_der%rhs, 2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = size(fdm_der%lhs, 1)                ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%mode_fdm = fdm_der%type
        fdmi%lambda = lambda
        fdmi%bc = ibc

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A), independent of lambda
        if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        allocate (fdmi%rhs, source=fdm_der%lhs)

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C), dependent on lambda
        !                    | 0 a_12 |
        !   C = B + lamnda h | 0 A_22 |     for BCS_MIN
        !
        ! and correspondingly for BCS_MAX
        if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        allocate (fdmi%lhs, source=fdm_der%rhs)

        select case (fdmi%bc)
        case (BCS_MIN)          ! first column of lambda h A is zero
            nmin = 2; nmax = nx
        case (BCS_MAX)          ! last column of lambda h A is zero
            nmin = 1; nmax = nx - 1
        end select

        fdmi%lhs(nmin:nmax, idr) = fdmi%lhs(nmin:nmax, idr) + lambda*fdm_der%lhs(nmin:nmax, idl)  ! center diagonal
        do ic = 1, idl - 1                                                                  ! off-center diagonals
            fdmi%lhs(nmin + ic:nx, idr - ic) = fdmi%lhs(nmin + ic:nx, idr - ic) + lambda*fdm_der%lhs(nmin + ic:nx, idl - ic)
            fdmi%lhs(1:nmax - ic, idr + ic) = fdmi%lhs(1:nmax - ic, idr + ic) + lambda*fdm_der%lhs(1:nmax - ic, idl + ic)
        end do

        ! -------------------------------------------------------------------
        ! Reduction; extending 2 diagonals the rhs at the boundaries
        if (allocated(fdmi%rhs_b1)) deallocate (fdmi%rhs_b1)
        if (allocated(fdmi%rhs_t1)) deallocate (fdmi%rhs_t1)
        allocate (fdmi%rhs_b1(max(idr, idl) + 1, 1:max(ndl, ndr) + 2), source=0.0_wp)
        allocate (fdmi%rhs_t1(max(idr, idl) + 1, 1:max(ndl, ndr) + 2), source=0.0_wp)

        ndr_b = size(fdmi%rhs_b1, 2)    ! they can have a different number of diagonals than rhs
        idr_b = ndr_b/2 + 1
        ndr_t = size(fdmi%rhs_t1, 2)
        idr_t = ndr_t/2 + 1

        nx_t = size(fdmi%rhs_t1, 1)

        if (allocated(locRhs_b)) deallocate (locRhs_b)
        allocate (locRhs_b(nx_t, ndr_b), source=0.0_wp)
        if (allocated(locRhs_t)) deallocate (locRhs_t)
        allocate (locRhs_t(nx_t, ndr_t), source=0.0_wp)
        call FDM_Bcs_Reduce(fdmi%bc, fdmi%rhs, fdmi%lhs, rhs_b=locRhs_b, rhs_t=locRhs_t)

        select case (fdmi%bc)
        case (BCS_MIN)
            fdmi%lhs(1:nx_t, 1:ndr) = locRhs_b(1:nx_t, idr_b - ndr/2:idr_b + ndr/2)

            fdmi%rhs_b1(1:nx_t, idr_b - ndl/2:idr_b + ndl/2) = fdmi%rhs(1:nx_t, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_b1(1 + ir, idr_b - ir) = -locRhs_b(1 + ir, idr_b - ir)
            end do

            ! reducing system in the opposite end to account for the case of extended stencils
            call FDM_Bcs_Reduce(BCS_MAX, fdmi%lhs, fdmi%rhs, rhs_t=fdmi%rhs_t1)

        case (BCS_MAX)
            fdmi%lhs(nx - nx_t + 1:nx, 1:ndr) = locRhs_t(1:nx_t, idr_t - ndr/2:idr_t + ndr/2)

            fdmi%rhs_t1(1:nx_t, idr_t - ndl/2:idr_t + ndl/2) = fdmi%rhs(nx - nx_t + 1:nx, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_t1(nx_t - ir, idr_t + ir) = -locRhs_t(nx_t - ir, idr_t + ir)
            end do

            call FDM_Bcs_Reduce(BCS_MIN, fdmi%lhs, fdmi%rhs, rhs_b=fdmi%rhs_b1)

        end select

        ! moving extended stencil in last element of old array to natural position
        ir = 1
        fdmi%rhs_b1(ir, idr_b + ndl/2 + 1) = fdmi%rhs_b1(ir, idr_b - ndl/2)
        fdmi%rhs_b1(ir, idr_b - ndl/2) = 0.0_wp
        ir = nx_t
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

end module FDM_Integral_1

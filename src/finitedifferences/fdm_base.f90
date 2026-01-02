#include "tlab_error.h"

!########################################################################
! Building blocks to construct FDMs
! Lagrange polynomials
! Calculation of RHS for different stencil lengths and bcs (periodic|biased)
!########################################################################
module FDM_Base
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public Pi                   ! Product function defined over interval given by idx(:), Pi(x-x_j) for all j in idx
    public Pi_p                 ! First-order derivative of Pi
    public Pi_pp_3              ! Second-order derivative when idx has only 3 points
    public Lag                  ! Lagrange polynomials on idx(:) around i
    public Lag_p                ! First-order derivative of Lag
    public Lag_pp_3             ! Second-order derivative when idx has only 3 points

    public coef_e1n2_biased     ! coefficients for the biased, 2. order approximation to 1. order derivative
    public coef_e1n3_biased     ! coefficients for the biased, 3. order approximation to 1. order derivative

    public FDM_Bcs_Reduce       ! System reduction at the boundary points

    public MultiplyByDiagonal   ! To multiply by Jacobians in non-uniform grids

contains
    !########################################################################
    !########################################################################
    function Pi(x, j, idx) result(f)    ! Product function on interval idx(:) evaluated at x_j
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            f = f*(x(j) - x(idx(k)))
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Pi_p(x, j, idx) result(f)
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        real(wp) dummy
        integer(wi) k, m

        f = 0.0_wp
        do k = 1, size(idx)
            dummy = 1.0_wp
            do m = 1, size(idx)
                if (m /= k) then
                    dummy = dummy*(x(j) - x(idx(m)))
                end if
            end do
            f = f + dummy
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Pi_pp_3(x, j, idx) result(f)       ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        f = 2.0_wp*(x(j) - x(idx(1)) + x(j) - x(idx(2)) + x(j) - x(idx(3)))

        return
    end function

!########################################################################
!########################################################################
    function Lag(x, j, i, idx) result(f)        ! Lagrange polynomials on idx(:) around i evaluated at x_j
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f*(x(j) - x(idx(k)))/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Lag_p(x, j, i, idx) result(f)        ! 1. derivative of Lagrange polynomials on idx(:) around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k, m
        real(wp) den, dummy

        den = 1.0_wp
        f = 0.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                dummy = 1.0_wp
                do m = 1, size(idx)
                    if (idx(m) /= i .and. m /= k) then
                        dummy = dummy*(x(j) - x(idx(m)))
                    end if
                end do
                f = f + dummy
                den = den*(x(i) - x(idx(k)))
            end if
        end do
        f = f/den

        return
    end function

    ! -------------------------------------------------------------------
    function Lag_pp_3(x, j, i, idx) result(f)    ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i, idx(:)
        real(wp) f

        integer(wi) k

        f = 2.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

!########################################################################
!########################################################################
    function coef_e1n3_biased(x, i, backwards) result(coef) ! first-order derivative, explicit, 3. order
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        logical, intent(in), optional :: backwards
        real(wp) coef(4)

        integer(wi) stencil(4), k

        if (present(backwards)) then
            stencil = [i, i - 1, i - 2, i - 3]
        else
            stencil = [i, i + 1, i + 2, i + 3]
        end if

        do k = 1, size(stencil)
            coef(k) = Lag_p(x, i, stencil(k), stencil(:))
        end do
        ! if uniform, [ -11/6 3 -3/2 1/3 ]/h

        return
    end function

    ! -------------------------------------------------------------------
    function coef_e1n2_biased(x, i, backwards) result(coef) ! first-order derivative, explicit, 2. order
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        logical, intent(in), optional :: backwards
        real(wp) coef(3)

        integer(wi) stencil(3), k

        if (present(backwards)) then
            stencil = [i, i - 1, i - 2]
        else
            stencil = [i, i + 1, i + 2]
        end if

        do k = 1, size(stencil)
            coef(k) = Lag_p(x, i, stencil(k), stencil(:))
        end do
        ! if uniform, [ -3/2 2 -1/2 ]/h

        return
    end function

! #######################################################################
! #######################################################################
    subroutine FDM_Bcs_Reduce(ibc, lhs, rhs, rhs_b, rhs_t)
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: lhs(:, :)
        real(wp), intent(in), optional :: rhs(:, :)
        ! real(wp), intent(out), optional :: rhs_b(:, 0:), rhs_t(0:, :)
        real(wp), intent(out), optional :: rhs_b(:, :), rhs_t(:, :)

        integer(wi) idl, ndl, idr, ndr, ir, ic, nx, nx_t
        integer(wi) idr_t, ndr_t, idr_b, ndr_b
        real(wp) dummy

        ! -------------------------------------------------------------------
        ndl = size(lhs, 2)
        idl = ndl/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = ndr/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points

        ! -------------------------------------------------------------------
        if (any([BCS_MIN, BCS_BOTH] == ibc)) then
            dummy = 1.0_wp/lhs(1, idl)      ! normalize by l11

            ! reduced array A^R_{22}
            lhs(1, 1:ndl) = -lhs(1, 1:ndl); lhs(1, idl) = -lhs(1, idl)
            do ir = 1, idl - 1              ! rows
                do ic = idl + 1, ndl        ! columns
                    lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lhs(1, ic)*dummy
                end do
                ic = ndl + 1                ! longer stencil at the boundary
                lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lhs(1, 1)*dummy
            end do

            ! reduced array B^R_{22}
            if (present(rhs_b)) then
                if (size(rhs_b, 1) < idl + 1 .or. size(rhs_b, 2) < ndr) then
                    call TLab_Write_ASCII(efile, __FILE__//'. rhs_b array is too small.')
                    call TLab_Stop(DNS_ERROR_UNDEVELOP)
                end if
                ndr_b = size(rhs_b, 2)          ! they can have a different number of diagonals than rhs
                idr_b = ndr_b/2 + 1

                nx_t = size(rhs_b, 1)

                rhs_b(1:nx_t, idr_b - ndr/2:idr_b + ndr/2) = rhs(1:nx_t, 1:ndr)

                do ir = 1, idl - 1              ! rows
                    do ic = 0, ndr/2            ! columns; ic = 0 corresponds to vector b^R_{21}
                        rhs_b(1 + ir, idr_b + ic - ir) = rhs_b(1 + ir, idr_b + ic - ir) - lhs(1 + ir, idl - ir)*rhs_b(1, idr_b + ic)*dummy
                    end do
                    ic = ndr/2 + 1                ! longer stencil at the boundary
                    rhs_b(1 + ir, idr_b + ic - ir) = rhs_b(1 + ir, idr_b + ic - ir) - lhs(1 + ir, idl - ir)*rhs_b(1, idr_b - ndr/2)*dummy
                end do
            end if

        end if

        if (any([BCS_MAX, BCS_BOTH] == ibc)) then
            dummy = 1.0_wp/lhs(nx, idl)     ! normalize by lnn

            ! reduced array A^R_{11}
            lhs(nx, 1:ndl) = -lhs(nx, 1:ndl); lhs(nx, idl) = -lhs(nx, idl)
            do ir = 1, idl - 1              ! rows
                ic = 0                      ! longer stencil at the boundary
                lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lhs(nx, ndl)*dummy
                do ic = 1, idl - 1          ! columns
                    lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lhs(nx, ic)*dummy
                end do
            end do

            ! reduced array B^R_{11}
            if (present(rhs_t)) then
                if (size(rhs_t, 1) < idl + 1 .or. size(rhs_t, 2) < ndr) then
                    call TLab_Write_ASCII(efile, __FILE__//'. rhs_t array is too small.')
                    call TLab_Stop(DNS_ERROR_UNDEVELOP)
                end if
                ndr_t = size(rhs_t, 2)
                idr_t = ndr_t/2 + 1
                nx_t = size(rhs_t, 1)

                rhs_t(1:nx_t, idr_t - ndr/2:idr_t + ndr/2) = rhs(nx - nx_t + 1:nx, 1:ndr)

                do ir = 1, idl - 1              ! rows
                    do ic = 0, ndr/2            ! ic = 0 corresponds to vector b^R_{1n}
                        rhs_t(nx_t - ir, idr_t - ic + ir) = rhs_t(nx_t - ir, idr_t - ic + ir) - lhs(nx - ir, idl + ir)*rhs_t(nx_t, idr_t - ic)*dummy
                    end do
                    ic = ndr/2 + 1                ! longer stencil at the boundary
                    rhs_t(nx_t - ir, idr_t - ic + ir) = rhs_t(nx_t - ir, idr_t - ic + ir) - lhs(nx - ir, idl + ir)*rhs_t(nx_t, idr_t + ndr/2)*dummy
                end do
            end if

        end if

        return
    end subroutine FDM_Bcs_Reduce

    ! #######################################################################
    ! Multiply A dx, where A is band diagonal and dx is diagonal
    subroutine MultiplyByDiagonal(A, dx)
        real(wp), intent(inout) :: A(:, :)          ! Diagonals in band matrix
        real(wp), intent(in) :: dx(:)               ! Diagonal in Jacobian matrix

        integer id, ic

        id = size(A, 2)/2 + 1

        A(:, id) = A(:, id)*dx(:)                    ! center diagonal
        do ic = 1, id - 1                               ! off-diagonals
            A(:, id - ic) = A(:, id - ic)*cshift(dx(:), -ic)
            A(:, id + ic) = A(:, id + ic)*cshift(dx(:), +ic)
        end do

        return
    end subroutine MultiplyByDiagonal

end module FDM_Base

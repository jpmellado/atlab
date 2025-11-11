module Preconditioning
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public Precon_Rhs

contains
    ! #######################################################################
    ! #######################################################################
    ! normalization such that 1. upper-diagonal in rhs is 1,
    ! at the bottom switch to the 1. lower diagonal
    subroutine Precon_Rhs(lhs, rhs, rhs_b, rhs_t, periodic)
        real(wp), intent(inout) :: lhs(:, :)
        real(wp), intent(inout) :: rhs(:, :)
        real(wp), intent(inout), optional :: rhs_b(:, :), rhs_t(:, :)
        logical, optional :: periodic

        ! -------------------------------------------------------------------
        integer(wi) nx, nmin, nmax, ir
        integer(wi) ic                          ! diagonal used to normalize
        integer ndl, idl, ndr, idr
        integer nx_b, ndr_b, idr_b, nx_t, ndr_t, idr_t

        real(wp) dummy

        logical periodic_loc

        ! ###################################################################
        if (present(periodic)) then
            periodic_loc = periodic
        else
            periodic_loc = .false.
        end if

        ndl = size(lhs, 2)
        idl = ndl/2 + 1
        ndr = size(rhs, 2)
        idr = ndr/2 + 1
        nx = size(lhs, 1)

        nmin = 1; nmax = nx

        ! -------------------------------------------------------------------
        ! top boundary
        if (present(rhs_b)) then
            nx_b = size(rhs_b, 1)
            ndr_b = size(rhs_b, 2)
            idr_b = ndr_b/2 + 1

            ic = 0                      ! central diagonal
            ! ic = 1                      ! 1. upper diagonal
            do ir = 1, nx_b
                ! dummy = 1.0_wp/rhs_b(ir, idr_b + ic)
                dummy = 1.0_wp/rhs(ir, idr + ic)
                ! print *, ir, dummy
                rhs_b(ir, :) = rhs_b(ir, :)*dummy
                lhs(ir, :) = lhs(ir, :)*dummy
                rhs(ir, :) = rhs(ir, :)*dummy
            end do

            nmin = nx_b + 1

        end if

        ! -------------------------------------------------------------------
        ! top boundary
        if (present(rhs_t)) then
            nx_t = size(rhs_t, 1)
            ndr_t = size(rhs_t, 2)
            idr_t = ndr_t/2 + 1

            ic = 0                      ! central diagonal
            ! ic = 1                      ! 1. upper diagonal
            ! ic = -1                     ! 1. lower diagonal
            do ir = 1, nx_t
                ! dummy = 1.0_wp/rhs_t(nx_t - ir + 1, idr_t + ic)
                dummy = 1.0_wp/rhs(nx - ir + 1, idr + ic)
                ! print *, 'up', ir, dummy
                rhs_t(nx_t - ir + 1, :) = rhs_t(nx_t - ir + 1, :)*dummy
                lhs(nx - ir + 1, :) = lhs(nx - ir + 1, :)*dummy
                rhs(nx - ir + 1, :) = rhs(nx - ir + 1, :)*dummy
            end do

            nmax = nx - nx_t

        end if

        ! -------------------------------------------------------------------
        ! interior points
        ! ic = 0                      ! central diagonal
        ic = 1                      ! 1. upper diagonal
        do ir = nmin, nmax
            if (.not. periodic_loc .and. ir + ic > nx) then      ! switch to the other side of the central diagonal
                ic = -ic
            end if
            dummy = 1.0_wp/rhs(ir, idr + ic)
            lhs(ir, 1:ndl) = lhs(ir, 1:ndl)*dummy
            rhs(ir, 1:ndr) = rhs(ir, 1:ndr)*dummy
        end do

        return
    end subroutine Precon_Rhs

end module Preconditioning

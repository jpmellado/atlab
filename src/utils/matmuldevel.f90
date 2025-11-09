! calculate f = A u, where A is narrow banded with diagonals given by rhs
! Allowing for different band size at the bottom and at the top
! The points at the boundary, can have a longer stencil by using the left of the rhs array

module MatMulDevel
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: MatMul_X              ! Generic procedures

    public :: MatMul_3              ! Particular procedures to accelerate
    public :: MatMul_3_antisym
    public :: MatMul_3_sym

    public :: MatMul_5
    public :: MatMul_5_antisym
    public :: MatMul_5_sym

    public :: MatMul_7_antisym
    public :: MatMul_7_sym

    public :: MatMul_X_ThomasL_Y

contains
    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_X(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr, ic

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ic)*u(:, ir - ic) + &
                           rhs(ir, idr + ic)*u(:, ir + ic)
            end do
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_X

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_3(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir - 1)*rhs(ir, 1) &
                       + u(:, ir)*rhs(ir, 2) &
                       + u(:, ir + 1)*rhs(ir, 3)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_3_antisym(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir + 1) - u(:, ir - 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_3_antisym

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_3_sym(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr
        real(wp) r2_loc             ! center diagonal

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r2_loc = rhs(nx_b + 1, 2)   ! Assume it is the same value for all points
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir)*r2_loc &
                       + u(:, ir + 1) + u(:, ir - 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_3_sym

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir - 2)*rhs(ir, 1) &
                       + u(:, ir - 1)*rhs(ir, 2) &
                       + u(:, ir)*rhs(ir, 3) &
                       + u(:, ir + 1)*rhs(ir, 4) &
                       + u(:, ir + 2)*rhs(ir, 5)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_5

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_antisym(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr
        real(wp) r5_loc             ! 2. upper-diagonal

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r5_loc = rhs(nx_b + 1, 5)   ! Assume it is the same value for all points
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir + 1) - u(:, ir - 1) &
                       + (u(:, ir + 2) - u(:, ir - 2))*r5_loc
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_5_antisym

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_sym(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr
        real(wp) r3_loc             ! center diagonal
        real(wp) r5_loc             ! 2. upper-diagonal

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r3_loc = rhs(nx_b + 1, 3)   ! Assume it is the same value for all points
        r5_loc = rhs(nx_b + 1, 5)
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir)*r3_loc &
                       + u(:, ir + 1) + u(:, ir - 1) &
                       + (u(:, ir + 2) + u(:, ir - 2))*r5_loc
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_5_sym

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_antisym(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr
        real(wp) r6_loc             ! 2. upper diagonal
        real(wp) r7_loc             ! 3. upper diagonal

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r6_loc = rhs(nx_b + 1, 6)   ! Assume it is the same value for all points
        r7_loc = rhs(nx_b + 1, 7)
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir + 1) - u(:, ir - 1) &
                       + (u(:, ir + 2) - u(:, ir - 2))*r6_loc &
                       + (u(:, ir + 3) - u(:, ir - 3))*r7_loc
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_7_antisym

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_sym(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr
        real(wp) r4_loc             ! center diagonal
        real(wp) r6_loc             ! 2. upper-diagonal
        real(wp) r7_loc             ! 3. upper-diagonal

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! -------------------------------------------------------------------
        if (present(bcs_b)) then
            call MatMul_X_LowerBoundary(rhs_b, u, f, bcs_b)
        else
            call MatMul_X_LowerBoundary(rhs_b, u, f)
        end if

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r4_loc = rhs(nx_b + 1, 4)   ! Assume it is the same value for all points
        r6_loc = rhs(nx_b + 1, 6)
        r7_loc = rhs(nx_b + 1, 7)
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir)*r4_loc &
                       + u(:, ir + 1) + u(:, ir - 1) &
                       + (u(:, ir + 2) + u(:, ir - 2))*r6_loc &
                       + (u(:, ir + 3) + u(:, ir - 3))*r7_loc
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if

        return
    end subroutine MatMul_7_sym

    ! ###################################################################
    ! ###################################################################
    ! Assumes that ndl is less or equal than nx_b, nx_t
    subroutine MatMul_X_ThomasL_Y(rhs, rhs_b, rhs_t, u, f, L)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir
        integer ndr, idr, ic
        integer ndl

        ! ###################################################################
        nx_b = size(rhs_b, 1)       ! size of system in lower boundary
        nx = size(rhs, 1)           ! size of the system
        nx_t = size(rhs_t, 1)       ! size of system in lower boundary

        ! array L
        ndl = size(L, 2)
        if (any([nx_b, nx_t] < ndl)) then
            print *, __FILE__//'Error'
        end if

        ! -------------------------------------------------------------------
        call MatMul_X_LowerBoundary(rhs_b, u, f)
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 1, nx_b
            do ic = 1, min(ir - 1, ndl)
                f(:, ir) = f(:, ir) + f(:, ir - ic)*L(ir, ndl - ic + 1)
            end do
        end do

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)      ! # of diagonals
        idr = ndr/2 + 1         ! index of centerline diagonal

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ic)*u(:, ir - ic) + &
                           rhs(ir, idr + ic)*u(:, ir + ic)
            end do
            ! Thomas step
            do ic = 1, ndl
                f(:, ir) = f(:, ir) + f(:, ir - ic)*L(ir, ndl - ic + 1)
            end do
        end do

        ! -------------------------------------------------------------------
        call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        do ir = nx_t - 1, 0, -1
            do ic = 1, ndl
                f(:, nx - ir) = f(:, nx - ir) + f(:, nx - ir - ic)*L(nx - ir, ndl - ic + 1)
            end do
        end do

        return
    end subroutine MatMul_X_ThomasL_Y

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_X_LowerBoundary(rhs, u, f, bcs)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs(:)

        integer(wi) ir
        integer nx, ndr, idr, ic

        ! ###################################################################
        nx = size(rhs, 1)           ! size of system
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        if (present(bcs)) then
            do ir = 2, min(idr, nx)
                f(:, ir) = rhs(ir, idr - ir + 1)*bcs(:)
                do ic = 2, ndr/2 + ir
                    f(:, ir) = f(:, ir) + &
                               rhs(ir, idr - ir + ic)*u(:, ic)
                end do
            end do
            ir = 1
            bcs(:) = rhs(ir, idr - ir + 1)*bcs(:)
            do ic = 2, ndr/2 + ir
                bcs(:) = bcs(:) + &
                         rhs(ir, idr - ir + ic)*u(:, ic)
            end do
            ! rhs(1,1) used for longer stencil at the boundary
            ic = ndr/2 + ir + 1
            bcs(:) = bcs(:) + &
                     rhs(ir, 1)*u(:, ic)

        else
            do ir = 1, min(idr, nx)
                f(:, ir) = rhs(ir, idr - ir + 1)*u(:, 1)
                do ic = 2, ndr/2 + ir
                    f(:, ir) = f(:, ir) + &
                               rhs(ir, idr - ir + ic)*u(:, ic)
                end do
            end do
            ! rhs(1,1) used for longer stencil at the boundary
            ir = 1; ic = ndr/2 + ir + 1
            f(:, ir) = f(:, ir) + &
                       rhs(ir, 1)*u(:, ic)

        end if

        do ir = min(idr, nx) + 1, nx
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, ndr/2
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr + ic)*u(:, ir + ic) + &
                           rhs(ir, idr - ic)*u(:, ir - ic)
            end do
        end do

        return
    end subroutine MatMul_X_LowerBoundary

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_X_UpperBoundary(rhs, u, f, bcs)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs(:)

        integer(wi) ir
        integer nx, ndr, idr, ic

        ! ###################################################################
        nx = size(rhs, 1)           ! size of system
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        do ir = nx - 1, min(idr, nx), -1
            f(:, nx - ir) = rhs(nx - ir, idr)*u(:, nx - ir)
            do ic = 1, ndr/2
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs(nx - ir, idr - ic)*u(:, nx - ir - ic) + &
                                rhs(nx - ir, idr + ic)*u(:, nx - ir + ic)
            end do
        end do

        if (present(bcs)) then
            do ir = min(idr, nx) - 1, 1, -1
                f(:, nx - ir) = rhs(nx - ir, idr + ir)*bcs(:)
                do ic = 1, ndr/2 + ir
                    f(:, nx - ir) = f(:, nx - ir) + &
                                    rhs(nx - ir, idr + ir - ic)*u(:, nx - ic)
                end do
            end do
            ir = 0
            bcs(:) = rhs(nx - ir, idr + ir)*bcs(:)
            do ic = 1, ndr/2 + ir
                bcs(:) = bcs(:) + &
                         rhs(nx - ir, idr + ir - ic)*u(:, nx - ic)
            end do
            ! rhs(nx,ndr) used for longer stencil at the boundary
            ic = ndr/2 + ir + 1
            bcs(:) = bcs(:) + &
                     rhs(nx, ndr)*u(:, nx - ic)

        else
            do ir = min(idr, nx) - 1, 0, -1
                f(:, nx - ir) = rhs(nx - ir, idr + ir)*u(:, nx)
                do ic = 1, ndr/2 + ir
                    f(:, nx - ir) = f(:, nx - ir) + &
                                    rhs(nx - ir, idr + ir - ic)*u(:, nx - ic)
                end do
            end do
            ! rhs(nx,ndr) used for longer stencil at the boundary
            ir = 0; ic = ndr/2 + ir + 1
            f(:, nx - ir) = f(:, nx - ir) + &
                            rhs(nx, ndr)*u(:, nx - ic)

        end if

        return
    end subroutine MatMul_X_UpperBoundary

end module MatMulDevel

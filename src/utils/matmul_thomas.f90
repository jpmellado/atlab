! Combining matrix multiplication with elimination step in Thomas algorithm
! Follows closely matmul procedures
! Assumes that ndl is less or equal than nx_b, nx_t

module MatMul_Thomas
    use TLab_Constants, only: wp, wi
    use MatMulDevel, only: MatMul_X_LowerBoundary, MatMul_X_UpperBoundary
    implicit none
    private

    public :: MatMul_X_ThomasL_Y                ! Generic procedures

    public :: MatMul_3_ThomasL_3                ! Particular procedures
    public :: MatMul_3_antisym_ThomasL_3
    public :: MatMul_3_sym_ThomasL_3            ! Check normalization below
    !                                             The first upper diagonal is normalized to 1
    public :: MatMul_5_ThomasL_3
    public :: MatMul_5_antisym_ThomasL_3
    public :: MatMul_5_sym_ThomasL_3
    public :: MatMul_5_ThomasL_5
    public :: MatMul_5_antisym_ThomasL_5
    public :: MatMul_5_sym_ThomasL_5
    public :: MatMul_5_sym_add_3_ThomasL_3
    public :: MatMul_5_sym_add_3_ThomasL_5

    public :: MatMul_7_antisym_ThomasL_3
    public :: MatMul_7_sym_ThomasL_3
    public :: MatMul_7_antisym_ThomasL_5
    public :: MatMul_7_sym_ThomasL_5
    public :: MatMul_7_sym_add_3_ThomasL_3
    public :: MatMul_7_sym_add_3_ThomasL_5

contains
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
    subroutine MatMul_3_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir - 1)*rhs(ir, 1) &
                       + u(:, ir)*rhs(ir, 2) &
                       + u(:, ir + 1)*rhs(ir, 3)
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_3_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_3_antisym_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir + 1) - u(:, ir - 1)
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            nx_thomas = nx
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_3_antisym_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_3_sym_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r2_loc = rhs(nx_b + 1, 2)   ! Assume it is the same value for all points
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir)*r2_loc &
                       + u(:, ir + 1) + u(:, ir - 1)
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_3_sym_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

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
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_antisym_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r5_loc = rhs(nx_b + 1, 5)   ! Assume it is the same value for all points
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir + 1) - u(:, ir - 1) &
                       + (u(:, ir + 2) - u(:, ir - 2))*r5_loc
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_antisym_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_sym_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

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
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_sym_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_ThomasL_5(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        ir = 2
        f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2)
        do ir = 3, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

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
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_ThomasL_5

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_antisym_ThomasL_5(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        ir = 2
        f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2)
        do ir = 3, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        ! interior points
        ndr = size(rhs, 2)          ! # of diagonals
        idr = ndr/2 + 1             ! index of center diagonal

        r5_loc = rhs(nx_b + 1, 5)   ! Assume it is the same value for all points
        do ir = nx_b + 1, nx - nx_t
            f(:, ir) = u(:, ir + 1) - u(:, ir - 1) &
                       + (u(:, ir + 2) - u(:, ir - 2))*r5_loc
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_antisym_ThomasL_5

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_sym_ThomasL_5(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        ir = 2
        f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2)
        do ir = 3, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

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
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_sym_ThomasL_5

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_sym_add_3_ThomasL_3(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: rhs_add(:, :)
        real(wp), intent(in) :: u_add(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Add second array
        ir = 1
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir)*rhs_add(ir, 2) &
                   + u_add(:, ir + 1)*rhs_add(ir, 3) &
                   + u_add(:, ir + 2)*rhs_add(ir, 1)
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

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
            ! Add second array
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Add second array
        do ir = nx - nx_t + 1, nx - 1
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ir = nx
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir - 2)*rhs_add(ir, 3) &
                   + u_add(:, ir - 1)*rhs_add(ir, 1) &
                   + u_add(:, ir)*rhs_add(ir, 2)
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_sym_add_3_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_5_sym_add_3_ThomasL_5(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: rhs_add(:, :)
        real(wp), intent(in) :: u_add(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Add second array
        ir = 1
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir)*rhs_add(ir, 2) &
                   + u_add(:, ir + 1)*rhs_add(ir, 3) &
                   + u_add(:, ir + 2)*rhs_add(ir, 1)
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        ir = 2
        f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2)
        do ir = 3, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

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
            ! Add second array
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Add second array
        do ir = nx - nx_t + 1, nx - 1
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ir = nx
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir - 2)*rhs_add(ir, 3) &
                   + u_add(:, ir - 1)*rhs_add(ir, 1) &
                   + u_add(:, ir)*rhs_add(ir, 2)
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        return
    end subroutine MatMul_5_sym_add_3_ThomasL_5

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_antisym_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

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
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_7_antisym_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_sym_ThomasL_3(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

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
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_7_sym_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_antisym_ThomasL_5(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        ir = 2
        f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2)
        do ir = 3, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

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
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        return
    end subroutine MatMul_7_antisym_ThomasL_5

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_sym_ThomasL_5(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        ir = 2
        f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2)
        do ir = 3, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

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
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        return
    end subroutine MatMul_7_sym_ThomasL_5

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_sym_add_3_ThomasL_3(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: rhs_add(:, :)
        real(wp), intent(in) :: u_add(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Add second array
        ir = 1
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir)*rhs_add(ir, 2) &
                   + u_add(:, ir + 1)*rhs_add(ir, 3) &
                   + u_add(:, ir + 2)*rhs_add(ir, 1)
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

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
            ! Add second array
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Add second array
        do ir = nx - nx_t + 1, nx - 1
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ir = nx
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir - 2)*rhs_add(ir, 3) &
                   + u_add(:, ir - 1)*rhs_add(ir, 1) &
                   + u_add(:, ir)*rhs_add(ir, 2)
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 1)
        end do

        return
    end subroutine MatMul_7_sym_add_3_ThomasL_3

    ! ###################################################################
    ! ###################################################################
    subroutine MatMul_7_sym_add_3_ThomasL_5(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, L, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: rhs_add(:, :)
        real(wp), intent(in) :: u_add(:, :)
        real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, nx_b, nx_t, ir, nx_thomas
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
        ! Add second array
        ir = 1
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir)*rhs_add(ir, 2) &
                   + u_add(:, ir + 1)*rhs_add(ir, 3) &
                   + u_add(:, ir + 2)*rhs_add(ir, 1)
        do ir = 2, nx_b
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ! Thomas step: nx_b is typically small, no need to interlace it with matmul loop
        ir = 2
        f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2)
        do ir = 3, nx_b
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

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
            ! Add second array
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
            ! Thomas step
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        ! -------------------------------------------------------------------
        if (present(bcs_t)) then
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx), bcs_t)
            nx_thomas = nx - 1
        else
            call MatMul_X_UpperBoundary(rhs_t, u(:, nx - nx_t + 1:nx), f(:, nx - nx_t + 1:nx))
            nx_thomas = nx
        end if
        ! Add second array
        do ir = nx - nx_t + 1, nx - 1
            f(:, ir) = f(:, ir) &
                       + u_add(:, ir - 1)*rhs_add(ir, 1) &
                       + u_add(:, ir)*rhs_add(ir, 2) &
                       + u_add(:, ir + 1)*rhs_add(ir, 3)
        end do
        ir = nx
        f(:, ir) = f(:, ir) &
                   + u_add(:, ir - 2)*rhs_add(ir, 3) &
                   + u_add(:, ir - 1)*rhs_add(ir, 1) &
                   + u_add(:, ir)*rhs_add(ir, 2)
        ! Thomas step: nx_t is typically small, no need to interlace it with matmul loop
        ! do not change last row if bcs are given
        do ir = nx - nx_t + 1, nx_thomas
            f(:, ir) = f(:, ir) + f(:, ir - 1)*L(ir, 2) + f(:, ir - 2)*L(ir, 1)
        end do

        return
    end subroutine MatMul_7_sym_add_3_ThomasL_5

end module MatMul_Thomas

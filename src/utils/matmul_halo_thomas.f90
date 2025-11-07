! Combining matrix multiplication with elimination step in Thomas algorithm
! Follows closely matmul procedures

module Matmul_Halo_Thomas
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: MatMul_Halo_X_ThomasL_Y               ! Generic procedures

    public :: MatMul_Halo_3_antisym_ThomasL_3       ! Particular procedures
    public :: MatMul_Halo_3_sym_ThomasL_3           ! Check normalization below
    !                                                 The first upper diagonal is normalized to 1
    public :: MatMul_Halo_5_antisym_ThomasL_3
    public :: MatMul_Halo_5_sym_ThomasL_3
    public :: MatMul_Halo_5_antisym_ThomasL_5
    public :: MatMul_Halo_5_sym_ThomasL_5

    public :: MatMul_Halo_7_antisym_ThomasL_3
    public :: MatMul_Halo_7_sym_ThomasL_3
    public :: MatMul_Halo_7_antisym_ThomasL_5
    public :: MatMul_Halo_7_sym_ThomasL_5

contains
    ! ###################################################################
    ! ###################################################################
    ! Assumes that ndl is less or equal than ndr/2
    subroutine MatMul_Halo_X_ThomasL_Y(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)
        real(wp), intent(in) :: L(:, :)

        integer(wi) nx, ir
        integer(wi) ndr, idr, ic
        integer(wi) ndl

        ! ###################################################################
        ndr = size(rhs, 2)      ! # of diagonals
        idr = ndr/2 + 1         ! index of centerline diagonal
        nx = size(rhs, 1)       ! size of the system

        ndl = size(L, 2)
        if (ndr/2 < ndl) then
            print *, __FILE__//'Error'
        end if

        ! -------------------------------------------------------------------
        ! lower boundary
        do ir = 1, ndr/2
            f(:, ir) = rhs(ir, idr - ir + 1)*u(:, 1)
            do ic = 2, ndr/2 + ir
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ir + ic)*u(:, ic)
            end do
            do ic = 0, ndr/2 - ir
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ir - ic)*u_halo_m(:, ndr/2 - ic)
            end do
            do ic = 1, min(ir - 1, ndl)  ! solve L
                f(:, ir) = f(:, ir) + f(:, ir - ic)*L(ir, ndl - ic + 1)
            end do
        end do

        ! -------------------------------------------------------------------
        ! interior points
        do ir = idr, nx - idr + 1
            f(:, ir) = rhs(ir, idr)*u(:, ir)
            do ic = 1, idr - 1
                f(:, ir) = f(:, ir) + &
                           rhs(ir, idr - ic)*u(:, ir - ic) + &
                           rhs(ir, idr + ic)*u(:, ir + ic)
            end do
            do ic = 1, ndl      ! solve L
                f(:, ir) = f(:, ir) + f(:, ir - ic)*L(ir, ndl - ic + 1)
            end do
        end do

        ! -------------------------------------------------------------------
        ! upper boundary
        do ir = ndr/2 - 1, 0, -1
            f(:, nx - ir) = rhs(nx - ir, idr + ir)*u(:, nx)
            do ic = 1, ndr/2 + ir
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs(nx - ir, idr + ir - ic)*u(:, nx - ic)
            end do
            do ic = 1, ndr/2 - ir
                f(:, nx - ir) = f(:, nx - ir) + &
                                rhs(nx - ir, idr + ir + ic)*u_halo_p(:, ic)
            end do
            do ic = 1, ndl      ! solve L
                f(:, nx - ir) = f(:, nx - ir) + f(:, nx - ir - ic)*L(nx - ir, ndl - ic + 1)
            end do
        end do

        return
    end subroutine MatMul_Halo_X_ThomasL_Y

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_3_antisym_ThomasL_3(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! #######################################################################
        nx = size(f, 2)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = u(:, n + 1) - u_halo_m(:, 1)

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 2, nx - 1
            f(:, n) = u(:, n + 1) - u(:, n - 1)
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx
        f(:, n) = u_halo_p(:, 1) - u(:, n - 1)
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        return
    end subroutine MatMul_Halo_3_antisym_ThomasL_3

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_3_sym_ThomasL_3(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r2_loc     ! center diagonal

        ! #######################################################################
        nx = size(f, 2)
        r2_loc = rhs(2)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = r2_loc*u(:, n) &
                  + u(:, n + 1) + u_halo_m(:, 1)

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 2, nx - 1
            f(:, n) = r2_loc*u(:, n) &
                      + u(:, n + 1) + u(:, n - 1)
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx
        f(:, n) = r2_loc*u(:, n) &
                  + u_halo_p(:, 1) + u(:, n - 1)
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        return
    end subroutine MatMul_Halo_3_sym_ThomasL_3

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_5_antisym_ThomasL_3(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r5_loc     ! 2. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r5_loc = rhs(5)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = u(:, n + 1) - u_halo_m(:, 2) &
                  + r5_loc*(u(:, n + 2) - u_halo_m(:, 1))

        n = 2
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r5_loc*(u(:, n + 2) - u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 3, nx - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) - u(:, n - 2))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 1
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 1) - u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = nx
        f(:, n) = u_halo_p(:, 1) - u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 2) - u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        return
    end subroutine MatMul_Halo_5_antisym_ThomasL_3

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_5_sym_ThomasL_3(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r5_loc = rhs(5)
        r3_loc = rhs(3)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = r3_loc*u(:, n) &
                  + u(:, n + 1) + u_halo_m(:, 2) &
                  + r5_loc*(u(:, n + 2) + u_halo_m(:, 1))

        n = 2
        f(:, n) = r3_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r5_loc*(u(:, n + 2) + u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 3, nx - 2
            f(:, n) = r3_loc*u(:, n) &
                      + u(:, n + 1) + u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) + u(:, n - 2))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 1
        f(:, n) = r3_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 1) + u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = nx
        f(:, n) = r3_loc*u(:, n) &
                  + u_halo_p(:, 1) + u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 2) + u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        return
    end subroutine MatMul_Halo_5_sym_ThomasL_3

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_5_antisym_ThomasL_5(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r5_loc     ! 2. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r5_loc = rhs(5)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = u(:, n + 1) - u_halo_m(:, 2) &
                  + r5_loc*(u(:, n + 2) - u_halo_m(:, 1))

        n = 2
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r5_loc*(u(:, n + 2) - u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2)  ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 3, nx - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) - u(:, n - 2))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 1
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 1) - u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        n = nx
        f(:, n) = u_halo_p(:, 1) - u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 2) - u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        return
    end subroutine MatMul_Halo_5_antisym_ThomasL_5

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_5_sym_ThomasL_5(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r5_loc = rhs(5)
        r3_loc = rhs(3)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = r3_loc*u(:, n) &
                  + u(:, n + 1) + u_halo_m(:, 2) &
                  + r5_loc*(u(:, n + 2) + u_halo_m(:, 1))

        n = 2
        f(:, n) = r3_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r5_loc*(u(:, n + 2) + u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2)  ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 3, nx - 2
            f(:, n) = r3_loc*u(:, n) &
                      + u(:, n + 1) + u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) + u(:, n - 2))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 1
        f(:, n) = r3_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 1) + u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        n = nx
        f(:, n) = r3_loc*u(:, n) &
                  + u_halo_p(:, 1) + u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 2) + u(:, n - 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        return
    end subroutine MatMul_Halo_5_sym_ThomasL_5

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_7_antisym_ThomasL_3(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r7_loc = rhs(7)
        r6_loc = rhs(6)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = u(:, n + 1) - u_halo_m(:, 3) &
                  + r6_loc*(u(:, n + 2) - u_halo_m(:, 2)) &
                  + r7_loc*(u(:, n + 3) - u_halo_m(:, 1))

        n = 2
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) - u_halo_m(:, 3)) &
                  + r7_loc*(u(:, n + 3) - u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = 3
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) - u(:, n - 2)) &
                  + r7_loc*(u(:, n + 3) - u_halo_m(:, 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 4, nx - 3
            f(:, n) = u(:, n + 1) - u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) - u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) - u(:, n - 3))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 2
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) - u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 1) - u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = nx - 1
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 1) - u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 2) - u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = nx
        f(:, n) = u_halo_p(:, 1) - u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 2) - u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 3) - u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        return
    end subroutine MatMul_Halo_7_antisym_ThomasL_3

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_7_sym_ThomasL_3(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r7_loc = rhs(7)
        r6_loc = rhs(6)
        r4_loc = rhs(4)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u_halo_m(:, 3) &
                  + r6_loc*(u(:, n + 2) + u_halo_m(:, 2)) &
                  + r7_loc*(u(:, n + 3) + u_halo_m(:, 1))

        n = 2
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u_halo_m(:, 3)) &
                  + r7_loc*(u(:, n + 3) + u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = 3
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                  + r7_loc*(u(:, n + 3) + u_halo_m(:, 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 4, nx - 3
            f(:, n) = r4_loc*u(:, n) &
                      + u(:, n + 1) + u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) + u(:, n - 3))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 2
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 1) + u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = nx - 1
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 1) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 2) + u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        n = nx
        f(:, n) = r4_loc*u(:, n) &
                  + u_halo_p(:, 1) + u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 2) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 3) + u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 1)  ! solve L

        return
    end subroutine MatMul_Halo_7_sym_ThomasL_3

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_7_antisym_ThomasL_5(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r7_loc = rhs(7)
        r6_loc = rhs(6)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = u(:, n + 1) - u_halo_m(:, 3) &
                  + r6_loc*(u(:, n + 2) - u_halo_m(:, 2)) &
                  + r7_loc*(u(:, n + 3) - u_halo_m(:, 1))

        n = 2
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) - u_halo_m(:, 3)) &
                  + r7_loc*(u(:, n + 3) - u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2)  ! solve L

        n = 3
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) - u(:, n - 2)) &
                  + r7_loc*(u(:, n + 3) - u_halo_m(:, 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 4, nx - 3
            f(:, n) = u(:, n + 1) - u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) - u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) - u(:, n - 3))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 2
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) - u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 1) - u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        n = nx - 1
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 1) - u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 2) - u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        n = nx
        f(:, n) = u_halo_p(:, 1) - u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 2) - u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 3) - u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        return
    end subroutine MatMul_Halo_7_antisym_ThomasL_5

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_7_sym_ThomasL_5(rhs, u, u_halo_m, u_halo_p, f, L)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u
        real(wp), intent(in) :: L(:, :)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal

        ! #######################################################################
        nx = size(f, 2)
        r7_loc = rhs(7)
        r6_loc = rhs(6)
        r4_loc = rhs(4)

        ! -------------------------------------------------------------------
        ! Halo left
        n = 1
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u_halo_m(:, 3) &
                  + r6_loc*(u(:, n + 2) + u_halo_m(:, 2)) &
                  + r7_loc*(u(:, n + 3) + u_halo_m(:, 1))

        n = 2
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u_halo_m(:, 3)) &
                  + r7_loc*(u(:, n + 3) + u_halo_m(:, 2))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2)  ! solve L

        n = 3
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                  + r7_loc*(u(:, n + 3) + u_halo_m(:, 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 4, nx - 3
            f(:, n) = r4_loc*u(:, n) &
                      + u(:, n + 1) + u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) + u(:, n - 3))
            f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 2
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 1) + u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        n = nx - 1
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 1) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 2) + u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        n = nx
        f(:, n) = r4_loc*u(:, n) &
                  + u_halo_p(:, 1) + u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 2) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 3) + u(:, n - 3))
        f(:, n) = f(:, n) + f(:, n - 1)*L(n, 2) + f(:, n - 2)*L(n, 1) ! solve L

        return
    end subroutine MatMul_Halo_7_sym_ThomasL_5

end module Matmul_Halo_Thomas

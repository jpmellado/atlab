module FDM_MatMul_Halo
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: MatMul_Halo_5d_antisym
    public :: MatMul_Halo_7d_sym

contains
    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_5d_antisym(rhs, u, u_halo_m, u_halo_p, f)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u

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

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 3, nx - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) - u(:, n - 2))
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 1
        f(:, n) = u(:, n + 1) - u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 1) - u(:, n - 2))

        n = nx
        f(:, n) = u_halo_p(:, 1) - u(:, n - 1) &
                  + r5_loc*(u_halo_p(:, 2) - u(:, n - 2))

        return
    end subroutine MatMul_Halo_5d_antisym

    !########################################################################
    !########################################################################
    subroutine MatMul_Halo_7d_sym(rhs, u, u_halo_m, u_halo_p, f)
        real(wp), intent(in) :: rhs(:)              ! diagonals of B
        real(wp), intent(in) :: u(:, :)             ! vector u
        real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        real(wp), intent(out) :: f(:, :)            ! vector f = B u

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

        n = 3
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                  + r7_loc*(u(:, n + 3) + u_halo_m(:, 3))

        ! -------------------------------------------------------------------
        ! Interior points
        do n = 4, nx - 3
            f(:, n) = r4_loc*u(:, n) &
                      + u(:, n + 1) + u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) + u(:, n - 3))
        end do

        ! -------------------------------------------------------------------
        ! Halo right
        n = nx - 2
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 1) + u(:, n - 3))
        n = nx - 1
        f(:, n) = r4_loc*u(:, n) &
                  + u(:, n + 1) + u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 1) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 2) + u(:, n - 3))

        n = nx
        f(:, n) = r4_loc*u(:, n) &
                  + u_halo_p(:, 1) + u(:, n - 1) &
                  + r6_loc*(u_halo_p(:, 2) + u(:, n - 2)) &
                  + r7_loc*(u_halo_p(:, 3) + u(:, n - 3))

        return
    end subroutine MatMul_Halo_7d_sym

end module FDM_MatMul_Halo

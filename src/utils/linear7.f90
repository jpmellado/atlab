! Vectorized Thomas algorithm for pentadiagonal systems
! LU factorization stage with unit diagonal in L
! First equation is divided by d1 to get 1 in the first diagonal elem.

module Thomas7
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: Thomas7_LU, Thomas7_Solve

contains
    ! #######################################################################
    ! #######################################################################
    subroutine Thomas7_LU(nmax, a, b, c, d, e, f, g)
        integer(wi) nmax
        real(wp), dimension(nmax), intent(INOUT) :: a, b, c, d, e, f, g

        ! -----------------------------------------------------------------------
        integer(wi) n

        ! #######################################################################
        g(1) = g(1)/d(1)
        f(1) = f(1)/d(1)
        e(1) = e(1)/d(1)
        c(1) = 1.0_wp/d(1) ! padding, and used in Thomas7_Solve to normalize 1st eqn.
        !   b(1) = 1.0_wp      ! padding
        !   a(1) = 1.0_wp      ! padding
        d(1) = 1.0_wp

        !   a(2) = 1.0_wp ! padding
        !   b(2) = 1.0_wp ! padding
        c(2) = c(2)/d(1)
        d(2) = d(2) - c(2)*e(1)
        e(2) = e(2) - c(2)*f(1)
        f(2) = f(2) - c(2)*g(1)

        !   a(3) = 1.0_wp ! padding
        b(3) = b(3)/d(1)
        c(3) = (c(3) - b(3)*e(1))/d(2)
        d(3) = d(3) - c(3)*e(2) - b(3)*f(1)
        e(3) = e(3) - c(3)*f(2) - b(3)*g(1)
        f(3) = f(3) - c(3)*g(2)

        do n = 4, nmax - 2
            a(n) = a(n)/d(n - 3)
            b(n) = (b(n) - a(n)*e(n - 3))/d(n - 2)
            c(n) = (c(n) - b(n)*e(n - 2) - a(n)*f(n - 3))/d(n - 1)
            d(n) = d(n) - c(n)*e(n - 1) - b(n)*f(n - 2) - a(n)*g(n - 3)
            e(n) = e(n) - c(n)*f(n - 1) - b(n)*g(n - 2)
            f(n) = f(n) - c(n)*g(n - 1)
        end do
        !   g(n-1) = 1.0_wp ! padding

        n = nmax - 1
        a(n) = a(n)/d(n - 3)
        b(n) = (b(n) - a(n)*e(n - 3))/d(n - 2)
        c(n) = (c(n) - b(n)*e(n - 2) - a(n)*f(n - 3))/d(n - 1)
        d(n) = d(n) - c(n)*e(n - 1) - b(n)*f(n - 2) - a(n)*g(n - 3)
        e(n) = e(n) - c(n)*f(n - 1) - b(n)*g(n - 2)
        !   f(n) = 1.0_wp ! padding
        !   g(n) = 1.0_wp ! padding

        n = nmax
        a(n) = a(n)/d(n - 3)
        b(n) = (b(n) - a(n)*e(n - 3))/d(n - 2)
        c(n) = (c(n) - b(n)*e(n - 2) - a(n)*f(n - 3))/d(n - 1)
        d(n) = d(n) - c(n)*e(n - 1) - b(n)*f(n - 2) - a(n)*g(n - 3)
        !   e(n) = 1.0_wp ! padding
        !   f(n) = 1.0_wp ! padding
        !   g(n) = 1.0_wp ! padding

        ! Final operations
        a(4:) = -a(4:)
        b(3:) = -b(3:)
        c(2:) = -c(2:)
        d = 1.0_wp/d
        e(:nmax - 1) = -e(:nmax - 1)
        f(:nmax - 2) = -f(:nmax - 2)
        g(:nmax - 3) = -g(:nmax - 3)

        return
    end subroutine Thomas7_LU

    ! #######################################################################
    ! Backward substitution step in the Thomas algorithm
    ! #######################################################################
    subroutine Thomas7_Solve(nmax, len, a, b, c, d, e, f, g, frc)
        integer(wi) nmax, len
        real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e, f, g
        real(wp), dimension(len, nmax), intent(INOUT) :: frc

        ! -----------------------------------------------------------------------
        integer(wi) n
        real(wp) :: dummy1, dummy2, dummy3

        ! #######################################################################
        ! -----------------------------------------------------------------------
        ! Solve Ly=frc, forward
        ! -----------------------------------------------------------------------
        frc(:, 1) = frc(:, 1)*c(1) ! Normalize first eqn. See Thomas7_LU
        frc(:, 2) = frc(:, 2) + frc(:, 1)*c(2)
        frc(:, 3) = frc(:, 3) + frc(:, 2)*c(3) + frc(:, 1)*b(3)

        do n = 4, nmax
            frc(:, n) = frc(:, n) + frc(:, n - 1)*c(n) + frc(:, n - 2)*b(n) + frc(:, n - 3)*a(n)
        end do

        ! -----------------------------------------------------------------------
        ! Solve Ux=y, backward
        ! -----------------------------------------------------------------------
        n = nmax
        frc(:, n) = frc(:, n)*d(n)

        n = nmax - 1
        dummy1 = d(n)*e(n)
        frc(:, n) = d(n)*frc(:, n) + dummy1*frc(:, n + 1)

        n = nmax - 2
        dummy1 = d(n)*e(n)
        dummy2 = d(n)*f(n)
        frc(:, n) = d(n)*frc(:, n) + dummy1*frc(:, n + 1) + dummy2*frc(:, n + 2)

        do n = nmax - 3, 1, -1
            dummy1 = d(n)*e(n)
            dummy2 = d(n)*f(n)
            dummy3 = d(n)*g(n)
            frc(:, n) = d(n)*frc(:, n) + dummy1*frc(:, n + 1) + dummy2*frc(:, n + 2) + dummy3*frc(:, n + 3)
        end do

        return
    end subroutine Thomas7_Solve

end module Thomas7

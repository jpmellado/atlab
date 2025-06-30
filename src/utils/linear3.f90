! Vectorized Thomas algorithm for tridiagonal systems

module Thomas3
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: Thomas3_LU, Thomas3_Solve
    public :: Thomas3P_LU, Thomas3P_Solve   ! circulant systems (periodic boundary conditions)

contains
    ! #######################################################################
    ! #######################################################################
    ! LU factorization stage; L with diagonal unity
    subroutine Thomas3_LU(nmax, a, b, c)
        integer(wi), intent(IN) :: nmax
        real(wp), dimension(nmax), intent(INOUT) :: a, b, c

        ! -----------------------------------------------------------------------
        integer(wi) n

        ! #######################################################################
        do n = 2, nmax
            a(n) = a(n)/b(n - 1)
            b(n) = b(n) - a(n)*c(n - 1)
        end do

        ! Final operations
        a = -a
        b = 1.0_wp/b
        c = -c*b

        return
    end subroutine Thomas3_LU

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas3_Solve(nmax, len, a, b, c, f)
        integer(wi), intent(IN) :: nmax                         ! dimension of tridiagonal systems
        integer(wi), intent(IN) :: len                          ! number of systems to solve
        real(wp), intent(IN) :: a(nmax), b(nmax), c(nmax)       ! factored LHS
        real(wp), intent(INOUT) :: f(len, nmax)                 ! RHS and solution

        ! -------------------------------------------------------------------
        integer(wi) :: n
        real(wp) :: dummy1, dummy2

        ! ###################################################################
        if (len <= 0) return

        ! -----------------------------------------------------------------------
        ! Forward sweep
        ! -----------------------------------------------------------------------
        do n = 2, nmax
            dummy1 = a(n)
            f(:, n) = f(:, n) + dummy1*f(:, n - 1)
        end do

        ! -----------------------------------------------------------------------
        ! Backward sweep
        ! -----------------------------------------------------------------------
        dummy1 = b(nmax)
        f(:, nmax) = f(:, nmax)*dummy1

        do n = nmax - 1, 1, -1
            dummy1 = c(n)
            dummy2 = b(n)
            f(:, n) = dummy2*f(:, n) + dummy1*f(:, n + 1)

        end do

        return
    end subroutine Thomas3_Solve

    !########################################################################
    !########################################################################
    subroutine Thomas3P_LU(nmax, a, b, c, d, e)
        integer(wi) nmax
        real(wp), dimension(nmax) :: a, b, c, d, e

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) sum

        ! ###################################################################
        ! Generate first elements of LU
        c(1) = c(1)/b(1)
        e(1) = a(1)/b(1)
        d(1) = c(nmax)

        ! Generate n=2 to n=n-2 elements of LU
        do n = 2, nmax - 2
            b(n) = b(n) - a(n)*c(n - 1)
            c(n) = c(n)/b(n)
            e(n) = -a(n)*e(n - 1)/b(n)
            d(n) = -d(n - 1)*c(n - 1)
        end do

        ! Generate n-1 elements
        b(nmax - 1) = b(nmax - 1) - a(nmax - 1)*c(nmax - 2)
        e(nmax - 1) = (c(nmax - 1) - a(nmax - 1)*e(nmax - 2))/b(nmax - 1)
        d(nmax - 1) = a(nmax) - d(nmax - 2)*c(nmax - 2)

        ! Generate the n-th element
        sum = 0.0_wp
        do n = 1, nmax - 1
            sum = sum + d(n)*e(n)
        end do
        b(nmax) = b(nmax) - sum

        ! Final operations
        do n = 1, nmax
            b(n) = 1.0_wp/b(n)
            a(n) = -a(n)*b(n)
            c(n) = -c(n)
            e(n) = -e(n)
        end do

        return
    end subroutine Thomas3P_LU

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas3P_Solve(nmax, len, a, b, c, d, e, f, wrk)
        integer(wi), intent(IN) :: nmax                         ! dimension of tridiagonal systems
        integer(wi), intent(IN) :: len                          ! number of systems to solve
        real(wp), intent(IN) :: a(nmax), b(nmax), c(nmax)       ! factored LHS
        real(wp), intent(IN) :: d(nmax), e(nmax)
        real(wp), intent(inout) :: f(len, nmax)                 ! RHS and solution
        real(wp), intent(inout) :: wrk(len)

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) :: dummy1, dummy2

        ! ###################################################################
        if (len <= 0) return

        ! -------------------------------------------------------------------
        ! Forward sweep
        ! -------------------------------------------------------------------
        dummy1 = b(1)
        f(:, 1) = f(:, 1)*dummy1

        do n = 2, nmax - 1
            dummy1 = a(n)
            dummy2 = b(n)
            f(:, n) = f(:, n)*dummy2 + dummy1*f(:, n - 1)
        end do

        wrk(:) = 0.0_wp

        do n = 1, nmax - 1
            dummy1 = d(n)
            wrk(:) = wrk(:) + dummy1*f(:, n)
        end do

        dummy1 = b(nmax)
        f(:, nmax) = (f(:, nmax) - wrk(:))*dummy1

        ! -------------------------------------------------------------------
        ! Backward sweep
        ! -------------------------------------------------------------------
        dummy1 = e(nmax - 1)
        f(:, nmax - 1) = dummy1*f(:, nmax) + f(:, nmax - 1)

        do n = nmax - 2, 1, -1
            dummy1 = c(n)
            dummy2 = e(n)
            f(:, n) = f(:, n) + dummy1*f(:, n + 1) + dummy2*f(:, nmax)
        end do

        return
    end subroutine Thomas3P_Solve

end module Thomas3

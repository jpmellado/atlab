#include "tlab_error.h"

! Vectorized Thomas algorithm for pentadiagonal systems
! LU factorization stage with unit diagonal in L

module Thomas5
    use TLab_Constants, only: wp, wi, small_wp, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: Thomas5_FactorLU, Thomas5_SolveLU
    public :: Thomas5C_SMW_LU, Thomas5C_SMW_Solve   ! circulant systems (periodic boundary conditions)

contains
    ! #######################################################################
    ! #######################################################################
    subroutine Thomas5_FactorLU(nmax, a, b, c, d, e)
        integer(wi) nmax
        real(wp), intent(inout) :: a(nmax), b(nmax), c(nmax), d(nmax), e(nmax)

        ! -----------------------------------------------------------------------
        integer(wi) n

        ! #######################################################################
        n = 2
        b(n) = b(n)/c(n - 1)
        c(n) = c(n) - b(n)*d(n - 1)
        d(n) = d(n) - b(n)*e(n - 1)

        do n = 3, nmax - 1
            a(n) = a(n)/c(n - 2)
            b(n) = (b(n) - a(n)*d(n - 2))/c(n - 1)
            c(n) = c(n) - b(n)*d(n - 1) - a(n)*e(n - 2)
            d(n) = d(n) - b(n)*e(n - 1)
        end do

        n = nmax
        a(n) = a(n)/c(n - 2)
        b(n) = (b(n) - a(n)*d(n - 2))/c(n - 1)
        c(n) = c(n) - b(n)*d(n - 1) - a(n)*e(n - 2)

        ! Final operations
        a(3:) = -a(3:)
        b(2:) = -b(2:)
        c = 1.0_wp/c
        d(:nmax - 1) = -d(:nmax - 1)*c(:nmax - 1)
        e(:nmax - 2) = -e(:nmax - 2)*c(:nmax - 2)

        return
    end subroutine Thomas5_FactorLU

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas5_SolveLU(nmax, len, a, b, c, d, e, f)
        integer(wi) nmax, len
        real(wp), intent(in) :: a(nmax), b(nmax), c(nmax), d(nmax), e(nmax)
        real(wp), intent(inout) :: f(len, nmax)

        ! -----------------------------------------------------------------------
        integer(wi) n
        real(wp) :: dummy1, dummy2

        ! #######################################################################
        ! Solve Ly=f, forward elimination
        n = 2
        f(:, n) = f(:, n) + f(:, n - 1)*b(2)

        do n = 3, nmax
            f(:, n) = f(:, n) + f(:, n - 1)*b(n) + f(:, n - 2)*a(n)
        end do

        ! Solve Ux=y, backward substitution
        n = nmax
        f(:, n) = f(:, n)*c(n)

        n = nmax - 1
        dummy1 = d(n)
        f(:, n) = c(n)*f(:, n) + dummy1*f(:, n + 1)

        do n = nmax - 2, 1, -1
            dummy1 = d(n)
            dummy2 = e(n)
            f(:, n) = c(n)*f(:, n) + dummy1*f(:, n + 1) + dummy2*f(:, n + 2)
        end do

        return
    end subroutine Thomas5_SolveLU

    !########################################################################
    !#
    !# 2021/12/09 - J. Kostelecky
    !#              Created
    !#
    !# Pentadiagonal Toeplitz solver for a circulant system (periodic bcs),
    !# based on the Sherman–Morrison–Woodbury Formula, described in: "A new
    !# algorithm for solving nearly penta-diagonal Toeplitz linear systems"
    !# https://doi.org/10.1016/j.camwa.2011.12.044
    !# Algorithm 2.2 is implemented here.
    !# Ensure that Matrix M is invertible!
    !#
    !########################################################################

    ! #######################################################################
    ! LU factorization stage
    ! #######################################################################
    subroutine Thomas5C_SMW_LU(nmax, a, b, c, d, e, f, g)
        integer(wi), intent(in) :: nmax
        real(wp), intent(inout) :: a(nmax), b(nmax), c(nmax), d(nmax), e(nmax) ! Diagonals
        real(wp), intent(out) :: f(nmax), g(nmax)       ! Additional (u1,u2)

        ! -----------------------------------------------------------------------
        real(wp) :: a0, b0, en, dn
        real(wp) :: m1, m2, m3, m4

        ! #######################################################################
        ! Build regular modified pentadiagonal matrix A1 (eq. 2.8)

        ! Save cyclic entries
        a0 = a(1)    ! upper-right corner
        b0 = b(1)
        en = e(nmax) ! lower-left  corner
        dn = d(nmax)

        ! Modified entries of A1
        b(2) = b(2) - d(nmax)
        c(1) = c(1) - e(nmax)
        c(2) = c(2) - e(nmax)

        c(nmax - 1) = c(nmax - 1) - a(1)
        c(nmax) = c(nmax) - a(1)
        d(nmax - 1) = d(nmax - 1) - b(1)

        ! Set off-digonal entries to zero for A1
        a(1) = 0.0_wp
        a(2) = 0.0_wp
        b(1) = 0.0_wp

        d(nmax) = 0.0_wp
        e(nmax) = 0.0_wp
        e(nmax - 1) = 0.0_wp

        ! Regular forward step for A1
        call Thomas5_FactorLU(nmax, a, b, c, d, e)

        ! Save cyclic entries again
        a(1) = a0 ! upper-right corner
        b(1) = b0
        e(nmax) = en ! lower-left corner
        d(nmax) = dn

        ! Define matrix u [here: u1,u2 stored in additional diagonals f, g] (eq. 2.8)
        f = 0.0_wp ! u1
        f(1) = 1.0_wp
        f(nmax - 1) = 1.0_wp
        g = 0.0_wp ! u2
        g(2) = 1.0_wp
        g(nmax) = 1.0_wp

        ! Regular backward step for u1, u2
        call Thomas5_SolveLU(nmax, 1, a, b, c, d, e, f)
        call Thomas5_SolveLU(nmax, 1, a, b, c, d, e, g)

        ! Compute entries of matrix M[2x2] once
        m1 = e(nmax)*f(1) + a(1)*f(nmax - 1) + b(1)*f(nmax) + 1.0_wp
        m2 = e(nmax)*g(1) + a(1)*g(nmax - 1) + b(1)*g(nmax)
        m3 = d(nmax)*f(1) + e(nmax)*f(2) + a(1)*f(nmax)
        m4 = d(nmax)*g(1) + e(nmax)*g(2) + a(1)*g(nmax) + 1.0_wp
        ! Check if M is invertible (eq. 2.9)
        if ((m1*m4 - m2*m3) < small_wp) then
            call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
            call TLab_Stop(DNS_ERROR_THOMAS)
        end if

        return
    end subroutine Thomas5C_SMW_LU

    ! #######################################################################
    ! Backward substitution step in the Thomas algorithm
    ! #######################################################################
    subroutine Thomas5C_SMW_Solve(nmax, len, a, b, c, d, e, f, g, frc)
        integer(wi), intent(IN) :: nmax, len
        real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e ! Diagonals
        real(wp), dimension(nmax), intent(IN) :: f, g       ! Additional (z,w)
        real(wp), dimension(len, nmax), intent(INOUT) :: frc       ! Rhs

        ! -----------------------------------------------------------------------
        integer(wi) :: n, l
        real(wp) :: m1, m2, m3, m4
        real(wp) :: di, d11, d12, d13, d14, d21, d22, d23, d24
        real(wp) :: dummy1, dummy2

        ! #######################################################################

        ! Regular backward step for rhs
        call Thomas5_SolveLU(nmax, len, a, b, c, d, e, frc)

        ! Compute entries of matrix m[2x2]
        m1 = e(nmax)*f(1) + a(1)*f(nmax - 1) + b(1)*f(nmax) + 1.0_wp
        m2 = e(nmax)*g(1) + a(1)*g(nmax - 1) + b(1)*g(nmax)
        m3 = d(nmax)*f(1) + e(nmax)*f(2) + a(1)*f(nmax)
        m4 = d(nmax)*g(1) + e(nmax)*g(2) + a(1)*g(nmax) + 1.0_wp

        ! Compute coefficients
        di = 1/(m1*m4 - m2*m3)
        d11 = di*(m4*e(nmax) - m2*d(nmax))
        d12 = di*(m4*b(1) - m2*a(1))
        d13 = di*m4*a(1)
        d14 = di*m2*e(nmax)
        d21 = di*(m1*d(nmax) - m3*e(nmax))
        d22 = di*(m1*a(1) - m3*b(1))
        d23 = di*m3*a(1)
        d24 = di*m1*e(nmax)

        ! Solve
        do n = 3, nmax - 3, 1 ! Main loop
            do l = 1, len, 1
                dummy1 = d11*frc(l, 1) + d12*frc(l, nmax) + d13*frc(l, nmax - 1) - d14*frc(l, 2)
                dummy2 = d21*frc(l, 1) + d22*frc(l, nmax) - d23*frc(l, nmax - 1) + d24*frc(l, 2)
                !
                frc(l, n) = frc(l, n) - dummy1*f(n) - dummy2*g(n)
            end do
        end do
        !
        do l = 1, len, 1     ! Boundaries
            dummy1 = d11*frc(l, 1) + d12*frc(l, nmax) + d13*frc(l, nmax - 1) - d14*frc(l, 2)
            dummy2 = d21*frc(l, 1) + d22*frc(l, nmax) - d23*frc(l, nmax - 1) + d24*frc(l, 2)
            do n = 1, 2, 1
                frc(l, n) = frc(l, n) - dummy1*f(n) - dummy2*g(n)
            end do
            do n = nmax - 2, nmax, 1
                frc(l, n) = frc(l, n) - dummy1*f(n) - dummy2*g(n)
            end do
        end do

        return
    end subroutine Thomas5C_SMW_Solve

end module Thomas5

! !# Reverse ordering.

! ! #######################################################################
! ! LU factorization stage
! ! #######################################################################
! subroutine Thomas5_FactorLU2(nmax, a, b, c, d, e)
!     use TLab_Constants, only: wp, wi

!     implicit none

!     integer(wi) nmax
!     real(wp), dimension(nmax), intent(INOUT) :: a, b, c, d, e

!     ! -----------------------------------------------------------------------
!     integer(wi) n

!     ! #######################################################################
!     n = nmax
!     e(n) = 1.0_wp ! padding
!     d(n) = 1.0_wp ! padding
!     c(n) = c(n)
!     b(n) = b(n)

!     n = nmax - 1
!     e(n) = 1.0_wp ! padding
!     d(n) = (d(n))/c(n + 1)
!     c(n) = c(n) - d(n)*b(n + 1)
!     b(n) = b(n) - d(n)*a(n + 1)

!     do n = nmax - 2, 3, -1
!         e(n) = e(n)/c(n + 2)
!         d(n) = (d(n) - e(n)*b(n + 2))/c(n + 1)
!         c(n) = c(n) - d(n)*b(n + 1) - e(n)*a(n + 2)
!         b(n) = b(n) - d(n)*a(n + 1)
!     end do

!     n = 2
!     e(n) = e(n)/c(n + 2)
!     d(n) = (d(n) - e(n)*b(n + 2))/c(n + 1)
!     c(n) = c(n) - d(n)*b(n + 1) - e(n)*a(n + 2)
!     b(n) = b(n) - d(n)*a(n + 1)
!     a(n) = 1.0_wp ! padding

!     n = 1
!     e(n) = e(n)/c(n + 2)
!     d(n) = (d(n) - e(n)*b(n + 2))/c(n + 1)
!     c(n) = c(n) - d(n)*b(n + 1) - e(n)*a(n + 2)
!     b(n) = 1.0_wp ! padding
!     a(n) = 1.0_wp ! padding

!     ! Final operations
!     c = 1.0_wp/c
!     a(3:) = -a(3:)*c(3:)
!     b(2:) = -b(2:)*c(2:)
!     d(:nmax - 1) = -d(:nmax - 1)
!     e(:nmax - 2) = -e(:nmax - 2)

!     return
! end subroutine Thomas5_FactorLU2

! ! #######################################################################
! ! Backward substitution step in the Thomas algorith
! ! #######################################################################
! subroutine Thomas5_SolveLU2(nmax, len, a, b, c, d, e, f)
!     use TLab_Constants, only: wp, wi

!     implicit none

!     integer(wi) nmax, len
!     real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
!     real(wp), dimension(len, nmax), intent(INOUT) :: f

!     ! -----------------------------------------------------------------------
!     integer(wi) n
!     real(wp) :: dummy1, dummy2

!     ! #######################################################################
!     ! -----------------------------------------------------------------------
!     ! Solve Ly=f, forward
!     ! -----------------------------------------------------------------------
!     n = nmax - 1
!     f(:, n) = f(:, n) + f(:, n + 1)*d(n)

!     do n = nmax - 2, 1, -1
!         f(:, n) = f(:, n) + f(:, n + 1)*d(n) + f(:, n + 2)*e(n)
!     end do

!     ! -----------------------------------------------------------------------
!     ! Solve Ux=y, backward
!     ! -----------------------------------------------------------------------
!     n = 1
!     f(:, n) = f(:, n)*c(n)

!     n = 2
!     dummy1 = b(n)
!     f(:, n) = c(n)*f(:, n) + dummy1*f(:, n - 1)

!     do n = 3, nmax
!         dummy1 = b(n)
!         dummy2 = a(n)
!         f(:, n) = c(n)*f(:, n) + dummy1*f(:, n - 1) + dummy2*f(:, n - 2)
!     end do

!     return
! end subroutine Thomas5_SolveLU2

#include "tlab_error.h"

! Vectorized Thomas algorithm for tridiagonal systems

module Thomas3
    use TLab_Constants, only: wp, wi, small_wp !, roundoff_wp
    use TLab_Constants, only: efile!, lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: Thomas3_FactorLU, Thomas3_SolveLU
    public :: Thomas3_FactorUL, Thomas3_SolveUL
    public :: Thomas3C_LU, Thomas3C_Solve           ! circulant systems (periodic boundary conditions)
    public :: Thomas3C_SMW_LU, Thomas3C_SMW_Solve

    ! integer(wi) :: n_smw_decay

contains
    ! #######################################################################
    ! #######################################################################
    ! LU factorization stage
    ! L has subdiagonal a and diagonal 1
    ! U has diagonal b and upperdiagonal c
    subroutine Thomas3_FactorLU(nmax, a, b, c)
        integer(wi), intent(IN) :: nmax
        real(wp), intent(INOUT) :: a(nmax), b(nmax), c(nmax)

        ! -----------------------------------------------------------------------
        integer(wi) n

        ! #######################################################################
        do n = 2, nmax
            a(n) = a(n)/b(n - 1)
            b(n) = b(n) - a(n)*c(n - 1)
        end do

        ! Final operations
        a(2:) = -a(2:)
        b = 1.0_wp/b
        c(:nmax - 1) = -c(:nmax - 1)*b(:nmax - 1)

        return
    end subroutine Thomas3_FactorLU

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas3_SolveLU(nmax, len, a, b, c, f)
        integer(wi), intent(IN) :: nmax                         ! dimension of tridiagonal systems
        integer(wi), intent(IN) :: len                          ! number of systems to solve
        real(wp), intent(IN) :: a(nmax), b(nmax), c(nmax)       ! factored LHS
        real(wp), intent(INOUT) :: f(len, nmax)                 ! RHS and solution

        ! -------------------------------------------------------------------
        integer(wi) :: n
        real(wp) :: dummy1, dummy2

        ! ###################################################################
        if (len <= 0) return

        ! Solve Ly=f, forward elimination
        do n = 2, nmax
            dummy1 = a(n)
            f(:, n) = f(:, n) + dummy1*f(:, n - 1)
        end do

        ! Solve Ux=y, backward substitution
        dummy1 = b(nmax)
        f(:, nmax) = f(:, nmax)*dummy1

        do n = nmax - 1, 1, -1
            dummy1 = c(n)
            dummy2 = b(n)
            f(:, n) = dummy2*f(:, n) + dummy1*f(:, n + 1)
        end do

        return
    end subroutine Thomas3_SolveLU

    ! #######################################################################
    ! #######################################################################
    ! UL  factorization stage
    ! L has subdiagonal a and diagonal b
    ! U has diagonal 1 and upperdiagonal c
    subroutine Thomas3_FactorUL(nmax, a, b, c)
        integer(wi), intent(in) :: nmax
        real(wp), intent(inout) :: a(nmax), b(nmax), c(nmax)

        ! -----------------------------------------------------------------------
        integer(wi) n

        ! #######################################################################
        do n = nmax - 1, 1, -1
            c(n) = c(n)/b(n + 1)
            b(n) = b(n) - c(n)*a(n + 1)
        end do

        ! Final operations
        c = -c
        b = 1.0_wp/b
        a = -a*b

        return
    end subroutine Thomas3_FactorUL

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas3_SolveUL(nmax, len, a, b, c, f)
        integer(wi), intent(in) :: nmax                         ! dimension of tridiagonal systems
        integer(wi), intent(in) :: len                          ! number of systems to solve
        real(wp), intent(in) :: a(nmax), b(nmax), c(nmax)       ! factored LHS
        real(wp), intent(inout) :: f(len, nmax)                 ! RHS and solution

        ! -------------------------------------------------------------------
        integer(wi) :: n
        real(wp) :: dummy1, dummy2

        ! ###################################################################
        if (len <= 0) return

        ! Solve Uy=f, backward elimination
        do n = nmax - 1, 1, -1
            dummy1 = c(n)
            f(:, n) = f(:, n) + dummy1*f(:, n + 1)
        end do

        ! Solve Lx=y, forward substitution
        dummy1 = b(1)
        f(:, 1) = f(:, 1)*dummy1

        do n = 2, nmax
            dummy1 = a(n)
            dummy2 = b(n)
            f(:, n) = dummy2*f(:, n) + dummy1*f(:, n - 1)
        end do

        return
    end subroutine Thomas3_SolveUL

    !########################################################################
    !########################################################################
    subroutine Thomas3C_LU(a, b, c, d, e)
        real(wp), intent(inout) :: a(:), b(:), c(:), d(:), e(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nmax
        real(wp) sum

        ! ###################################################################
        nmax = size(a)

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
    end subroutine Thomas3C_LU

    ! #######################################################################
    ! #######################################################################
    subroutine Thomas3C_Solve(a, b, c, d, e, f, wrk)
        real(wp), intent(in) :: a(:), b(:), c(:)
        real(wp), intent(in) :: d(:), e(:)
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution
        real(wp), intent(inout) :: wrk(:)

        ! -------------------------------------------------------------------
        integer(wi) nmax, len, n
        real(wp) :: dummy1, dummy2, dummy3

        ! ###################################################################
        len = size(f, 1)
        nmax = size(f, 2)

        if (len <= 0) return

        ! Forward elimination
        dummy1 = b(1)
        f(:, 1) = f(:, 1)*dummy1

        dummy3 = d(1)
        wrk(:) = dummy3*f(:, 1)

        do n = 2, nmax - 1
            dummy1 = a(n)
            dummy2 = b(n)
            f(:, n) = f(:, n)*dummy2 + dummy1*f(:, n - 1)

            dummy3 = d(n)
            wrk(:) = wrk(:) + dummy3*f(:, n)

        end do

        dummy1 = b(nmax)
        f(:, nmax) = (f(:, nmax) - wrk(:))*dummy1

        ! Backward substitution
        dummy1 = e(nmax - 1)
        f(:, nmax - 1) = dummy1*f(:, nmax) + f(:, nmax - 1)

        do n = nmax - 2, 1, -1
            dummy1 = c(n)
            dummy2 = e(n)
            f(:, n) = f(:, n) + dummy1*f(:, n + 1) + dummy2*f(:, nmax)
        end do

        return
    end subroutine Thomas3C_Solve

    !########################################################################
    !########################################################################
    ! Using Sherman-Morrison-Woodbury formula
    ! Adapted from 10.1016/j.camwa.2011.12.044
    ! Marginally slower because one more call to memory for array f, but clearer

    subroutine Thomas3C_SMW_LU(a, b, c, z)
        real(wp), intent(inout) :: a(:), b(:), c(:)
        real(wp), intent(out) :: z(:)

        ! -------------------------------------------------------------------
        integer(wi) nmax
        real(wp) a1, cn, m
        ! character(len=32) str

        ! ###################################################################
        nmax = size(a)

        a1 = a(1)
        cn = c(nmax)

        ! -------------------------------------------------------------------
        ! Generate matrix A1
        b(1) = b(1) - cn
        b(nmax) = b(nmax) - a1
        call Thomas3_FactorLU(nmax, a, b, c)

        ! -------------------------------------------------------------------
        ! Generate vector z
        z(:) = 0.0_wp
        z(1) = 1.0_wp
        z(nmax) = 1.0_wp
        call Thomas3_SolveLU(nmax, 1, a, b, c, z)

        ! -------------------------------------------------------------------
        ! Calculate normalized coefficients a1 and cn
        m = 1.0_wp + cn*z(1) + a1*z(nmax)
        if (abs(m) < small_wp) then
            call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
            call TLab_Stop(DNS_ERROR_THOMAS)
        end if

        c(nmax) = -cn/m
        a(1) = -a1/m

        ! ! -------------------------------------------------------------------
        ! ! Calculate decay index
        ! do n_smw_decay = 2, nmax
        !     if (abs(z(n_smw_decay)/z(1)) < roundoff_wp) exit
        !     ! print *, abs(z(n_smw_decay)/z(1)
        ! end do
        ! write (str, *) n_smw_decay
        ! call TLab_Write_ASCII(lfile, 'Decay to round-off in SMW algorithm in '//trim(adjustl(str))//' indexes.')

        return
    end subroutine Thomas3C_SMW_LU

    !########################################################################
    !########################################################################
    subroutine Thomas3C_SMW_Solve(a, b, c, z, f, wrk)
        real(wp), intent(in) :: a(:), b(:), c(:)
        real(wp), intent(in) :: z(:)
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution
        real(wp), intent(inout) :: wrk(:)

        integer(wi) nmax, len, n

        len = size(f, 1)
        nmax = size(f, 2)

        call Thomas3_SolveLU(nmax, len, a, b, c, f)

        wrk(:) = c(nmax)*f(:, 1) + a(1)*f(:, nmax)
        do n = 1, nmax
            f(:, n) = f(:, n) + wrk(:)*z(n)
        end do

        ! This would save time in the serial case, but we are interested in the parallel case
        ! n_smw_decay = 64
        ! do n = 1, min(nmax/2, n_smw_decay)
        !     f(:, n) = f(:, n) + wrk(:)*z(n)
        !     f(:, nmax - n + 1) = f(:, nmax - n + 1) + wrk(:)*z(nmax - n + 1)
        ! end do

        return
    end subroutine Thomas3C_SMW_Solve

end module Thomas3

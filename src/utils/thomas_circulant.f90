#include "tlab_error.h"

! Using Sherman-Morrison-Woodbury formula
! Adapted from 10.1016/j.camwa.2011.12.044
! Marginally slower because one more call to memory for array f, but clearer

! Written in terms of generic splitting algorithm

module Thomas_Circulant
    use TLab_Constants, only: wp, wi
    use Thomas_Split, only: thomas_split_dt
    implicit none
    private

    public :: thomas_circulant_dt

    ! -----------------------------------------------------------------------
    type, extends(thomas_split_dt) :: thomas_circulant_dt
    contains
        procedure :: initialize => thomas_initialize_dt
        ! the next is needed because pentadiagonal case not yet in generic form
        procedure :: reduce => thomas_reduce_dt
    end type

contains
    ! #######################################################################
    ! #######################################################################
    subroutine thomas_initialize_dt(self, lhs)
        class(thomas_circulant_dt), intent(out) :: self
        real(wp), intent(in) :: lhs(:, :)

        ! split at the end of the array
        call self%initialize_split(lhs, index=size(lhs, 1))

        ! this case is not yet in split format, only here for circulant cases
        select case (size(lhs, 2))
        case (5)
            call ThomasCirculant_5_Initialize(self%L, self%U, self%z)
        end select

        return
    end subroutine

    subroutine thomas_reduce_dt(self, f)
        use TLab_Arrays, only: wrk2d
        class(thomas_circulant_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%reduce_split(f)

        ! this case is not yet in split format, only here for circulant cases
        select case (size(self%L, 2))
        case (2)
            call ThomasCirculant_5_Reduce(self%L, &
                                          self%U, &
                                          self%z, &
                                          f)!, wrk2d)
        end select

        return
    end subroutine

    ! #######################################################################
    ! #######################################################################
#define a(i) L(i,1)
#define b(i) L(i,2)
#define c(i) U(i,1)
#define d(i) U(i,2)
#define e(i) U(i,3)

#define z1(i) z_mem(1,i)
#define z2(i) z_mem(2,i)

    subroutine ThomasCirculant_5_Initialize(L, U, z_mem)
        use TLab_Constants, only: small_wp
        use Thomas, only: Thomas5_FactorLU_InPlace, Thomas5_SolveL, Thomas5_SolveU
        use TLab_Constants, only: efile
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        real(wp), intent(inout) :: L(:, :), U(:, :)
        real(wp), intent(inout) :: z_mem(2, size(L, 1))

        ! -----------------------------------------------------------------------
        integer nmax
        real(wp) :: m1, m2, m3, m4

        ! #######################################################################
        nmax = size(L, 1)

        ! -------------------------------------------------------------------
        ! Generate matrix A1
        b(2) = b(2) - d(nmax)
        c(1) = c(1) - e(nmax - 1)
        c(2) = c(2) - e(nmax)

        c(nmax - 1) = c(nmax - 1) - a(1)
        c(nmax) = c(nmax) - a(2)
        d(nmax - 1) = d(nmax - 1) - b(1)

        call Thomas5_FactorLU_InPlace(L, U)

        ! -------------------------------------------------------------------
        ! Generate vector z
        z1(:) = 0.0_wp ! u1
        z1(1) = 1.0_wp
        z1(nmax - 1) = 1.0_wp

        z2(:) = 0.0_wp ! u2
        z2(2) = 1.0_wp
        z2(nmax) = 1.0_wp

        call Thomas5_SolveL(L, z_mem)
        call Thomas5_SolveU(U, z_mem)

        ! Compute entries of matrix M[2x2] once
        m1 = e(nmax - 1)*z1(1) + a(1)*z1(nmax - 1) + b(1)*z1(nmax) + 1.0_wp
        m2 = e(nmax - 1)*z2(1) + a(1)*z2(nmax - 1) + b(1)*z2(nmax)
        m3 = d(nmax)*z1(1) + e(nmax)*z1(2) + a(2)*z1(nmax)
        m4 = d(nmax)*z2(1) + e(nmax)*z2(2) + a(2)*z2(nmax) + 1.0_wp
        ! Check if M is invertible (eq. 2.9)
        if ((m1*m4 - m2*m3) < small_wp) then
            call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
            call TLab_Stop(DNS_ERROR_THOMAS)
        end if

        return
    end subroutine ThomasCirculant_5_Initialize

    ! #######################################################################
    ! #######################################################################
    subroutine ThomasCirculant_5_Reduce(L, U, z_mem, f)
        real(wp), intent(in) :: L(:, :), U(:, :)
        real(wp), intent(in) :: z_mem(2, size(L, 1))
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution

        ! -----------------------------------------------------------------------
        integer(wi) :: nmax, n, ll
        real(wp) :: m1, m2, m3, m4
        real(wp) :: di, d11, d12, d13, d14, d21, d22, d23, d24
        real(wp) :: dummy1, dummy2

        ! #######################################################################
        if (size(f, 1) <= 0) return

        nmax = size(f, 2)

        ! Compute entries of matrix m[2x2]
        m1 = e(nmax - 1)*z1(1) + a(1)*z1(nmax - 1) + b(1)*z1(nmax) + 1.0_wp
        m2 = e(nmax - 1)*z2(1) + a(1)*z2(nmax - 1) + b(1)*z2(nmax)
        m3 = d(nmax)*z1(1) + e(nmax)*z1(2) + a(2)*z1(nmax)
        m4 = d(nmax)*z2(1) + e(nmax)*z2(2) + a(2)*z2(nmax) + 1.0_wp

        ! Compute coefficients
        di = 1.0_wp/(m1*m4 - m2*m3)
        d11 = di*(m4*e(nmax - 1) - m2*d(nmax))
        d12 = di*(m4*b(1) - m2*a(2))
        d13 = di*m4*a(1)
        d14 = di*m2*e(nmax)
        d21 = di*(m1*d(nmax) - m3*e(nmax - 1))
        d22 = di*(m1*a(2) - m3*b(1))
        d23 = di*m3*a(1)
        d24 = di*m1*e(nmax)

        ! Solve
        do n = 3, nmax - 2, 1 ! Main loop
            do ll = 1, size(f, 1)
                dummy1 = d11*f(ll, 1) + d12*f(ll, nmax) + d13*f(ll, nmax - 1) - d14*f(ll, 2)
                dummy2 = d21*f(ll, 1) + d22*f(ll, nmax) - d23*f(ll, nmax - 1) + d24*f(ll, 2)
                !
                f(ll, n) = f(ll, n) - dummy1*z1(n) - dummy2*z2(n)
            end do
        end do
        !
        do ll = 1, size(f, 1)    ! Boundaries
            dummy1 = d11*f(ll, 1) + d12*f(ll, nmax) + d13*f(ll, nmax - 1) - d14*f(ll, 2)
            dummy2 = d21*f(ll, 1) + d22*f(ll, nmax) - d23*f(ll, nmax - 1) + d24*f(ll, 2)
            do n = 1, 2
                f(ll, n) = f(ll, n) - dummy1*z1(n) - dummy2*z2(n)
            end do
            do n = nmax - 1, nmax
                f(ll, n) = f(ll, n) - dummy1*z1(n) - dummy2*z2(n)
            end do
        end do

        return
    end subroutine ThomasCirculant_5_Reduce

end module Thomas_Circulant

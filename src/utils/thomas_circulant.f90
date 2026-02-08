#include "tlab_error.h"

! Using Sherman-Morrison-Woodbury formula
! Adapted from 10.1016/j.camwa.2011.12.044
! Marginally slower because one more call to memory for array f, but clearer

module Thomas_Circulant
    use TLab_Constants, only: wp, wi, small_wp !, roundoff_wp
    use TLab_Constants, only: efile!, lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas
    implicit none
    private

    public :: ThomasCirculant_3_Initialize
    public :: ThomasCirculant_3_Reduce

    public :: ThomasCirculant_5_Initialize
    public :: ThomasCirculant_5_Reduce

    public :: thomas_circulant_dt

    ! -----------------------------------------------------------------------
    type, extends(thomas_dt) :: thomas_circulant_dt
        real(wp), allocatable :: z(:, :)
    contains
        procedure :: initialize => thomas_initialize_dt
        procedure :: solveL => thomas_solveL_dt
        procedure :: solveU => thomas_solveU_dt
        procedure :: reduce => thomas_reduce_dt
    end type

contains
    ! #######################################################################
    ! #######################################################################
    subroutine thomas_initialize_dt(self, lhs)
        class(thomas_circulant_dt), intent(out) :: self
        real(wp), intent(in) :: lhs(:, :)

        integer ndl

        ndl = size(lhs, 2)
        allocate (self%L, source=lhs(:, 1:ndl/2))
        allocate (self%U, source=lhs(:, ndl/2 + 1:ndl))
        allocate (self%z(ndl/2, size(lhs, 1)))
        select case (ndl)
        case (3)
            call ThomasCirculant_3_Initialize(self%L, self%U, self%z)
        case (5)
            call ThomasCirculant_5_Initialize(self%L, self%U, self%z)
        end select

        select case (ndl)
        case (3)
            self%ptr_solveL => Thomas3_SolveL
            self%ptr_solveU => Thomas3_SolveU
        case (5)
            self%ptr_solveL => Thomas5_SolveL
            self%ptr_solveU => Thomas5_SolveU
        case (7)
            self%ptr_solveL => Thomas7_SolveL
            self%ptr_solveU => Thomas7_SolveU
        end select

        return
    end subroutine

    subroutine thomas_solveL_dt(self, f)
        class(thomas_circulant_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%ptr_solveL(self%L, f)

        return
    end subroutine

    ! subroutine thomas_solveU_dt(self, f)
    !     use TLab_Arrays, only: wrk2d
    !     class(thomas_circulant_dt), intent(in) :: self
    !     real(wp), intent(inout) :: f(:, :)

    !     call self%ptr_solveU(self%U, f)

    !     select case (size(self%L, 2))
    !     case (1)
    !         call ThomasCirculant_3_Reduce(self%L, &
    !                                       self%U, &
    !                                       self%z(1, :), &
    !                                       f, wrk2d(:, 1))
    !     case (2)
    !         call ThomasCirculant_5_Reduce(self%L, &
    !                                       self%U, &
    !                                       self%z, &
    !                                       f)!, wrk2d)
    !     end select

    !     return
    ! end subroutine

    subroutine thomas_solveU_dt(self, f)
        class(thomas_circulant_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%ptr_solveU(self%U, f)

        return
    end subroutine

    subroutine thomas_reduce_dt(self, f)
        use TLab_Arrays, only: wrk2d
        class(thomas_circulant_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        select case (size(self%L, 2))
        case (1)
            call ThomasCirculant_3_Reduce(self%L, &
                                          self%U, &
                                          self%z(1, :), &
                                          f, wrk2d(:, 1))
        case (2)
            call ThomasCirculant_5_Reduce(self%L, &
                                          self%U, &
                                          self%z, &
                                          f)!, wrk2d)
        end select

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine ThomasCirculant_3_Initialize(L, U, z_mem)
        real(wp), intent(inout) :: L(:, :), U(:, :)
        real(wp), intent(inout) :: z_mem(1, size(L, 1))

        ! -----------------------------------------------------------------------
        integer nmax
        real(wp) a1, cn, m

        ! #######################################################################
        nmax = size(L, 1)

        a1 = L(1, 1)
        cn = U(nmax, 2)

        ! -------------------------------------------------------------------
        ! Generate matrix A1
        U(1, 1) = U(1, 1) - cn

        U(nmax, 1) = U(nmax, 1) - a1

        call Thomas3_FactorLU_InPlace(L, U)

        ! -------------------------------------------------------------------
        ! Generate vector z
#define z(i) z_mem(1,i)

        z(:) = 0.0_wp
        z(1) = 1.0_wp
        z(nmax) = 1.0_wp

        call Thomas3_SolveL(L, z_mem)
        call Thomas3_SolveU(U, z_mem)

        ! -------------------------------------------------------------------
        ! Calculate normalized coefficients a1 and cn
        m = 1.0_wp + cn*z(1) + a1*z(nmax)
        if (abs(m) < small_wp) then
            call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
            call TLab_Stop(DNS_ERROR_THOMAS)
        end if

        U(nmax, 2) = -cn/m
        L(1, 1) = -a1/m

        ! ! -------------------------------------------------------------------
        ! ! Calculate decay index
        ! do n_smw_decay = 2, nmax
        !     if (abs(z(n_smw_decay)/z(1)) < roundoff_wp) exit
        !     ! print *, abs(z(n_smw_decay)/z(1)
        ! end do
        ! write (str, *) n_smw_decay
        ! call TLab_Write_ASCII(lfile, 'Decay to round-off in SMW algorithm in '//trim(adjustl(str))//' indexes.')

#undef z

        return
    end subroutine ThomasCirculant_3_Initialize

    !########################################################################
    !########################################################################
    subroutine ThomasCirculant_3_Reduce(L, U, z, f, wrk)
        real(wp), intent(in) :: L(:, :), U(:, :), z(:)
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution
        real(wp), intent(inout) :: wrk(size(f, 1))

        ! -------------------------------------------------------------------
        integer(wi) nmax, n

        ! ###################################################################
        if (size(f, 1) <= 0) return

        nmax = size(f, 2)
        wrk(:) = U(nmax, 2)*f(:, 1) + L(1, 1)*f(:, nmax)
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
    end subroutine ThomasCirculant_3_Reduce

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

#include "tlab_error.h"

! Splitting algorithm of linear solver

module Thomas_Split
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: efile, lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas, only: thomas_base_dt
    use Thomas, only: Thomas_FactorLU_InPlace, Thomas3_SolveL, Thomas3_SolveU
    implicit none
    private

    public :: thomas_split_dt
    public :: Thomas_3_Split_InPlace
    public :: Thomas_3_Split_Reduce
    public :: decay_index

    ! -----------------------------------------------------------------------
    type, extends(thomas_base_dt) :: thomas_split_dt
        real(wp), allocatable :: z(:, :)
        real(wp) :: alpha(2) = [0.0_wp, 0.0_wp]
    contains
        procedure :: initialize_split => thomas_initialize_dt
        procedure :: reduce_split => thomas_reduce_dt
        procedure :: solve => thomas_solve_dt
    end type

contains
    ! #######################################################################
    ! #######################################################################
    subroutine thomas_initialize_dt(self, lhs, index)
        class(thomas_split_dt), intent(out) :: self
        real(wp), intent(in) :: lhs(:, :)
        integer, intent(in) :: index

        integer ndl

        call self%initialize_base(lhs)

        ! corrections for splitting case
        ndl = size(lhs, 2)
        self%L(:, :) = lhs(:, 1:ndl/2)
        self%U(:, :) = lhs(:, ndl/2 + 1:ndl)
        allocate (self%z(ndl/2, size(lhs, 1)))
        select case (ndl)
        case (3)
            call Thomas_3_Split_InPlace(self%L, self%U, self%z, index)
            ! case (5)
            !     call ThomasCirculant_5_Initialize(self%L, self%U, self%z)
        end select

        return
    end subroutine

    subroutine thomas_reduce_dt(self, f)
        use TLab_Arrays, only: wrk2d
        class(thomas_split_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        select case (size(self%L, 2))
        case (1)
            call Thomas_3_Split_Reduce(self%L, &
                                       self%U, &
                                       self%z(1, :), &
                                       f, wrk2d(:, 1))
            ! case (2)
            !     call ThomasCirculant_5_Reduce(self%L, &
            !                                   self%U, &
            !                                   self%z, &
            !                                   f)!, wrk2d)
        end select

        return
    end subroutine

    subroutine thomas_solve_dt(self, f)
        class(thomas_split_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%solveL(f)
        call self%solveU(f)
        call self%reduce_split(f)

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Thomas_3_Split_InPlace(L, U, z, index)
        use TLab_Constants, only: small_wp
        real(wp), intent(inout) :: L(:, :), U(:, :)
        real(wp), intent(out) :: z(1, size(L, 1))
        integer, intent(in) :: index

        ! -----------------------------------------------------------------------
        integer nmax, p, p_plus_1
        real(wp) alpha(2), delta

        ! #######################################################################
        nmax = size(L, 1)

        p = mod(index - 1, nmax) + 1                ! in circulant cases, this n
        p_plus_1 = mod(p + 1 - 1, nmax) + 1         ! in circulant cases, this 1

#define a(i) L(i,1)
#define b(i) U(i,1)
#define c(i) U(i,2)

        ! Start definition of alpha; in circulant cases, this is a1 and cn
        alpha(1) = a(p_plus_1)
        alpha(2) = c(p)

        ! Generate matrix A1
        b(p) = b(p) - a(p_plus_1)
        a(p_plus_1) = 0.0_wp
        b(p_plus_1) = b(p_plus_1) - c(p)
        c(p) = 0.0_wp

        ! call Thomas3_FactorLU_InPlace(L, U)
        call Thomas_FactorLU_InPlace(L, U)

        ! Generate vector z1
        z(1, :) = 0.0_wp
        z(1, p) = 1.0_wp
        z(1, p_plus_1) = 1.0_wp

        call Thomas3_SolveL(L, z)
        call Thomas3_SolveU(U, z)

        ! Calculate normalized alpha coefficients
        delta = 1.0_wp + alpha(1)*z(1, p) + alpha(2)*z(1, p_plus_1)
        if (abs(delta) < small_wp) then
            call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
            call TLab_Stop(DNS_ERROR_THOMAS)
        end if

        a(p_plus_1) = -alpha(1)/delta
        c(p) = -alpha(2)/delta

        ! -------------------------------------------------------------------
        ! Calculate decay index
        ! call decay_index(z(1, p_plus_1:))!, n_smw_decay)

#undef a
#undef b
#undef c

        return
    end subroutine Thomas_3_Split_InPlace

    !########################################################################
    !########################################################################
    subroutine Thomas_3_Split_Reduce(L, U, z, f, wrk)
        real(wp), intent(in) :: L(:, :), U(:, :), z(:)
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution
        real(wp), intent(inout) :: wrk(size(f, 1))

        ! -------------------------------------------------------------------
        integer(wi) nmax, n

        ! ###################################################################
        if (size(f, 1) <= 0) return

#define cn U(nmax, 2)
#define a1 L(1, 1)
        nmax = size(f, 2)
        wrk(:) = cn*f(:, 1) + a1*f(:, nmax)
        do n = 1, nmax
            f(:, n) = f(:, n) + wrk(:)*z(n)
        end do

#undef cn
#undef a1

        ! This would save time in the serial case, but we are interested in the parallel case
        ! n_smw_decay = 64
        ! do n = 1, min(nmax/2, n_smw_decay)
        !     f(:, n) = f(:, n) + wrk(:)*z(n)
        !     f(:, nmax - n + 1) = f(:, nmax - n + 1) + wrk(:)*z(nmax - n + 1)
        ! end do

        return
    end subroutine Thomas_3_Split_Reduce

    !########################################################################
    !########################################################################
    subroutine decay_index(z, index)
        use TLab_Constants, only: roundoff_wp
        real(wp), intent(in) :: z(:)
        integer, intent(out), optional :: index

        integer n
        character(len=32) str

        do n = 2, size(z)
            if (abs(z(n)/z(1)) < roundoff_wp) exit
            ! print *, abs(z(n)/z(1)
        end do
        write (str, *) n
        call TLab_Write_ASCII(lfile, 'Decay to round-off splitting algorithm in '//trim(adjustl(str))//' indexes.')

        if (present(index)) index = n

        return
    end subroutine

end module Thomas_Split

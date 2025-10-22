#include "tlab_error.h"

module Thomas_Circulant
    use TLab_Constants, only: wp, wi, small_wp !, roundoff_wp
    use TLab_Constants, only: efile!, lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas
    implicit none
    private

    public :: Thomas3_C_SMW_Initialize
    public :: Thomas3_C_SMW_Solve

contains
    !########################################################################
    !########################################################################
    ! Using Sherman-Morrison-Woodbury formula
    ! Adapted from 10.1016/j.camwa.2011.12.044
    ! Marginally slower because one more call to memory for array f, but clearer

    subroutine Thomas3_C_SMW_Initialize(L, U, z_mem)
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
#define z(i) z_mem(i,1)

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
    end subroutine Thomas3_C_SMW_Initialize

    !########################################################################
    !########################################################################
    subroutine Thomas3_C_SMW_Solve(L, U, z, f, wrk)
        real(wp), intent(in) :: L(:, :), U(:, :), z(:)
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution
        real(wp), intent(inout) :: wrk(:)

        ! -------------------------------------------------------------------
        integer(wi) nmax, n

        ! ###################################################################
        if (size(f, 1) <= 0) return

        call Thomas3_SolveL(L, f)
        call Thomas3_SolveU(U, f)

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
    end subroutine Thomas3_C_SMW_Solve

end module Thomas_Circulant

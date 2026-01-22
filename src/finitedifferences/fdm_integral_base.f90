module FDM_Integral_Base
    use TLab_Constants, only: wp
    implicit none
    ! everything is public, so no private statement

    type, public :: fdm_integral_dt
        sequence
        integer mode_fdm                                    ! original finite-difference method; only informative
        real(wp) :: lambda                                  ! constant of the equation
        integer :: bc                                       ! type of boundary condition, [ BCS_MIN, BCS_MAX ]
        real(wp), allocatable :: lhs(:, :)                  ! Often overwritten to LU decomposition.
        real(wp), allocatable :: rhs(:, :)
        real(wp), allocatable :: rhs_b1(:, :), rhs_t1(:, :) ! boundary conditions
        procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul_thomas => null()
        procedure(thomas_ice), pointer, nopass :: thomasL => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
    end type fdm_integral_dt
    ! This type is used in elliptic operators for different eigenvalues. This can lead to fragmented memory.
    ! One could use pointers instead of allocatable for lhs and rhs, and point the pointers to the
    ! corresponding memory space.

    ! -----------------------------------------------------------------------
    abstract interface
        subroutine thomas_ice(A, f)
            import wp
            real(wp), intent(in) :: A(:, :)
            real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        end subroutine
    end interface

    abstract interface
        subroutine matmul_ice(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
            import wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine
    end interface

    abstract interface
        subroutine matmul_thomas_ice(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
            import wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(in) :: L(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine
    end interface

end module FDM_Integral_Base

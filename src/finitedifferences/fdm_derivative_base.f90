module FDM_Derivative_Base
    use TLab_Constants, only: wp, wi
    use Thomas, only: thomas_dt
    use Thomas_Circulant, only: thomas_circulant_dt
    implicit none
    private

    public :: der_dt
    public :: der_periodic
    public :: der_biased

    ! public :: matmul_halo_ice
    public :: matmul_halo_thomas_ice
    public :: matmul_halo_thomas_combined_ice
    ! public :: matmul_ice
    public :: matmul_thomas_ice

    ! -----------------------------------------------------------------------
    type, abstract :: der_dt
        integer type                                ! finite-difference method
        real(wp), allocatable :: lhs(:, :)          ! A diagonals of system A u' = B u
        real(wp), allocatable :: rhs(:, :)          ! B diagonals of system A u' = B u
    contains
        procedure(initialize_ice), deferred :: initialize
        procedure(compute_ice), deferred :: compute
    end type
    abstract interface
        subroutine initialize_ice(self, x, fdm_type)
            import der_dt, wp
            class(der_dt), intent(out) :: self
            real(wp), intent(in) :: x(:)
            integer, intent(in) :: fdm_type
        end subroutine
        subroutine compute_ice(self, nlines, u, result)
            import der_dt, wp, wi
            class(der_dt), intent(in) :: self
            integer(wi), intent(in) :: nlines
            real(wp), intent(in) :: u(nlines, size(self%lhs, 1))
            real(wp), intent(out) :: result(nlines, size(self%lhs, 1))
        end subroutine
    end interface

    type, extends(der_dt), abstract :: der_periodic
        ! procedure(matmul_halo_ice), pointer, nopass :: matmul => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul => null()
        type(thomas_circulant_dt) :: thomas
    contains
        procedure :: compute => der_periodic_compute
    end type

    type, extends(der_dt), abstract :: der_biased
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul => null()
        type(thomas_dt) :: thomas
    contains
        procedure :: compute => der_biased_compute
    end type

    ! -----------------------------------------------------------------------
    abstract interface
        ! subroutine matmul_halo_ice(rhs, u, u_halo_m, u_halo_p, f)
        !     use TLab_Constants, only: wp
        !     real(wp), intent(in) :: rhs(:)              ! diagonals of B
        !     real(wp), intent(in) :: u(:, :)             ! vector u
        !     real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        !     real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        !     real(wp), intent(out) :: f(:, :)            ! vector f = B u
        ! end subroutine

        subroutine matmul_halo_thomas_ice(rhs, u, u_halo_m, u_halo_p, f, L)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:)              ! diagonals of B
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B u
            real(wp), intent(in) :: L(:, :)
        end subroutine

        subroutine matmul_halo_thomas_combined_ice(rhs1, rhs2, u, u_halo_m, u_halo_p, f, L1, g, L2)
            import wp
            real(wp), intent(in) :: rhs1(:)             ! diagonals of B1
            real(wp), intent(in) :: rhs2(:)             ! diagonals of B2
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B1 u
            real(wp), intent(in) :: L1(:, :)
            real(wp), intent(out) :: g(:, :)            ! vector g = B2 u
            real(wp), intent(in) :: L2(:, :)
        end subroutine

        ! subroutine matmul_ice(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        !     use TLab_Constants, only: wp
        !     real(wp), intent(in) :: rhs(:, :)
        !     real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        !     real(wp), intent(in) :: u(:, :)
        !     real(wp), intent(out) :: f(:, :)
        !     real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        ! end subroutine

        subroutine matmul_thomas_ice(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(in) :: L(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine

        ! subroutine thomas_ice(A, f)
        !     use TLab_Constants, only: wp
        !     real(wp), intent(in) :: A(:, :)
        !     real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        ! end subroutine

    end interface

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der_periodic_compute(self, nlines, u, result)
        class(der_periodic), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%lhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%lhs, 1))

        integer nx, ndr

        ! ###################################################################
        nx = size(self%lhs, 1)
        ndr = size(self%rhs, 2)

        ! Calculate RHS in system of equations A u' = B u
        ! call self%matmul(rhs=self%rhs(1, 1:ndr), &
        !                    u=u, &
        !                    u_halo_m=u(:, nx - ndr/2 + 1:nx), &
        !                    u_halo_p=u(:, 1:ndr/2), &
        !                    f=result)
        call self%matmul(rhs=self%rhs(1, 1:ndr), &
                         u=u, &
                         u_halo_m=u(:, nx - ndr/2 + 1:nx), &
                         u_halo_p=u(:, 1:ndr/2), &
                         f=result, &
                         L=self%thomas%L)

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomas%solveL(result)
        call self%thomas%solveU(result)
        call self%thomas%reduce(result)

        return
    end subroutine der_periodic_compute

    ! ###################################################################
    ! ###################################################################
    subroutine der_biased_compute(self, nlines, u, result)
        class(der_biased), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%rhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%rhs, 1))

        integer nx, ndr

        ! ###################################################################
        nx = size(self%rhs, 1)
        ndr = size(self%rhs, 2)

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs(1:ndr/2, 1:ndr), &
                         rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
                         u=u, &
                         f=result, &
                         L=self%thomas%L)

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomas%solveL(result)
        call self%thomas%solveU(result)

        return
    end subroutine der_biased_compute

end module FDM_Derivative_Base

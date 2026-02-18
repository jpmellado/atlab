! Coefficients in explicit formula to calculate boundary value in Neumann conditions
module FDM_derivative_Neumann
    use TLab_Constants, only: wp, wi, roundoff_wp
    use TLab_Constants, only: BCS_ND, BCS_DN, BCS_NN
    use FDM_Derivative_1order, only: der1_biased_extended
    implicit none
    private

    public :: FDM_Der1_Neumann_Initialize           ! General, full version
    public :: FDM_Der1_NeumannMin_Initialize        ! Truncated versions
    public :: FDM_Der1_NeumannMax_Initialize

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_Neumann_Initialize(ibc, der, c_b, c_t, u, z)
        integer, intent(in) :: ibc                          ! Boundary condition [BCS_ND, BCS_DN, BCS_NN]
        type(der1_biased_extended), intent(in) :: der
        real(wp), intent(out) :: c_b(size(der%lhs, 1))      ! coefficients for bottom value and top value
        real(wp), intent(out) :: c_t(size(der%lhs, 1))
        real(wp), intent(inout) :: u(1, size(der%lhs, 1))   ! Working arrays
        real(wp), intent(inout) :: z(1, size(der%lhs, 1))

        ! -------------------------------------------------------------------
        integer n, nx
        real(wp) bcs_b(1), bcs_t(1)

        ! ###################################################################
        nx = size(der%lhs, 1)
        do n = 1, nx
            u(1, :) = 0.0_wp                ! Create delta-function forcing term
            u(1, n) = 1.0_wp

            select case (ibc)
            case (BCS_ND)
                bcs_b(1:1) = u(1, 1)
                call der%bcsND%compute(1, u, z, bcs_b)
                c_b(n) = bcs_b(1)

            case (BCS_DN)
                bcs_t(1:1) = u(1, nx)
                call der%bcsDN%compute(1, u, z, bcs_t)
                c_t(n) = bcs_t(1)

            case (BCS_NN)
                bcs_b(1:1) = u(1, 1)
                bcs_t(1:1) = u(1, nx)
                call der%bcsNN%compute(1, u, z, bcs_b(:), bcs_t)
                c_b(n) = bcs_b(1)
                c_t(n) = bcs_t(1)

            end select

        end do

        return
    end subroutine FDM_Der1_Neumann_Initialize

    ! ###################################################################
    ! ###################################################################
    ! Truncated version for BCS_ND
    subroutine FDM_Der1_NeumannMin_Initialize(der, c_b, u, z, n_bcs)
        type(der1_biased_extended), intent(in) :: der
        real(wp), intent(out) :: c_b(size(der%lhs, 1))      ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, size(der%lhs, 1))   ! Working arrays
        real(wp), intent(inout) :: z(1, size(der%lhs, 1))
        integer, intent(out) :: n_bcs                       ! Index of truncation

        ! -------------------------------------------------------------------
        integer n, nx
        real(wp) bcs_b(1)

        ! ###################################################################
        nx = size(der%lhs, 1)

        c_b(:) = 0.0_wp
        do n = 1, nx
            u(1, :) = 0.0_wp                ! Create delta-function forcing term
            u(1, n) = 1.0_wp

            bcs_b(1) = u(1, 1)
            call der%bcsND%compute(1, u, z, bcs_b)
            c_b(n) = bcs_b(1)

            if (abs(c_b(n)/c_b(1)) < roundoff_wp) exit

        end do
        n_bcs = min(n, nx)

        return
    end subroutine FDM_Der1_NeumannMin_Initialize

    ! ###################################################################
    ! ###################################################################
    ! Truncated version for BCS_DN
    subroutine FDM_Der1_NeumannMax_Initialize(der, c_t, u, z, n_bcs)
        type(der1_biased_extended), intent(in) :: der
        real(wp), intent(out) :: c_t(size(der%lhs, 1))      ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, size(der%lhs, 1))   ! Working arrays
        real(wp), intent(inout) :: z(1, size(der%lhs, 1))
        integer, intent(out) :: n_bcs                       ! Index of truncation

        ! -------------------------------------------------------------------
        integer n, nx
        real(wp) bcs_t(1)

        ! ###################################################################
        nx = size(der%lhs, 1)

        c_t(:) = 0.0_wp
        do n = nx, 1, -1
            u(1, :) = 0.0_wp                ! Create delta-function forcing term
            u(1, n) = 1.0_wp

            bcs_t(1) = u(1, nx)
            call der%bcsDN%compute(1, u, z, bcs_t)
            c_t(n) = bcs_t(1)
            if (abs(c_t(n)/c_t(nx)) < roundoff_wp) exit

        end do
        n_bcs = min(nx - n + 1, nx)

        return
    end subroutine FDM_Der1_NeumannMax_Initialize

end module FDM_derivative_Neumann

! Coefficients in explicit formula to calculate boundary value in Neumann conditions
module FDM_derivative_Neumann
    use TLab_Constants, only: wp, wi, roundoff_wp
    use TLab_Constants, only: BCS_ND, BCS_DN, BCS_NN
    use FDM_Derivative
    use Thomas
    ! use Thomas3
    ! use Thomas5
    implicit none
    private

    public :: FDM_Der1_Neumann_Initialize           ! General, full version
    public :: FDM_Der1_NeumannMin_Initialize        ! Truncated versions
    public :: FDM_Der1_NeumannMax_Initialize

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_Neumann_Initialize(ibc, g, c_b, c_t, u, z)
        integer, intent(in) :: ibc                          ! Boundary condition [BCS_ND, BCS_DN, BCS_NN]
        type(fdm_derivative_dt), intent(in) :: g
        real(wp), intent(out) :: c_b(g%size), c_t(g%size)   ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, g%size)             ! Working arrays
        real(wp), intent(inout) :: z(1, g%size)

        ! -------------------------------------------------------------------
        integer ndl, idl, ndr, idr, ic
        integer(wi) nmin, nmax, nsize, n, ip
        real(wp) bcs_hb(1), bcs_ht(1)

        ! ###################################################################
        ndl = g%nb_diag(1)                  ! number of diagonals in lhs
        idl = ndl/2 + 1
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1

        do n = 1, g%size
            u(1, :) = 0.0_wp                ! Create delta-function forcing term
            u(1, n) = 1.0_wp

            nmin = 1; nmax = g%size
            if (any([BCS_ND, BCS_NN] == ibc)) then
                z(1, 1) = u(1, 1)
                nmin = nmin + 1
            end if
            if (any([BCS_DN, BCS_NN] == ibc)) then
                z(1, g%size) = u(1, g%size)
                nmax = nmax - 1
            end if
            nsize = nmax - nmin + 1

            ! -------------------------------------------------------------------
            ! Calculate RHS in system of equations A u' = B u
            call g%matmul(g%rhs, u, z, ibc, g%rhs_b, g%rhs_t, bcs_hb, bcs_ht)

            ! -------------------------------------------------------------------
            ! Solve for u' in system of equations A u' = B u
            ip = ibc*5
            select case (ndl)
            case (3)
                call Thomas3_SolveL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
                call Thomas3_SolveU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))
                ! call Thomas3_SolveLU(nsize, 1, &
                !                      g%lu(nmin:nmax, ip + 1), &
                !                      g%lu(nmin:nmax, ip + 2), &
                !                      g%lu(nmin:nmax, ip + 3), &
                !                      z(:, nmin:nmax))
            case (5)
                call Thomas5_SolveL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
                call Thomas5_SolveU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))
                ! call Thomas5_SolveLU(nsize, 1, &
                !                      g%lu(nmin:nmax, ip + 1), &
                !                      g%lu(nmin:nmax, ip + 2), &
                !                      g%lu(nmin:nmax, ip + 3), &
                !                      g%lu(nmin:nmax, ip + 4), &
                !                      g%lu(nmin:nmax, ip + 5), &
                !                      z(:, nmin:nmax))
            end select

            if (any([BCS_ND, BCS_NN] == ibc)) then
                do ic = 1, idl - 1
                    bcs_hb(1) = bcs_hb(1) + g%lu(1, ip + idl + ic)*z(1, 1 + ic)
                end do
                c_b(n) = bcs_hb(1)/g%rhs(1, idr)
            end if
            if (any([BCS_DN, BCS_NN] == ibc)) then
                do ic = 1, idl - 1
                    bcs_ht(1) = bcs_ht(1) + g%lu(g%size, ip + idl - ic)*z(1, g%size - ic)
                end do
                c_t(n) = bcs_ht(1)/g%rhs(g%size, idr)
            end if

        end do

        return
    end subroutine FDM_Der1_Neumann_Initialize

! ###################################################################
! ###################################################################
! Truncated version for BCS_ND
    subroutine FDM_Der1_NeumannMin_Initialize(g, c_b, u, z, n_bcs)
        type(fdm_derivative_dt), intent(in) :: g
        real(wp), intent(out) :: c_b(g%size)                ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, g%size)             ! Working arrays
        real(wp), intent(inout) :: z(1, g%size)
        integer, intent(out) :: n_bcs                       ! Index of truncation

        ! -------------------------------------------------------------------
        integer ndl, idl, ndr, idr, ic
        integer ibc
        integer(wi) nmin, nmax, nsize, n, ip
        real(wp) bcs_hb(1), bcs_ht(1)

        ! ###################################################################
        ibc = BCS_ND

        c_b(:) = 0.0_wp

        nmin = 2; nmax = g%size
        nsize = nmax - nmin + 1

        ndl = g%nb_diag(1)                  ! number of diagonals in lhs
        idl = ndl/2 + 1
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1

        do n = 1, g%size
            u(1, :) = 0.0_wp                ! Create delta-function forcing term
            u(1, n) = 1.0_wp

            z(1, 1) = u(1, 1)

            ! -------------------------------------------------------------------
            ! Calculate RHS in system of equations A u' = B u
            call g%matmul(g%rhs, u, z, ibc, g%rhs_b, g%rhs_t, bcs_hb, bcs_ht)

            ! -------------------------------------------------------------------
            ! Solve for u' in system of equations A u' = B u
            ip = ibc*5
            select case (ndl)
            case (3)
                call Thomas3_SolveL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
                call Thomas3_SolveU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))
                ! call Thomas3_SolveLU(nsize, 1, &
                !                      g%lu(nmin:nmax, ip + 1), &
                !                      g%lu(nmin:nmax, ip + 2), &
                !                      g%lu(nmin:nmax, ip + 3), &
                !                      z(:, nmin:nmax))
            case (5)
                call Thomas5_SolveL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
                call Thomas5_SolveU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))
                ! call Thomas5_SolveLU(nsize, 1, &
                !                      g%lu(nmin:nmax, ip + 1), &
                !                      g%lu(nmin:nmax, ip + 2), &
                !                      g%lu(nmin:nmax, ip + 3), &
                !                      g%lu(nmin:nmax, ip + 4), &
                !                      g%lu(nmin:nmax, ip + 5), &
                !                      z(:, nmin:nmax))
            end select

            do ic = 1, idl - 1
                bcs_hb(1) = bcs_hb(1) + g%lu(1, ip + idl + ic)*z(1, 1 + ic)
            end do
            c_b(n) = bcs_hb(1)/g%rhs(1, idr)

            if (abs(c_b(n)/c_b(1)) < roundoff_wp) exit

        end do
        n_bcs = min(n, g%size)

        return
    end subroutine FDM_Der1_NeumannMin_Initialize

! ###################################################################
! ###################################################################
! Truncated version for BCS_DN
    subroutine FDM_Der1_NeumannMax_Initialize(g, c_t, u, z, n_bcs)
        type(fdm_derivative_dt), intent(in) :: g
        real(wp), intent(out) :: c_t(g%size)                ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, g%size)             ! Working arrays
        real(wp), intent(inout) :: z(1, g%size)
        integer, intent(out) :: n_bcs                       ! Index of truncation

        ! -------------------------------------------------------------------
        integer ndl, idl, ndr, idr, ic
        integer ibc
        integer(wi) nmin, nmax, nsize, n, ip
        real(wp) bcs_hb(1), bcs_ht(1)

        ! ###################################################################
        ibc = BCS_DN

        c_t(:) = 0.0_wp

        nmin = 1; nmax = g%size - 1
        nsize = nmax - nmin + 1

        ndl = g%nb_diag(1)                  ! number of diagonals in lhs
        idl = ndl/2 + 1
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1

        do n = g%size, 1, -1
            u(1, :) = 0.0_wp                ! Create delta-function forcing term
            u(1, n) = 1.0_wp

            z(1, g%size) = u(1, g%size)

            ! -------------------------------------------------------------------
            ! Calculate RHS in system of equations A u' = B u
            call g%matmul(g%rhs, u, z, ibc, g%rhs_b, g%rhs_t, bcs_hb, bcs_ht)

            ! -------------------------------------------------------------------
            ! Solve for u' in system of equations A u' = B u
            ip = ibc*5
            select case (ndl)
            case (3)
                call Thomas3_SolveL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
                call Thomas3_SolveU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))
                ! call Thomas3_SolveLU(nsize, 1, &
                !                      g%lu(nmin:nmax, ip + 1), &
                !                      g%lu(nmin:nmax, ip + 2), &
                !                      g%lu(nmin:nmax, ip + 3), &
                !                      z(:, nmin:nmax))
            case (5)
                call Thomas5_SolveL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
                call Thomas5_SolveU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))
                ! call Thomas5_SolveLU(nsize, 1, &
                !                      g%lu(nmin:nmax, ip + 1), &
                !                      g%lu(nmin:nmax, ip + 2), &
                !                      g%lu(nmin:nmax, ip + 3), &
                !                      g%lu(nmin:nmax, ip + 4), &
                !                      g%lu(nmin:nmax, ip + 5), &
                !                      z(:, nmin:nmax))
            end select

            do ic = 1, idl - 1
                bcs_ht(1) = bcs_ht(1) + g%lu(g%size, ip + idl - ic)*z(1, g%size - ic)
            end do
            c_t(n) = bcs_ht(1)/g%rhs(g%size, idr)

            if (abs(c_t(n)/c_t(g%size)) < roundoff_wp) exit

        end do
        n_bcs = min(g%size - n + 1, g%size)

        return
    end subroutine FDM_Der1_NeumannMax_Initialize

end module FDM_derivative_Neumann

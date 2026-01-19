! Coefficients in explicit formula to calculate boundary value in Neumann conditions
module FDM_derivative_Neumann
    use TLab_Constants, only: wp, wi, roundoff_wp
    use TLab_Constants, only: BCS_ND, BCS_DN, BCS_NN
    use FDM_Derivative_1order_X, only: der1_biased
    ! use FDM_Derivative
    implicit none
    private

    public :: FDM_Der1_Neumann_Initialize           ! General, full version
    public :: FDM_Der1_NeumannMin_Initialize        ! Truncated versions
    public :: FDM_Der1_NeumannMax_Initialize

contains
    ! ###################################################################
    ! ###################################################################
    ! subroutine FDM_Der1_Neumann_Initialize(ibc, g, c_b, c_t, u, z)
    !     use FDM, only: fdm_der1_Z
    !     use FDM_Derivative_1order_X, only: der1_biased
    !     integer, intent(in) :: ibc                          ! Boundary condition [BCS_ND, BCS_DN, BCS_NN]
    !     type(fdm_derivative_dt), intent(in) :: g
    !     real(wp), intent(out) :: c_b(g%size), c_t(g%size)   ! coefficients for bottom value and top value
    !     real(wp), intent(inout) :: u(1, g%size)             ! Working arrays
    !     real(wp), intent(inout) :: z(1, g%size)

    !     ! -------------------------------------------------------------------
    !     integer ndl, idl, ndr, idr, ic
    !     integer(wi) nmin, nmax, nsize, n, ip
    !     real(wp) bcs_hb(1), bcs_ht(1)

    !     ! ###################################################################
    !     ndl = g%nb_diag(1)                  ! number of diagonals in lhs
    !     idl = ndl/2 + 1
    !     ndr = g%nb_diag(2)
    !     idr = ndr/2 + 1

    !     ip = ibc*5

    !     do n = 1, g%size
    !         u(1, :) = 0.0_wp                ! Create delta-function forcing term
    !         u(1, n) = 1.0_wp

    !         nmin = 1; nmax = g%size
    !         if (any([BCS_ND, BCS_NN] == ibc)) then
    !             bcs_hb(1) = u(1, 1)
    !             nmin = nmin + 1
    !         end if
    !         if (any([BCS_DN, BCS_NN] == ibc)) then
    !             bcs_ht(1) = u(1, g%size)
    !             nmax = nmax - 1
    !         end if
    !         nsize = nmax - nmin + 1

    !         ! -------------------------------------------------------------------
    !         ! Calculate RHS in system of equations A u' = B u
    !         select case (ibc)
    !         case (BCS_ND)
    !             ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
    !             !                    rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
    !             !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
    !             !                    u=u, &
    !             !                    f=z, bcs_b=bcs_hb(:))
    !             call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
    !                                  rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
    !                                  rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
    !                                  u=u, &
    !                                  L=g%lu(:, ip + 1:ip + ndl/2), &
    !                                  f=z, bcs_b=bcs_hb(:))
    !         case (BCS_DN)
    !             ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
    !             !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
    !             !                    rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
    !             !                    u=u, &
    !             !                    f=z, bcs_t=bcs_ht(:))
    !             call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
    !                                  rhs_b=g%rhs(1:ndr/2, 1:ndr), &
    !                                  rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
    !                                  u=u, &
    !                                  L=g%lu(:, ip + 1:ip + ndl/2), &
    !                                  f=z, bcs_t=bcs_ht(:))
    !         case (BCS_NN)
    !             ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
    !             !                    rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
    !             !                    rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
    !             !                    u=u, &
    !             !                    f=z, bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))
    !             call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
    !                                  rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
    !                                  rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
    !                                  u=u, &
    !                                  L=g%lu(:, ip + 1:ip + ndl/2), &
    !                                  f=z, bcs_b=bcs_hb(:), bcs_t=bcs_ht(:))

    !         end select
    !         ! -------------------------------------------------------------------
    !         ! Solve for u' in system of equations A u' = B u
    !         ! call g%thomasL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
    !         call g%thomasU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))

    !         if (any([BCS_ND, BCS_NN] == ibc)) then
    !             do ic = 1, idl - 1
    !                 bcs_hb(1) = bcs_hb(1) + g%lu(1, ip + idl + ic)*z(1, 1 + ic)
    !             end do
    !             c_b(n) = bcs_hb(1)/g%rhs(1, idr)
    !         end if
    !         if (any([BCS_DN, BCS_NN] == ibc)) then
    !             do ic = 1, idl - 1
    !                 bcs_ht(1) = bcs_ht(1) + g%lu(g%size, ip + idl - ic)*z(1, g%size - ic)
    !             end do
    !             c_t(n) = bcs_ht(1)/g%rhs(g%size, idr)
    !         end if

    !     end do

    !     return
    ! end subroutine FDM_Der1_Neumann_Initialize

    subroutine FDM_Der1_Neumann_Initialize(ibc, der, c_b, c_t, u, z)
        integer, intent(in) :: ibc                          ! Boundary condition [BCS_ND, BCS_DN, BCS_NN]
        type(der1_biased), intent(in) :: der
        real(wp), intent(out) :: c_b(size(der%lhs, 1)), c_t(size(der%lhs, 1))   ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, size(der%lhs, 1))             ! Working arrays
        real(wp), intent(inout) :: z(1, size(der%lhs, 1))

        ! -------------------------------------------------------------------
        integer n, nx
        real(wp) bcs_b(1), bcs_t(1)

        ! ###################################################################
        nx = size(der%lhs, 1)
        do n = 1, nx
            u(1, :) = 0.0_wp                ! Create delta-function forcing term
            u(1, n) = 1.0_wp

            ! select type (der)
            ! type is (der1_biased)
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
            ! end select

        end do

        return
    end subroutine FDM_Der1_Neumann_Initialize

! ###################################################################
! ###################################################################
! Truncated version for BCS_ND
    ! subroutine FDM_Der1_NeumannMin_Initialize(g, c_b, u, z, n_bcs)
    !     type(fdm_derivative_dt), intent(in) :: g
    !     real(wp), intent(out) :: c_b(g%size)                ! coefficients for bottom value and top value
    !     real(wp), intent(inout) :: u(1, g%size)             ! Working arrays
    !     real(wp), intent(inout) :: z(1, g%size)
    !     integer, intent(out) :: n_bcs                       ! Index of truncation

    !     ! -------------------------------------------------------------------
    !     integer ndl, idl, ndr, idr, ic
    !     integer ibc
    !     integer(wi) nmin, nmax, nsize, n, ip
    !     real(wp) bcs_hb(1)!, bcs_ht(1)

    !     ! ###################################################################
    !     ibc = BCS_ND

    !     c_b(:) = 0.0_wp

    !     nmin = 2; nmax = g%size
    !     nsize = nmax - nmin + 1

    !     ndl = g%nb_diag(1)                  ! number of diagonals in lhs
    !     idl = ndl/2 + 1
    !     ndr = g%nb_diag(2)
    !     idr = ndr/2 + 1

    !     ip = ibc*5

    !     do n = 1, g%size
    !         u(1, :) = 0.0_wp                ! Create delta-function forcing term
    !         u(1, n) = 1.0_wp

    !         bcs_hb(1) = u(1, 1)

    !         ! -------------------------------------------------------------------
    !         ! Calculate RHS in system of equations A u' = B u
    !         ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
    !         !                    rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
    !         !                    rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
    !         !                    u=u, &
    !         !                    f=z, bcs_b=bcs_hb(:))
    !         call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
    !                              rhs_b=g%rhs_b1(1:max(idl, idr + 1), 1:ndr + 2), &
    !                              rhs_t=g%rhs(g%size - ndr/2 + 1:g%size, 1:ndr), &
    !                              u=u, &
    !                              L=g%lu(:, ip + 1:ip + ndl/2), &
    !                              f=z, bcs_b=bcs_hb(:))

    !         ! -------------------------------------------------------------------
    !         ! Solve for u' in system of equations A u' = B u
    !         ! call g%thomasL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
    !         call g%thomasU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))

    !         do ic = 1, idl - 1
    !             bcs_hb(1) = bcs_hb(1) + g%lu(1, ip + idl + ic)*z(1, 1 + ic)
    !         end do
    !         c_b(n) = bcs_hb(1)/g%rhs(1, idr)

    !         if (abs(c_b(n)/c_b(1)) < roundoff_wp) exit

    !     end do
    !     n_bcs = min(n, g%size)

    !     return
    ! end subroutine FDM_Der1_NeumannMin_Initialize

    subroutine FDM_Der1_NeumannMin_Initialize(der, c_b, u, z, n_bcs)
        type(der1_biased), intent(in) :: der
        real(wp), intent(out) :: c_b(size(der%lhs, 1))                ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, size(der%lhs, 1))             ! Working arrays
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

            ! select type (der)
            ! type is (der1_biased)
            bcs_b(1) = u(1, 1)
            call der%bcsND%compute(1, u, z, bcs_b)
            c_b(n) = bcs_b(1)
            ! end select
            if (abs(c_b(n)/c_b(1)) < roundoff_wp) exit

        end do
        n_bcs = min(n, nx)

        return
    end subroutine FDM_Der1_NeumannMin_Initialize

! ###################################################################
! ###################################################################
! Truncated version for BCS_DN
    ! subroutine FDM_Der1_NeumannMax_Initialize(g, c_t, u, z, n_bcs)
    !     type(fdm_derivative_dt), intent(in) :: g
    !     real(wp), intent(out) :: c_t(g%size)                ! coefficients for bottom value and top value
    !     real(wp), intent(inout) :: u(1, g%size)             ! Working arrays
    !     real(wp), intent(inout) :: z(1, g%size)
    !     integer, intent(out) :: n_bcs                       ! Index of truncation

    !     ! -------------------------------------------------------------------
    !     integer ndl, idl, ndr, idr, ic
    !     integer ibc
    !     integer(wi) nmin, nmax, nsize, n, ip
    !     real(wp) bcs_ht(1)

    !     ! ###################################################################
    !     ibc = BCS_DN

    !     c_t(:) = 0.0_wp

    !     nmin = 1; nmax = g%size - 1
    !     nsize = nmax - nmin + 1

    !     ndl = g%nb_diag(1)                  ! number of diagonals in lhs
    !     idl = ndl/2 + 1
    !     ndr = g%nb_diag(2)
    !     idr = ndr/2 + 1

    !     ip = ibc*5

    !     do n = g%size, 1, -1
    !         u(1, :) = 0.0_wp                ! Create delta-function forcing term
    !         u(1, n) = 1.0_wp

    !         bcs_ht(1) = u(1, g%size)

    !         ! -------------------------------------------------------------------
    !         ! Calculate RHS in system of equations A u' = B u
    !         ! call g%matmul(rhs=g%rhs(:, 1:ndr), &
    !         !                    rhs_b=g%rhs(1:ndr/2, 1:ndr), &
    !         !                    rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
    !         !                    u=u, &
    !         !                    f=z, bcs_t=bcs_ht(:))
    !         call g%matmul_thomas(rhs=g%rhs(:, 1:ndr), &
    !                              rhs_b=g%rhs(1:ndr/2, 1:ndr), &
    !                              rhs_t=g%rhs_t1(1:max(idl, idr + 1), 1:ndr + 2), &
    !                              u=u, &
    !                              L=g%lu(:, ip + 1:ip + ndl/2), &
    !                              f=z, bcs_t=bcs_ht(:))

    !         ! -------------------------------------------------------------------
    !         ! Solve for u' in system of equations A u' = B u
    !         ! call g%thomasL(g%lu(nmin:nmax, ip + 1:ip + ndl/2), z(:, nmin:nmax))
    !         call g%thomasU(g%lu(nmin:nmax, ip + ndl/2 + 1:ip + ndl), z(:, nmin:nmax))

    !         do ic = 1, idl - 1
    !             bcs_ht(1) = bcs_ht(1) + g%lu(g%size, ip + idl - ic)*z(1, g%size - ic)
    !         end do
    !         c_t(n) = bcs_ht(1)/g%rhs(g%size, idr)

    !         if (abs(c_t(n)/c_t(g%size)) < roundoff_wp) exit

    !     end do
    !     n_bcs = min(g%size - n + 1, g%size)

    !     return
    ! end subroutine FDM_Der1_NeumannMax_Initialize

    subroutine FDM_Der1_NeumannMax_Initialize(der, c_t, u, z, n_bcs)
        type(der1_biased), intent(in) :: der
        real(wp), intent(out) :: c_t(size(der%lhs, 1))                ! coefficients for bottom value and top value
        real(wp), intent(inout) :: u(1, size(der%lhs, 1))             ! Working arrays
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

            ! select type (der)
            ! type is (der1_biased)
            bcs_t(1) = u(1, nx)
            call der%bcsDN%compute(1, u, z, bcs_t)
            c_t(n) = bcs_t(1)
            ! end select
            if (abs(c_t(n)/c_t(nx)) < roundoff_wp) exit

        end do
        n_bcs = min(nx - n + 1, nx)

        return
    end subroutine FDM_Der1_NeumannMax_Initialize

end module FDM_derivative_Neumann
